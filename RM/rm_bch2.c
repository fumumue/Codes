
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define N 16
#define K 5
#define M 32
#define MAX_ITER 30

typedef uint16_t Word;

static Word RM[M];
static Word BCH[M];

static Word H_RM[11];
static Word H_BCH[11];

static int parity(Word x)
{
    int p = 0;
    while (x) {
        p ^= 1;
        x &= (Word)(x - 1);
    }
    return p;
}

static int popcount16(Word x)
{
    int n = 0;
    while (x) {
        n++;
        x &= (Word)(x - 1);
    }
    return n;
}

/* --------------------------------------------------
   RM(1,4) [16,5,8]
   f(x0,x1,x2,x3) = a0 + a1*x0 + ... + a4*x3
   -------------------------------------------------- */

static Word rm_encode(int msg)
{
    Word w = 0;

    for (int x = 0; x < 16; x++) {
        int bit = (msg & 1);

        for (int j = 0; j < 4; j++)
            bit ^= ((msg >> (j + 1)) & 1)
                   & ((x >> j) & 1);

        if (bit)
            w |= (Word)(1u << x);
    }

    return w;
}

/* --------------------------------------------------
   Extended binary BCH [16,5,8]

   BCH(15,5,7):
   g(x) = x^10+x^8+x^5+x^4+x^2+x+1

   g(x) = 0x537
   Position 15: overall even parity
   -------------------------------------------------- */

static Word bch_encode(int msg)
{
    Word w = 0;

    for (int i = 0; i < 5; i++) {
        if ((msg >> i) & 1)
            w ^= (Word)(0x537u << i);
    }

    /* Extend to even parity */
    if (parity(w))
        w |= (Word)(1u << 15);

    return w;
}

/* --------------------------------------------------
   Generate all 32 codewords
   -------------------------------------------------- */

static void init_codes(void)
{
    for (int i = 0; i < M; i++) {
        RM[i] = rm_encode(i);
        BCH[i] = bch_encode(i);
    }
}

/* --------------------------------------------------
   Construct parity-check basis.

   Find 11 independent vectors orthogonal to
   the five generator codewords.
   -------------------------------------------------- */

static int make_parity_check(
    Word code[M],
    Word H[11])
{
    Word basis[16] = {0};
    int count = 0;

    for (int v = 1; v < 65536 && count < 11; v++) {

        Word h = (Word)v;
        int valid = 1;

        for (int i = 0; i < K; i++) {
            if (parity(h & code[1 << i])) {
                valid = 0;
                break;
            }
        }

        if (!valid)
            continue;

        Word t = h;

        for (int b = 15; b >= 0; b--) {
            if (!(t & (1u << b)))
                continue;

            if (basis[b])
                t ^= basis[b];
            else {
                basis[b] = t;
                H[count++] = h;
                break;
            }
        }
    }

    return count;
}

/* --------------------------------------------------
   Syndrome calculation: 11 bits
   -------------------------------------------------- */

static unsigned int syndrome(
    Word received,
    Word H[11])
{
    unsigned int s = 0;

    for (int i = 0; i < 11; i++) {
        if (parity(received & H[i]))
            s |= (1u << i);
    }

    return s;
}

/* --------------------------------------------------
   Exhaustive nearest-codeword decoder
   -------------------------------------------------- */

static Word decode_line(
    Word received,
    Word code[M],
    int *distance)
{
    int best = 17;
    Word result = 0;

    for (int i = 0; i < M; i++) {
        int d = popcount16(received ^ code[i]);

        if (d < best) {
            best = d;
            result = code[i];
        }

        if (best == 0)
            break;
    }

    *distance = best;
    return result;
}

/* --------------------------------------------------
   Matrix display
   -------------------------------------------------- */

static void print_matrix(
    const char *name,
    int a[N][N])
{
    printf("\n%s\n", name);

    for (int i = 0; i < N; i++) {
        printf("%2d: ", i);

        for (int j = 0; j < N; j++)
            printf("%d", a[i][j]);

        printf("\n");
    }
}

/* --------------------------------------------------
   Convert a row into a 16-bit word
   -------------------------------------------------- */

static Word get_row(int a[N][N], int r)
{
    Word w = 0;

    for (int j = 0; j < N; j++)
        if (a[r][j])
            w |= (Word)(1u << j);

    return w;
}

/* --------------------------------------------------
   Convert a column into a 16-bit word
   -------------------------------------------------- */

static Word get_col(int a[N][N], int c)
{
    Word w = 0;

    for (int i = 0; i < N; i++)
        if (a[i][c])
            w |= (Word)(1u << i);

    return w;
}

static void set_row(
    int a[N][N], int r, Word w)
{
    for (int j = 0; j < N; j++)
        a[r][j] = (w >> j) & 1;
}

static void set_col(
    int a[N][N], int c, Word w)
{
    for (int i = 0; i < N; i++)
        a[i][c] = (w >> i) & 1;
}

/* --------------------------------------------------
   Product-code encoder

   1. Encode five information rows with RM.
   2. Encode each column with extended BCH.
   -------------------------------------------------- */

static void encode_product(
    int info[5][5],
    int codeword[N][N])
{
    int temp[5][N] = {{0}};

    for (int i = 0; i < 5; i++) {
        int msg = 0;

        for (int j = 0; j < 5; j++)
            msg |= info[i][j] << j;

        Word w = RM[msg];

        for (int j = 0; j < N; j++)
            temp[i][j] = (w >> j) & 1;
    }

    for (int j = 0; j < N; j++) {
        int msg = 0;

        for (int i = 0; i < 5; i++)
            msg |= temp[i][j] << i;

        Word w = BCH[msg];
        set_col(codeword, j, w);
    }
}

/* --------------------------------------------------
   Add random errors
   -------------------------------------------------- */

static void add_errors(
    int received[N][N],
    int error[N][N],
    int ne)
{
    int positions[256];

    for (int i = 0; i < 256; i++)
        positions[i] = i;

    for (int i = 255; i > 0; i--) {
        int j = rand() % (i + 1);
        int t = positions[i];
        positions[i] = positions[j];
        positions[j] = t;
    }

    for (int i = 0; i < ne; i++) {
        int p = positions[i];
        int r = p / 16;
        int c = p % 16;

        error[r][c] = 1;
        received[r][c] ^= 1;
    }
}

/* --------------------------------------------------
   Display all row and column syndromes
   -------------------------------------------------- */

static void print_syndromes(
    const char *name,
    int a[N][N])
{
    printf("\n%s\n", name);

    printf("RM row syndromes:\n");

    for (int i = 0; i < N; i++) {
        unsigned int s = syndrome(get_row(a, i), H_RM);
        printf("%2d: %03X  ", i, s);

        if (i % 4 == 3)
            printf("\n");
    }

    printf("\nBCH column syndromes:\n");

    for (int j = 0; j < N; j++) {
        unsigned int s = syndrome(get_col(a, j), H_BCH);
        printf("%2d: %03X  ", j, s);

        if (j % 4 == 3)
            printf("\n");
    }

    printf("\n");
}

/* --------------------------------------------------
   Alternating product-code decoder
   -------------------------------------------------- */

static int decode_product(
    int received[N][N],
    int decoded[N][N])
{
    memcpy(decoded, received, sizeof(int) * 256);

    int iterations = 0;

    for (int iter = 0; iter < MAX_ITER; iter++) {

        int changes = 0;

        /* RM decoding: row direction */
        for (int i = 0; i < N; i++) {
            Word old = get_row(decoded, i);
            int d;

            Word nw = decode_line(old, RM, &d);

            if (nw != old) {
                set_row(decoded, i, nw);
                changes++;
            }
        }

        /* BCH decoding: column direction */
        for (int j = 0; j < N; j++) {
            Word old = get_col(decoded, j);
            int d;

            Word nw = decode_line(old, BCH, &d);

            if (nw != old) {
                set_col(decoded, j, nw);
                changes++;
            }
        }

        iterations = iter + 1;

        printf("Iteration %2d: changes = %d\n",
               iterations, changes);

        if (changes == 0)
            break;
    }

    return iterations;
}

/* --------------------------------------------------
   Main
   -------------------------------------------------- */

int main(void)
{
    int info[5][5];
    int codeword[N][N] = {{0}};
    int received[N][N] = {{0}};
    int error[N][N] = {{0}};
    int decoded[N][N] = {{0}};

    int ne;
    unsigned int seed;

    init_codes();

    if (make_parity_check(RM, H_RM) != 11 ||
        make_parity_check(BCH, H_BCH) != 11) {
        printf("Parity-check construction failed.\n");
        return 1;
    }

    printf("RM(1,4) x Extended BCH(16,5,8)\n");
    printf("Product code: [256,25,64]\n");

    printf("Number of errors (0-256): ");
    if (scanf("%d", &ne) != 1 || ne < 0 || ne > 256)
        return 1;

    printf("Random seed: ");
    if (scanf("%u", &seed) != 1)
        return 1;

    srand(seed);

    /* Generate 25 information bits */
    printf("\nInformation matrix (5x5):\n");

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            info[i][j] = rand() & 1;
            printf("%d", info[i][j]);
        }
        printf("\n");
    }

    /* Encoding */
    encode_product(info, codeword);

    print_matrix("Original codeword", codeword);

    memcpy(received, codeword, sizeof(codeword));

    /* Add errors */
    add_errors(received, error, ne);

    print_matrix("Error matrix", error);
    print_matrix("Received word", received);

    /* Initial syndrome */
    print_syndromes("Initial product syndromes", received);

    /* Decoding */
    printf("\n--- Alternating decoding ---\n");

    int iterations = decode_product(received, decoded);

    print_matrix("Decoded word", decoded);

    /* Final syndrome */
    print_syndromes("Final product syndromes", decoded);

    /* Count errors and residual errors */
    int total_errors = 0;
    int residual = 0;
    int mismatch = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            total_errors += error[i][j];

            if (decoded[i][j] != codeword[i][j])
                residual++;

            if (decoded[i][j] != received[i][j])
                mismatch++;
        }
    }

    printf("\n========== RESULT ==========\n");
    printf("Injected errors : %d\n", total_errors);
    printf("Changed bits    : %d\n", mismatch);
    printf("Residual errors : %d\n", residual);
    printf("Iterations      : %d\n", iterations);

    if (residual == 0)
        printf("RESULT: PASS\n");
    else
        printf("RESULT: FAIL\n");

    printf("============================\n");

    return 0;
}
