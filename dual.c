
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define N 16
#define K 10
#define R 6
#define T 3

typedef uint8_t GF;

static GF gf_exp[30];
static GF gf_log[16];
static GF eval[N];
static GF weight[N];
static GF H[R][N];

/* ------------------------------------------------
   GF(16): primitive polynomial x^4+x+1 = 0x13
   ------------------------------------------------ */

static void gf_init(void)
{
    GF x = 1;

    for (int i = 0; i < 15; i++) {
        gf_exp[i] = x;
        gf_log[x] = i;

        x <<= 1;

        if (x & 16)
            x ^= 0x13;
    }

    for (int i = 15; i < 30; i++)
        gf_exp[i] = gf_exp[i - 15];

    gf_log[0] = 0;
}

static GF mul(GF a, GF b)
{
    if (!a || !b)
        return 0;

    return gf_exp[gf_log[a] + gf_log[b]];
}

static GF inv(GF a)
{
    if (!a)
        return 0;

    return gf_exp[15 - gf_log[a]];
}

static GF power(GF a, int n)
{
    if (n == 0)
        return 1;

    if (!a)
        return 0;

    return gf_exp[(gf_log[a] * n) % 15];
}

static GF div_gf(GF a, GF b)
{
    if (!a)
        return 0;

    if (!b)
        return 0;

    return mul(a, inv(b));
}

/* ------------------------------------------------
   Linear system solver over GF(16)
   ------------------------------------------------ */

static int solve(
    GF a[10][11],
    int n,
    GF *answer)
{
    for (int col = 0; col < n; col++) {

        int pivot = col;

        while (pivot < n && !a[pivot][col])
            pivot++;

        if (pivot == n)
            return 0;

        if (pivot != col) {
            for (int j = col; j <= n; j++) {
                GF t = a[col][j];
                a[col][j] = a[pivot][j];
                a[pivot][j] = t;
            }
        }

        GF scale = inv(a[col][col]);

        for (int j = col; j <= n; j++)
            a[col][j] = mul(a[col][j], scale);

        for (int i = 0; i < n; i++) {
            if (i == col)
                continue;

            GF factor = a[i][col];

            for (int j = col; j <= n; j++)
                a[i][j] ^= mul(factor, a[col][j]);
        }
    }

    for (int i = 0; i < n; i++)
        answer[i] = a[i][n];

    return 1;
}

/* ------------------------------------------------
   Extended Reed-Solomon [16,10,7]

   Evaluation code:
       c[i] = f(eval[i])
       deg(f) < 10

   All 16 GF(16) elements are evaluation points.

   Parity checks:
       sum weight[i]*c[i]*eval[i]^j = 0
       j = 0,...,5

   This is an extended GRS code, not a
   binary BCH code.
   ------------------------------------------------ */

static void init_code(void)
{
    for (int i = 0; i < N; i++)
        eval[i] = (GF)i;

    for (int i = 0; i < N; i++) {

        GF product = 1;

        for (int j = 0; j < N; j++) {
            if (i != j)
                product = mul(
                    product,
                    eval[i] ^ eval[j]
                );
        }

        weight[i] = inv(product);
    }

    for (int j = 0; j < R; j++) {
        for (int i = 0; i < N; i++) {
            H[j][i] = mul(
                weight[i],
                power(eval[i], j)
            );
        }
    }
}

/* ------------------------------------------------
   Encode 10 information symbols into 16 symbols.

   Interpolate polynomial through first 10 points,
   then evaluate it at all 16 points.
   ------------------------------------------------ */

static void encode(
    const GF msg[K],
    GF code[N])
{
    GF matrix[10][11] = {{0}};
    GF coeff[K];

    for (int i = 0; i < K; i++) {

        GF p = 1;

        for (int j = 0; j < K; j++) {
            matrix[i][j] = p;
            p = mul(p, eval[i]);
        }

        matrix[i][K] = msg[i];
    }

    if (!solve(matrix, K, coeff)) {
        printf("Encoding matrix error\n");
        exit(1);
    }

    for (int i = 0; i < N; i++) {

        GF y = 0;

        for (int j = K - 1; j >= 0; j--)
            y = mul(y, eval[i]) ^ coeff[j];

        code[i] = y;
    }
}

/* ------------------------------------------------
   Syndrome
   ------------------------------------------------ */

static void syndrome(
    const GF received[N],
    GF s[R])
{
    for (int j = 0; j < R; j++) {

        s[j] = 0;

        for (int i = 0; i < N; i++)
            s[j] ^= mul(H[j][i], received[i]);
    }
}

static int zero_syndrome(const GF s[R])
{
    for (int i = 0; i < R; i++)
        if (s[i])
            return 0;

    return 1;
}

/* ------------------------------------------------
   Decode by searching error locations.

   Maximum 3 errors.
   Solve error magnitudes from first e syndromes,
   then verify all six syndromes.

   This is a reference decoder, not BMS.
   ------------------------------------------------ */

static int try_positions(
    const GF s[R],
    const int pos[T],
    int e,
    GF errors[T])
{
    GF a[10][11] = {{0}};

    for (int row = 0; row < e; row++) {

        for (int col = 0; col < e; col++)
            a[row][col] = H[row][pos[col]];

        a[row][e] = s[row];
    }

    if (!solve(a, e, errors))
        return 0;

    for (int j = 0; j < R; j++) {

        GF check = 0;

        for (int i = 0; i < e; i++)
            check ^= mul(H[j][pos[i]], errors[i]);

        if (check != s[j])
            return 0;
    }

    return 1;
}

static int decode(
    GF received[N],
    int error_pos[T],
    GF error_mag[T])
{
    GF s[R];

    syndrome(received, s);

    if (zero_syndrome(s))
        return 0;

    int pos[T];
    GF errors[T];

    /* One error */
    for (int a = 0; a < N; a++) {

        pos[0] = a;

        if (try_positions(s, pos, 1, errors)) {

            error_pos[0] = a;
            error_mag[0] = errors[0];

            received[a] ^= errors[0];

            return 1;
        }
    }

    /* Two errors */
    for (int a = 0; a < N; a++) {
        for (int b = a + 1; b < N; b++) {

            pos[0] = a;
            pos[1] = b;

            if (try_positions(s, pos, 2, errors)) {

                for (int i = 0; i < 2; i++) {
                    error_pos[i] = pos[i];
                    error_mag[i] = errors[i];
                    received[pos[i]] ^= errors[i];
                }

                return 2;
            }
        }
    }

    /* Three errors */
    for (int a = 0; a < N; a++) {
        for (int b = a + 1; b < N; b++) {
            for (int c = b + 1; c < N; c++) {

                pos[0] = a;
                pos[1] = b;
                pos[2] = c;

                if (try_positions(s, pos, 3, errors)) {

                    for (int i = 0; i < 3; i++) {
                        error_pos[i] = pos[i];
                        error_mag[i] = errors[i];
                        received[pos[i]] ^= errors[i];
                    }

                    return 3;
                }
            }
        }
    }

    return -1;
}

/* ------------------------------------------------
   Product code encoder
   ------------------------------------------------ */

static void encode_product(
    GF info[K][K],
    GF code[N][N])
{
    GF temp[K][N];

    for (int i = 0; i < K; i++)
        encode(info[i], temp[i]);

    for (int j = 0; j < N; j++) {

        GF column[K];
        GF result[N];

        for (int i = 0; i < K; i++)
            column[i] = temp[i][j];

        encode(column, result);

        for (int i = 0; i < N; i++)
            code[i][j] = result[i];
    }
}

/* ------------------------------------------------
   Product code decoder
   ------------------------------------------------ */

static int decode_product(
    GF received[N][N],
    int verbose)
{
    int total = 0;

    for (int iteration = 0; iteration < 20; iteration++) {

        int changes = 0;

        /* Row decoding */
        for (int i = 0; i < N; i++) {

            int pos[T];
            GF mag[T];

            int n = decode(received[i], pos, mag);

            if (n > 0) {
                changes += n;
                total += n;

                if (verbose) {
                    printf("Row %2d: corrected %d symbols\n",
                           i, n);
                }
            }
        }

        /* Column decoding */
        for (int j = 0; j < N; j++) {

            GF column[N];

            for (int i = 0; i < N; i++)
                column[i] = received[i][j];

            int pos[T];
            GF mag[T];

            int n = decode(column, pos, mag);

            if (n > 0) {

                for (int i = 0; i < N; i++)
                    received[i][j] = column[i];

                changes += n;
                total += n;

                if (verbose) {
                    printf("Col  %2d: corrected %d symbols\n",
                           j, n);
                }
            }
        }

        if (verbose)
            printf("Iteration %d: %d corrections\n",
                   iteration + 1, changes);

        if (changes == 0)
            break;
    }

    return total;
}

/* ------------------------------------------------
   Display
   ------------------------------------------------ */

static void print_matrix(
    const char *name,
    GF a[N][N])
{
    printf("\n%s\n", name);

    for (int i = 0; i < N; i++) {
        printf("%2d: ", i);

        for (int j = 0; j < N; j++)
            printf("%X ", a[i][j]);

        printf("\n");
    }
}

/* ------------------------------------------------
   Main
   ------------------------------------------------ */

int main(void)
{
    GF info[K][K];
    GF original[N][N] = {{0}};
    GF received[N][N] = {{0}};
    GF error[N][N] = {{0}};

    int ne;
    unsigned int seed;

    gf_init();
    init_code();

    printf("GF(16) Extended GRS Product Code\n");
    printf("[256,100,49] over GF(16)\n");
    printf("Information: 400 bits\n");
    printf("Codeword: 1024 bits\n");
    printf("Correction: up to 3 symbols per row/column\n");

    printf("\nNumber of errors (0-256): ");

    if (scanf("%d", &ne) != 1 || ne < 0 || ne > 256)
        return 1;

    printf("Random seed: ");

    if (scanf("%u", &seed) != 1)
        return 1;

    srand(seed);

    /* Information */
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
            info[i][j] = rand() & 15;

    encode_product(info, original);

    memcpy(received, original, sizeof(original));

    /* Random error positions */
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
        int r = p / N;
        int c = p % N;

        GF e = (GF)(1 + rand() % 15);

        received[r][c] ^= e;
        error[r][c] = e;
    }

    print_matrix("Original codeword", original);
    print_matrix("Error values", error);
    print_matrix("Received word", received);

    printf("\n--- Decoding ---\n");

    int corrected = decode_product(received, 1);

    print_matrix("Decoded word", received);

    int residual = 0;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (received[i][j] != original[i][j])
                residual++;

    printf("\n========== RESULT ==========\n");
    printf("Injected errors : %d\n", ne);
    printf("Corrections     : %d\n", corrected);
    printf("Residual errors : %d\n", residual);
    printf("Result          : %s\n",
           residual == 0 ? "PASS" : "FAIL");
    printf("============================\n");

    return 0;
}
