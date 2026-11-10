
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 16
#define K 5
#define R 11
#define WORDS 32
#define MAX_ITER 20

typedef int Bit;

static Bit GR[K][N], GB[K][N];
static Bit HR[R][N], HB[R][N];

static Bit rm_code[WORDS][N];
static Bit bch_code[WORDS][N];

static Bit info[K][K];
static Bit codeword[N][N];
static Bit error[N][N];
static Bit received[N][N];
static Bit work[N][N];

static int popcount(unsigned int x)
{
    int n = 0;
    while (x) {
        n += x & 1u;
        x >>= 1;
    }
    return n;
}

/* ==========================================
   RM(1,4) encoder [16,5,8]
   ========================================== */

static void rm_encode(int msg, Bit out[N])
{
    for (int x = 0; x < N; x++) {
        int v = msg & 1;

        for (int j = 0; j < 4; j++)
            v ^= ((msg >> (j+1)) & 1) &
                 ((x >> j) & 1);

        out[x] = v;
    }
}

/* ==========================================
   BCH(15,5,7) + overall parity
   Extended BCH [16,5,8]

   g(x) = x^10+x^8+x^6+x^5+x^4+x^2+x+1
   ========================================== */

static void bch_encode(int msg, Bit out[N])
{
    unsigned int c = 0;
    const unsigned int g = 0x537;

    for (int i = 0; i < K; i++)
        if ((msg >> i) & 1)
            c ^= g << i;

    c |= (popcount(c) & 1u) << 15;

    for (int i = 0; i < N; i++)
        out[i] = (c >> i) & 1u;
}

/* ==========================================
   Generate generator matrices
   ========================================== */

static void make_generators(void)
{
    for (int i = 0; i < K; i++) {
        rm_encode(1 << i, GR[i]);
        bch_encode(1 << i, GB[i]);
    }

    for (int m = 0; m < WORDS; m++) {
        rm_encode(m, rm_code[m]);
        bch_encode(m, bch_code[m]);
    }
}

/* ==========================================
   Generate parity check matrix from G

   H is a basis of the nullspace of G.
   ========================================== */

static int make_parity_check(
    Bit G[K][N],
    Bit H[R][N])
{
    Bit A[K][N];
    int pivot[N];
    int rank = 0;

    memcpy(A, G, sizeof(A));

    for (int col = 0; col < N && rank < K; col++) {

        int p = -1;

        for (int i = rank; i < K; i++) {
            if (A[i][col]) {
                p = i;
                break;
            }
        }

        if (p < 0)
            continue;

        for (int j = 0; j < N; j++) {
            Bit t = A[rank][j];
            A[rank][j] = A[p][j];
            A[p][j] = t;
        }

        for (int i = 0; i < K; i++) {
            if (i != rank && A[i][col]) {
                for (int j = 0; j < N; j++)
                    A[i][j] ^= A[rank][j];
            }
        }

        pivot[rank] = col;
        rank++;
    }

    if (rank != K)
        return -1;

    int hrow = 0;

    for (int f = 0; f < N; f++) {

        int is_pivot = 0;

        for (int i = 0; i < rank; i++)
            if (pivot[i] == f)
                is_pivot = 1;

        if (is_pivot)
            continue;

        for (int j = 0; j < N; j++)
            H[hrow][j] = 0;

        H[hrow][f] = 1;

        for (int i = 0; i < rank; i++)
            H[hrow][pivot[i]] = A[i][f];

        hrow++;
    }

    return hrow;
}

/* ==========================================
   Matrix display
   ========================================== */

static void print_matrix(
    const char *name,
    Bit *a,
    int rows,
    int cols)
{
    printf("\n%s (%d x %d)\n", name, rows, cols);

    for (int i = 0; i < rows; i++) {
        printf("%2d: ", i);

        for (int j = 0; j < cols; j++)
            printf("%d", a[i * cols + j]);

        printf("\n");
    }
}

static void print_16x16(
    const char *name,
    Bit a[N][N])
{
    printf("\n%s\n", name);

    for (int i = 0; i < N; i++) {
        printf("%2d: ", i);

        for (int j = 0; j < N; j++)
            printf("%d", a[i][j]);

        printf("\n");
    }
}

static void print_info(
    const char *name,
    Bit a[K][K])
{
    printf("\n%s\n", name);

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++)
            printf("%d", a[i][j]);

        printf("\n");
    }
}

/* ==========================================
   Verify G * H^T = 0
   ========================================== */

static int verify_orthogonality(
    Bit G[K][N],
    Bit H[R][N])
{
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < R; j++) {

            int s = 0;

            for (int k = 0; k < N; k++)
                s ^= G[i][k] & H[j][k];

            if (s)
                return 0;
        }
    }

    return 1;
}

/* ==========================================
   Minimum distance
   ========================================== */

static int min_distance(Bit code[WORDS][N])
{
    int dmin = N;

    for (int i = 0; i < WORDS; i++) {
        for (int j = i+1; j < WORDS; j++) {

            int d = 0;

            for (int k = 0; k < N; k++)
                d += code[i][k] != code[j][k];

            if (d < dmin)
                dmin = d;
        }
    }

    return dmin;
}

/* ==========================================
   Syndrome calculation
   ========================================== */

static void syndrome(
    Bit H[R][N],
    Bit v[N],
    Bit s[R])
{
    for (int i = 0; i < R; i++) {
        s[i] = 0;

        for (int j = 0; j < N; j++)
            s[i] ^= H[i][j] & v[j];
    }
}

/* ==========================================
   Product code syndromes

   Rows    : RM
   Columns : BCH
   ========================================== */

static void print_product_syndromes(
    const char *title,
    Bit a[N][N])
{
    printf("\n====================================\n");
    printf("%s\n", title);
    printf("====================================\n");

    printf("\nRM row syndromes:\n");

    for (int i = 0; i < N; i++) {
        Bit s[R];

        syndrome(HR, a[i], s);

        printf("Row %2d: ", i);

        for (int j = 0; j < R; j++)
            printf("%d", s[j]);

        printf("\n");
    }

    printf("\nBCH column syndromes:\n");

    for (int j = 0; j < N; j++) {
        Bit col[N], s[R];

        for (int i = 0; i < N; i++)
            col[i] = a[i][j];

        syndrome(HB, col, s);

        printf("Col %2d: ", j);

        for (int i = 0; i < R; i++)
            printf("%d", s[i]);

        printf("\n");
    }
}

/* ==========================================
   Nearest codeword decoder
   ========================================== */

static int nearest_decode(
    Bit input[N],
    Bit code[WORDS][N],
    Bit output[N],
    int *distance)
{
    int best = N + 1;
    int best_msg = 0;

    for (int m = 0; m < WORDS; m++) {

        int d = 0;

        for (int i = 0; i < N; i++)
            d += input[i] != code[m][i];

        if (d < best) {
            best = d;
            best_msg = m;
        }
    }

    memcpy(output, code[best_msg], sizeof(Bit) * N);

    *distance = best;

    return best_msg;
}

/* ==========================================
   Product encoder

   RM rows, then BCH columns.
   ========================================== */

static void encode_product(void)
{
    Bit temp[K][N];

    for (int i = 0; i < K; i++) {

        int msg = 0;

        for (int j = 0; j < K; j++)
            msg |= info[i][j] << j;

        rm_encode(msg, temp[i]);
    }

    for (int j = 0; j < N; j++) {

        int msg = 0;

        for (int i = 0; i < K; i++)
            msg |= temp[i][j] << i;

        Bit col[N];

        bch_encode(msg, col);

        for (int i = 0; i < N; i++)
            codeword[i][j] = col[i];
    }
}

/* ==========================================
   Iterative product decoder
   ========================================== */

static int decode_product(void)
{
    memcpy(work, received, sizeof(work));

    int iterations = 0;

    for (int iter = 0; iter < MAX_ITER; iter++) {

        int changed = 0;

        /* RM row decoding */

        for (int i = 0; i < N; i++) {

            Bit input[N], output[N];
            int d;

            for (int j = 0; j < N; j++)
                input[j] = work[i][j];

            nearest_decode(input, rm_code, output, &d);

            for (int j = 0; j < N; j++) {
                if (work[i][j] != output[j])
                    changed++;

                work[i][j] = output[j];
            }
        }

        /* BCH column decoding */

        for (int j = 0; j < N; j++) {

            Bit input[N], output[N];
            int d;

            for (int i = 0; i < N; i++)
                input[i] = work[i][j];

            nearest_decode(input, bch_code, output, &d);

            for (int i = 0; i < N; i++) {
                if (work[i][j] != output[i])
                    changed++;

                work[i][j] = output[i];
            }
        }

        iterations = iter + 1;

        if (changed == 0)
            break;
    }

    return iterations;
}

/* ==========================================
   Recover information bits
   ========================================== */

static void recover_information(Bit result[K][K])
{
    Bit temp[K][N];

    for (int j = 0; j < N; j++) {

        Bit col[N], output[N];
        int d;

        for (int i = 0; i < N; i++)
            col[i] = work[i][j];

        int msg = nearest_decode(
            col, bch_code, output, &d
        );

        for (int i = 0; i < K; i++)
            temp[i][j] = (msg >> i) & 1;
    }

    for (int i = 0; i < K; i++) {

        Bit output[N];
        int d;

        int msg = nearest_decode(
            temp[i], rm_code, output, &d
        );

        for (int j = 0; j < K; j++)
            result[i][j] = (msg >> j) & 1;
    }
}

/* ==========================================
   Main
   ========================================== */

int main(int argc, char **argv)
{
    int errors = 20;
    unsigned int seed = (unsigned int)time(NULL);

    if (argc > 1)
        errors = atoi(argv[1]);

    if (argc > 2)
        seed = (unsigned int)strtoul(argv[2], NULL, 10);

    if (errors < 0 || errors > 256) {
        printf("Errors must be 0..256\n");
        return 1;
    }

    srand(seed);

    make_generators();

    if (make_parity_check(GR, HR) != R ||
        make_parity_check(GB, HB) != R) {
        printf("Parity check construction failed\n");
        return 1;
    }

    printf("====================================\n");
    printf(" RM x Extended BCH Product Code\n");
    printf("====================================\n");

    printf("RM minimum distance  : %d\n",
           min_distance(rm_code));

    printf("BCH minimum distance : %d\n",
           min_distance(bch_code));

    printf("Product parameters   : [256,25,64]\n");

    printf("RM G*H^T = 0         : %s\n",
           verify_orthogonality(GR, HR) ? "PASS" : "FAIL");

    printf("BCH G*H^T = 0        : %s\n",
           verify_orthogonality(GB, HB) ? "PASS" : "FAIL");

    /* Display generator and parity-check matrices */

    print_matrix("RM Generator G", &GR[0][0], K, N);
    print_matrix("RM Parity Check H", &HR[0][0], R, N);

    print_matrix("BCH Generator G", &GB[0][0], K, N);
    print_matrix("BCH Parity Check H", &HB[0][0], R, N);

    /* Random information */

    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
            info[i][j] = rand() & 1;

    print_info("Original information", info);

    /* Encode */

    encode_product();

    print_16x16("Encoded codeword C", codeword);

    /* Generate random distinct errors */

    int positions[256];

    for (int i = 0; i < 256; i++)
        positions[i] = i;

    for (int i = 0; i < 256; i++) {
        int j = i + rand() % (256 - i);

        int t = positions[i];
        positions[i] = positions[j];
        positions[j] = t;
    }

    memset(error, 0, sizeof(error));

    for (int i = 0; i < errors; i++) {
        int p = positions[i];
        error[p / N][p % N] = 1;
    }

    /* R = C + E */

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            received[i][j] =
                codeword[i][j] ^ error[i][j];

    print_16x16("Error pattern E", error);
    print_16x16("Received word R = C + E", received);

    printf("\nInjected errors: %d\n", errors);

    /* Before decoding */

    print_product_syndromes(
        "SYNDROMES BEFORE DECODING",
        received
    );

    /* Decode */

    int iterations = decode_product();

    print_16x16("Decoded codeword", work);

    /* After decoding */

    print_product_syndromes(
        "SYNDROMES AFTER DECODING",
        work
    );

    /* Recover information */

    Bit recovered[K][K];

    recover_information(recovered);

    print_info("Recovered information", recovered);

    /* Compare */

    int info_ok = 1;
    int code_ok = 1;
    int remaining = 0;

    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
            if (info[i][j] != recovered[i][j])
                info_ok = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (codeword[i][j] != work[i][j])
                code_ok = 0;

            if (codeword[i][j] != received[i][j])
                remaining++;
        }
    }

    /* Note: remaining above counts errors in RECEIVED,
       not errors remaining after decoding. */

    int decoded_errors = 0;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            decoded_errors += codeword[i][j] != work[i][j];

    printf("\n====================================\n");
    printf("RESULT\n");
    printf("====================================\n");

    printf("Iterations              : %d\n", iterations);
    printf("Injected errors         : %d\n", errors);
    printf("Errors in received word : %d\n", remaining);
    printf("Errors after decoding   : %d\n", decoded_errors);

    printf("Information recovery    : %s\n",
           info_ok ? "PASS" : "FAIL");

    printf("Codeword recovery       : %s\n",
           code_ok ? "PASS" : "FAIL");

    printf("====================================\n");

    return (info_ok && code_ok) ? 0 : 1;
}

