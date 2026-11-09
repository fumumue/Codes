// Experimental Niederreiter-style round trip using the existing BM/Patterson code.
// This is a correctness experiment, not production cryptography: rand() is not
// cryptographically secure, and this demo does not encode arbitrary messages.
#define VMO_NO_MAIN
#include "vmo.c"

static short public_h[K][M];

static int generate_secret_key(vec *g)
{
    fflush(stdout);
    int saved_stdout = dup(STDOUT_FILENO);
    if (saved_stdout < 0)
    {
        perror("dup");
        return 0;
    }
    int sink = open("/dev/null", O_WRONLY);
    if (sink < 0)
    {
        perror("open /dev/null");
        close(saved_stdout);
        return 0;
    }
    if (dup2(sink, STDOUT_FILENO) < 0)
    {
        perror("dup2");
        close(sink);
        close(saved_stdout);
        return 0;
    }
    close(sink);

    mkmf();
    makefg();
    de();
    memset(mat, 0, sizeof(mat));
    *g = vv(K);

    fflush(stdout);
    if (dup2(saved_stdout, STDOUT_FILENO) < 0)
    {
        perror("restore stdout");
        close(saved_stdout);
        return 0;
    }
    close(saved_stdout);
    return 1;
}

static void random_invertible_matrix(short a[K][K], short a_inv[K][K])
{
    do
    {
        for (int i = 0; i < K; i++)
            for (int j = 0; j < K; j++)
                a[i][j] = (short)(rand() % O);
    } while (!q_matinv(a, a_inv));
}

static void make_permutation(int permutation[M])
{
    for (int i = 0; i < M; i++)
        permutation[i] = i;
    for (int i = M - 1; i > 0; i--)
    {
        int j = rand() % (i + 1);
        int tmp = permutation[i];
        permutation[i] = permutation[j];
        permutation[j] = tmp;
    }
}

static void make_public_matrix(const short a[K][K], const int permutation[M])
{
    for (int i = 0; i < K; i++)
        for (int j = 0; j < M; j++)
        {
            short sum = 0;
            for (int r = 0; r < K; r++)
                sum = q_add(sum, q_mul(a[i][r], mat[r][permutation[j]]));
            public_h[i][j] = sum;
        }
}

static void make_error(short error[M], int weight)
{
    int positions[M];
    memset(error, 0, sizeof(short) * M);
    for (int i = 0; i < M; i++)
        positions[i] = i;
    for (int i = 0; i < weight; i++)
    {
        int j = i + rand() % (M - i);
        int tmp = positions[i];
        positions[i] = positions[j];
        positions[j] = tmp;
        error[positions[i]] = 1;
    }
}

static void encrypt_syndrome(const short error[M], short ciphertext[K])
{
    for (int row = 0; row < K; row++)
    {
        short sum = 0;
        for (int col = 0; col < M; col++)
            if (error[col])
                sum = q_add(sum, q_mul(public_h[row][col], error[col]));
        ciphertext[K - 1 - row] = sum;
    }
}

static void remove_row_scramble(const short a_inv[K][K],
                                const short ciphertext[K],
                                short bm_syndrome[K])
{
    for (int row = 0; row < K; row++)
    {
        short sum = 0;
        for (int r = 0; r < K; r++)
            sum = q_add(sum, q_mul(a_inv[row][r], ciphertext[K - 1 - r]));
        bm_syndrome[K - 1 - row] = sum;
    }
}

static void print_syndrome(const char *label, const short *syndrome)
{
    printf("%s:", label);
    for (int i = 0; i < K; i++)
        printf(" %d", syndrome[i]);
    printf("\n");
}


static int bm_decode_errors(const short *syndrome, short decoded[M], int *locator_degree)
{
    short syndrome[K];
    for (int n = 0; n < K; n++)
        syndrome[n] = syn[K - 1 - n];   // T_n = syn[K-1-n]

    short locator[K + 1] = {1}, previous[K + 1] = {1};
    memset(decoded, 0, sizeof(short) * M);
    int degree = 0, shift = 1;
    short last_discrepancy = 1;

    for (int n = 0; n < K; n++)
    {
        short discrepancy = syndrome[n];
        for (int i = 1; i <= degree; i++)
            discrepancy = q_add(discrepancy,
                                q_mul(locator[i], syndrome[n - i]));
        if (discrepancy == 0)
        {
            shift++;
            continue;
        }

        short old_locator[K + 1];
        memcpy(old_locator, locator, sizeof(old_locator));
        short factor = q_mul(discrepancy, q_inv(last_discrepancy));
        for (int i = 0; i + shift <= K; i++)
            locator[i + shift] = q_sub(locator[i + shift],
                                       q_mul(factor, previous[i]));

        if (2 * degree <= n)
        {
            degree = n + 1 - degree;
            memcpy(previous, old_locator, sizeof(previous));
            last_discrepancy = discrepancy;
            shift = 1;
        }
        else
        {
            shift++;
        }
    }

    *locator_degree = degree;
    int roots = 0;
    for (int position = 0; position < M; position++)
    {
        short value = 0;
        for (int i = degree; i >= 0; i--)
            value = q_add(q_mul(value, (short)position), locator[i]);
        if (value == 0)
        {
            decoded[position] = 1;
            roots++;
        }
    }
    return roots;
}

int main(int argc, char **argv)
{
    static short a[K][K], a_inv[K][K];
    static short error[M], ciphertext[K], bm_syndrome[K], bm_decoded[M];
    static short decoded[Q_MAXC][M];
    static int permutation[M];

    int weight = argc > 1 ? atoi(argv[1]) : K;
    if (weight < 1 || weight > K)
    {
        fprintf(stderr, "error weight must be in [1, %d]\n", K);
        return 2;
    }
    int raw = argc > 2 && strcmp(argv[2], "raw") == 0;
    if (argc > 2 && !raw)
    {
        fprintf(stderr, "mode must be 'raw'\n");
        return 2;
    }

    srand((unsigned)time(NULL));
    vec g;
    if (!generate_secret_key(&g))
        return 1;
    if (raw)
    {
        for (int i = 0; i < K; i++)
            a[i][i] = a_inv[i][i] = 1;
        for (int i = 0; i < M; i++)
            permutation[i] = i;
    }
    else
    {
        random_invertible_matrix(a, a_inv);
        make_permutation(permutation);
    }
    make_public_matrix(a, permutation);
    make_error(error, weight);
    encrypt_syndrome(error, ciphertext);
    remove_row_scramble(a_inv, ciphertext, bm_syndrome);

    for (int row = 0; row < K; row++)
    {
        short expected = 0;
        for (int col = 0; col < M; col++)
            if (error[col])
                expected = q_add(expected, q_mul(mat[row][permutation[col]], error[col]));
        if (bm_syndrome[K - 1 - row] != expected)
        {
            fprintf(stderr, "row-scramble inversion failed at row %d\n", row);
            return 1;
        }
    }

    if (!q_setup(g.x))
    {
        fprintf(stderr, "generated g is not square-free\n");
        return 1;
    }
    q_res expected_patterson = q_zero();
    for (int col = 0; col < M; col++)
        if (error[col])
        {
            q_res factor = q_zero(), inverse;
            factor.c[0] = q_neg((short)permutation[col]);
            factor.c[1] = 1;
            if (!q_sinv(factor, &inverse))
            {
                fprintf(stderr, "error location is a root of g\n");
                return 1;
            }
            expected_patterson = q_radd(
                expected_patterson, q_rscale(inverse, error[col]));
        }
    q_res converted_patterson = q_bm_to_patterson_syndrome(bm_syndrome);
    for (int i = 0; i < K; i++)
        if (converted_patterson.c[i] != expected_patterson.c[i])
        {
            fprintf(stderr, "BM-to-Patterson conversion mismatch at coefficient %d\n", i);
            return 1;
        }

    printf("\nNiederreiter-style demo: K=%d, error weight=%d, mode=%s\n",
           K, weight, raw ? "raw (no A/P scramble)" : "scrambled");
    print_syndrome(raw ? "ciphertext syndrome (raw H)" : "ciphertext syndrome (scrambled)",
                   ciphertext);
    print_syndrome(raw ? "BM syndrome (raw H)" : "BM syndrome (after removing row scramble)",
                   bm_syndrome);
    printf("Patterson syndrome s(x) coefficients [x^i]:");
    for (int i = 0; i < K; i++)
        printf(" %d", converted_patterson.c[i]);
    printf("\n");

    int bm_locator_degree = 0;
    int bm_roots = bm_decode_errors(bm_syndrome, bm_decoded, &bm_locator_degree);
    int bm_recovered = bm_roots == weight;
    int bm_syndrome_matches = 1;
    for (int row = 0; row < K; row++)
    {
        short sum = 0;
        for (int position = 0; position < M; position++)
            if (bm_decoded[position])
                sum = q_add(sum, q_mul(mat[row][position], bm_decoded[position]));
        if (sum != bm_syndrome[K - 1 - row])
            bm_syndrome_matches = 0;
    }
    for (int col = 0; col < M; col++)
        if (bm_decoded[permutation[col]] != error[col])
            bm_recovered = 0;

    printf("\nBM decoder locator degree: %d (unique-decoding limit K/2=%d)\n",
           bm_locator_degree, K / 2);
    printf("BM locator roots (not necessarily valid errors): %d\n", bm_roots);
    printf("BM decoded support matches syndrome: %s\n",
           bm_syndrome_matches ? "yes" : "no");
    printf("BM original error recovered: %s\n", bm_recovered ? "yes" : "no");

    int candidates = patterson_q(bm_syndrome, (const short *)g.x, 1, Pr - 1, decoded);
    int recovered = 0;
    printf("\n| BM column | BM error | BM locator root | Pat location | Pat decoded |\n");
    printf("|-----------:|---------:|---------------:|-------------:|------------:|\n");
    if (candidates == 0)
    {
        for (int j = 0; j < M; j++)
            if (error[j])
                printf("| %9d | %8d | %14d | %12d | %11s |\n",
                       j, error[j], bm_decoded[permutation[j]], permutation[j], "なし");
    }
    for (int c = 0; c < candidates; c++)
    {
        int matches = 1;
        for (int j = 0; j < M; j++)
            if (decoded[c][permutation[j]] != error[j])
            {
                matches = 0;
                break;
            }
        recovered |= matches;

        for (int j = 0; j < M; j++)
        {
            short decoded_value = decoded[c][permutation[j]];
            if (error[j] || bm_decoded[permutation[j]] || decoded_value)
                printf("| %9d | %8d | %14d | %12d | %11d |%s\n",
                       j, error[j], bm_decoded[permutation[j]], permutation[j], decoded_value,
                       error[j] == decoded_value ? "" : " mismatch");
        }
        printf("candidate %d matches injected errors: %s\n",
               c + 1, matches ? "yes" : "no");
    }

    printf("Patterson candidates: %d\n", candidates);
    printf("original error recovered: %s\n", recovered ? "yes" : "no");

    return recovered ? 0 : 1;
}
