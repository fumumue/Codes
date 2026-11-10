
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define RM_N 16
#define RM_K 5
#define G_N 24
#define G_K 12
#define N (RM_N * G_N)
#define SYNDROMES 4096
#define MAX_ITER 30

typedef uint8_t bit;

/* 拡張ゴーレイ符号 */
static bit G[12][24];
static bit H[12][24];
static int pivot[12];
static uint32_t syndrome_table[4096];
static int table_valid[4096];

/* ---------- 基本関数 ---------- */

static int parity32(uint32_t x)
{
    return __builtin_parity(x);
}

static int parity_bits(const bit *v, int n)
{
    int p = 0;
    for (int i = 0; i < n; i++)
        p ^= v[i];
    return p;
}

/* ---------- ゴーレイ符号の生成 ---------- */

/*
 * 23ビットGolay符号の生成多項式
 *
 * g(x) = x^11+x^9+x^7+x^6+x^5+x+1
 *
 * これに全体パリティを追加して
 * 拡張Golay [24,12,8] を作る。
 */

static void init_golay(void)
{
    memset(G, 0, sizeof(G));

    const int exponent[] = {0,1,5,6,7,9,11};

    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 7; j++) {
            G[i][i + exponent[j]] = 1;
        }

        G[i][23] = parity_bits(G[i], 23);
    }

    /* GF(2)で生成行列を簡約化 */
    int rank = 0;

    for (int col = 0; col < 24 && rank < 12; col++) {

        int selected = -1;

        for (int r = rank; r < 12; r++) {
            if (G[r][col]) {
                selected = r;
                break;
            }
        }

        if (selected < 0)
            continue;

        for (int j = 0; j < 24; j++) {
            bit t = G[rank][j];
            G[rank][j] = G[selected][j];
            G[selected][j] = t;
        }

        for (int r = 0; r < 12; r++) {
            if (r == rank || !G[r][col])
                continue;

            for (int j = 0; j < 24; j++)
                G[r][j] ^= G[rank][j];
        }

        pivot[rank] = col;
        rank++;
    }

    if (rank != 12) {
        fprintf(stderr, "Golay matrix rank error\n");
        exit(1);
    }

    /* Hを構築する */
    int free_col[12];
    int nf = 0;

    for (int j = 0; j < 24; j++) {
        int found = 0;

        for (int i = 0; i < 12; i++)
            if (pivot[i] == j)
                found = 1;

        if (!found)
            free_col[nf++] = j;
    }

    if (nf != 12) {
        fprintf(stderr, "Golay parity matrix error\n");
        exit(1);
    }

    memset(H, 0, sizeof(H));

    for (int r = 0; r < 12; r++) {
        H[r][free_col[r]] = 1;

        for (int i = 0; i < 12; i++)
            H[r][pivot[i]] = G[i][free_col[r]];
    }
}

/* ---------- ゴーレイ符号のシンドローム ---------- */

static int golay_syndrome(const bit *v)
{
    int s = 0;

    for (int i = 0; i < 12; i++) {
        int p = 0;

        for (int j = 0; j < 24; j++)
            p ^= H[i][j] & v[j];

        s |= p << i;
    }

    return s;
}

static int error_syndrome(uint32_t e)
{
    bit v[24];

    for (int i = 0; i < 24; i++)
        v[i] = (e >> i) & 1;

    return golay_syndrome(v);
}

/*
 * 重み3以下の誤りパターンを登録する。
 * 全符号語の探索ではなく、訂正可能な
 * 2325個の誤りパターンだけを登録する。
 */

static void register_error(uint32_t e)
{
    int s = error_syndrome(e);

    if (table_valid[s]) {
        fprintf(stderr,
                "Golay syndrome collision: %d\n", s);
        exit(1);
    }

    syndrome_table[s] = e;
    table_valid[s] = 1;
}

static void init_syndrome_table(void)
{
    memset(table_valid, 0, sizeof(table_valid));

    register_error(0);

    for (int i = 0; i < 24; i++)
        register_error(1u << i);

    for (int i = 0; i < 24; i++)
        for (int j = i + 1; j < 24; j++)
            register_error((1u << i) | (1u << j));

    for (int i = 0; i < 24; i++)
        for (int j = i + 1; j < 24; j++)
            for (int k = j + 1; k < 24; k++)
                register_error(
                    (1u << i) | (1u << j) | (1u << k));
}

/* ゴーレイ符号化 */

static void golay_encode(const bit *msg, bit *code)
{
    memset(code, 0, 24);

    for (int i = 0; i < 12; i++) {
        if (!msg[i])
            continue;

        for (int j = 0; j < 24; j++)
            code[j] ^= G[i][j];
    }
}

/*
 * Golay復号
 * 成功: 1
 * 訂正不能: 0
 */

static int golay_decode(bit *v)
{
    int s = golay_syndrome(v);

    if (!table_valid[s])
        return 0;

    uint32_t e = syndrome_table[s];

    for (int i = 0; i < 24; i++)
        v[i] ^= (e >> i) & 1;

    return 1;
}

/* ---------- RM(1,4) ---------- */

/*
 * Walsh-Hadamard変換
 *
 * 16点の変換をO(16 log 16)で計算。
 */

static void hadamard(int *a, int n)
{
    for (int len = 1; len < n; len <<= 1) {

        for (int i = 0; i < n; i += 2 * len) {

            for (int j = 0; j < len; j++) {

                int x = a[i + j];
                int y = a[i + j + len];

                a[i + j] = x + y;
                a[i + j + len] = x - y;
            }
        }
    }
}

/* RM符号化 */

static void rm_encode(const bit *msg, bit *code)
{
    for (int x = 0; x < 16; x++) {

        int p = msg[0];

        for (int j = 0; j < 4; j++)
            p ^= msg[j + 1] & ((x >> j) & 1);

        code[x] = p;
    }
}

/*
 * RM(1,4)の最尤復号
 * Walsh-Hadamard変換で16個の相関値を求める。
 *
 * 戻り値:
 * 1 = 受信語がRM符号語だった
 * 0 = 誤りを含んでいた
 */

static int rm_decode(bit *v, bit *msg)
{
    int a[16];

    int original = 1;

    for (int i = 0; i < 16; i++) {
        a[i] = v[i] ? -1 : 1;
    }

    hadamard(a, 16);

    int best = 0;
    int best_abs = -1;

    for (int i = 0; i < 16; i++) {

        int z = abs(a[i]);

        if (z > best_abs) {
            best_abs = z;
            best = i;
        }
    }

    msg[0] = (a[best] < 0);

    for (int j = 0; j < 4; j++)
        msg[j + 1] = (best >> j) & 1;

    bit corrected[16];
    rm_encode(msg, corrected);

    for (int i = 0; i < 16; i++) {
        if (corrected[i] != v[i])
            original = 0;

        v[i] = corrected[i];
    }

    return original;
}

/* ---------- 積符号 ---------- */

/*
 * 情報行 5 x 12
 * Golay符号化後 5 x 24
 * RM符号化後 16 x 24
 */

static void product_encode(
    bit msg[5][12],
    bit code[16][24])
{
    bit temp[5][24];

    for (int i = 0; i < 5; i++)
        golay_encode(msg[i], temp[i]);

    for (int j = 0; j < 24; j++) {

        bit m[5], c[16];

        for (int i = 0; i < 5; i++)
            m[i] = temp[i][j];

        rm_encode(m, c);

        for (int i = 0; i < 16; i++)
            code[i][j] = c[i];
    }
}

/* ---------- 積符号の反復復号 ---------- */

static int is_product_codeword(bit v[16][24])
{
    /* 全行がGolay符号語か確認 */

    for (int i = 0; i < 16; i++) {
        if (golay_syndrome(v[i]) != 0)
            return 0;
    }

    /* 全列がRM符号語か確認 */

    for (int j = 0; j < 24; j++) {

        bit col[16], msg[5];

        for (int i = 0; i < 16; i++)
            col[i] = v[i][j];

        if (!rm_decode(col, msg))
            return 0;
    }

    return 1;
}

/*
 * 行: Golay
 * 列: RM
 *
 * 交互復号を繰り返す。
 */

static int product_decode(
    bit v[16][24],
    int *iterations)
{
    for (int iter = 0; iter < MAX_ITER; iter++) {

        /* Golay行復号 */

        for (int i = 0; i < 16; i++)
            golay_decode(v[i]);

        /* RM列復号 */

        for (int j = 0; j < 24; j++) {

            bit col[16], msg[5];

            for (int i = 0; i < 16; i++)
                col[i] = v[i][j];

            rm_decode(col, msg);

            for (int i = 0; i < 16; i++)
                v[i][j] = col[i];
        }

        *iterations = iter + 1;

        if (is_product_codeword(v))
            return 1;
    }

    return 0;
}

/* ---------- 情報語の復元 ---------- */

static void product_extract(
    bit code[16][24],
    bit msg[5][12])
{
    bit temp[5][24];

    for (int j = 0; j < 24; j++) {

        bit col[16], m[5];

        for (int i = 0; i < 16; i++)
            col[i] = code[i][j];

        rm_decode(col, m);

        for (int i = 0; i < 5; i++)
            temp[i][j] = m[i];
    }

    for (int i = 0; i < 5; i++) {

        bit row[24];

        memcpy(row, temp[i], 24);

        if (!golay_decode(row)) {
            fprintf(stderr, "Extraction error\n");
            exit(1);
        }

        for (int j = 0; j < 12; j++)
            msg[i][j] = row[pivot[j]];
    }
}

/* ---------- 表示 ---------- */

static void print_matrix(const char *title, bit v[16][24])
{
    printf("\n%s\n", title);

    for (int i = 0; i < 16; i++) {

        printf("%2d: ", i);

        for (int j = 0; j < 24; j++)
            printf("%d", v[i][j]);

        printf("\n");
    }
}

static int matrix_equal(
    bit a[16][24],
    bit b[16][24])
{
    return memcmp(a, b, 16 * 24) == 0;
}

/* ---------- メイン ---------- */

int main(int argc, char **argv)
{
    int errors = 20;
    unsigned seed = 1;

    if (argc >= 2)
        errors = atoi(argv[1]);

    if (argc >= 3)
        seed = (unsigned)atoi(argv[2]);

    if (errors < 0 || errors > 384) {
        printf("Error count must be 0..384\n");
        return 1;
    }

    srand(seed);

    init_golay();
    init_syndrome_table();

    printf("RM(1,4) x Extended Golay [24,12,8]\n");
    printf("Product code: [384,60,64]\n");
    printf("Golay syndrome patterns: 2325\n");
    printf("Random errors: %d\n", errors);
    printf("Seed: %u\n", seed);

    /* 情報語 */

    bit message[5][12];

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 12; j++)
            message[i][j] = rand() & 1;

    /* 符号化 */

    bit original[16][24];

    product_encode(message, original);

    print_matrix("Original codeword", original);

    /* 誤りを付加 */

    bit received[16][24];

    memcpy(received, original, sizeof(original));

    int positions[384];

    for (int i = 0; i < 384; i++)
        positions[i] = i;

    /* Fisher-Yates shuffle */

    for (int i = 383; i > 0; i--) {
        int j = rand() % (i + 1);

        int t = positions[i];
        positions[i] = positions[j];
        positions[j] = t;
    }

    printf("\nError positions:\n");

    for (int i = 0; i < errors; i++) {

        int p = positions[i];
        int r = p / 24;
        int c = p % 24;

        received[r][c] ^= 1;

        printf("(%d,%d) ", r, c);
    }

    printf("\n");

    print_matrix("Received word", received);

    /* 復号 */

    int iterations = 0;

    clock_t start = clock();

    int success = product_decode(received, &iterations);

    clock_t end = clock();

    double elapsed =
        1000.0 * (double)(end - start) / CLOCKS_PER_SEC;

    print_matrix("Decoded codeword", received);

    /* 判定 */

    printf("\n========== RESULT ==========\n");

    printf("Decoder status: %s\n",
           success ? "VALID PRODUCT CODEWORD" : "FAILED");

    printf("Iterations: %d\n", iterations);

    printf("Time: %.3f ms\n", elapsed);

    printf("Original codeword match: %s\n",
           matrix_equal(original, received) ? "YES" : "NO");

    if (success && matrix_equal(original, received)) {

        bit recovered[5][12];

        product_extract(received, recovered);

        printf("\nRecovered information:\n");

        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 12; j++)
                printf("%d", recovered[i][j]);

            printf("\n");
        }

        printf("Information match: %s\n",
            memcmp(message, recovered, sizeof(message)) == 0
            ? "YES" : "NO");
    }

    return 0;
}
