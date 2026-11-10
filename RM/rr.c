
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define N 16
#define K 5
#define G_K (K*K)
#define CN 256
#define HN (N-K)

typedef uint8_t bit;

static bit G[K][N];
static bit H[HN][N];

static bit parity(unsigned x)
{
    return __builtin_parity(x);
}

/* RM(1,4) の生成行列 */
void make_G(void)
{
    for (int j = 0; j < N; j++) {
        G[0][j] = 1;

        for (int b = 0; b < 4; b++)
            G[b+1][j] = (j >> b) & 1;
    }
}

/* RM(1,4) の検査行列を生成 */
void make_H(void)
{
    int count = 0;
    uint16_t basis[N] = {0};

    for (unsigned v = 1; v < (1u << N); v++) {
        int valid = 1;

        for (int a = 0; a < K; a++) {
            unsigned dot = 0;

            for (int j = 0; j < N; j++)
                dot ^= G[a][j] & ((v >> j) & 1);

            if (dot) {
                valid = 0;
                break;
            }
        }

        if (!valid) continue;

        unsigned x = v;

        for (int b = N-1; b >= 0; b--) {
            if (!(x & (1u << b))) continue;

            if (basis[b])
                x ^= basis[b];
            else {
                basis[b] = x;
                H[count][0] = 0;

                for (int j = 0; j < N; j++)
                    H[count][j] = (x >> j) & 1;

                count++;
                break;
            }
        }

        if (count == HN) break;
    }

    if (count != HN) {
        fprintf(stderr, "H generation failed\n");
        exit(1);
    }
}

/* 生成行列表示 */
void print_G(void)
{
    printf("\n=== RM(1,4) Generator G (5 x 16) ===\n");

    for (int a = 0; a < K; a++) {
        for (int j = 0; j < N; j++)
            printf("%d", G[a][j]);
        printf("\n");
    }
}

/* 検査行列表示 */
void print_H(void)
{
    printf("\n=== RM(1,4) Parity Check H (11 x 16) ===\n");

    for (int a = 0; a < HN; a++) {
        for (int j = 0; j < N; j++)
            printf("%d", H[a][j]);
        printf("\n");
    }
}

/* 2次元積符号の生成 */
void encode(bit msg[K][K], bit C[N][N])
{
    memset(C, 0, N*N);

    for (int a = 0; a < K; a++) {
        for (int b = 0; b < K; b++) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    C[i][j] ^=
                        msg[a][b] & G[a][i] & G[b][j];
                }
            }
        }
    }
}

/* 行方向のシンドローム */
void syndrome_rows(bit C[N][N], bit S[N][HN])
{
    for (int i = 0; i < N; i++) {
        for (int a = 0; a < HN; a++) {
            S[i][a] = 0;

            for (int j = 0; j < N; j++)
                S[i][a] ^= H[a][j] & C[i][j];
        }
    }
}

/* 列方向のシンドローム */
void syndrome_cols(bit C[N][N], bit S[N][HN])
{
    for (int j = 0; j < N; j++) {
        for (int a = 0; a < HN; a++) {
            S[j][a] = 0;

            for (int i = 0; i < N; i++)
                S[j][a] ^= H[a][i] & C[i][j];
        }
    }
}

void print_matrix(const char *name, bit C[N][N])
{
    printf("\n=== %s ===\n", name);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            printf("%d", C[i][j]);
        printf("\n");
    }
}

void print_syndrome(const char *name, bit S[N][HN])
{
    printf("\n=== %s ===\n", name);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < HN; j++)
            printf("%d", S[i][j]);
        printf("\n");
    }
}

int main(void)
{
    bit msg[K][K];
    bit C[N][N];
    bit R[N][N];
    bit E[N][N];

    bit SR[N][HN];
    bit SC[N][HN];

    make_G();
    make_H();

    print_G();
    print_H();

    /* メッセージ */
    srand(12345);

    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
            msg[i][j] = rand() & 1;

    printf("\n=== Message (5 x 5) ===\n");

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++)
            printf("%d", msg[i][j]);
        printf("\n");
    }

    /* 符号化 */
    encode(msg, C);

    print_matrix("Encoded Codeword (16 x 16)", C);

    /* エラーを追加 */
    memset(E, 0, sizeof(E));

    E[2][3] = 1;
    E[7][8] = 1;
    E[10][4] = 1;

    print_matrix("Error Matrix", E);

    /* 受信語 */
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            R[i][j] = C[i][j] ^ E[i][j];

    print_matrix("Received Matrix", R);

    /* シンドローム */
    syndrome_rows(R, SR);
    syndrome_cols(R, SC);

    print_syndrome("Row Syndrome (16 x 11)", SR);
    print_syndrome("Column Syndrome (16 x 11)", SC);

    /* 符号語の検査 */
    syndrome_rows(C, SR);
    syndrome_cols(C, SC);

    int valid = 1;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < HN; j++)
            if (SR[i][j] || SC[i][j])
                valid = 0;

    printf("\nCodeword verification: %s\n",
           valid ? "PASS" : "FAIL");

    return 0;
}
