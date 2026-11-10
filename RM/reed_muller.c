
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define M 7
#define N (1 << M)
#define K 8

// GF(2)内積の偶奇
static int parity(unsigned int x)
{
    return __builtin_parity(x);
}

// RM(1,7) 符号化
// message = 8 bit
void rm_encode(uint8_t message, uint8_t *code)
{
    unsigned int a = message & 0x7f;
    int b = (message >> 7) & 1;

    for (int i = 0; i < N; i++) {
        code[i] = parity(a & i) ^ b;
    }
}

// 複製符号化
void duplicated_encode(
    uint8_t message,
    int mult,
    uint8_t *code)
{
    uint8_t base[N];

    rm_encode(message, base);

    for (int d = 0; d < mult; d++) {
        memcpy(code + d * N, base, N);
    }
}

// 高速ウォルシュ・アダマール変換
void fwht(int32_t *a, int n)
{
    for (int len = 1; len < n; len <<= 1) {

        for (int i = 0; i < n; i += len * 2) {

            for (int j = 0; j < len; j++) {

                int32_t u = a[i + j];
                int32_t v = a[i + j + len];

                a[i + j]       = u + v;
                a[i + j + len] = u - v;
            }
        }
    }
}

// HQC型 複製RM復号
uint8_t duplicated_decode(
    const uint8_t *received,
    int mult,
    int *best_corr)
{
    int32_t F[N] = {0};

    // 1. 複製ビットのソフト判定加算
    for (int i = 0; i < N; i++) {

        for (int d = 0; d < mult; d++) {

            uint8_t bit = received[d * N + i];

            F[i] += bit ? -1 : 1;
        }
    }

    // 2. 高速アダマール変換
    fwht(F, N);

    // 3. 最大絶対値の探索
    int best = 0;
    int32_t max_abs = -1;

    for (int i = 0; i < N; i++) {

        int32_t v = F[i];
        int32_t av = (v < 0) ? -v : v;

        if (av > max_abs) {
            max_abs = av;
            best = i;
        }
    }


    void hadamard_matrix(int8_t H[128][128])
{
    H[0][0] = 1;

    for (int n = 1; n < 128; n *= 2) {

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {

                H[i][j+n] = H[i][j];
                H[i+n][j] = H[i][j];
                H[i+n][j+n] = -H[i][j];
            }
        }
    }
}


    // 4. 情報ビットの復元
    // 下位7ビット: Walsh係数の位置
    // 最上位ビット: Walsh係数の符号
    uint8_t message = (uint8_t)best;

    if (F[best] < 0)
        message |= 0x80;

    *best_corr = (int)max_abs;

    return message;
}

int main(void)
{
    int mult;

    printf("Duplication factor: ");

    if (scanf("%d", &mult) != 1 ||
        mult < 1 || mult > 100) {
        printf("Invalid multiplicity\n");
        return 1;
    }

    int len = N * mult;

    uint8_t *code = malloc(len);
    uint8_t *received = malloc(len);

    if (!code || !received) {
        printf("Memory allocation error\n");
        free(code);
        free(received);
        return 1;
    }

    unsigned int input;

    printf("Input message (0-255): ");

    if (scanf("%u", &input) != 1 || input > 255) {
        printf("Invalid message\n");
        free(code);
        free(received);
        return 1;
    }

    uint8_t message = (uint8_t)input;

    // 符号化
    duplicated_encode(message, mult, code);

    memcpy(received, code, len);

    // 誤り注入
    int errors;

    printf("Number of errors: ");

    if (scanf("%d", &errors) != 1 ||
        errors < 0 || errors > len) {
        printf("Invalid error count\n");
        free(code);
        free(received);
        return 1;
    }

    printf("Error positions (0-%d):\n", len - 1);

    for (int i = 0; i < errors; i++) {

        int pos;

        if (scanf("%d", &pos) != 1 ||
            pos < 0 || pos >= len) {
            printf("Invalid position\n");
            free(code);
            free(received);
            return 1;
        }

        received[pos] ^= 1;
    }

    // 復号
    int correlation;

    uint8_t decoded =
        duplicated_decode(received, mult, &correlation);

    printf("\nOriginal message : %02X\n", message);
    printf("Decoded message  : %02X\n", decoded);
    printf("Correlation      : %d\n", correlation);

    printf("Result           : %s\n",
           message == decoded ? "PASS" : "FAIL");

    free(code);
    free(received);

    return 0;
}
