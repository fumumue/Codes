
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define M 3
#define N 3
#define LEN (M*N)
#define TOTAL (1u << LEN)

// 2次元巡回畳み込み
// g(x,y) = 1 + x + y

uint16_t encode(uint16_t u)
{
    uint16_t c = 0;

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {

            int p = i * N + j;

            if ((u >> p) & 1) {

                // 1
                c ^= (uint16_t)(1u << p);

                // x
                int px = ((i + 1) % M) * N + j;
                c ^= (uint16_t)(1u << px);

                // y
                int py = i * N + (j + 1) % N;
                c ^= (uint16_t)(1u << py);
            }
        }
    }

    return c;
}

// ハミング重み
int weight(uint16_t x)
{
    return __builtin_popcount((unsigned int)x);
}

// ハミング距離
int distance(uint16_t a, uint16_t b)
{
    return weight(a ^ b);
}

// 最小距離の計算
int minimum_distance(void)
{
    int dmin = LEN + 1;

    for (unsigned int u = 1; u < TOTAL; u++) {

        uint16_t c = encode((uint16_t)u);

        int w = weight(c);

        if (w && w < dmin)
            dmin = w;
    }

    return dmin;
}

// 最尤復号（全符号語探索）
int decode(uint16_t received,
           uint16_t *decoded,
           int *ties)
{
    int best_distance = LEN + 1;
    uint16_t best_code = 0;

    *ties = 0;

    for (unsigned int u = 0; u < TOTAL; u++) {

        uint16_t c = encode((uint16_t)u);

        int d = distance(received, c);

        if (d < best_distance) {

            best_distance = d;
            best_code = c;
            *ties = 1;

        } else if (d == best_distance) {

            (*ties)++;
        }
    }

    *decoded = best_code;

    return best_distance;
}

// 3x3行列表示
void print_matrix(uint16_t x)
{
    for (int i = 0; i < M; i++) {

        for (int j = 0; j < N; j++) {

            int p = i * N + j;

            printf("%d ", (x >> p) & 1);
        }

        printf("\n");
    }
}

int main(void)
{
    unsigned int input;
    int errors;

    printf("2D Cyclic Code\n");
    printf("Generator: g(x,y) = 1+x+y\n");

    printf("Input message (0-511): ");

    if (scanf("%u", &input) != 1 ||
        input >= TOTAL) {
        printf("Invalid input\n");
        return 1;
    }

    uint16_t message = (uint16_t)input;

    // 符号化
    uint16_t code = encode(message);

    printf("\nInformation:\n");
    print_matrix(message);

    printf("\nEncoded codeword:\n");
    print_matrix(code);

    // 誤り注入
    uint16_t received = code;

    printf("\nNumber of errors: ");

    if (scanf("%d", &errors) != 1 ||
        errors < 0 || errors > LEN) {
        printf("Invalid error count\n");
        return 1;
    }

    for (int i = 0; i < errors; i++) {

        int pos;

        printf("Error position (0-8): ");

        if (scanf("%d", &pos) != 1 ||
            pos < 0 || pos >= LEN) {
            printf("Invalid position\n");
            return 1;
        }

        received ^= (uint16_t)(1u << pos);
    }

    printf("\nReceived:\n");
    print_matrix(received);

    // 復号
    uint16_t decoded;
    int ties;

    int d = decode(received, &decoded, &ties);

    printf("\nDecoded:\n");
    print_matrix(decoded);

    // 最小距離
    int dmin = minimum_distance();

    printf("\nMinimum distance = %d\n", dmin);
    printf("Correction capacity = %d\n",
           (dmin - 1) / 2);

    printf("Decoding distance = %d\n", d);
    printf("Nearest codeword candidates = %d\n", ties);

    if (ties == 1 && decoded == code)
        printf("Result: PASS\n");
    else if (ties > 1)
        printf("Result: AMBIGUOUS\n");
    else
        printf("Result: FAIL\n");

    return 0;
}
