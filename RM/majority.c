
/* rm_majority.c
 *
 * RM(1,8) majority logic decoder
 *
 * n = 256
 * k = 9
 * d = 128
 * guaranteed correction: 63 errors
 *
 * gcc -O3 -Wall -Wextra -o rm_majority rm_majority.c
 * ./rm_majority 63
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define M 8
#define N (1 << M)
#define K (M + 1)

static int parity(unsigned int x)
{
    return __builtin_parity(x);
}

/* RM(1,8) encoding */
void encode(const int msg[K], int code[N])
{
    for (int x = 0; x < N; x++) {
        int v = msg[0];

        for (int i = 0; i < M; i++)
            v ^= msg[i + 1] & ((x >> i) & 1);

        code[x] = v;
    }
}

/* Majority logic decoding */
int decode(const int received[N], int msg[K])
{
    int ones, zeros;

    /* Estimate a_i by 128 pair checks */
    for (int i = 0; i < M; i++) {
        ones = 0;
        zeros = 0;

        int mask = 1 << i;

        for (int x = 0; x < N; x++) {
            if (x & mask)
                continue;

            int s = received[x] ^ received[x | mask];

            if (s)
                ones++;
            else
                zeros++;
        }

        printf("a%d: 0=%3d 1=%3d ", i + 1,
               zeros, ones);

        if (ones > zeros)
            msg[i + 1] = 1;
        else if (zeros > ones)
            msg[i + 1] = 0;
        else {
            printf("TIE\n");
            return 0;
        }

        printf("=> %d\n", msg[i + 1]);
    }

    /* Estimate constant term b */
    ones = 0;
    zeros = 0;

    for (int x = 0; x < N; x++) {
        int v = received[x];

        for (int i = 0; i < M; i++)
            v ^= msg[i + 1] & ((x >> i) & 1);

        if (v)
            ones++;
        else
            zeros++;
    }

    printf("b : 0=%3d 1=%3d ", zeros, ones);

    if (ones > zeros)
        msg[0] = 1;
    else if (zeros > ones)
        msg[0] = 0;
    else {
        printf("TIE\n");
        return 0;
    }

    printf("=> %d\n", msg[0]);

    return 1;
}

int main(int argc, char **argv)
{
    int errors = 30;

    if (argc >= 2)
        errors = atoi(argv[1]);

    if (errors < 0 || errors > N) {
        printf("Error count must be 0..256\n");
        return 1;
    }

    srand((unsigned int)time(NULL));

    int msg[K], decoded[K];
    int code[N], received[N];
    int error[N] = {0};
    int permutation[N];

    /* Random information */
    for (int i = 0; i < K; i++)
        msg[i] = rand() & 1;

    /* Encode */
    encode(msg, code);

    for (int i = 0; i < N; i++) {
        received[i] = code[i];
        permutation[i] = i;
    }

    /* Select distinct error positions */
    for (int i = N - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int tmp = permutation[i];
        permutation[i] = permutation[j];
        permutation[j] = tmp;
    }

    for (int i = 0; i < errors; i++) {
        int p = permutation[i];
        error[p] = 1;
        received[p] ^= 1;
    }

    printf("===== RM(1,8) =====\n");
    printf("Length       : %d\n", N);
    printf("Dimension    : %d\n", K);
    printf("Min distance : 128\n");
    printf("Errors       : %d\n\n", errors);

    printf("Original message:\n");
    for (int i = 0; i < K; i++)
        printf("%d", msg[i]);
    printf("\n\n");

    printf("Error positions:\n");
    for (int i = 0; i < N; i++) {
        if (error[i])
            printf("%d ", i);
    }
    printf("\n\n");

    printf("===== Majority decoding =====\n");

    int success = decode(received, decoded);

    printf("\n===== Result =====\n");

    if (!success) {
        printf("Decoding failed (tie).\n");
        return 0;
    }

    printf("Decoded message:\n");
    for (int i = 0; i < K; i++)
        printf("%d", decoded[i]);
    printf("\n");

    int correct = 1;

    for (int i = 0; i < K; i++) {
        if (msg[i] != decoded[i])
            correct = 0;
    }

    if (correct)
        printf("PASS: Original message recovered!\n");
    else
        printf("FAIL: Decoding error.\n");

    /* Re-encode decoded message and count residual errors */
    encode(decoded, code);

    int residual = 0;
    for (int i = 0; i < N; i++)
        residual += code[i] ^ received[i];

    printf("Residual error weight: %d\n", residual);

    return 0;
}
