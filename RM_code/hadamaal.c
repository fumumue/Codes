
/*
 * rm_hadamard.c
 *
 * RM(1,m) maximum likelihood decoder
 * using Fast Walsh-Hadamard Transform.
 *
 * gcc -O3 -Wall -Wextra -o rm_hadamard rm_hadamard.c
 *
 * Usage:
 * ./rm_hadamard m errors
 *
 * Examples:
 * ./rm_hadamard 8 30
 * ./rm_hadamard 8 63
 * ./rm_hadamard 9 100
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define MMAX 15
#define NMAX (1 << MMAX)

static void encode(int m, int b, int a, int code[])
{
    int n = 1 << m;

    for (int x = 0; x < n; x++) {
        code[x] = b ^ __builtin_parity((unsigned)(a & x));
    }
}

/* Fast Walsh-Hadamard Transform */
static void fwht(int64_t v[], int n)
{
    for (int len = 1; len < n; len <<= 1) {

        for (int i = 0; i < n; i += 2 * len) {

            for (int j = 0; j < len; j++) {

                int64_t u = v[i + j];
                int64_t w = v[i + j + len];

                v[i + j]       = u + w;
                v[i + j + len] = u - w;
            }
        }
    }
}

/* Maximum likelihood decoder */
static int decode(int m, int received[],
                  int *decoded_b, int *decoded_a,
                  int64_t *best_correlation,
                  int *tie)
{
    int n = 1 << m;
    int64_t v[NMAX];

    /* Convert bits to bipolar representation */
    for (int x = 0; x < n; x++)
        v[x] = received[x] ? -1 : 1;

    /* Transform */
    fwht(v, n);

    /* Find maximum absolute correlation */
    int best = 0;
    int64_t maxabs = 0;
    *tie = 0;

    for (int a = 0; a < n; a++) {

        int64_t c = v[a];
        int64_t ac = (c < 0) ? -c : c;

        if (ac > maxabs) {
            maxabs = ac;
            best = a;
            *tie = 0;
        }
        else if (ac == maxabs) {
            *tie = 1;
        }
    }

    *decoded_a = best;
    *decoded_b = (v[best] < 0);
    *best_correlation = v[best];

    return 1;
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        printf("Usage: ./rm_hadamard m errors\n");
        return 1;
    }

    int m = atoi(argv[1]);
    int errors = atoi(argv[2]);

    if (m < 1 || m > MMAX) {
        printf("m must be between 1 and %d\n", MMAX);
        return 1;
    }

    int n = 1 << m;

    if (errors < 0 || errors > n) {
        printf("Invalid error count\n");
        return 1;
    }

    srand((unsigned)time(NULL));

    /* Random message */
    int b = rand() & 1;
    int a = rand() % n;

    int code[NMAX];
    int received[NMAX];
    int error[NMAX] = {0};
    int permutation[NMAX];

    encode(m, b, a, code);

    for (int i = 0; i < n; i++) {
        received[i] = code[i];
        permutation[i] = i;
    }

    /* Select distinct random error positions */
    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int tmp = permutation[i];
        permutation[i] = permutation[j];
        permutation[j] = tmp;
    }

    for (int i = 0; i < errors; i++) {
        int p = permutation[i];
        received[p] ^= 1;
        error[p] = 1;
    }

    printf("\n===== RM(1,%d) Hadamard Decoder =====\n", m);
    printf("Length       : %d\n", n);
    printf("Dimension    : %d\n", m + 1);
    printf("Distance     : %d\n", n / 2);
    printf("Guaranteed t : %d\n", (n / 2 - 1) / 2);
    printf("Errors       : %d\n", errors);

    printf("\nOriginal message:\n");
    printf("b=%d a=%d\n", b, a);

    printf("\nError positions:\n");
    for (int i = 0; i < n; i++)
        if (error[i])
            printf("%d ", i);

    printf("\n\n===== WHT decoding =====\n");

    int db, da, tie;
    int64_t correlation;

    decode(m, received, &db, &da,
           &correlation, &tie);

    printf("Decoded b          : %d\n", db);
    printf("Decoded a          : %d\n", da);
    printf("Correlation        : %lld\n",
           (long long)correlation);
    printf("Maximum tie        : %s\n",
           tie ? "YES" : "NO");

    int decoded[NMAX];
    encode(m, db, da, decoded);

    int residual = 0;

    for (int i = 0; i < n; i++)
        residual += decoded[i] ^ received[i];

    printf("Residual error     : %d\n", residual);

    if (db == b && da == a)
        printf("PASS: Original message recovered!\n");
    else
        printf("FAIL: Decoding error.\n");

    return 0;
}
