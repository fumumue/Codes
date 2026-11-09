
/*
 * rm_hadamard.c
 *
 * General RM(r,m) decoder
 *
 * RM(1,m) : Walsh-Hadamard transform
 * RM(r,m) : Recursive Plotkin decoding
 *
 * gcc -O3 -Wall -Wextra -o rm_hadamard rm_hadamard.c
 *
 * Usage:
 * ./rm_hadamard r m errors
 *
 * Examples:
 * ./rm_hadamard 1 8 30
 * ./rm_hadamard 2 8 20
 * ./rm_hadamard 3 8 10
 * ./rm_hadamard 4 8 5
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define MMAX 13
#define NMAX (1 << MMAX)

/* --------------------------------
   Walsh-Hadamard transform
   -------------------------------- */

static void fwht(int64_t v[], int n)
{
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += 2 * len) {
            for (int j = 0; j < len; j++) {
                int64_t a = v[i+j];
                int64_t b = v[i+j+len];

                v[i+j] = a+b;
                v[i+j+len] = a-b;
            }
        }
    }
}

/* RM(1,m) WHT decoder */

static void decode_wht(int m, const int received[],
                       int decoded[])
{
    int n = 1 << m;
    int64_t v[NMAX];

    for (int i = 0; i < n; i++)
        v[i] = received[i] ? -1 : 1;

    fwht(v, n);

    int best = 0;
    int64_t maxval = -1;

    for (int a = 0; a < n; a++) {
        int64_t c = v[a];
        int64_t ac = c < 0 ? -c : c;

        if (ac > maxval) {
            maxval = ac;
            best = a;
        }
    }

    int b = v[best] < 0;

    for (int x = 0; x < n; x++) {
        decoded[x] =
            b ^ __builtin_parity((unsigned)(best & x));
    }
}

/* --------------------------------
   Recursive Plotkin decoder
   -------------------------------- */

static void decode_rm(int r, int m,
                      const int received[],
                      int decoded[])
{
    int n = 1 << m;

    /* RM(r,m) with r < 0 is the zero code */
    if (r < 0) {
        for (int i = 0; i < n; i++)
            decoded[i] = 0;
        return;
    }

    /* RM(m,m) is the entire binary space */
    if (r >= m) {
        for (int i = 0; i < n; i++)
            decoded[i] = received[i];
        return;
    }

    /* RM(0,m): repetition code */
    if (r == 0) {
        int ones = 0;

        for (int i = 0; i < n; i++)
            ones += received[i];

        int bit = (ones > n / 2);

        for (int i = 0; i < n; i++)
            decoded[i] = bit;

        return;
    }

    /* RM(1,m): use WHT */
    if (r == 1) {
        decode_wht(m, received, decoded);
        return;
    }

    /*
       Plotkin decomposition:

       received = (y0,y1)

       z = y0 XOR y1
       decode z -> v

       u observations:
       y0 and y1 XOR v

       Majority combination -> u

       decoded = (u,u XOR v)
    */

    int half = n / 2;

    int z[NMAX / 2];
    int v[NMAX / 2];
    int uobs[NMAX / 2];
    int u[NMAX / 2];

    /* Decode v in RM(r-1,m-1) */
    for (int i = 0; i < half; i++)
        z[i] = received[i] ^ received[i+half];

    decode_rm(r-1, m-1, z, v);

    /* Combine observations of u */
    for (int i = 0; i < half; i++) {
        int a = received[i];
        int b = received[i+half] ^ v[i];

        /*
         * Majority of two observations.
         * On a tie, choose the first observation.
         */
        if (a == b)
            uobs[i] = a;
        else
            uobs[i] = a;
    }

    /* Decode u in RM(r,m-1) */
    decode_rm(r, m-1, uobs, u);

    /* Reconstruct codeword */
    for (int i = 0; i < half; i++) {
        decoded[i] = u[i];
        decoded[i+half] = u[i] ^ v[i];
    }
}

/* --------------------------------
   General RM encoder
   -------------------------------- */

static void encode(int r, int m,
                   const int msg[],
                   int code[])
{
    int n = 1 << m;

    for (int x = 0; x < n; x++) {
        int bit = 0;

        for (int s = 0; s < n; s++) {
            if (__builtin_popcount((unsigned)s) <= r &&
                (x & s) == s)
                bit ^= msg[s];
        }

        code[x] = bit;
    }
}


/* --------------------------------
   Majority Logic Decoding
   General RM(r,m)
   -------------------------------- */

static void decode_majority_rm(
    int r, int m,
    const int received[],
    int decoded[])
{
    int n = 1 << m;

    if (r < 0) {
        for (int i = 0; i < n; i++)
            decoded[i] = 0;
        return;
    }

    if (r >= m) {
        for (int i = 0; i < n; i++)
            decoded[i] = received[i];
        return;
    }

    /* RM(0,m): repetition code */
    if (r == 0) {
        int ones = 0;

        for (int i = 0; i < n; i++)
            ones += received[i];

        int b = (ones > n / 2);

        for (int i = 0; i < n; i++)
            decoded[i] = b;

        return;
    }

    int top[NMAX] = {0};
    int residual[NMAX];

    /*
       Estimate all degree-r coefficients.

       For each monomial S of degree r,
       use disjoint 2^r-point blocks.

       The parity of each block estimates
       the coefficient of monomial S.
    */

    for (int S = 1; S < n; S++) {

        if (__builtin_popcount((unsigned)S) != r)
            continue;

        int ones = 0;
        int total = 0;

        /* Enumerate all blocks */
        for (int x = 0; x < n; x++) {

            if (x & S)
                continue;

            int p = 0;

            /* XOR over all subsets of S */
            int U = S;

            while (1) {
                p ^= received[x | U];

                if (U == 0)
                    break;

                U = (U - 1) & S;
            }

            ones += p;
            total++;
        }

        /* Majority vote */
        int coefficient = (ones > total / 2);

        /* Reconstruct degree-r contribution */
        if (coefficient) {
            for (int x = 0; x < n; x++) {
                if ((x & S) == S)
                    top[x] ^= 1;
            }
        }
    }

    /* Remove estimated highest-degree terms */
    for (int i = 0; i < n; i++)
        residual[i] = received[i] ^ top[i];

    /* Recursively decode lower-degree terms */
    int lower[NMAX];

    decode_majority_rm(
        r - 1, m,
        residual,
        lower
    );

    /* Reconstruct complete codeword */
    for (int i = 0; i < n; i++)
        decoded[i] = top[i] ^ lower[i];
}

/* --------------------------------
   Main
   -------------------------------- */

int main(int argc, char **argv)
{
    if (argc != 4) {
        printf("Usage: ./rm_hadamard r m errors\n");
        return 1;
    }

    int r = atoi(argv[1]);
    int m = atoi(argv[2]);
    int errors = atoi(argv[3]);

    if (m < 1 || m > MMAX || r < 0 || r > m) {
        printf("Invalid parameters\n");
        return 1;
    }

    int n = 1 << m;

    if (errors < 0 || errors > n) {
        printf("Invalid error count\n");
        return 1;
    }

    srand((unsigned)time(NULL));

    int msg[NMAX] = {0};
    int code[NMAX];
    int received[NMAX];
    int decoded[NMAX];
    int permutation[NMAX];
    int error[NMAX] = {0};
    int decoded_majority[NMAX];

    int k = 0;

    /* Random information coefficients */
    for (int s = 0; s < n; s++) {
        if (__builtin_popcount((unsigned)s) <= r) {
            msg[s] = rand() & 1;
            k++;
        }
    }

    encode(r, m, msg, code);

    for (int i = 0; i < n; i++) {
        received[i] = code[i];
        permutation[i] = i;
    }

    /* Random distinct error positions */
    for (int i = n-1; i > 0; i--) {
        int j = rand() % (i+1);

        int tmp = permutation[i];
        permutation[i] = permutation[j];
        permutation[j] = tmp;
    }

    for (int i = 0; i < errors; i++) {
        int p = permutation[i];

        received[p] ^= 1;
        error[p] = 1;
    }

    int distance = 1 << (m-r);
    int guaranteed = (distance-1)/2;

    printf("\n===== RM(%d,%d) =====\n", r, m);
    printf("Length       : %d\n", n);
    printf("Dimension    : %d\n", k);
    printf("Distance     : %d\n", distance);
    printf("Guaranteed t : %d\n", guaranteed);
    printf("Errors       : %d\n", errors);

    printf("\nError positions:\n");

    for (int i = 0; i < n; i++)
        if (error[i])
            printf("%d ", i);

    printf("\n\n===== Decoding =====\n");

    clock_t start = clock();

    decode_rm(r, m, received, decoded);

    clock_t end = clock();

    double elapsed =
        1000.0 * (double)(end-start) / CLOCKS_PER_SEC;

    printf("Decoder      : ");

    if (r == 1)
        printf("Walsh-Hadamard\n");
    else
        printf("Recursive Plotkin\n");

    printf("Time         : %.3f ms\n", elapsed);

    /* Compare original and decoded codewords */
    int message_ok = 1;

    for (int i = 0; i < n; i++) {
        if (code[i] != decoded[i])
            message_ok = 0;
    }

    int residual = 0;

    for (int i = 0; i < n; i++)
        residual += received[i] ^ decoded[i];

    printf("\n===== Result =====\n");
    printf("Residual error weight : %d\n", residual);

    if (message_ok)
        printf("PASS: Original codeword recovered!\n");
    else
        printf("FAIL: Decoding error.\n");

        printf("\n===== Majority Logic Decoding =====\n");

clock_t ms = clock();

decode_majority_rm(
    r, m,
    received,
    decoded_majority
);

clock_t me = clock();

double majority_time =
    1000.0 * (double)(me-ms) / CLOCKS_PER_SEC;

int majority_residual = 0;
int majority_ok = 1;

for (int i = 0; i < n; i++) {
    majority_residual +=
        received[i] ^ decoded_majority[i];

    if (code[i] != decoded_majority[i])
        majority_ok = 0;
}

printf("Time           : %.3f ms\n", majority_time);
printf("Residual error : %d\n", majority_residual);

printf("Result         : %s\n",
       majority_ok ? "PASS" : "FAIL");

       
    return 0;
}
