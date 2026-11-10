
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MAX_M 12
#define MAX_N (1 << MAX_M)

// --------------------------------------------------
// RM(1,m) encoder
// message: m+1 bits
// output: +1 or -1
// --------------------------------------------------

void rm_encode(int message, int m, int8_t *out)
{
    int n = 1 << m;

    for (int x = 0; x < n; x++) {
        int bit = (message & 1)
                ^ __builtin_parity((message >> 1) & x);

        out[x] = bit ? -1 : +1;
    }
}

// --------------------------------------------------
// Duplicated RM encoder
// --------------------------------------------------

void rm_dup_encode(
    int message,
    int m,
    int multiplicity,
    int8_t *out)
{
    int n = 1 << m;

    int8_t base[MAX_N];

    rm_encode(message, m, base);

    for (int d = 0; d < multiplicity; d++) {
        memcpy(out + d*n, base, n);
    }
}

// --------------------------------------------------
// Fast Walsh-Hadamard Transform
// --------------------------------------------------

void fwht(int32_t *a, int n)
{
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += 2*len) {
            for (int j = 0; j < len; j++) {
                int32_t x = a[i+j];
                int32_t y = a[i+j+len];

                a[i+j] = x+y;
                a[i+j+len] = x-y;
            }
        }
    }
}

// --------------------------------------------------
// Duplicated RM decoder
// input: signed values (+1/-1)
// --------------------------------------------------

int rm_dup_decode(
    const int8_t *received,
    int m,
    int multiplicity,
    int32_t *metric)
{
    int n = 1 << m;

    int32_t F[MAX_N] = {0};

    // Add corresponding signed positions
    for (int x = 0; x < n; x++) {
        for (int d = 0; d < multiplicity; d++) {
            F[x] += received[d*n+x];
        }
    }

    // Hadamard transform
    fwht(F, n);

    // Maximum absolute correlation
    int best = 0;
    int32_t maxval = -1;

    for (int i = 0; i < n; i++) {
        int32_t v = abs(F[i]);

        if (v > maxval) {
            maxval = v;
            best = i;
        }
    }

    if (metric)
        *metric = maxval;

    int message = best << 1;

    if (F[best] < 0)
        message |= 1;

    return message;
}

// --------------------------------------------------
// Main
// --------------------------------------------------

int main(void)
{
    int information_bits;
    int multiplicity;
    int errors;

    printf("Information bits (11 or 13): ");
    scanf("%d", &information_bits);

    int m;

    if (information_bits == 11)
        m = 10;
    else if (information_bits == 13)
        m = 12;
    else
        return 1;

    printf("Multiplicity: ");
    scanf("%d", &multiplicity);

    if (multiplicity < 1 || multiplicity > 100)
        return 1;

    printf("Number of errors: ");
    scanf("%d", &errors);

    int n = 1 << m;
    int total = n * multiplicity;

    int8_t *code = calloc(total, sizeof(int8_t));

    if (!code)
        return 1;

    int message = rand() & ((1 << (m+1))-1);

    printf("\nOriginal message: %04X\n", message);

    rm_dup_encode(message, m, multiplicity, code);

    printf("Code length: %d\n", total);

    // Random bit errors
    for (int i = 0; i < errors; i++) {
        int pos = rand() % total;
        code[pos] = -code[pos];
    }

    int32_t metric;

    int decoded =
        rm_dup_decode(code, m, multiplicity, &metric);

    printf("Decoded message:  %04X\n", decoded);
    printf("Metric: %d\n", metric);

    printf("Result: %s\n",
           message == decoded ? "PASS" : "FAIL");

    free(code);

    return 0;
}
