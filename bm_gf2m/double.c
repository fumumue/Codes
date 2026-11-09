#include <stdio.h>
#include <stdint.h>

#define M 3
#define Q (1 << M)
#define N 8
#define T 2

typedef uint8_t gf;

/* GF(8), reduction polynomial x^3 + x + 1 */
static gf gf_mul(gf a, gf b)
{
    gf result = 0;

    while (b) {
        if (b & 1)
            result ^= a;

        b >>= 1;
        a <<= 1;

        if (a & Q)
            a ^= 0xB; /* x^3 + x + 1 */
    }

    return result & (Q - 1);
}

static gf gf_pow(gf a, int e)
{
    gf result = 1;

    while (e > 0) {
        if (e & 1)
            result = gf_mul(result, a);

        a = gf_mul(a, a);
        e >>= 1;
    }

    return result;
}

static gf gf_inv(gf a)
{
    if (a == 0)
        return 0;

    return gf_pow(a, Q - 2);
}

/* g(x) = x^2 + x + 1 */
static gf g_eval(gf x)
{
    return gf_mul(x, x) ^ x ^ 1;
}

static void doubled_syndrome(
    const gf *s,
    const gf *alpha,
    gf *S)
{
    for (int j = 0; j < 2 * T; j++)
        S[j] = 0;

    /*
     * v = (s | 0)
     * Here s is a binary vector of length N.
     */
    for (int i = 0; i < N; i++) {
        if (s[i] != 1)
            continue;

        gf gx = g_eval(alpha[i]);
        gf weight = gf_inv(gf_mul(gx, gx));

        gf power = 1;

        for (int j = 0; j < 2 * T; j++) {
            S[j] ^= gf_mul(power, weight);
            power = gf_mul(power, alpha[i]);
        }
    }
}

int main(void)
{
    gf alpha[N] = {0, 1, 2, 3, 4, 5, 6, 7};

    /*
     * Example ciphertext syndrome s.
     * v = (s | 0)
     */
    gf s[N] = {1, 0, 1, 0, 1, 0, 0, 0};

    gf S[2 * T];

    doubled_syndrome(s, alpha, S);

    printf("Doubled syndrome:\n");

    for (int j = 0; j < 2 * T; j++)
        printf("S[%d] = %u\n", j, S[j]);

    return 0;
}
