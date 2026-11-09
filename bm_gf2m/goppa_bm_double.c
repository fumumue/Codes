/*
 * Binary Goppa-style syndrome decoding with Berlekamp-Massey over GF(2^13).
 *
 * - Field: GF(2^13) = GF(2)[x]/(x^13+x^4+x^3+x+1)
 * - Support: all 8192 field elements, represented by 13-bit polynomial-basis values.
 * - Goppa polynomial g(x): random monic degree-13 polynomial, checked to be
 *   square-free and nonzero at every support point.
 * - Ordinary moments:  S[j]  = sum_{i in E} a_i^j / g(a_i)
 * - Doubled moments:   S2[j] = sum_{i in E} a_i^j / g(a_i)^2
 *
 * The doubled moments are the 2t parity-check moments for the alternant code
 * defined by g(x)^2. For square-free g, the binary Goppa codes Gamma(L,g)
 * and Gamma(L,g^2) coincide. BM is run independently on the ordinary and
 * doubled sequences; their locators must agree. This is a cross-check, not
 * an invalid concatenation of two different sequences.
 *
 * The coordinate a=0 is handled separately: BM is run on moments S[1..2t],
 * which removes its contribution, and S[0] is then used to decide whether
 * coordinate zero was in error.
 *
 * Build:  cc -O2 -std=c11 -Wall -Wextra -o goppa_bm_double goppa_bm_double.c
 * Run:    ./goppa_bm_double [error_weight]
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define M       13
#define Q       (1u << M)
#define N       Q
#define T       32
#define NSYN    (2 * T)
#define POLY    0x201bu /* x^13 + x^4 + x^3 + x + 1 */
#define MAXP    (T + 1)

static uint32_t rng_state = 0x6d2b79f5u;
static uint32_t rnd32(void)
{
    uint32_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_state = x;
    return x;
}

static uint16_t gf_mul(uint16_t a, uint16_t b)
{
    uint16_t r = 0;
    while (b) {
        if (b & 1u) r ^= a;
        b >>= 1;
        a <<= 1;
        if (a & Q) a ^= POLY;
    }
    return r & (Q - 1u);
}

static uint16_t gf_pow(uint16_t a, unsigned e)
{
    uint16_t r = 1;
    while (e) {
        if (e & 1u) r = gf_mul(r, a);
        a = gf_mul(a, a);
        e >>= 1;
    }
    return r;
}

static uint16_t gf_inv(uint16_t a)
{
    if (!a) {
        fprintf(stderr, "GF inverse of zero\n");
        exit(EXIT_FAILURE);
    }
    return gf_pow(a, Q - 2u);
}

static uint16_t poly_eval(const uint16_t *p, int d, uint16_t x)
{
    uint16_t y = p[d];
    for (int i = d - 1; i >= 0; --i)
        y = gf_mul(y, x) ^ p[i];
    return y;
}

typedef struct {
    int d;
    uint16_t c[2 * T + 2];
} Poly;

static void poly_trim(Poly *p)
{
    while (p->d > 0 && p->c[p->d] == 0) --p->d;
}

static Poly poly_remainder(Poly a, const Poly *b)
{
    if (b->d == 0 && b->c[0] == 0) {
        fprintf(stderr, "polynomial division by zero\n");
        exit(EXIT_FAILURE);
    }
    uint16_t invlead = gf_inv(b->c[b->d]);
    while (a.d >= b->d && !(a.d == 0 && a.c[0] == 0)) {
        int shift = a.d - b->d;
        uint16_t factor = gf_mul(a.c[a.d], invlead);
        for (int j = 0; j <= b->d; ++j)
            a.c[j + shift] ^= gf_mul(factor, b->c[j]);
        poly_trim(&a);
    }
    return a;
}

static int g_is_square_free(const uint16_t *g)
{
    Poly a = { .d = T }, b = { .d = T - 1 };
    for (int i = 0; i <= T; ++i) a.c[i] = g[i];
    memset(b.c, 0, sizeof(b.c));
    /* Formal derivative in characteristic two. */
    for (int i = 1; i <= T; i += 2) b.c[i - 1] = g[i];
    poly_trim(&b);
    while (!(b.d == 0 && b.c[0] == 0)) {
        Poly r = poly_remainder(a, &b);
        a = b;
        b = r;
    }
    /* gcd(g,g') is a nonzero constant iff g is square-free. */
    return a.d == 0 && a.c[0] != 0;
}

static void make_goppa_polynomial(uint16_t *g)
{
    for (unsigned attempt = 0; attempt < 10000; ++attempt) {
        for (int i = 0; i < T; ++i) g[i] = (uint16_t)(rnd32() & (Q - 1u));
        g[T] = 1;
        if (!g_is_square_free(g)) continue;
        int good = 1;
        for (unsigned a = 0; a < Q; ++a) {
            if (poly_eval(g, T, (uint16_t)a) == 0) {
                good = 0;
                break;
            }
        }
        if (good) return;
    }
    fprintf(stderr, "could not find a suitable square-free g(x)\n");
    exit(EXIT_FAILURE);
}

/* Berlekamp-Massey. s[0..len-1] is a power-sum sequence. */
static int berlekamp_massey(const uint16_t *s, int len, uint16_t *C)
{
    uint16_t B[NSYN + 1] = {0}, oldC[NSYN + 1] = {0};
    memset(C, 0, (NSYN + 1) * sizeof(*C));
    C[0] = B[0] = 1;
    int L = 0, m = 1;
    uint16_t b = 1;

    for (int n = 0; n < len; ++n) {
        uint16_t d = s[n];
        for (int i = 1; i <= L && i <= n; ++i)
            d ^= gf_mul(C[i], s[n - i]);
        if (d == 0) {
            ++m;
            continue;
        }
        memcpy(oldC, C, sizeof(oldC));
        uint16_t coef = gf_mul(d, gf_inv(b));
        for (int i = 0; i + m <= len; ++i)
            C[i + m] ^= gf_mul(coef, B[i]);
        if (2 * L <= n) {
            L = n + 1 - L;
            memcpy(B, oldC, sizeof(B));
            b = d;
            m = 1;
        } else {
            ++m;
        }
    }
    return L;
}

static int eval_locator(const uint16_t *C, int L, uint16_t z)
{
    uint16_t y = C[L];
    for (int i = L - 1; i >= 0; --i) y = gf_mul(y, z) ^ C[i];
    return y == 0;
}

/* Decode nonzero support points. Locator roots are z = 1/a. */
static int chien_search(const uint16_t *C, int L, int *positions)
{
    int count = 0;
    for (unsigned a = 1; a < Q; ++a) {
        if (eval_locator(C, L, gf_inv((uint16_t)a))) {
            if (count < T) positions[count] = (int)a;
            ++count;
        }
    }
    return count;
}

static void make_syndromes(const uint8_t *error, const uint16_t *g,
                           uint16_t *S, uint16_t *S2)
{
    memset(S, 0, (NSYN + 1) * sizeof(*S));
    memset(S2, 0, (NSYN + 1) * sizeof(*S2));
    for (unsigned a = 0; a < Q; ++a) {
        if (!error[a]) continue;
        uint16_t ga = poly_eval(g, T, (uint16_t)a);
        uint16_t w = gf_inv(ga);
        uint16_t w2 = gf_mul(w, w);
        uint16_t p = 1;
        for (int j = 0; j <= NSYN; ++j) {
            S[j] ^= gf_mul(p, w);
            S2[j] ^= gf_mul(p, w2);
            p = gf_mul(p, (uint16_t)a);
        }
    }
}

static int same_positions(const int *a, int na, const int *b, int nb)
{
    if (na != nb || na < 0 || na > T) return 0;
    for (int i = 0; i < na; ++i) {
        int found = 0;
        for (int j = 0; j < nb; ++j) if (a[i] == b[j]) found = 1;
        if (!found) return 0;
    }
    return 1;
}

static int recover_zero_coordinate(uint16_t S0, uint16_t S20,
                                   const int *pos, int count,
                                   const uint16_t *g)
{
    uint16_t predicted = 0, predicted2 = 0;
    for (int i = 0; i < count; ++i) {
        uint16_t w = gf_inv(poly_eval(g, T, (uint16_t)pos[i]));
        predicted ^= w;
        predicted2 ^= gf_mul(w, w);
    }
    uint16_t r = S0 ^ predicted;
    uint16_t r2 = S20 ^ predicted2;
    uint16_t w0 = gf_inv(g[0]);
    uint16_t w02 = gf_mul(w0, w0);
    int ordinary = (r == 0) ? 0 : (r == w0) ? 1 : -1;
    int doubled  = (r2 == 0) ? 0 : (r2 == w02) ? 1 : -1;
    if (ordinary < 0 || doubled < 0 || ordinary != doubled) return -1;
    return ordinary;
}

int main(int argc, char **argv)
{
    int weight = T;
    if (argc >= 2) {
        char *end = NULL;
        long v = strtol(argv[1], &end, 10);
        if (!end || *end || v < 0 || v > T) {
            fprintf(stderr, "usage: %s [error_weight 0..%d]\n", argv[0], T);
            return EXIT_FAILURE;
        }
        weight = (int)v;
    }

    /* Verify that x=2 is a primitive element of the selected field. */
    if (gf_pow(2, Q - 1u) != 1 || gf_pow(2, 1) == 1) {
        fprintf(stderr, "field polynomial / primitive-element check failed\n");
        return EXIT_FAILURE;
    }

    uint16_t g[T + 1];
    make_goppa_polynomial(g);
    printf("GF(2^%d), N=%u, t=%d, error weight=%d\n", M, N, T, weight);
    printf("Field polynomial: x^13+x^4+x^3+x+1\n");
    printf("g(x) coefficients (low degree first):\n");
    for (int i = 0; i <= T; ++i)
        printf("  g[%2d] = 0x%04x%s", i, g[i], (i % 4 == 3 || i == T) ? "\n" : "   ");

    uint8_t error[N] = {0}, decoded[N] = {0};
    int made = 0;
    /* Include coordinate zero when weight>0 to exercise its special case. */
    if (weight > 0) { error[0] = 1; made = 1; }
    while (made < weight) {
        unsigned a = 1u + (rnd32() % (Q - 1u));
        if (!error[a]) { error[a] = 1; ++made; }
    }

    uint16_t S[NSYN + 1], S2[NSYN + 1];
    make_syndromes(error, g, S, S2);

    /* S[1..2t] and S2[1..2t] omit the a=0 coordinate. */
    uint16_t seq[NSYN], seq2[NSYN], C[NSYN + 1], C2[NSYN + 1];
    for (int j = 0; j < NSYN; ++j) {
        seq[j] = S[j + 1];
        seq2[j] = S2[j + 1];
    }
    int L = berlekamp_massey(seq, NSYN, C);
    int L2 = berlekamp_massey(seq2, NSYN, C2);
    int pos[T + 1], pos2[T + 1];
    int count = chien_search(C, L, pos);
    int count2 = chien_search(C2, L2, pos2);

    printf("\nBM ordinary: locator degree=%d, roots=%d\n", L, count);
    printf("BM doubled : locator degree=%d, roots=%d\n", L2, count2);

    int ok = 1;
    if (L > T || L2 > T || count > T || count2 > T ||
        L != count || L2 != count2 || !same_positions(pos, count, pos2, count2)) {
        printf("FAIL: the two BM decoders did not produce the same nonzero positions\n");
        ok = 0;
    }

    int zero_error = recover_zero_coordinate(S[0], S2[0], pos, count, g);
    if (zero_error < 0) {
        printf("FAIL: zero-coordinate check is inconsistent\n");
        ok = 0;
    } else if (zero_error) {
        decoded[0] = 1;
    }

    for (int i = 0; i < count && i <= T; ++i) {
        if (pos[i] > 0 && pos[i] < (int)Q) decoded[pos[i]] = 1;
    }

    int original_weight = 0, decoded_weight = 0, mismatches = 0;
    printf("\nOriginal error positions:");
    for (unsigned i = 0; i < Q; ++i) {
        if (error[i]) { printf(" %u", i); ++original_weight; }
        if (decoded[i]) ++decoded_weight;
        if (error[i] != decoded[i]) ++mismatches;
    }
    printf("\nDecoded error positions :");
    for (unsigned i = 0; i < Q; ++i) if (decoded[i]) printf(" %u", i);
    printf("\nWeights: original=%d decoded=%d, mismatches=%d\n",
           original_weight, decoded_weight, mismatches);

    if (ok && mismatches == 0 && original_weight == weight) {
        printf("PASS: ordinary and doubled BM recover the exact binary error vector.\n");
        return EXIT_SUCCESS;
    }
    printf("FAIL\n");
    return EXIT_FAILURE;
}
