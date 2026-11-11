/*
 * intersecting_goppa.c
 *
 * Small binary Goppa product-code experiment:
 *   G(x) = g1(x) * g2(x)
 *   C    = Gamma(L, G)
 *
 * Builds binary parity-check matrices for g1, g2 and G, constructs
 * codewords from ker(H_G), then independently tries bounded-weight
 * syndrome decoding with H_g1 and H_g2.
 *
 * This is an EXHAUSTIVE small-parameter decoder, not Patterson.
 * It is intended to test the construction, not for cryptographic use.
 *
 * Compile:
 *   gcc -O2 -std=c11 -Wall -Wextra intersecting_goppa.c -o intersecting_goppa
 * Run:
 *   ./intersecting_goppa
 *
 * Defaults:
 *   GF(2^4), modulus x^4+x+1 (0x13)
 *   g1 = x^2+x+1
 *   g2 = x+2  (the field element 2)
 *
 * Limits: n <= 63, m*deg(g) <= 20 for the syndrome lookup table.
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 63
#define MAX_M 8
#define MAX_DEG 8
#define MAX_ROWS 20
#define MAX_BASIS 63

typedef uint64_t u64;

static int M = 4;
static int IRR = 0x13;
static int Q;
static int g1[MAX_DEG + 1] = {1, 1, 1};
static int d1 = 2;
static int g2[MAX_DEG + 1] = {2, 1, 1};
static int d2 = 2;

static int gf_mul(int a, int b) {
    int out = 0;
    int top = 1 << M;
    while (b) {
        if (b & 1) out ^= a;
        b >>= 1;
        a <<= 1;
        if (a & top) a ^= IRR;
    }
    return out;
}

static int gf_pow(int a, int e) {
    int out = 1;
    while (e > 0) {
        if (e & 1) out = gf_mul(out, a);
        a = gf_mul(a, a);
        e >>= 1;
    }
    return out;
}

static int gf_inv(int a) {
    if (a == 0) {
        fprintf(stderr, "ERROR: inverse of zero\n");
        exit(1);
    }
    return gf_pow(a, Q - 2);
}

static int poly_eval(const int *p, int degree, int x) {
    int v = 0;
    for (int i = degree; i >= 0; --i)
        v = gf_mul(v, x) ^ p[i];
    return v;
}

static void poly_mul(const int *a, int da, const int *b, int db,
                     int *out, int *dout) {
    for (int i = 0; i <= da + db; ++i) out[i] = 0;
    for (int i = 0; i <= da; ++i)
        for (int j = 0; j <= db; ++j)
            out[i + j] ^= gf_mul(a[i], b[j]);
    *dout = da + db;
    while (*dout > 0 && out[*dout] == 0) --(*dout);
}

static void print_poly(const char *name, const int *p, int d) {
    printf("%s(x) = ", name);
    int first = 1;
    for (int i = d; i >= 0; --i) {
        if (!p[i]) continue;
        if (!first) printf(" + ");
        if (i == 0) printf("%d", p[i]);
        else if (i == 1) printf("%d*x", p[i]);
        else printf("%d*x^%d", p[i], i);
        first = 0;
    }
    if (first) printf("0");
    printf("\n");
}

/* Rows are binary masks; column i corresponds to support[i]. */
static int build_check(const int *g, int deg, const int *support, int n,
                       u64 *rows) {
    int nr = M * deg;
    if (nr > MAX_ROWS) {
        fprintf(stderr, "Too many check rows (%d > %d)\n", nr, MAX_ROWS);
        exit(1);
    }
    for (int r = 0; r < nr; ++r) rows[r] = 0;

    for (int col = 0; col < n; ++col) {
        int a = support[col];
        int gv = poly_eval(g, deg, a);
        if (gv == 0) {
            fprintf(stderr, "Internal error: support contains a root\n");
            exit(1);
        }
        int invg = gf_inv(gv);
        int power = 1;
        for (int j = 0; j < deg; ++j) {
            int z = gf_mul(power, invg);
            for (int b = 0; b < M; ++b)
                if ((z >> b) & 1)
                    rows[j * M + b] |= (u64)1 << col;
            power = gf_mul(power, a);
        }
    }
    return nr;
}

static int parity(u64 x) {
    return __builtin_parityll(x);
}

static u64 syndrome(u64 v, const u64 *rows, int nr) {
    u64 s = 0;
    for (int i = 0; i < nr; ++i)
        if (parity(rows[i] & v)) s |= (u64)1 << i;
    return s;
}

/* Row-echelon reduction; pivots are highest set bits. */
static int row_reduce(const u64 *input, int nr, int n, u64 *basis_rows,
                      int *pivot_cols) {
    int rank = 0;
    for (int i = 0; i < nr; ++i) {
        u64 v = input[i];
        for (;;) {
            if (!v) break;
            int p = 63 - __builtin_clzll(v);
            int found = -1;
            for (int j = 0; j < rank; ++j)
                if (pivot_cols[j] == p) { found = j; break; }
            if (found >= 0) v ^= basis_rows[found];
            else {
                basis_rows[rank] = v;
                pivot_cols[rank] = p;
                ++rank;
                break;
            }
        }
    }
    (void)n;
    return rank;
}

static int make_nullspace(const u64 *rows, int nr, int n,
                          u64 *vectors, int *dimension) {
    u64 rr[MAX_ROWS];
    int piv[MAX_ROWS];
    int rank = row_reduce(rows, nr, n, rr, piv);
    int is_pivot[MAX_N] = {0};
    for (int i = 0; i < rank; ++i) is_pivot[piv[i]] = 1;

    int k = 0;
    for (int f = 0; f < n; ++f) {
        if (is_pivot[f]) continue;
        u64 v = (u64)1 << f;
        /* Highest-pivot equations depend only on lower columns:
           solve pivots in ascending order. */
        for (int p = 0; p < n; ++p) {
            for (int r = 0; r < rank; ++r) {
                if (piv[r] == p && parity(rr[r] & v))
                    v |= (u64)1 << p;
            }
        }
        vectors[k++] = v;
    }
    *dimension = k;
    return rank;
}

/* Store one error mask for every syndrome, if found. */
typedef struct {
    u64 *error;
    unsigned char *found;
    int syndrome_count;
    int collisions;
} Decoder;

static void add_error(Decoder *dec, u64 e, u64 s) {
    if (!dec->found[s]) {
        dec->found[s] = 1;
        dec->error[s] = e;
    } else if (dec->error[s] != e) {
        dec->collisions++;
    }
}

static void enumerate_errors_rec(Decoder *dec, const u64 *rows, int nr,
                                 int n, int start, int left, u64 e) {
    add_error(dec, e, syndrome(e, rows, nr));
    if (left == 0) return;
    for (int i = start; i <= n - left; ++i)
        enumerate_errors_rec(dec, rows, nr, n, i + 1, left - 1,
                             e | ((u64)1 << i));
}

static Decoder make_decoder(const u64 *rows, int nr, int n, int radius) {
    Decoder dec;
    dec.syndrome_count = 1 << nr;
    dec.error = (u64 *)calloc((size_t)dec.syndrome_count, sizeof(u64));
    dec.found = (unsigned char *)calloc((size_t)dec.syndrome_count, 1);
    dec.collisions = 0;
    if (!dec.error || !dec.found) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    for (int w = 0; w <= radius; ++w)
        enumerate_errors_rec(&dec, rows, nr, n, 0, w, 0);
    return dec;
}

static int decode(const Decoder *dec, u64 received,
                  const u64 *rows, int nr, u64 *candidate) {
    u64 s = syndrome(received, rows, nr);
    if (!dec->found[s]) return 0;
    *candidate = received ^ dec->error[s];
    return 1;
}

static int is_codeword(u64 v, const u64 *rows, int nr) {
    return syndrome(v, rows, nr) == 0;
}

static u64 random_codeword(const u64 *basis, int k) {
    u64 c = 0;
    for (int i = 0; i < k; ++i)
        if (rand() & 1) c ^= basis[i];
    return c;
}

static u64 random_error(int n, int weight) {
    u64 e = 0;
    while (__builtin_popcountll(e) < weight) {
        int p = rand() % n;
        e |= (u64)1 << p;
    }
    return e;
}

static void print_bits(u64 v, int n) {
    for (int i = n - 1; i >= 0; --i)
        putchar((v >> i) & 1 ? '1' : '0');
}

int main(int argc, char **argv) {
    int trials = 30;
    int error_weight = -1;
    unsigned seed = 7;

    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--trials") && i + 1 < argc)
            trials = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--errors") && i + 1 < argc)
            error_weight = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc)
            seed = (unsigned)strtoul(argv[++i], NULL, 10);
        else {
            fprintf(stderr, "Usage: %s [--trials N] [--errors W] [--seed N]\n",
                    argv[0]);
            return 1;
        }
    }

    Q = 1 << M;
    srand(seed);

    int G[MAX_DEG + 1];
    int dG;
    poly_mul(g1, d1, g2, d2, G, &dG);

    int support[MAX_N];
    int n = 0;
    for (int a = 0; a < Q; ++a)
        if (poly_eval(G, dG, a) != 0)
            support[n++] = a;

    if (n == 0 || n > MAX_N) {
        fprintf(stderr, "Unsupported support length n=%d\n", n);
        return 1;
    }

    u64 H1[MAX_ROWS], H2[MAX_ROWS], HG[MAX_ROWS];
    int r1 = build_check(g1, d1, support, n, H1);
    int r2 = build_check(g2, d2, support, n, H2);
    int rG = build_check(G, dG, support, n, HG);

    u64 basis[MAX_BASIS];
    int k = 0;
    int rankG = make_nullspace(HG, rG, n, basis, &k);
    u64 tmp[MAX_ROWS];
    int pivtmp[MAX_ROWS];
    int rank1 = row_reduce(H1, r1, n, tmp, pivtmp);
    int rank2 = row_reduce(H2, r2, n, tmp, pivtmp);

    printf("=== Product Goppa code ===\n");
    printf("GF(2^%d), modulus = 0x%x, q=%d\n", M, IRR, Q);
    print_poly("g1", g1, d1);
    print_poly("g2", g2, d2);
    print_poly("G = g1*g2", G, dG);
    printf("support length n = %d\n", n);
    printf("rank(H1)=%d, rank(H2)=%d, rank(HG)=%d\n",
           rank1, rank2, rankG);
    printf("dimension k=%d, rate=%.4f\n", k, (double)k / n);

    int inclusion_ok = 1;
    for (int i = 0; i < k; ++i)
        if (!is_codeword(basis[i], H1, r1) ||
            !is_codeword(basis[i], H2, r2))
            inclusion_ok = 0;
    printf("Every basis vector of ker(HG) is in ker(H1) and ker(H2): %s\n",
           inclusion_ok ? "YES" : "NO");

    int t1 = d1;
    int t2 = d2;
    if (error_weight < 0) error_weight = (t1 > t2 ? t1 : t2);
    if (error_weight > n || error_weight < 0) {
        fprintf(stderr, "Invalid error weight\n");
        return 1;
    }
    if (r1 > MAX_ROWS || r2 > MAX_ROWS || rG > MAX_ROWS) return 1;

    printf("\nBuilding exhaustive syndrome tables (radii %d and %d)...\n",
           t1, t2);
    Decoder dec1 = make_decoder(H1, r1, n, t1);
    Decoder dec2 = make_decoder(H2, r2, n, t2);
    printf("g1 table collisions=%d; g2 table collisions=%d\n",
           dec1.collisions, dec2.collisions);

    int ok1 = 0, ok2 = 0, either = 0;
    printf("\n=== Trials: %d, error weight: %d ===\n", trials, error_weight);
    for (int tr = 0; tr < trials; ++tr) {
        u64 c = random_codeword(basis, k);
        u64 e = random_error(n, error_weight);
        u64 received = c ^ e;
        u64 c1 = 0, c2 = 0;
        int found1 = decode(&dec1, received, H1, r1, &c1);
        int found2 = decode(&dec2, received, H2, r2, &c2);
        int good1 = found1 && c1 == c && is_codeword(c1, HG, rG);
        int good2 = found2 && c2 == c && is_codeword(c2, HG, rG);
        ok1 += good1;
        ok2 += good2;
        either += (good1 || good2);
        printf("%02d: ", tr + 1);
        print_bits(e, n);
        printf("  g1=%s g2=%s either=%s\n",
               good1 ? "OK" : "FAIL",
               good2 ? "OK" : "FAIL",
               (good1 || good2) ? "OK" : "FAIL");
    }

    printf("\n=== Summary ===\n");
    printf("g1 decoder: %d/%d\n", ok1, trials);
    printf("g2 decoder: %d/%d\n", ok2, trials);
    printf("either decoder: %d/%d\n", either, trials);
    printf("\nNOTE: exhaustive bounded-weight decoding, not Patterson.\n");
    printf("A component decoder may return a word outside the product code;\n");
    printf("the program checks membership in C=Gamma(L,g1*g2).\n");

    free(dec1.error); free(dec1.found);
    free(dec2.error); free(dec2.found);
    return 0;
}

