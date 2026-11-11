/*
 * 多点 Miura(C_ab) 符号 + Skorobogatov-Vladut(SV) 復号  (C言語版, 単一ファイル)
 *
 * 曲線:  y^a = f(x),  f(x) = prod_{i<b} (x - alpha_i),  gcd(a,b)=1,  p は a を割らない
 *        無限遠点 P_inf は1点、種数 g = (a-1)(b-1)/2
 *        Q_i = (alpha_i, 0) は有理点
 * 因子:  G = m0*P_inf + sum_i m_i*Q_i   (3点符号なら m1, m2 の2つを指定)
 * 符号:  C = { c : sum_k c_k h(P_k) = 0  (全ての h in L(G)) }   (= C_Omega(D,G))
 * 性能:  最小距離 >= deg G - 2g + 2,  SV法で t <= (deg G - 3g + 1)/2 個の誤りを訂正
 *
 * L(G) の基底:  y^j x^k / prod_i (x-alpha_i)^{n_i}
 *   n_i = floor((m_i + j)/a),  0 <= k <= sum n_i + floor((m0 - b j)/a),  0 <= j < a
 *
 * ビルド:  gcc -O2 -o miura miura_multipoint.c
 * 実行:    ./miura --p 2 --m 6 --a 3 --b 4 --m0 20 --ms 6,6 --trials 30
 *          ./miura --help でオプション一覧を表示
 * 注意:    加算表を q*q 個持つので、q <= 4096 程度までを想定しています。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXR 8

/* ------------------------------------------------------------ 有限体 GF(p^m) */
static int P, M, Q;
static int *EXPT, *LOGT, *NEGT, *ADDT;

static int gf_init(int p, int m)
{
    int q = 1, i, a, b;
    for (i = 0; i < m; i++) q *= p;
    P = p; M = m; Q = q;
    int *seq = malloc(sizeof(int) * (q > 1 ? q - 1 : 1));
    int *seen = malloc(sizeof(int) * q);
    int *c = calloc(m, sizeof(int));
    int *v = calloc(m, sizeof(int));
    int found = 0;
    long total = 1;
    for (i = 0; i < m; i++) total *= p;
    for (long idx = 0; idx < total && !found; idx++) {
        long t = idx;
        for (i = 0; i < m; i++) { c[i] = (int)(t % p); t /= p; }
        if (c[0] == 0) continue;
        memset(seen, 0, sizeof(int) * q);
        for (i = 0; i < m; i++) v[i] = 0;
        v[0] = 1;
        int ok = 1;
        for (int step = 0; step < q - 1; step++) {
            int val = 0, pw = 1;
            for (i = 0; i < m; i++) { val += v[i] * pw; pw *= p; }
            if (seen[val]) { ok = 0; break; }
            seen[val] = 1;
            seq[step] = val;
            int top = v[m - 1];
            for (i = m - 1; i > 0; i--) v[i] = v[i - 1];
            v[0] = 0;
            for (i = 0; i < m; i++) v[i] = (((v[i] - top * c[i]) % p) + p) % p;
        }
        if (ok) {
            int one = (v[0] == 1);
            for (i = 1; i < m; i++) if (v[i]) one = 0;
            if (one) found = 1;
        }
    }
    if (!found) { fprintf(stderr, "原始多項式が見つかりません\n"); return 0; }
    EXPT = malloc(sizeof(int) * 2 * (q - 1));
    LOGT = malloc(sizeof(int) * q);
    NEGT = malloc(sizeof(int) * q);
    for (i = 0; i < q - 1; i++) { EXPT[i] = EXPT[i + q - 1] = seq[i]; LOGT[seq[i]] = i; }
    LOGT[0] = -1;
    ADDT = NULL;
    for (a = 0; a < q; a++) {
        int val = 0, pw = 1, t = a;
        for (i = 0; i < m; i++) { int d = t % p; t /= p; val += ((p - d) % p) * pw; pw *= p; }
        NEGT[a] = val;
    }
    if (p != 2) {
        ADDT = malloc(sizeof(int) * (size_t)q * q);
        for (a = 0; a < q; a++)
            for (b = 0; b < q; b++) {
                int x = a, y = b, val = 0, pw = 1;
                for (i = 0; i < m; i++) { val += (((x % p) + (y % p)) % p) * pw; x /= p; y /= p; pw *= p; }
                ADDT[(size_t)a * q + b] = val;
            }
    }
    free(seq); free(seen); free(c); free(v);
    return 1;
}

static inline int gf_add(int a, int b) { return P == 2 ? (a ^ b) : ADDT[(size_t)a * Q + b]; }
static inline int gf_sub(int a, int b) { return gf_add(a, NEGT[b]); }
static inline int gf_neg(int a) { return NEGT[a]; }
static inline int gf_mul(int a, int b) { return (a == 0 || b == 0) ? 0 : EXPT[LOGT[a] + LOGT[b]]; }
static inline int gf_inv(int a) { return EXPT[(Q - 1) - LOGT[a]]; }
static int gf_pow(int a, int e)
{
    if (e == 0) return 1;
    if (a == 0) return 0;               /* e < 0 は呼ばない前提 */
    long k = ((long)LOGT[a] * e) % (Q - 1);
    if (k < 0) k += Q - 1;
    return EXPT[k];
}

static int floordiv(int a, int b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }
static int gcd_int(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

/* ------------------------------------------------------------ 線形代数 (行優先の平坦配列) */
/* A (rows x cols) を既約階段形に直す(その場で)。ピボット列を piv に入れ、階数を返す。 */
static int rref(int *A, int rows, int cols, int *piv)
{
    int r = 0;
    for (int c = 0; c < cols && r < rows; c++) {
        int pr = -1;
        for (int i = r; i < rows; i++) if (A[(size_t)i * cols + c]) { pr = i; break; }
        if (pr < 0) continue;
        if (pr != r)
            for (int j = 0; j < cols; j++) {
                int t = A[(size_t)r * cols + j];
                A[(size_t)r * cols + j] = A[(size_t)pr * cols + j];
                A[(size_t)pr * cols + j] = t;
            }
        int iv = gf_inv(A[(size_t)r * cols + c]);
        for (int j = 0; j < cols; j++) A[(size_t)r * cols + j] = gf_mul(iv, A[(size_t)r * cols + j]);
        for (int i = 0; i < rows; i++) {
            if (i == r) continue;
            int f = A[(size_t)i * cols + c];
            if (!f) continue;
            for (int j = 0; j < cols; j++)
                A[(size_t)i * cols + j] = gf_sub(A[(size_t)i * cols + j], gf_mul(f, A[(size_t)r * cols + j]));
        }
        piv[r++] = c;
    }
    return r;
}

/* A x = 0 の基底(nvec x cols)を返す。呼び出し側で free。 */
static int *nullspace(const int *A, int rows, int cols, int *nvec)
{
    int *R = malloc(sizeof(int) * (size_t)(rows ? rows : 1) * cols);
    int *piv = malloc(sizeof(int) * (cols + 1));
    int *ispiv = calloc(cols, sizeof(int));
    if (rows) memcpy(R, A, sizeof(int) * (size_t)rows * cols);
    int rk = rows ? rref(R, rows, cols, piv) : 0;
    for (int i = 0; i < rk; i++) ispiv[piv[i]] = 1;
    int nf = cols - rk;
    int *out = malloc(sizeof(int) * (size_t)(nf ? nf : 1) * cols);
    int vi = 0;
    for (int fc = 0; fc < cols; fc++) {
        if (ispiv[fc]) continue;
        int *v = out + (size_t)vi * cols;
        memset(v, 0, sizeof(int) * cols);
        v[fc] = 1;
        for (int i = 0; i < rk; i++) v[piv[i]] = gf_neg(R[(size_t)i * cols + fc]);
        vi++;
    }
    free(R); free(piv); free(ispiv);
    *nvec = nf;
    return out;
}

/* A (rows x n) x = s を解く。解がない/一意でないなら 0。解は x に入れる。 */
static int solve_unique(const int *A, int rows, int n, const int *s, int *x)
{
    int cols = n + 1;
    int *aug = malloc(sizeof(int) * (size_t)rows * cols);
    int *piv = malloc(sizeof(int) * (cols + 1));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < n; j++) aug[(size_t)i * cols + j] = A[(size_t)i * n + j];
        aug[(size_t)i * cols + n] = s[i];
    }
    int rk = rref(aug, rows, cols, piv);
    int ok = 1;
    for (int i = 0; i < rk; i++) if (piv[i] == n) ok = 0;
    if (ok && rk < n) ok = 0;
    if (ok) {
        for (int j = 0; j < n; j++) x[j] = 0;
        for (int i = 0; i < rk; i++) x[piv[i]] = aug[(size_t)i * cols + n];
    }
    free(aug); free(piv);
    return ok;
}

/* ------------------------------------------------------------ Miura 多点符号 */
typedef struct { int j, k; int ns[MAXR]; } Mono;

typedef struct {
    int a, b, m0, r, ms[MAXR], roots[256];
    int g, degG, tmax;
    int n, *px, *py;
    int nH, *H;          /* 検査行列 (nH x n) */
    int k;
} Code;

static Mono *basis_L(const Code *C, int m0, const int *ms, int *count)
{
    int a = C->a, b = C->b, r = C->r, cnt = 0;
    int j, i, k;
    for (j = 0; j < a; j++) {
        int N = 0;
        for (i = 0; i < r; i++) N += floordiv(ms[i] + j, a);
        int D = N + floordiv(m0 - b * j, a);
        if (D >= 0) cnt += D + 1;
    }
    Mono *B = malloc(sizeof(Mono) * (cnt ? cnt : 1));
    int t = 0;
    for (j = 0; j < a; j++) {
        int ns[MAXR], N = 0;
        for (i = 0; i < r; i++) { ns[i] = floordiv(ms[i] + j, a); N += ns[i]; }
        int D = N + floordiv(m0 - b * j, a);
        for (k = 0; k <= D; k++) {
            B[t].j = j; B[t].k = k;
            for (i = 0; i < r; i++) B[t].ns[i] = ns[i];
            t++;
        }
    }
    *count = cnt;
    return B;
}

static int *evaluate(const Code *C, const Mono *B, int cnt)
{
    int *E = malloc(sizeof(int) * (size_t)(cnt ? cnt : 1) * C->n);
    for (int t = 0; t < cnt; t++)
        for (int p = 0; p < C->n; p++) {
            int x = C->px[p], y = C->py[p];
            int v = gf_mul(gf_pow(y, B[t].j), gf_pow(x, B[t].k));
            for (int i = 0; i < C->r; i++)
                if (B[t].ns[i]) v = gf_mul(v, gf_pow(gf_sub(x, C->roots[i]), -B[t].ns[i]));
            E[(size_t)t * C->n + p] = v;
        }
    return E;
}

static unsigned long long rng_state = 88172645463325252ULL;
static unsigned rnd(void)
{
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return (unsigned)(rng_state >> 11);
}

static int code_init(Code *C, int a, int b, int m0, int r, const int *ms,
                     const int *roots_in, unsigned seed)
{
    int i;
    memset(C, 0, sizeof(*C));
    if (gcd_int(a, b) != 1) { fprintf(stderr, "gcd(a,b)=1 が必要です\n"); return 0; }
    if (a % P == 0) { fprintf(stderr, "標数 p が a を割らないようにしてください\n"); return 0; }
    if (b > Q || b > 256) { fprintf(stderr, "b <= q (かつ 256 以下) が必要です\n"); return 0; }
    if (r > b || r > MAXR) { fprintf(stderr, "ms の個数が多すぎます\n"); return 0; }
    for (i = 0; i < r; i++) if (ms[i] < 0) { fprintf(stderr, "ms は非負\n"); return 0; }
    C->a = a; C->b = b; C->m0 = m0; C->r = r;
    for (i = 0; i < r; i++) C->ms[i] = ms[i];
    rng_state ^= (unsigned long long)seed * 2654435761ULL + 1;
    /* f の根: 指定があればそれを使い、なければ GF(q) から b 個を重複なく選ぶ */
    if (roots_in) {
        for (i = 0; i < b; i++) {
            if (roots_in[i] < 0 || roots_in[i] >= Q) { fprintf(stderr, "根は 0..q-1 の範囲で指定してください\n"); return 0; }
            for (int j = 0; j < i; j++)
                if (roots_in[j] == roots_in[i]) { fprintf(stderr, "根が重複しています\n"); return 0; }
            C->roots[i] = roots_in[i];
        }
    } else {
        int *perm = malloc(sizeof(int) * Q);
        for (i = 0; i < Q; i++) perm[i] = i;
        for (i = Q - 1; i > 0; i--) { int j = rnd() % (i + 1); int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
        for (i = 0; i < b; i++) C->roots[i] = perm[i];
        free(perm);
    }
    C->g = (a - 1) * (b - 1) / 2;
    C->degG = m0;
    for (i = 0; i < r; i++) C->degG += ms[i];
    C->tmax = (C->degG - 3 * C->g + 1) / 2;
    if (C->degG - 3 * C->g + 1 < 0) C->tmax = 0;
    /* 有理点 */
    C->px = malloc(sizeof(int) * (size_t)Q * a);
    C->py = malloc(sizeof(int) * (size_t)Q * a);
    int *pa = malloc(sizeof(int) * Q);
    for (int y = 0; y < Q; y++) pa[y] = gf_pow(y, a);
    pa[0] = 0;
    for (int x = 0; x < Q; x++) {
        int fx = 1;
        for (i = 0; i < b; i++) fx = gf_mul(fx, gf_sub(x, C->roots[i]));
        for (int y = 0; y < Q; y++) {
            if (pa[y] != fx) continue;
            int excl = 0;
            for (i = 0; i < r; i++) if (ms[i] > 0 && x == C->roots[i] && y == 0) excl = 1;
            if (excl) continue;
            C->px[C->n] = x; C->py[C->n] = y; C->n++;
        }
    }
    free(pa);
    /* 検査行列: L(G) の基底を評価し、既約階段形にして独立な行だけ残す */
    int cnt;
    Mono *B = basis_L(C, m0, ms, &cnt);
    int *E = evaluate(C, B, cnt);
    int *piv = malloc(sizeof(int) * (C->n + 1));
    int rk = cnt ? rref(E, cnt, C->n, piv) : 0;
    C->nH = rk;
    C->H = malloc(sizeof(int) * (size_t)(rk ? rk : 1) * C->n);
    if (rk) memcpy(C->H, E, sizeof(int) * (size_t)rk * C->n);
    C->k = C->n - rk;
    free(B); free(E); free(piv);
    return 1;
}

static void code_free(Code *C) { free(C->px); free(C->py); free(C->H); }

static void syndrome(const Code *C, const int *r, int *s)
{
    for (int i = 0; i < C->nH; i++) {
        int acc = 0;
        for (int k = 0; k < C->n; k++) {
            int h = C->H[(size_t)i * C->n + k];
            if (h && r[k]) acc = gf_add(acc, gf_mul(h, r[k]));
        }
        s[i] = acc;
    }
}

static void random_codeword(const Code *C, int *c)
{
    int nv, i, k;
    int *V = nullspace(C->H, C->nH, C->n, &nv);
    for (k = 0; k < C->n; k++) c[k] = 0;
    for (i = 0; i < nv; i++) {
        int co = rnd() % Q;
        if (!co) continue;
        for (k = 0; k < C->n; k++) c[k] = gf_add(c[k], gf_mul(co, V[(size_t)i * C->n + k]));
    }
    free(V);
}

/* SV復号。成功なら 1 を返し、out に復号語、*nerr に誤り個数を入れる。 */
static int decode(const Code *C, const int *r, int *out, int *nerr)
{
    int n = C->n, g = C->g, i, k;
    int *s = malloc(sizeof(int) * (C->nH ? C->nH : 1));
    syndrome(C, r, s);
    int zero = 1;
    for (i = 0; i < C->nH; i++) if (s[i]) zero = 0;
    if (zero) { memcpy(out, r, sizeof(int) * n); *nerr = 0; free(s); return 1; }

    int zeros[MAXR] = {0};
    for (int t = 1; t <= C->tmax; t++) {
        int degF = t + g, nF, nFp;
        Mono *bF = basis_L(C, degF, zeros, &nF);
        Mono *bFp = basis_L(C, C->m0 - degF, C->ms, &nFp);
        if (nF == 0 || nFp == 0) { free(bF); free(bFp); continue; }
        int *EF = evaluate(C, bF, nF);
        int *EFp = evaluate(C, bFp, nFp);
        /* M[j][i] = sum_k r_k EFp[j][k] EF[i][k] */
        int *Mx = malloc(sizeof(int) * (size_t)nFp * nF);
        int *w = malloc(sizeof(int) * n);
        for (int j = 0; j < nFp; j++) {
            for (k = 0; k < n; k++) w[k] = gf_mul(r[k], EFp[(size_t)j * n + k]);
            for (i = 0; i < nF; i++) {
                int acc = 0;
                for (k = 0; k < n; k++) {
                    int wk = w[k], fk = EF[(size_t)i * n + k];
                    if (wk && fk) acc = gf_add(acc, gf_mul(wk, fk));
                }
                Mx[(size_t)j * nF + i] = acc;
            }
        }
        int nk;
        int *ker = nullspace(Mx, nFp, nF, &nk);
        int success = 0;
        if (nk > 0) {
            /* 核の各関数の値を計算し、全てが 0 になる位置の集合 Z を求める */
            char *inZ = malloc(n);
            memset(inZ, 1, n);
            for (int v = 0; v < nk; v++)
                for (k = 0; k < n; k++) {
                    int val = 0;
                    for (i = 0; i < nF; i++) {
                        int ci = ker[(size_t)v * nF + i];
                        if (ci) val = gf_add(val, gf_mul(ci, EF[(size_t)i * n + k]));
                    }
                    if (val) inZ[k] = 0;
                }
            int *Z = malloc(sizeof(int) * n), nz = 0;
            for (k = 0; k < n; k++) if (inZ[k]) Z[nz++] = k;
            if (nz > 0 && nz <= t + g) {
                int *Hz = malloc(sizeof(int) * (size_t)C->nH * nz);
                for (i = 0; i < C->nH; i++)
                    for (int z = 0; z < nz; z++) Hz[(size_t)i * nz + z] = C->H[(size_t)i * n + Z[z]];
                int *e = malloc(sizeof(int) * nz);
                if (solve_unique(Hz, C->nH, nz, s, e)) {
                    int wt = 0;
                    for (int z = 0; z < nz; z++) if (e[z]) wt++;
                    if (wt <= C->tmax) {
                        memcpy(out, r, sizeof(int) * n);
                        for (int z = 0; z < nz; z++) if (e[z]) out[Z[z]] = gf_sub(out[Z[z]], e[z]);
                        *nerr = wt;
                        success = 1;
                    }
                }
                free(Hz); free(e);
            }
            free(inZ); free(Z);
        }
        free(ker); free(Mx); free(w); free(EF); free(EFp); free(bF); free(bFp);
        if (success) { free(s); return 1; }
    }
    free(s);
    return 0;
}

/* ------------------------------------------------------------ コマンドライン */
static void usage(const char *prog)
{
    printf(
"使い方: %s [オプション]\n"
"\n"
"  --p N         体の標数 p                       (既定 2)\n"
"  --m N         拡大次数。体は GF(p^m)           (既定 6)\n"
"  --a N         曲線 y^a = f(x) の a             (既定 3)\n"
"  --b N         f の次数 b (gcd(a,b)=1)          (既定 4)\n"
"  --m0 N        無限遠点 P_inf の重み            (既定 20)\n"
"  --ms L        Q_i の重み。カンマ区切り         (既定 6,6)  例: --ms 4,5,3\n"
"                  2個なら3点符号(P_inf,Q1,Q2)、1個なら2点、空(--ms 0)なら1点符号\n"
"  --roots L     f の根 alpha_i を b 個、カンマ区切りで指定 (体の元は 0..q-1 の整数)\n"
"                  省略すると乱数で選ぶ\n"
"  --errors T    誤り個数を T に固定してテスト    (省略: 0..t_max+2 を順に試す)\n"
"  --trials N    各誤り個数での試行回数           (既定 30)\n"
"  --seed N      乱数の種                         (既定 1)\n"
"  --help        この説明を表示\n"
"\n"
"値は「--p 2」でも「--p=2」でも指定できます。\n"
"例: %s --p 2 --m 6 --a 3 --b 4 --m0 20 --ms 6,6 --trials 30\n",
        prog, prog);
}

/* "1,2,3" を整数の配列に読む。個数を返す(失敗は -1)。 */
static int parse_list(const char *str, int *out, int maxn)
{
    int n = 0;
    const char *p = str;
    while (*p) {
        char *end;
        long v = strtol(p, &end, 10);
        if (end == p) return -1;
        if (n >= maxn) return -1;
        out[n++] = (int)v;
        p = end;
        if (*p == ',') p++;
        else if (*p) return -1;
    }
    return n;
}

static int parse_int(const char *name, const char *str, int *out)
{
    char *end;
    long v = strtol(str, &end, 10);
    if (*str == '\0' || *end != '\0') { fprintf(stderr, "%s には整数を指定してください: %s\n", name, str); return 0; }
    *out = (int)v;
    return 1;
}

/* ------------------------------------------------------------ 動作テスト */
int main(int argc, char **argv)
{
    int p = 2, m = 6, a = 3, b = 4, m0 = 20, trials = 30, errors = -1, seed = 1;
    int ms[MAXR] = { 6, 6 }, r = 2;
    int roots[256], nroots = 0;

    for (int i = 1; i < argc; i++) {
        char *arg = argv[i], *val = NULL;
        if (strncmp(arg, "--", 2) != 0) { fprintf(stderr, "不明な引数: %s\n", arg); usage(argv[0]); return 1; }
        char *eq = strchr(arg, '=');
        char name[32];
        size_t len = eq ? (size_t)(eq - arg) : strlen(arg);
        if (len >= sizeof(name)) { fprintf(stderr, "不明なオプション: %s\n", arg); return 1; }
        memcpy(name, arg, len); name[len] = '\0';
        if (strcmp(name, "--help") == 0 || strcmp(name, "-h") == 0) { usage(argv[0]); return 0; }
        if (eq) val = eq + 1;
        else if (i + 1 < argc) val = argv[++i];
        else { fprintf(stderr, "%s には値が必要です\n", name); return 1; }

        if (!strcmp(name, "--p")) { if (!parse_int(name, val, &p)) return 1; }
        else if (!strcmp(name, "--m")) { if (!parse_int(name, val, &m)) return 1; }
        else if (!strcmp(name, "--a")) { if (!parse_int(name, val, &a)) return 1; }
        else if (!strcmp(name, "--b")) { if (!parse_int(name, val, &b)) return 1; }
        else if (!strcmp(name, "--m0")) { if (!parse_int(name, val, &m0)) return 1; }
        else if (!strcmp(name, "--trials")) { if (!parse_int(name, val, &trials)) return 1; }
        else if (!strcmp(name, "--errors")) { if (!parse_int(name, val, &errors)) return 1; }
        else if (!strcmp(name, "--seed")) { if (!parse_int(name, val, &seed)) return 1; }
        else if (!strcmp(name, "--ms")) {
            r = parse_list(val, ms, MAXR);
            if (r < 0) { fprintf(stderr, "--ms は整数のカンマ区切りで、%d 個以下にしてください\n", MAXR); return 1; }
            if (r == 1 && ms[0] == 0) r = 0;       /* "--ms 0" は1点符号 */
        }
        else if (!strcmp(name, "--roots")) {
            nroots = parse_list(val, roots, 256);
            if (nroots < 0) { fprintf(stderr, "--roots は整数のカンマ区切りで指定してください\n"); return 1; }
        }
        else { fprintf(stderr, "不明なオプション: %s\n", name); usage(argv[0]); return 1; }
    }

    if (p < 2 || m < 1 || a < 2 || b < 1 || trials < 1) {
        fprintf(stderr, "パラメーターが不正です(p>=2, m>=1, a>=2, b>=1, trials>=1)\n");
        return 1;
    }
    for (int d = 2; d * d <= p; d++)
        if (p % d == 0) { fprintf(stderr, "p は素数にしてください\n"); return 1; }
    if (nroots && nroots != b) {
        fprintf(stderr, "--roots は b=%d 個必要です(指定は %d 個)\n", b, nroots);
        return 1;
    }
    {
        double qd = 1; for (int i = 0; i < m; i++) qd *= p;
        if (qd > 4096) { fprintf(stderr, "q=p^m が大きすぎます(4096 以下にしてください)\n"); return 1; }
    }

    if (!gf_init(p, m)) return 1;
    Code C;
    if (!code_init(&C, a, b, m0, r, ms, nroots ? roots : NULL, (unsigned)seed)) return 1;

    printf("GF(%d^%d) 上の y^%d=f(x), deg f=%d, 種数 g=%d\n", p, m, a, b, C.g);
    printf("f の根:");
    for (int i = 0; i < b; i++) printf(" %d", C.roots[i]);
    printf("\nG = %d*P_inf", m0);
    for (int i = 0; i < r; i++) printf(" + %d*Q%d", ms[i], i + 1);
    printf("  (deg G=%d, %d点符号)\n", C.degG, r + 1);
    printf("符号長 n=%d, 次元 k=%d, 設計最小距離 >= %d, 訂正可能 t <= %d\n",
           C.n, C.k, C.degG - 2 * C.g + 2, C.tmax);
    if (C.k <= 0) { printf("次元が0以下です。m0, ms を小さくしてください。\n"); code_free(&C); return 0; }

    int n = C.n;
    int *c = malloc(sizeof(int) * n), *rr = malloc(sizeof(int) * n);
    int *out = malloc(sizeof(int) * n), *idx = malloc(sizeof(int) * n);
    int tlo = 0, thi = C.tmax + 2;
    if (errors >= 0) tlo = thi = errors;
    for (int t = tlo; t <= thi && t <= n; t++) {
        int good = 0;
        for (int tr = 0; tr < trials; tr++) {
            random_codeword(&C, c);
            memcpy(rr, c, sizeof(int) * n);
            for (int i = 0; i < n; i++) idx[i] = i;
            for (int i = 0; i < t; i++) {                       /* t 個の位置を選ぶ */
                int j = i + rnd() % (n - i);
                int tmp = idx[i]; idx[i] = idx[j]; idx[j] = tmp;
                rr[idx[i]] = gf_add(rr[idx[i]], 1 + rnd() % (Q - 1));
            }
            int ne;
            if (decode(&C, rr, out, &ne) && memcmp(out, c, sizeof(int) * n) == 0) good++;
        }
        printf("  誤り %2d 個: %d/%d 成功%s\n", t, good, trials, t > C.tmax ? "  (保証範囲外)" : "");
    }
    free(c); free(rr); free(out); free(idx);
    code_free(&C);
    return 0;
}

