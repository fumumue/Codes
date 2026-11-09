// wild Goppa 符号 Γ_p(L, g^{p-1}) (= Γ_p(L, g^p)) の復号実験。L = GF(p^m) 全体 (n = q = p^m)。
//  - 符号の構成: 法 g^{p-1} の検査行列を F_p に展開して零空間から符号語 c をランダムに取る
//  - 受信語 r = c + e (e は F_p 上、重み w)。復号側は「法 g^p のシンドローム」を計算し、
//    鍵方程式 sigma*S = omega (mod g^p) を拡張ユークリッドで解く (deg g^p / 2 = floor(p*t/2) 個まで)
//  gcc -O2 -o wild_goppa wild_goppa.c
//  ./wild_goppa p m t w trials [old] [verbose] [errors]
//    old     : 1 なら法 g^{p-1} の従来方式で復号 (0 なら wild: 法 g^p)
//    verbose : 最初の verbose 回分について、入れた誤りと復号結果を表示
//    errors  : 誤りを手で指定 例 "3:1,10:2,50:1" (位置は元の番号 0..q-1、値は F_p)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXD 1024

static int p, m, q;
static int *EXP_, *LOG_, *NEGT;
static uint16_t *ADDT;
static int PRIMC = 0;   // 体の定義多項式 (係数を p 進数で並べた番号)

static uint64_t rs = 88172645463325252ULL;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

/* ---------- GF(p^m) ---------- */
static void todig(int a, int *d) { for (int j = 0; j < m; j++) { d[j] = a % p; a /= p; } }
static int fromdig(const int *d) { int a = 0; for (int j = m - 1; j >= 0; j--) a = a * p + d[j]; return a; }

static inline int fadd(int a, int b) { return ADDT[(size_t)a * q + b]; }
static inline int fsub(int a, int b) { return fadd(a, NEGT[b]); }
static inline int fmul(int a, int b) {
    if (!a || !b) return 0;
    int s = LOG_[a] + LOG_[b]; if (s >= q - 1) s -= q - 1;
    return EXP_[s];
}
static inline int finv(int a) { return EXP_[(q - 1 - LOG_[a]) % (q - 1)]; }

static int try_prim(int c) {
    int f[16], d[16] = {0}, nd[16];
    todig(c, f);
    if (f[0] == 0) return 0;
    d[0] = 1;
    for (int i = 0; i < q - 1; i++) {
        EXP_[i] = fromdig(d);
        if (i > 0 && EXP_[i] == 1) return 0;
        int top = d[m - 1];
        for (int k = m - 1; k >= 1; k--) nd[k] = (((d[k - 1] - top * f[k]) % p) + p) % p;
        nd[0] = (((-top * f[0]) % p) + p) % p;
        memcpy(d, nd, sizeof(int) * m);
    }
    return fromdig(d) == 1;
}

static void build_field(void) {
    EXP_ = malloc(sizeof(int) * q); LOG_ = malloc(sizeof(int) * q);
    NEGT = malloc(sizeof(int) * q); ADDT = malloc(sizeof(uint16_t) * (size_t)q * q);
    int ok = 0;
    for (int c = 1; c < q && !ok; c++) { ok = try_prim(c); if (ok) PRIMC = c; }
    if (!ok) { fprintf(stderr, "no primitive poly\n"); exit(1); }
    LOG_[0] = 0;
    for (int i = 0; i < q - 1; i++) LOG_[EXP_[i]] = i;
    int da[16], db[16], dc[16];
    for (int a = 0; a < q; a++) {
        todig(a, da);
        for (int j = 0; j < m; j++) dc[j] = (p - da[j]) % p;
        NEGT[a] = fromdig(dc);
        for (int b = 0; b < q; b++) {
            todig(b, db);
            for (int j = 0; j < m; j++) dc[j] = (da[j] + db[j]) % p;
            ADDT[(size_t)a * q + b] = (uint16_t)fromdig(dc);
        }
    }
}

/* ---------- 多項式 over GF(q) ---------- */
static int pdeg(const int *a, int n) { while (n > 0 && a[n - 1] == 0) n--; return n - 1; }

static void pmul(const int *a, int na, const int *b, int nb, int *o) {
    for (int i = 0; i < na + nb - 1; i++) o[i] = 0;
    for (int i = 0; i < na; i++) if (a[i])
        for (int j = 0; j < nb; j++) if (b[j]) o[i + j] = fadd(o[i + j], fmul(a[i], b[j]));
}

static void pdivmod(int *a, int na, const int *b, int db, int *qt) {
    int ib = finv(b[db]);
    for (int i = 0; i < na; i++) qt[i] = 0;
    for (int da = pdeg(a, na); da >= db; da = pdeg(a, na)) {
        int c = fmul(a[da], ib);
        qt[da - db] = c;
        for (int k = 0; k <= db; k++) a[da - db + k] = fsub(a[da - db + k], fmul(c, b[k]));
    }
}

static int peval(const int *a, int n, int x) {
    int r = 0;
    for (int k = n - 1; k >= 0; k--) r = fadd(fmul(r, x), a[k]);
    return r;
}

static int squarefree(const int *g, int t) {
    int a[MAXD], b[MAXD], tmp[MAXD], qt[MAXD];
    memset(a, 0, sizeof a); memset(b, 0, sizeof b);
    for (int i = 0; i <= t; i++) a[i] = g[i];
    for (int i = 1; i <= t; i++) b[i - 1] = fmul(i % p, g[i]);
    int n = t + 2;
    while (pdeg(b, n) >= 0) {
        pdivmod(a, n, b, pdeg(b, n), qt);
        memcpy(tmp, a, sizeof(int) * n); memcpy(a, b, sizeof(int) * n); memcpy(b, tmp, sizeof(int) * n);
    }
    return pdeg(a, n) == 0;
}

/* ---------- 符号 ---------- */
static int T, G[MAXD];
static int D1, D2, GQ1[MAXD], GQ2[MAXD];   // g^{p-1}, g^p
static int *INV1, *INV2;                   // (x-L_i)^{-1} mod GQ1 / GQ2

static int *build_inv(const int *M, int D) {
    int *out = malloc(sizeof(int) * (size_t)q * D);
    for (int i = 0; i < q; i++) {
        int b[MAXD];
        b[D - 1] = M[D];
        for (int k = D - 1; k >= 1; k--) b[k - 1] = fadd(M[k], fmul(i, b[k]));
        int rem = fadd(M[0], fmul(i, b[0]));
        if (rem == 0) { fprintf(stderr, "root in L\n"); exit(1); }
        int c = NEGT[finv(rem)];
        for (int k = 0; k < D; k++) out[(size_t)i * D + k] = fmul(c, b[k]);
    }
    return out;
}

static void syndrome(const int *v, const int *INV, int D, int *s) {
    for (int k = 0; k < D; k++) s[k] = 0;
    for (int i = 0; i < q; i++) if (v[i])
        for (int k = 0; k < D; k++) s[k] = fadd(s[k], fmul(v[i], INV[(size_t)i * D + k]));
}

static void make_g(int t) {
    T = t;
    for (;;) {
        for (int i = 0; i < t; i++) G[i] = (int)(rnd() % q);
        G[t] = 1;
        int ok = 1;
        for (int x = 0; x < q && ok; x++) if (peval(G, t + 1, x) == 0) ok = 0;
        if (ok && squarefree(G, t)) break;
    }
    // g^{p-1}, g^p
    int cur[MAXD], tmp[MAXD], len = t + 1;
    memcpy(cur, G, sizeof(int) * len);
    for (int e = 2; e <= p; e++) {
        pmul(cur, len, G, t + 1, tmp);
        len += t;
        memcpy(cur, tmp, sizeof(int) * len);
        if (e == p - 1) { D1 = (p - 1) * t; memcpy(GQ1, cur, sizeof(int) * len); }
    }
    if (p == 2) { D1 = t; memcpy(GQ1, G, sizeof(int) * (t + 1)); }   // g^{p-1} = g
    D2 = p * t; memcpy(GQ2, cur, sizeof(int) * (D2 + 1));
}

/* ---------- F_p 上の零空間 ---------- */
static int pinv(int a) { for (int x = 1; x < p; x++) if (a * x % p == 1) return x; return 0; }

static int nullvec(const int *INV, int D, int *c) {
    int R = m * D, n = q;
    int *H = calloc((size_t)R * n, sizeof(int));
    for (int i = 0; i < n; i++)
        for (int k = 0; k < D; k++) {
            int d[16]; todig(INV[(size_t)i * D + k], d);
            for (int j = 0; j < m; j++) H[(size_t)(k * m + j) * n + i] = d[j];
        }
    int *piv = malloc(sizeof(int) * R), rank = 0;
    for (int col = 0; col < n && rank < R; col++) {
        int r = -1;
        for (int i = rank; i < R; i++) if (H[(size_t)i * n + col]) { r = i; break; }
        if (r < 0) continue;
        if (r != rank) for (int k = 0; k < n; k++) { int tmp = H[(size_t)r * n + k]; H[(size_t)r * n + k] = H[(size_t)rank * n + k]; H[(size_t)rank * n + k] = tmp; }
        int iv = pinv(H[(size_t)rank * n + col]);
        for (int k = 0; k < n; k++) H[(size_t)rank * n + k] = H[(size_t)rank * n + k] * iv % p;
        for (int i = 0; i < R; i++) if (i != rank && H[(size_t)i * n + col]) {
            int f = H[(size_t)i * n + col];
            for (int k = 0; k < n; k++) H[(size_t)i * n + k] = ((H[(size_t)i * n + k] - f * H[(size_t)rank * n + k]) % p + p) % p;
        }
        piv[rank++] = col;
    }
    int *isp = calloc(n, sizeof(int));
    for (int i = 0; i < rank; i++) isp[piv[i]] = 1;
    for (int i = 0; i < n; i++) c[i] = isp[i] ? 0 : (int)(rnd() % p);
    for (int i = 0; i < rank; i++) {
        long s = 0;
        for (int k = 0; k < n; k++) if (!isp[k]) s += (long)H[(size_t)i * n + k] * c[k];
        c[piv[i]] = (int)(((-s) % p + p) % p);
    }
    free(H); free(piv); free(isp);
    return n - rank;   // 次元
}

/* ---------- 復号: 鍵方程式 sigma*S = omega (mod g^p) を拡張ユークリッド ---------- */
static int decode(const int *S, const int *M, int D, int *e_out) {
    int wmax = D / 2;
    int rp[MAXD], rc[MAXD], vp[MAXD], vc[MAXD], qt[MAXD], tmp[MAXD], prod[MAXD];
    memset(rp, 0, sizeof rp); memset(rc, 0, sizeof rc); memset(vp, 0, sizeof vp); memset(vc, 0, sizeof vc);
    memcpy(rp, M, sizeof(int) * (D + 1));
    memcpy(rc, S, sizeof(int) * D);
    vc[0] = 1;
    while (pdeg(rc, D + 1) >= wmax) {
        int dr = pdeg(rc, D + 1);
        int a[MAXD]; memcpy(a, rp, sizeof a);
        pdivmod(a, D + 1, rc, dr, qt);                      // a = rp mod rc
        // vnext = vp - qt*vc
        int dq = pdeg(qt, D + 1), dv = pdeg(vc, D + 1);
        for (int i = 0; i < MAXD; i++) tmp[i] = vp[i];
        if (dq >= 0 && dv >= 0) {
            pmul(qt, dq + 1, vc, dv + 1, prod);
            for (int i = 0; i <= dq + dv; i++) tmp[i] = fsub(tmp[i], prod[i]);
        }
        memcpy(rp, rc, sizeof rp); memcpy(rc, a, sizeof rc);
        memcpy(vp, vc, sizeof vp); memcpy(vc, tmp, sizeof vc);
    }
    int ds = pdeg(vc, MAXD);
    if (ds < 0) return 0;
    int ilc = finv(vc[ds]);
    for (int i = 0; i <= ds; i++) vc[i] = fmul(vc[i], ilc);
    for (int i = 0; i <= D; i++) rc[i] = fmul(rc[i], ilc);
    // 根探索 (L = GF(q) 全体)
    memset(e_out, 0, sizeof(int) * q);
    int nroots = 0;
    int dsig[MAXD]; memset(dsig, 0, sizeof dsig);
    for (int k = 1; k <= ds; k++) dsig[k - 1] = fmul(k % p, vc[k]);
    for (int x = 0; x < q; x++) {
        if (peval(vc, ds + 1, x) != 0) continue;
        nroots++;
        int d = peval(dsig, ds, x);
        if (d == 0) return -1;                              // 重根
        int ev = fmul(peval(rc, D + 1, x), finv(d));
        if (ev >= p) return -1;                             // F_p の元でない
        e_out[x] = ev;
    }
    if (nroots != ds) return -1;                            // L 上で完全分解しない
    return ds;
}


/* ---------- 表示: 体とゴッパ多項式 ---------- */
static void print_poly(const char *name, const int *a, int deg) {
    printf("  %s =", name);
    int first = 1;
    for (int k = deg; k >= 0; k--) {
        if (!a[k]) continue;
        printf(first ? " " : " + ");
        first = 0;
        if (a[k] == 1) { if (k == 0) printf("1"); }
        else printf("a^%d", LOG_[a[k]]);
        if (k >= 2) printf(a[k] == 1 ? "x^%d" : "*x^%d", k);
        else if (k == 1) printf(a[k] == 1 ? "x" : "*x");
    }
    printf("\n");
}

static void print_setup(void) {
    int f[16]; todig(PRIMC, f);
    printf("体 GF(%d^%d) (元の番号 = sum d_j p^j ↔ sum d_j a^j): a の最小多項式 x^%d", p, m, m);
    for (int k = m - 1; k >= 0; k--) if (f[k]) { if (f[k] == 1) printf(" + x^%d", k); else printf(" + %d*x^%d", f[k], k); }
    printf(" = 0\n");
    printf("ゴッパ多項式 (次数 %d, monic, a^k は a のk乗):\n", T);
    print_poly("g(x)    ", G, T);
    printf("  係数の元の番号 (x^k: 番号):");
    for (int k = T; k >= 0; k--) printf(" x^%d:%d", k, G[k]);
    printf("\n");
    print_poly("g(x)^p  ", GQ2, D2);
    printf("  (g^{p-1} は次数 %d, g^p は次数 %d。符号は Γ(L,g^{p-1}) = Γ(L,g^p))\n", D1, D2);
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s p m t w trials\n", argv[0]); return 1; }
    p = atoi(argv[1]); m = atoi(argv[2]);
    int t = atoi(argv[3]), w = atoi(argv[4]), trials = atoi(argv[5]);
    q = 1; for (int i = 0; i < m; i++) q *= p;
    // パラメータ検査 (配列サイズ・符号が自明でないこと)
    if (m < 1 || m > 15 || p < 2) { fprintf(stderr, "p>=2, 1<=m<=15 にしてください\n"); return 1; }
    if (q > 4096) { fprintf(stderr, "q=p^m=%d は大きすぎます (q<=4096。加算表が q^2 になるため)\n", q); return 1; }
    if (p * t + 1 > MAXD) {
        fprintf(stderr, "g^p の次数 p*t=%d が配列サイズ MAXD=%d を超えます。p が大きい場合は wild 方式は現実的ではありません\n", p * t, MAXD);
        return 1;
    }
    if ((long)m * (p - 1) * t >= q) {
        fprintf(stderr, "冗長度 m(p-1)t=%d が符号長 n=q=%d 以上で、符号の次元が0以下です。p を小さく、m を大きく、t を小さくしてください\n", m * (p - 1) * t, q);
        return 1;
    }
    build_field();
    make_g(t);
    INV1 = build_inv(GQ1, D1);
    INV2 = build_inv(GQ2, D2);

    int verbose = argc > 7 ? atoi(argv[7]) : 0;      // 最初の verbose 回分の誤りと復号結果を表示
    const char *spec = argc > 8 ? argv[8] : NULL;    // 誤りの手動指定 "pos:val,pos:val,..."
    int useold = argc > 6 && atoi(argv[6]) == 1;   // 1: 法 g^{p-1} で復号 (従来方式)
    if (verbose > 0) print_setup();
    int ok = 0, detected = 0, miscorr = 0, thmfail = 0, dim = 0;
    for (int tr = 0; tr < trials; tr++) {
        int c[4096], e[4096], r[4096], S[MAXD], eh[4096];
        dim = nullvec(INV1, D1, c);
        // 理論の確認: 符号語の法 g^p シンドロームは 0 のはず
        syndrome(c, INV2, D2, S);
        int z = 1; for (int k = 0; k < D2; k++) if (S[k]) z = 0;
        if (!z) thmfail++;
        memset(e, 0, sizeof(int) * q);
        if (spec && spec[0]) {                       // 誤りを手で指定: "位置:値,位置:値,..."
            char buf[4096]; strncpy(buf, spec, 4095); buf[4095] = 0;
            for (char *tok = strtok(buf, ","); tok; tok = strtok(NULL, ",")) {
                int pos, val;
                if (sscanf(tok, "%d:%d", &pos, &val) == 2 && pos >= 0 && pos < q && val % p != 0)
                    e[pos] = ((val % p) + p) % p;
            }
            w = 0; for (int i = 0; i < q; i++) if (e[i]) w++;
        } else {                                      // ランダム: 重み w, 値は F_p の非零元
            for (int cnt = 0; cnt < w;) {
                int pos = (int)(rnd() % q);
                if (e[pos]) continue;
                e[pos] = 1 + (int)(rnd() % (p - 1));
                cnt++;
            }
        }
        for (int i = 0; i < q; i++) r[i] = (c[i] + e[i]) % p;
        syndrome(r, useold ? INV1 : INV2, useold ? D1 : D2, S);
        int res = decode(S, useold ? GQ1 : GQ2, useold ? D1 : D2, eh);
        if (res < 0) detected++;
        else if (memcmp(e, eh, sizeof(int) * q) == 0) ok++;
        else miscorr++;
        if (tr < verbose) {
            printf("--- trial %d (重み %d) ---\n  入れた誤り (位置:値):", tr + 1, w);
            for (int i = 0; i < q; i++) if (e[i]) printf(" %d:%d", i, e[i]);
            printf("\n  復号結果   (位置:値):");
            if (res < 0) printf(" 復号失敗(訂正不能と判定)");
            else for (int i = 0; i < q; i++) if (eh[i]) printf(" %d:%d", i, eh[i]);
            printf("\n  判定: %s\n", res < 0 ? "失敗を検出" : (memcmp(e, eh, sizeof(int) * q) == 0 ? "一致 (正しく訂正)" : "不一致 (誤訂正)"));
        }
    }
    printf("[復号の法 g^%d] p=%d m=%d q=%d t=%d | dim=%d (>= %d) | wild bound floor(pt/2)=%d, old bound floor((p-1)t/2)=%d | w=%d: ok=%d fail_detected=%d miscorrected=%d theory_check_fail=%d (trials=%d)\n",
           useold ? p - 1 : p, p, m, q, t, dim, q - m * D1, (p * t) / 2, ((p - 1) * t) / 2, w, ok, detected, miscorr, thmfail, trials);
    return 0;
}