// 一般化Patterson復号(格子なし版): 鍵方程式 phi*sigma' = sigma*s (mod g) を
// sigma の係数についての連立一次方程式として解く。素体 F_p (m=1) 上の参照実装。
//   gcc -O2 -o gpat_linear gpat_linear.c
//   ./gpat_linear p t w trials mode [verbose] [errors]
//     mode    : 0 = 誤りの値が全部同じ, 1 = 値がばらばら
//     verbose : 1以上で、ゴッパ多項式と、最初の verbose 回分の誤りベクトル・復号候補を表示
//     errors  : 誤りを手で指定 例 "3:5,100:5,999:5" (位置 i は L_i = i、値は F_p の非零元)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXT 64
#define MAXP 4096
#define MAXLIST 16

static int P;
static inline int md(long a) { a %= P; if (a < 0) a += P; return (int)a; }
static int inv(int a) {
    long r = 1, b = a, e = P - 2;
    while (e) { if (e & 1) r = r * b % P; b = b * b % P; e >>= 1; }
    return (int)r;
}

static uint64_t rs = 88172645463325252ULL;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

/* ---------- 符号 ---------- */
static int T;                 // deg g
static int G[MAXT + 1];       // monic g
static int INVL[MAXP][MAXT];  // (x - L_i)^{-1} mod g,  L_i = i, n = P

static int pdeg(const int *a, int n) { while (n > 0 && a[n - 1] == 0) n--; return n - 1; }

static void pmod_inplace(int *a, int na, const int *b, int nb) {
    int db = pdeg(b, nb), ib = inv(b[db]);
    for (int da = pdeg(a, na); da >= db; da = pdeg(a, na)) {
        int c = md((long)a[da] * ib);
        for (int k = 0; k <= db; k++) a[da - db + k] = md(a[da - db + k] - (long)c * b[k]);
    }
}

static int squarefree(const int *g, int t) {
    int a[MAXT + 2], b[MAXT + 2];
    memset(a, 0, sizeof a); memset(b, 0, sizeof b);
    for (int i = 0; i <= t; i++) a[i] = g[i];
    for (int i = 1; i <= t; i++) b[i - 1] = md((long)i * g[i]);
    while (pdeg(b, t + 2) >= 0) {
        pmod_inplace(a, t + 2, b, t + 2);
        int tmp[MAXT + 2];
        memcpy(tmp, a, sizeof tmp); memcpy(a, b, sizeof a); memcpy(b, tmp, sizeof b);
    }
    return pdeg(a, t + 2) == 0;
}

static int geval(int x) {
    long r = 0;
    for (int k = T; k >= 0; k--) r = (r * x + G[k]) % P;
    return (int)r;
}

static void make_code(int t) {
    T = t;
    for (;;) {
        for (int i = 0; i < t; i++) G[i] = (int)(rnd() % P);
        G[t] = 1;
        int ok = 1;
        for (int x = 0; x < P && ok; x++) if (geval(x) == 0) ok = 0;   // Lに根なし
        if (ok && squarefree(G, t)) break;
    }
    for (int i = 0; i < P; i++) {
        int b[MAXT];
        b[T - 1] = G[T];
        for (int k = T - 1; k >= 1; k--) b[k - 1] = md(G[k] + (long)i * b[k]);
        int rem = md(G[0] + (long)i * b[0]);
        int c = md(-(long)inv(rem));
        for (int k = 0; k < T; k++) INVL[i][k] = md((long)c * b[k]);
    }
}

static void syndrome(const int *e, int *s) {
    memset(s, 0, sizeof(int) * T);
    for (int i = 0; i < P; i++) if (e[i])
        for (int k = 0; k < T; k++) s[k] = md(s[k] + (long)e[i] * INVL[i][k]);
}

/* ---------- 復号: 格子なし ---------- */
static int cand[MAXLIST][MAXP];

static int same(const int *a, const int *b) { return memcmp(a, b, sizeof(int) * P) == 0; }

static int decode(const int *s, int *ncand) {
    int found = 0;
    for(int phi = 1; phi <=1; phi++) {
        // 列 col_d = (phi*d*x^{d-1} - x^d * s) mod g,  d = 0..T
        int X[MAXT];
        memcpy(X, s, sizeof(int) * T);
        int bv[MAXT + 1][MAXT], bc[MAXT + 1][MAXT + 1], bp[MAXT + 1], nb = 0;
        int sigma[MAXT + 1], sd = -1;
        for (int d = 0; d <= T; d++) {
            int v[MAXT], c[MAXT + 1];
            for (int k = 0; k < T; k++) v[k] = md(-X[k]);
            if (d >= 1) v[d - 1] = md(v[d - 1] + (long)phi * d);
            memset(c, 0, sizeof c); c[d] = 1;
            for (int j = 0; j < nb; j++) {
                int f = v[bp[j]];
                if (!f) continue;
                for (int k = 0; k < T; k++) v[k] = md(v[k] - (long)f * bv[j][k]);
                for (int k = 0; k <= T; k++) c[k] = md(c[k] - (long)f * bc[j][k]);
            }
            int pv = -1;
            for (int k = 0; k < T; k++) if (v[k]) { pv = k; break; }
            if (pv < 0) { memcpy(sigma, c, sizeof c); sd = d; break; }   // 核ベクトル(monic, 最小次数)
            int iv = inv(v[pv]);
            for (int k = 0; k < T; k++) bv[nb][k] = md((long)v[k] * iv);
            for (int k = 0; k <= T; k++) bc[nb][k] = md((long)c[k] * iv);
            bp[nb++] = pv;
            // X <- x*X mod g
            int top = X[T - 1];
            for (int k = T - 1; k >= 1; k--) X[k] = md(X[k - 1] - (long)top * G[k]);
            X[0] = md(-(long)top * G[0]);
        }
        if (sd < 0) continue;
        // 根と重複度 -> e_i = phi * mu_i
        int e[MAXP]; memset(e, 0, sizeof(int) * P);
        int tot = 0;
        for (int i = 0; i < P; i++) {
            int q[MAXT + 1], cur[MAXT + 1], dg = sd;
            memcpy(cur, sigma, sizeof cur);
            int mu = 0;
            while (dg > 0) {
                q[dg - 1] = cur[dg];
                for (int k = dg - 1; k >= 1; k--) q[k - 1] = md(cur[k] + (long)i * q[k]);
                int rem = md(cur[0] + (long)i * q[0]);
                if (rem) break;
                mu++; dg--;
                memcpy(cur, q, sizeof(int) * (dg + 1));
            }
            tot += mu;
            e[i] = md((long)phi * mu);
        }
        if (tot != sd) continue;                // L上で完全に分解しない
        int s2[MAXT]; syndrome(e, s2);
        if (memcmp(s, s2, sizeof(int) * T) != 0) continue;
        int dup = 0;
        for (int j = 0; j < found; j++) if (same(cand[j], e)) dup = 1;
        if (!dup && found < MAXLIST) memcpy(cand[found++], e, sizeof(int) * P);
    }
    *ncand = found;
    return found;
}

static void print_g(void) {
    printf("ゴッパ多項式 g(x) (次数 %d, monic, F_%d 上):\n  g(x) =", T, P);
    int first = 1;
    for (int k = T; k >= 0; k--) {
        if (!G[k]) continue;
        printf(first ? " " : " + "); first = 0;
        if (k == 0) printf("%d", G[k]);
        else { if (G[k] != 1) printf("%d*", G[k]); if (k == 1) printf("x"); else printf("x^%d", k); }
    }
    printf("\n  係数 (x^0 から x^%d の順):", T);
    for (int k = 0; k <= T; k++) printf(" %d", G[k]);
    printf("\n  (サポート L = F_%d 全体, 位置 i は L_i = i, n = %d)\n", P, P);
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s p t w trials mode\n", argv[0]); return 1; }
    P = atoi(argv[1]); int t = atoi(argv[2]), w = atoi(argv[3]), trials = atoi(argv[4]), mode = atoi(argv[5]);
    if (P > MAXP || t > MAXT || t >= P || w > P) { fprintf(stderr, "bad params\n"); return 1; }
    int verbose = argc > 6 ? atoi(argv[6]) : 0;
    const char *spec = argc > 7 ? argv[7] : NULL;
    make_code(t);
    if (verbose > 0) print_g();
    int ok = 0, only_wrong = 0, none = 0, lightest = 0; long lsum = 0;
    for (int tr = 0; tr < trials; tr++) {
        int e[MAXP]; memset(e, 0, sizeof(int) * P);
        int phi0 = 1; //+ (int)(rnd() % (P - 1));
        if (spec && spec[0]) {                       // 誤りを手で指定
            char buf[4096]; strncpy(buf, spec, 4095); buf[4095] = 0;
            for (char *tok = strtok(buf, ","); tok; tok = strtok(NULL, ",")) {
                int pos, val;
                if (sscanf(tok, "%d:%d", &pos, &val) == 2 && pos >= 0 && pos < P && md(val) != 0)
                    e[pos] = md(val);
            }
            w = 0; for (int i = 0; i < P; i++) if (e[i]) w++;
        } else {
            for (int c = 0; c < w;) {
                int pos = (int)(rnd() % P);
                if (e[pos]) continue;
                e[pos] = mode == 0 ? phi0 : 1 + (int)(rnd() % (P - 1));
                c++;
            }
        }
        int s[MAXT]; syndrome(e, s);
        int n; decode(s, &n);
        lsum += n;
        int hit = 0;
        for (int j = 0; j < n; j++) if (same(cand[j], e)) hit = 1;
        if (hit) ok++; else if (n) only_wrong++; else none++;
        if (hit) {                                   // 正解が候補の中で唯一の最軽量(重み最小)か
            int wt = 0; for (int i = 0; i < P; i++) if (e[i]) wt++;
            int uniq = 1;
            for (int j = 0; j < n; j++) if (!same(cand[j], e)) {
                int wj = 0; for (int i = 0; i < P; i++) if (cand[j][i]) wj++;
                if (wj <= wt) uniq = 0;
            }
            lightest += uniq;
        }
        if (tr < verbose) {
            printf("--- trial %d (重み %d) ---\n  誤りベクトル (位置:値):", tr + 1, w);
            for (int i = 0; i < P; i++) if (e[i]) printf(" %d:%d", i, e[i]);
            printf("\n  復号候補 %d 個\n", n);
            for (int j = 0; j < n; j++) {
                printf("   [%d]", j + 1);
                for (int i = 0; i < P; i++) if (cand[j][i]) printf(" %d:%d", i, cand[j][i]);
                printf("%s\n", same(cand[j], e) ? "   <- 入れた誤りと一致" : "");
            }
            if (n == 0) printf("   (候補なし = 復号失敗)\n");
        }
    }
    printf("p=%d t=%d w=%d mode=%d trials=%d : truth_found=%d  only_wrong=%d  no_candidate=%d  avg_list=%.2f  truth_is_unique_lightest=%d\n",
           P, t, w, mode, trials, ok, only_wrong, none, (double)lsum / trials, lightest);
    return 0;
}
