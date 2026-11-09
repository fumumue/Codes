// 一般化Patterson復号(格子なし版): 鍵方程式 phi*sigma' = sigma*s (mod g) を
// sigma の係数についての連立一次方程式として解く。素体 F_p (m=1) 上の参照実装。
//   gcc -O2 -o gpat_linear gpat_linear.c
//   ./gpat_linear p t w trials mode [verbose] [errors]
//     mode    : 0 = 誤りの値が全部同じ(値は不明。復号側は全 phi を試す)
//               1 = 値がばらばら(従来の候補生成)
//               2 = 誤りの値が全部 1 (復号側も phi=1 と知っている。2元Pattersonに相当)
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
static int PHI_LO = 1, PHI_HI = 1;   // 試す値 phi の範囲 (main で設定)
static int EQ_ONLY = 0;              // 1: 誤りの値が全部同じ候補(重根なし)だけを残す

static int same(const int *a, const int *b) { return memcmp(a, b, sizeof(int) * P) == 0; }

static int decode(const int *s, int *ncand) {
    int found = 0;
    for (int phi = PHI_LO; phi <= PHI_HI; phi++) {
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
        int tot = 0, bad = 0;
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
            if (EQ_ONLY && mu > 1) bad = 1;     // 値が揃っていない候補は棄却
            e[i] = md((long)phi * mu);
        }
        if (bad) continue;
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

/* ---------- BM 法 (べき乗和シンドローム) : 比較用 ----------
 * 同じ符号 Γ(L,g) を、通常のGoppa符号の復号として解く。
 *  s(x) = sum e_i/(x-L_i) mod g  から  T_j = sum_i (e_i/g(L_i)) L_i^j (j=0..t-1) を復元し、
 *  Berlekamp-Massey で位置多項式 prod(1 - L_i x) を求める。誤りの値は任意、訂正数は floor(t/2) まで。 */
static int powmod_(int a, int e) { long r = 1; for (int i = 0; i < e; i++) r = r * a % P; return (int)r; }

static int bm_decode(const int *s, int *e_out) {
    int N = T, Tt[MAXT + 1];
    for (int a = T - 1; a >= 0; a--) {            // s_a = - sum_b g_{a+1+b} T_b を T について解く
        int idx = T - 1 - a;
        long v = -(long)s[a];
        for (int b = 0; b < idx; b++) v -= (long)G[a + 1 + b] * Tt[b];
        Tt[idx] = md(v);
    }
    int C[MAXT + 2], B[MAXT + 2], tmp[MAXT + 2];
    memset(C, 0, sizeof C); memset(B, 0, sizeof B);
    C[0] = 1; B[0] = 1;
    int L = 0, mm = 1, bb = 1;
    for (int n = 0; n < N; n++) {
        long d = Tt[n];
        for (int i = 1; i <= L; i++) d += (long)C[i] * Tt[n - i];
        d = md(d);
        if (d == 0) { mm++; continue; }
        int coef = md((long)d * inv(bb));
        if (2 * L <= n) {
            memcpy(tmp, C, sizeof C);
            for (int i = 0; i + mm <= N; i++) C[i + mm] = md(C[i + mm] - (long)coef * B[i]);
            L = n + 1 - L; memcpy(B, tmp, sizeof B); bb = (int)d; mm = 1;
        } else {
            for (int i = 0; i + mm <= N; i++) C[i + mm] = md(C[i + mm] - (long)coef * B[i]);
            mm++;
        }
    }
    if (L > N / 2) return -1;                       // 訂正能力 floor(t/2) を超えている
    int dc = pdeg(C, MAXT + 2);
    if (dc < 0 || L - dc > 1 || dc > L) return -1;
    int pos[MAXT], nv = 0;
    for (int i = 1; i < P; i++) {                   // L_i = i != 0 の根: Lambda(1/L_i) = 0
        int xi = inv(i); long r = 0;
        for (int k = dc; k >= 0; k--) r = (r * xi + C[k]) % P;
        if (r == 0) pos[nv++] = i;
    }
    if (nv != dc) return -1;
    if (L - dc == 1) pos[nv++] = 0;                 // 位置 L=0 の誤り
    if (nv != L) return -1;
    // 誤りの値: Vandermonde 系  sum_k Y_k L_k^j = T_j (j<L)
    int M[MAXT][MAXT + 1];
    for (int j = 0; j < L; j++) {
        for (int k = 0; k < L; k++) M[j][k] = pos[k] == 0 ? (j == 0) : powmod_(pos[k], j);
        M[j][L] = Tt[j];
    }
    for (int c = 0; c < L; c++) {
        int r = -1;
        for (int i = c; i < L; i++) if (M[i][c]) { r = i; break; }
        if (r < 0) return -1;
        if (r != c) for (int k = 0; k <= L; k++) { int x = M[r][k]; M[r][k] = M[c][k]; M[c][k] = x; }
        int iv = inv(M[c][c]);
        for (int k = 0; k <= L; k++) M[c][k] = md((long)M[c][k] * iv);
        for (int i = 0; i < L; i++) if (i != c && M[i][c]) {
            int f = M[i][c];
            for (int k = 0; k <= L; k++) M[i][k] = md(M[i][k] - (long)f * M[c][k]);
        }
    }
    memset(e_out, 0, sizeof(int) * P);
    for (int k = 0; k < L; k++) e_out[pos[k]] = md((long)M[k][L] * geval(pos[k]));
    int s2[MAXT]; syndrome(e_out, s2);
    if (memcmp(s, s2, sizeof(int) * T) != 0) return -1;
    return L;
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
    EQ_ONLY = (mode != 1);
    PHI_LO = (mode == 2) ? 1 : 1; PHI_HI = (mode == 2) ? 1 : P - 1;
    make_code(t);
    if (verbose > 0) print_g();
    int ok = 0, only_wrong = 0, none = 0, lightest = 0, uok = 0, utie = 0, uwrong = 0, bm_ok = 0, bm_fail = 0, bm_wrong = 0; long lsum = 0;
    for (int tr = 0; tr < trials; tr++) {
        int e[MAXP]; memset(e, 0, sizeof(int) * P);
        int phi0 = mode == 2 ? 1 : 1 + (int)(rnd() % (P - 1));
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
                e[pos] = mode != 1 ? phi0 : 1 + (int)(rnd() % (P - 1));
                c++;
            }
        }
        int s[MAXT]; syndrome(e, s);
        int n; decode(s, &n);
        int ebm[MAXP]; int rb = bm_decode(s, ebm);
        if (rb < 0) bm_fail++; else if (memcmp(ebm, e, sizeof(int) * P) == 0) bm_ok++; else bm_wrong++;
        lsum += n;
        int hit = 0;
        for (int j = 0; j < n; j++) if (same(cand[j], e)) hit = 1;
        if (hit) ok++; else if (n) only_wrong++; else none++;
        if (n > 0) {                                 // 一意復号: 重み最小の候補を返す
            int best = 1 << 30, cnt = 0, bj = -1;
            for (int j = 0; j < n; j++) {
                int wj = 0; for (int i = 0; i < P; i++) if (cand[j][i]) wj++;
                if (wj < best) { best = wj; cnt = 1; bj = j; } else if (wj == best) cnt++;
            }
            if (cnt > 1) utie++; else if (same(cand[bj], e)) uok++; else uwrong++;
        }
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
            printf("  BM法の結果:");
            if (rb < 0) printf(" 復号失敗");
            else { for (int i = 0; i < P; i++) if (ebm[i]) printf(" %d:%d", i, ebm[i]); printf("%s", memcmp(ebm, e, sizeof(int) * P) == 0 ? "   <- 一致" : "   <- 不一致(誤訂正)"); }
            printf("\n");
        }
    }
    printf("p=%d t=%d w=%d mode=%d trials=%d : truth_found=%d  only_wrong=%d  no_candidate=%d  avg_list=%.2f  truth_is_unique_lightest=%d  | 一意復号: ok=%d tie=%d wrong=%d | BM法: ok=%d fail=%d wrong=%d\n",
           P, t, w, mode, trials, ok, only_wrong, none, (double)lsum / trials, lightest, uok, utie, uwrong, bm_ok, bm_fail, bm_wrong);
    return 0;
}