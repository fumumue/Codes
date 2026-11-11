// Alternant 符号の構成と復号 (標数 P 一般, 係数体 F_q, q = P^m, 符号は F_P 上)
// 自己完結版: gcc -O2 -o alternant alternant.c
//
// ---- 定義 ----
// 相異なる α_1..α_n ∈ F_q^* と、乗数 y_1..y_n ∈ F_q^* (すべて非零) を与える。
// 検査行列 H = [ y_i α_i^j ]  (j = 0..R-1, i = 1..n)   (R x n, F_q 上)
// Alternant 符号 = { c ∈ F_P^n | H c = 0 }  (F_q の元を F_P 上の m 桁に展開すれば F_P 上の検査行列 mR x n)
// 最小距離 >= R+1、次元 >= n - mR。⌊R/2⌋ 個までの誤りを訂正できる。
//
// ---- 乗数 y_i の選び方 (mode) ----
//   0: ランダムな非零元                      -> 一般の alternant 符号
//   1: y_i = 1/g(α_i)^e  (g はランダムな t 次多項式, t = R/e)
//        -> Goppa 符号 Γ(L, g^e)  (行数 R = e t)。e=1 なら古典 Goppa、e=P なら wild Goppa Γ(L,g^P) (g 無平方なら Γ(L,g^(P-1)) と一致)
//        乗数が 1/G(α_i) (G = 任意の R 次多項式) の形なら、G が和 g1+g2 でも積 g1*g2 でもこの mode と同じ
//   2: y_i = 1/g1(α_i) + 1/g2(α_i)  (g1, g2 はランダムな R 次多項式)  -> 一般には Goppa にならない alternant 符号
//   3: y_i = 1                                -> GRS の部分体部分符号 (BCH 型)
//
// ---- 復号 ----
// 受信語 r = c + e のシンドローム s_j = Σ_i r_i y_i α_i^j = Σ_{k in 誤り位置} Y_k X_k^j  (X_k = α_k, Y_k = e_k y_k)
// 1. Berlekamp-Massey 法で位置多項式 Λ(x) = Π (1 - X_k x) を求める
// 2. Λ(α_i^{-1}) = 0 となる i を探す (Chien 探索)
// 3. Ω(x) = S(x) Λ(x) mod x^R、Forney の式 Y_k = - X_k Ω(X_k^{-1}) / Λ'(X_k^{-1}) で値を求め、e_k = Y_k / y_k
//    (e_k が F_P の元 (非零) でなければ失敗。最後にシンドロームを再計算して照合)
//
// 使い方: ./alternant P m R w trials mode [e] [seed] [demo]
//   P      標数 (素数)
//   m      拡大次数 (q = P^m <= 16384)
//   R      検査行の数 (設計距離 R+1)。⌊R/2⌋ 個まで確実に訂正
//   w      入れる誤りの個数
//   trials 試行回数 (試行ごとに乗数/多項式を作り直す)
//   mode   0..3 (上記)
//   e      (mode 1 のみ) g の指数。既定 1。R は e の倍数
//   seed   乱数の種。既定 12345
//   demo   1 にすると、符号を F_P 上の検査行列から明示的に構成し (RREF)、ランダムな符号語を作って
//          誤りを加え、復号して送信語を復元する様子を 1 回表示する (trials は無視)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXQ 16384
#define MAXR 1024

static int P, Mm, Q, Q1;
static int pw[16];
static int EXPT[2 * MAXQ], LOGT[MAXQ], ZECH[MAXQ], NEG[MAXQ];

// ---------- F_q (q = P^m), 元は P 進 m 桁の整数 (多項式基底)。F_P の元は 0..P-1 ----------
static int dadd(int a, int b)
{
    int r = 0;
    for (int i = 0; i < Mm; i++)
        r += (((a / pw[i]) % P + (b / pw[i]) % P) % P) * pw[i];
    return r;
}
static int dneg(int a)
{
    int r = 0;
    for (int i = 0; i < Mm; i++)
        r += ((P - (a / pw[i]) % P) % P) * pw[i];
    return r;
}
static int dscale(int v, int c)
{
    int r = 0;
    for (int i = 0; i < Mm; i++)
        r += ((((v / pw[i]) % P) * c) % P) * pw[i];
    return r;
}
static inline int mul(int a, int b) { return (a && b) ? EXPT[LOGT[a] + LOGT[b]] : 0; }
static inline int inv(int a) { return EXPT[Q1 - LOGT[a]]; }
static inline int fadd(int a, int b) // Zech 対数表による加算 (表の大きさは q)
{
    if (!a)
        return b;
    if (!b)
        return a;
    int d = LOGT[b] - LOGT[a];
    if (d < 0)
        d += Q1;
    int z = ZECH[d];
    return z < 0 ? 0 : EXPT[LOGT[a] + z];
}
static inline int fsub(int a, int b) { return fadd(a, NEG[b]); }

static int init_field(int p, int m)
{
    P = p;
    Mm = m;
    pw[0] = 1;
    for (int i = 1; i <= m; i++)
        pw[i] = pw[i - 1] * P;
    Q = pw[m];
    if (Q > MAXQ)
        return 0;
    Q1 = Q - 1;
    int ok = 0;
    for (int pc = 1; pc < Q && !ok; pc++)
    {
        int negpc = dneg(pc), a = 1, ord = 0;
        do
        {
            int top = a / pw[m - 1], sh = (a % pw[m - 1]) * P;
            a = dadd(sh, dscale(negpc, top));
            ord++;
        } while (a != 1 && ord < Q);
        if (a == 1 && ord == Q1)
        {
            a = 1;
            for (int k = 0; k < Q1; k++)
            {
                EXPT[k] = EXPT[k + Q1] = a;
                LOGT[a] = k;
                int top = a / pw[m - 1], sh = (a % pw[m - 1]) * P;
                a = dadd(sh, dscale(negpc, top));
            }
            ok = 1;
        }
    }
    if (!ok)
        return 0;
    for (int d = 0; d < Q1; d++)
    {
        int sm = dadd(1, EXPT[d]);
        ZECH[d] = sm ? LOGT[sm] : -1;
    }
    for (int a = 0; a < Q; a++)
        NEG[a] = dneg(a);
    return 1;
}

// F_q 係数多項式 (c[0..d]) の Horner 評価
static int peval(const int *c, int d, int x)
{
    int r = 0;
    for (int i = d; i >= 0; i--)
        r = fadd(mul(r, x), c[i]);
    return r;
}

// ---------- 符号の構成 ----------
static int n;
static int Al[MAXQ], Yv[MAXQ]; // α_i, y_i

static int rnd_nonzero(void) { return 1 + rand() % Q1; }

static void build_code(int mode, int R, int ee)
{
    int g[MAXR + 2], g2[MAXR + 2], t = 0;
    if (mode == 1)
    {
        t = R / ee;
        for (int i = 0; i < t; i++)
            g[i] = rand() % Q;
        g[t] = 1;
    }
    if (mode == 2)
        for (int i = 0; i < R; i++)
        {
            g[i] = rand() % Q;
            g2[i] = rand() % Q;
            g[R] = g2[R] = 1;
        }
    n = 0;
    for (int k = 0; k < Q1; k++)
    {
        int a = EXPT[k], y;
        if (mode == 0)
            y = rnd_nonzero();
        else if (mode == 1)
        {
            int gv = peval(g, t, a);
            if (!gv)
                continue; // g(α)=0 の点は台から除く
            int pe = 1;
            for (int i = 0; i < ee; i++)
                pe = mul(pe, gv);
            y = inv(pe);
        }
        else if (mode == 2)
        {
            int v1 = peval(g, R, a), v2 = peval(g2, R, a);
            if (!v1 || !v2)
                continue;
            y = fadd(inv(v1), inv(v2));
            if (!y)
                continue;
        }
        else
            y = 1;
        Al[n] = a;
        Yv[n] = y;
        n++;
    }
}

// 受信語 (または誤りベクトル) v ∈ F_P^n のシンドローム s_j = Σ v_i y_i α_i^j (j < R)
static void syndromes(const int *v, int R, int *s)
{
    for (int j = 0; j < R; j++)
        s[j] = 0;
    for (int i = 0; i < n; i++)
        if (v[i])
        {
            int Yk = mul(v[i], Yv[i]), pwr = 1;
            for (int j = 0; j < R; j++)
            {
                s[j] = fadd(s[j], mul(Yk, pwr));
                pwr = mul(pwr, Al[i]);
            }
        }
}

// ---------- 復号 ----------
// Berlekamp-Massey: s[0..N-1] から C(x) (C[0]=1) と長さ L を返す。s_k + Σ_{i=1}^L C_i s_{k-i} = 0
static int berlekamp_massey(const int *s, int N, int *C)
{
    static int B[MAXR + 4], T[MAXR + 4];
    for (int i = 0; i <= N + 1; i++)
        C[i] = B[i] = 0;
    C[0] = B[0] = 1;
    int L = 0, mm = 1, b = 1;
    for (int k = 0; k < N; k++)
    {
        int d = s[k];
        for (int i = 1; i <= L; i++)
            d = fadd(d, mul(C[i], s[k - i]));
        if (d == 0)
        {
            mm++;
            continue;
        }
        int coef = mul(d, inv(b));
        if (2 * L <= k)
        {
            memcpy(T, C, sizeof(int) * (N + 2));
            for (int i = 0; i + mm <= N + 1; i++)
                C[i + mm] = fsub(C[i + mm], mul(coef, B[i]));
            L = k + 1 - L;
            memcpy(B, T, sizeof(int) * (N + 2));
            b = d;
            mm = 1;
        }
        else
        {
            for (int i = 0; i + mm <= N + 1; i++)
                C[i + mm] = fsub(C[i + mm], mul(coef, B[i]));
            mm++;
        }
    }
    return L;
}

// s[0..R-1] から誤りベクトル e (長さ n, F_P の元) を求める。成功 1 / 失敗 0
static int decode_alt(const int *s, int R, int *e)
{
    int Lam[MAXR + 4], Om[MAXR + 4], dL[MAXR + 4], pos[MAXR + 4], chk[MAXR + 4];
    memset(e, 0, sizeof(int) * n);
    int allz = 1;
    for (int j = 0; j < R; j++)
        if (s[j])
            allz = 0;
    if (allz)
        return 1;
    int L = berlekamp_massey(s, R, Lam);
    if (L == 0 || L > R / 2 || Lam[L] == 0)
        return 0;
    int cnt = 0;
    for (int i = 0; i < n; i++)
        if (peval(Lam, L, inv(Al[i])) == 0)
        {
            if (cnt >= L)
                return 0;
            pos[cnt++] = i;
        }
    if (cnt != L)
        return 0;
    for (int k = 0; k < L; k++) // Ω = S Λ mod x^R (次数 <= L-1)
    {
        int v = 0;
        for (int i = 0; i <= k; i++)
            v = fadd(v, mul(Lam[i], s[k - i]));
        Om[k] = v;
    }
    for (int i = 1; i <= L; i++) // Λ'(x): 係数 i*Λ_i (i は F_P の元として)
        dL[i - 1] = mul(i % P, Lam[i]);
    for (int k = 0; k < L; k++)
    {
        int i = pos[k], X = Al[i], xi = inv(X);
        int dv = peval(dL, L - 1, xi);
        if (dv == 0)
            return 0;
        int Y = NEG[mul(X, mul(peval(Om, L - 1, xi), inv(dv)))];
        int ev = mul(Y, inv(Yv[i]));
        if (ev == 0 || ev >= P)
            return 0; // F_P の元でない
        e[i] = ev;
    }
    syndromes(e, R, chk); // 最終照合
    for (int j = 0; j < R; j++)
        if (chk[j] != s[j])
            return 0;
    return 1;
}

// ---------- テスト ----------
static int is_prime(int p)
{
    if (p < 2)
        return 0;
    for (int d = 2; d * d <= p; d++)
        if (p % d == 0)
            return 0;
    return 1;
}
static void make_error(int *e, int w)
{
    memset(e, 0, sizeof(int) * n);
    int cnt = 0;
    while (cnt < w)
    {
        int i = rand() % n;
        if (e[i])
            continue;
        e[i] = 1 + rand() % (P - 1);
        cnt++;
    }
}
static int modinv_p(int a)
{
    long r = 1, b = a % P;
    for (int ex = P - 2; ex; ex >>= 1)
    {
        if (ex & 1)
            r = r * b % P;
        b = b * b % P;
    }
    return (int)r;
}

// F_P 上の検査行列 (mR x n) から符号を明示的に構成し、符号語を作って誤りを加え、復号して復元する
static void demo(int R, int w)
{
    size_t rows = (size_t)Mm * R;
    unsigned short *H = calloc(rows * n, sizeof(unsigned short));
    int *piv = malloc(rows * sizeof(int));
    char *isp = calloc(n, 1);
    int *c = calloc(n, sizeof(int)), *e = calloc(n, sizeof(int)), *r = calloc(n, sizeof(int)), *er = calloc(n, sizeof(int));
    for (int i = 0; i < n; i++)
    {
        int v = Yv[i]; // y_i α_i^j
        for (int j = 0; j < R; j++)
        {
            for (int d = 0; d < Mm; d++)
                H[((size_t)j * Mm + d) * n + i] = (unsigned short)((v / pw[d]) % P);
            v = mul(v, Al[i]);
        }
    }
    int rank = 0;
    for (int col = 0; col < n && rank < (int)rows; col++)
    {
        int p = -1;
        for (int j = rank; j < (int)rows; j++)
            if (H[(size_t)j * n + col])
            {
                p = j;
                break;
            }
        if (p < 0)
            continue;
        if (p != rank)
            for (int cc = 0; cc < n; cc++)
            {
                unsigned short tmp = H[(size_t)p * n + cc];
                H[(size_t)p * n + cc] = H[(size_t)rank * n + cc];
                H[(size_t)rank * n + cc] = tmp;
            }
        int iv = modinv_p(H[(size_t)rank * n + col]);
        for (int cc = col; cc < n; cc++)
            H[(size_t)rank * n + cc] = (unsigned short)((long)H[(size_t)rank * n + cc] * iv % P);
        for (int j = 0; j < (int)rows; j++)
            if (j != rank && H[(size_t)j * n + col])
            {
                int f = H[(size_t)j * n + col];
                for (int cc = col; cc < n; cc++)
                    H[(size_t)j * n + cc] = (unsigned short)((H[(size_t)j * n + cc] + (long)(P - f) * H[(size_t)rank * n + cc]) % P);
            }
        piv[rank] = col;
        isp[col] = 1;
        rank++;
    }
    int k = n - rank;
    printf("[構成] 符号長 n=%d, F_P 上の検査行列 %d x %d の階数=%d -> 次元 k=%d (下限 n-mR=%d), 最小距離 >= %d, 訂正可能 %d 個\n",
           n, (int)rows, n, rank, k, n - Mm * R, R + 1, R / 2);
    if (k <= 0)
    {
        printf("次元 0 の符号です (R が大きすぎます)\n");
        goto done;
    }
    // 符号語: 自由変数 (非ピボット列) をランダムに選び、ピボット変数を決める
    for (int i = 0; i < n; i++)
        if (!isp[i])
            c[i] = rand() % P;
    for (int rr = 0; rr < rank; rr++)
    {
        long sum = 0;
        for (int i = 0; i < n; i++)
            if (!isp[i])
                sum += (long)H[(size_t)rr * n + i] * c[i];
        c[piv[rr]] = (int)((P - sum % P) % P);
    }
    int *sc = malloc(sizeof(int) * R), *sr = malloc(sizeof(int) * R);
    syndromes(c, R, sc);
    int zero = 1;
    for (int j = 0; j < R; j++)
        if (sc[j])
            zero = 0;
    printf("[符号化] ランダムな符号語 c を作成。H c = 0 の確認: %s\n", zero ? "OK" : "NG");
    make_error(e, w);
    for (int i = 0; i < n; i++)
        r[i] = (c[i] + e[i]) % P;
    syndromes(r, R, sr); // 受信語から計算 (誤りだけのシンドロームと一致するはず)
    int ok = decode_alt(sr, R, er);
    if (!ok)
        printf("[復号] 失敗 (w=%d, 訂正可能 %d 個)\n", w, R / 2);
    else
    {
        int same_e = (memcmp(e, er, sizeof(int) * n) == 0), same_c = 1;
        for (int i = 0; i < n; i++)
            if (((r[i] - er[i]) % P + P) % P != c[i])
                same_c = 0;
        printf("[復号] w=%d 個の誤りを加えて復号 -> 誤りの復元: %s, 送信語の復元: %s\n", w, same_e ? "一致" : "不一致", same_c ? "一致" : "不一致");
    }
    free(sc);
    free(sr);
done:
    free(H);
    free(piv);
    free(isp);
    free(c);
    free(e);
    free(r);
    free(er);
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "使い方: %s P m R w trials mode [e] [seed] [demo]\n"
            "  P      標数 (素数)\n"
            "  m      拡大次数 (q=P^m <= %d)\n"
            "  R      検査行の数 (設計距離 R+1, 訂正可能 R/2 個)\n"
            "  w      入れる誤りの個数\n"
            "  trials 試行回数\n"
            "  mode   0: ランダム乗数 / 1: Goppa Γ(L,g^e) (t=R/e) / 2: 乗数 1/g1+1/g2 / 3: 乗数 1 (BCH 型)\n"
            "  e      (mode 1) g の指数。既定 1。wild Goppa なら e=P\n"
            "  seed   乱数の種。既定 12345\n"
            "  demo   1: 符号を明示的に構成して符号語を作り、誤りを加えて復号する様子を 1 回表示\n"
            "例: %s 3 4 24 12 200 1 3      (F_81, wild Goppa Γ(L,g^3), deg g=8, 誤り 12 個)\n"
            "    %s 3 4 24 10 1 1 3 1 1    (同じ符号を明示的に構成して復号する例)\n",
            prog, MAXQ, prog, prog);
}

int main(int argc, char **argv)
{
    if (argc < 7)
    {
        usage(argv[0]);
        return 1;
    }
    int p = atoi(argv[1]), m = atoi(argv[2]), R = atoi(argv[3]), w = atoi(argv[4]);
    int trials = atoi(argv[5]), mode = atoi(argv[6]);
    int ee = argc > 7 ? atoi(argv[7]) : 1;
    unsigned seed = argc > 8 ? (unsigned)atoi(argv[8]) : 12345u;
    int dm = argc > 9 ? atoi(argv[9]) : 0;
    if (!is_prime(p) || m < 1 || R < 1 || R > MAXR || w < 0 || mode < 0 || mode > 3 || ee < 1 || (mode == 1 && R % ee != 0))
    {
        usage(argv[0]);
        return 1;
    }
    if (!init_field(p, m))
    {
        fprintf(stderr, "q=P^m が大きすぎます (上限 %d)\n", MAXQ);
        return 1;
    }
    srand(seed);

    int *e = malloc(sizeof(int) * MAXQ), *er = malloc(sizeof(int) * MAXQ), s[MAXR + 2];
    if (dm)
    {
        build_code(mode, R, ee);
        if (w > n)
        {
            fprintf(stderr, "w が符号長 n=%d を超えています\n", n);
            return 1;
        }
        printf("P=%d m=%d (q=%d) R=%d mode=%d e=%d\n", P, Mm, Q, R, mode, ee);
        demo(R, w);
        return 0;
    }
    int okc = 0, wrong = 0, fail = 0;
    double t_dec = 0;
    for (int tr = 0; tr < trials; tr++)
    {
        build_code(mode, R, ee);
        if (w > n)
        {
            fprintf(stderr, "w が符号長 n=%d を超えています\n", n);
            return 1;
        }
        make_error(e, w);
        syndromes(e, R, s);
        clock_t c0 = clock();
        int ok = decode_alt(s, R, er);
        t_dec += (double)(clock() - c0) / CLOCKS_PER_SEC;
        if (!ok)
            fail++;
        else if (memcmp(e, er, sizeof(int) * n) == 0)
            okc++;
        else
            wrong++;
    }
    printf("P=%d m=%d (q=%d, 符号長 n=%d) R=%d mode=%d e=%d w=%d (訂正可能 %d): 成功=%d 誤訂正=%d 復号失敗(検出)=%d (試行 %d)\n",
           P, Mm, Q, n, R, mode, ee, w, R / 2, okc, wrong, fail, trials);
    printf("  復号の平均時間: %.3f ms/試行\n", 1000.0 * t_dec / trials);
    return 0;
}