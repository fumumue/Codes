// patterson_q.c : 一般化 Patterson 復号 (Barreto-Misoczki-Lindner, arXiv:1103.3296)
//
// 対象: L = GF(Pr^E) の全元(位置 i = 体の元 i)、g = 次数 K の無平方(既約でよい)多項式、
//       誤りの値は部分体 F_Pr の元 (整数 1..Pr-1)。
//
// 使い方: OVP のソースで vv(), mat, gf, fg, mlt, mkmf(), makefg() が定義された後、
//         main() の前に  #include "patterson_q.c"  する。
//
// 鍵方程式: φ σ' = σ s  (mod g),  s(x) = Σ e_i/(x-L_i) mod g,  σ = Π (x-L_i)^μ_i,  e_i = φ μ_i
//   σ = Σ_j x^j a_j(x)^p と書くと、 Σ_j a_j(x) v_j(x) = 0 (mod g),
//   v_j = (x^j - φ j x^(j-1)/s)^(1/p)  (mod g)。
//   a_j の係数を未知数とする連立一次方程式 (K 本、未知数 K+1 個) を解き、
//   次数が最小の解を σ とする。
//
// 注意: 体の演算は gf/fg/mlt だけに依存し、足し算は Pr 進の桁ごとに行う(plus と同じ結果)。
// 名前の衝突を避けるため、関数名・変数名は q_ で始めている。

#define Q_MAXC 4

typedef struct
{
    short c[K];
} q_res; // g を法とする剰余 (係数 K 個)

static short q_g[K + 1];      // モニックにした g
static short q_g_lead;        // BM 検査行列が使う元の g の最高次係数
static short q_Finv[K][K];    // Frobenius 行列 (列 i = x^(i*Pr) mod g) の逆行列
static long q_rootexp;        // Pr^(E-1) : 体の元の p 乗根を取る指数

// ---------- 体 GF(Pr^E) の演算 ----------
static short q_add(short a, short b)
{
    int r = 0, m = 1;
    for (int k = 0; k < E; k++)
    {
        r += (((a % Pr) + (b % Pr)) % Pr) * m;
        a /= Pr;
        b /= Pr;
        m *= Pr;
    }
    return (short)r;
}
static short q_neg(short a)
{
    int r = 0, m = 1;
    for (int k = 0; k < E; k++)
    {
        r += ((Pr - (a % Pr)) % Pr) * m;
        a /= Pr;
        m *= Pr;
    }
    return (short)r;
}
static short q_sub(short a, short b) { return q_add(a, q_neg(b)); }
static short q_mul(short a, short b)
{
    if (a == 0 || b == 0)
        return 0;
    return gf[mlt(fg[a], fg[b])];
}
static short q_pow(short a, long e)
{
    short r = 1, b = a;
    while (e > 0)
    {
        if (e & 1)
            r = q_mul(r, b);
        b = q_mul(b, b);
        e >>= 1;
    }
    return r;
}
static short q_inv(short a) { return q_pow(a, O - 2); }

// ---------- 剰余環 GF(q)[x]/g の演算 ----------
static q_res q_zero(void)
{
    q_res r;
    memset(&r, 0, sizeof r);
    return r;
}
static q_res q_radd(q_res a, q_res b)
{
    q_res r;
    for (int i = 0; i < K; i++)
        r.c[i] = q_add(a.c[i], b.c[i]);
    return r;
}
static q_res q_rsub(q_res a, q_res b)
{
    q_res r;
    for (int i = 0; i < K; i++)
        r.c[i] = q_sub(a.c[i], b.c[i]);
    return r;
}
static q_res q_rscale(q_res a, short c)
{
    q_res r;
    for (int i = 0; i < K; i++)
        r.c[i] = q_mul(a.c[i], c);
    return r;
}
// x を掛けて g で還元 (g はモニック)
static q_res q_mulx(q_res a)
{
    q_res r;
    short top = a.c[K - 1];
    for (int i = K - 1; i > 0; i--)
        r.c[i] = a.c[i - 1];
    r.c[0] = 0;
    if (top)
        for (int i = 0; i < K; i++)
            r.c[i] = q_sub(r.c[i], q_mul(top, q_g[i]));
    return r;
}
static q_res q_rmul(q_res a, q_res b)
{
    q_res r = q_zero();
    for (int i = K - 1; i >= 0; i--)
    {
        r = q_mulx(r);
        if (a.c[i])
            r = q_radd(r, q_rscale(b, a.c[i]));
    }
    return r;
}

// K x K 行列の逆行列。特異なら 0 を返す
static int q_matinv(short A[K][K], short Out[K][K])
{
    static short W[K][2 * K];
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
        {
            W[i][j] = A[i][j];
            W[i][K + j] = (i == j);
        }
    for (int c = 0; c < K; c++)
    {
        int pr = -1;
        for (int i = c; i < K; i++)
            if (W[i][c])
            {
                pr = i;
                break;
            }
        if (pr < 0)
            return 0;
        if (pr != c)
            for (int j = 0; j < 2 * K; j++)
            {
                short tmp = W[c][j];
                W[c][j] = W[pr][j];
                W[pr][j] = tmp;
            }
        short iv = q_inv(W[c][c]);
        for (int j = 0; j < 2 * K; j++)
            W[c][j] = q_mul(W[c][j], iv);
        for (int i = 0; i < K; i++)
            if (i != c && W[i][c])
            {
                short fct = W[i][c];
                for (int j = 0; j < 2 * K; j++)
                    W[i][j] = q_sub(W[i][j], q_mul(fct, W[c][j]));
            }
    }
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
            Out[i][j] = W[i][K + j];
    return 1;
}

// g を設定し、p 乗根のための Frobenius 逆行列を作る。g が無平方でなければ 0
static int q_setup(const short *gcoef)
{
    short lead = gcoef[K];
    if (lead == 0)
        return 0;
    q_g_lead = lead;
    short il = q_inv(lead);
    for (int i = 0; i <= K; i++)
        q_g[i] = q_mul(gcoef[i], il);

    q_rootexp = 1;
    for (int i = 0; i < E - 1; i++)
        q_rootexp *= Pr;

    q_res xp = q_zero();
    xp.c[0] = 1;
    for (int i = 0; i < Pr; i++)
        xp = q_mulx(xp); // x^Pr mod g

    static short A[K][K];
    q_res col = q_zero();
    col.c[0] = 1;
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < K; j++)
            A[j][i] = col.c[j];
        col = q_rmul(col, xp); // x^(i*Pr)
    }
    return q_matinv(A, q_Finv);
}

// z の p 乗根 (v^p = z mod g となる v)。 (Σ y_i x^i)^p = Σ y_i^p x^(i p)
static q_res q_root(q_res z)
{
    q_res r;
    for (int i = 0; i < K; i++)
    {
        short w = 0;
        for (int j = 0; j < K; j++)
            if (z.c[j])
                w = q_add(w, q_mul(q_Finv[i][j], z.c[j]));
        r.c[i] = q_pow(w, q_rootexp); // 係数の p 乗根は w^(Pr^(E-1))
    }
    return r;
}

// synd() の vec.x に格納された逆順シンドロームを、g の係数から Patterson 用に変換する。
static q_res q_bm_to_patterson_syndrome(const short *syndrome)
{
    static short transform[K][K];
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
            transform[i][j] = (i + j < K)
                                  ? q_mul(q_g_lead, q_g[i + j + 1])
                                  : 0;

    q_res s = q_zero();
    for (int i = 0; i < K; i++)
    {
        short acc = 0;
        for (int j = 0; j + i < K; j++)
            acc = q_add(acc, q_mul(transform[i][j], syndrome[K - 1 - j]));
        s.c[i] = q_neg(acc);
    }
    return s;
}

// s の g を法とする逆元。s が可逆でなければ 0
static int q_sinv(q_res s, q_res *out)
{
    static short A[K][K], Ai[K][K];
    q_res col = s;
    for (int k = 0; k < K; k++)
    {
        for (int j = 0; j < K; j++)
            A[j][k] = col.c[j];
        col = q_mulx(col);
    }
    if (!q_matinv(A, Ai))
        return 0;
    for (int j = 0; j < K; j++)
        out->c[j] = Ai[j][0];
    return 1;
}

static short q_eval(const short *c, int d, short x)
{
    short r = 0;
    for (int i = d; i >= 0; i--)
        r = q_add(q_mul(r, x), c[i]);
    return r;
}

// ---------- 復号本体 ----------
// syndrome[j] = synd() が返す vec.x[j] (通常順シンドローム)。
// gcoef[i] = g の x^i の係数 (i = 0..K)。
// φ = phimin..phimax を試す。誤りの値が全部 1 なら phimin=phimax=1 でよい。
// 戻り値: 見つかった候補数 (>=0)。-1: g が無平方でない/s が可逆でない、-2: s = 0 (誤りなし)。
// errs[c][i] = 候補 c の位置 i (体の元 i) の誤りの値 (0 は誤りなし)。
int patterson_q(const short *syndrome, const short *gcoef, int phimin, int phimax, short errs[Q_MAXC][M])
{
    if (!q_setup(gcoef))
        return -1;
    //for (int i = 0; i < K; i++)
    //    printf("input syndrome[%d]=%d\n", i, syndrome[i]);

    q_res s = q_bm_to_patterson_syndrome(syndrome);

    for (int i = 0; i < K; i++)
        printf("syndrome[%d] = bm[%4d]  pat[%4d]\n", i, syndrome[i],s.c[i]);
        //exit(1);

    int nz = 0;
    for (int i = 0; i < K; i++)
        if (s.c[i])
            nz = 1;
    if (!nz)
        return -2;
    q_res sinv;
    if (!q_sinv(s, &sinv))
        return -1;

    // φ に依存しない部分を前計算: v_j = root(x^j) - (φ j) * root(x^(j-1)/s)
    int jmax = (Pr - 1 < K) ? Pr - 1 : K;
    q_res xpow[K + 1], rootA[K + 1], rootB[K + 1];
    q_res cur = q_zero();
    cur.c[0] = 1;
    for (int j = 0; j <= jmax; j++)
    {
        xpow[j] = cur;
        cur = q_mulx(cur);
    }
    for (int j = 0; j <= jmax; j++)
    {
        rootA[j] = q_root(xpow[j]);
        if (j > 0)
            rootB[j] = q_root(q_rmul(sinv, xpow[j - 1]));
    }

    int found = 0;
    static short B[K][K + 1];
    for (int phi = phimin; phi <= phimax; phi++)
    {
        // 係数行列: 列 e = Pr*i + j  は  x^i * v_j mod g
        memset(B, 0, sizeof B);
        for (int j = 0; j <= jmax; j++)
        {
            short c = (short)((phi * j) % Pr);
            q_res v = rootA[j];
            if (j > 0 && c)
                v = q_rsub(v, q_rscale(rootB[j], c));
            for (int i = 0; Pr * i + j <= K; i++)
            {
                int e = Pr * i + j;
                for (int r = 0; r < K; r++)
                    B[r][e] = v.c[r];
                v = q_mulx(v);
            }
        }

        // 列を次数の小さい順に掃き出し、最初のピボットなしの列 f を自由変数にする。
        // これが次数最小の核ベクトル (= σ の次数 f)。
        int rr = 0, f = -1;
        for (int c = 0; c <= K; c++)
        {
            int pr = -1;
            for (int i = rr; i < K; i++)
                if (B[i][c])
                {
                    pr = i;
                    break;
                }
            if (pr < 0)
            {
                f = c;
                break;
            }
            if (pr != rr)
                for (int col = 0; col <= K; col++)
                {
                    short tmp = B[rr][col];
                    B[rr][col] = B[pr][col];
                    B[pr][col] = tmp;
                }
            short iv = q_inv(B[rr][c]);
            for (int col = c; col <= K; col++)
                B[rr][col] = q_mul(B[rr][col], iv);
            for (int i = 0; i < K; i++)
                if (i != rr && B[i][c])
                {
                    short fct = B[i][c];
                    for (int col = c; col <= K; col++)
                        B[i][col] = q_sub(B[i][col], q_mul(fct, B[rr][col]));
                }
            rr++;
        }
        if (f <= 0)
            continue; // f=0 なら σ=1 (誤りなし)で s≠0 と矛盾

        // σ_e = (x_e)^p  (σ = Σ_j x^j a_j(x)^p なので係数は p 乗)
        short sg[K + 1];
        memset(sg, 0, sizeof sg);
        sg[f] = 1;
        for (int c = 0; c < f; c++)
            sg[c] = q_pow(q_neg(B[c][f]), Pr);

        // σ の根と重複度 → e_i = φ μ_i
        short e[M];
        memset(e, 0, sizeof e);
        int total = 0, bad = 0;
        for (int a = 0; a < M && !bad; a++)
        {
            if (q_eval(sg, f, (short)a) != 0)
                continue;
            short cp[K + 1];
            memcpy(cp, sg, sizeof cp);
            int d = f, mu = 0;
            while (d > 0 && q_eval(cp, d, (short)a) == 0)
            {
                short qq[K + 1];
                qq[d - 1] = cp[d];
                for (int j = d - 1; j >= 1; j--)
                    qq[j - 1] = q_add(cp[j], q_mul((short)a, qq[j]));
                d--;
                for (int j = 0; j <= d; j++)
                    cp[j] = qq[j];
                mu++;
            }
            if (mu >= Pr)
            {
                bad = 1;
                break;
            }
            e[a] = (short)((phi * mu) % Pr);
            total += mu;
            if (total > f)
                bad = 1;
        }
        if (bad || total != f)
            continue;

        // synd() の出力と同じ通常順で検証する。
        int okv = 1;
        for (int j = 0; j < K && okv; j++)
        {
            short acc = 0;
            for (int i = 0; i < M; i++)
                if (e[i])
                    acc = q_add(acc, q_mul(e[i], mat[j][i]));
            if (acc != syndrome[K - 1 - j])
                okv = 0;
        }
        if (!okv)
            continue;

        int dup = 0;
        for (int c = 0; c < found; c++)
            if (memcmp(errs[c], e, sizeof e) == 0)
                dup = 1;
        if (!dup && found < Q_MAXC)
        {
            memcpy(errs[found], e, sizeof e);
            found++;
        }
    }
    return found;
}

// ---------- 実験 ----------
// wmin..wmax : 誤りの個数の範囲
// ng         : g を取り直す回数 (vv(K) を ng 回呼ぶ。1 回が重いなら小さく)
// trials     : g 1 個あたり、w ごとの試行回数
// maxval     : 誤りの値を 1..maxval から一様に選ぶ (1 なら全部 1、Pr-1 なら一様)
// phimax     : φ を 1..phimax で試す (maxval=1 なら 1 でよい)
void patterson_experiment(int wmin, int wmax, int ng, int trials, int maxval, int phimax)
{
    static int okc[K + 2], ambc[K + 2], failc[K + 2], skipc[K + 2];
    static short errs[Q_MAXC][M];
    static short ev[M];
    memset(okc, 0, sizeof okc);
    memset(ambc, 0, sizeof ambc);
    memset(failc, 0, sizeof failc);
    memset(skipc, 0, sizeof skipc);
    if (wmax > K)
        wmax = K;

    for (int gi = 0; gi < ng; gi++)
    {
        vec gp = vv(K); // mat も作り直される
        short gc[K + 1];
        for (int i = 0; i <= K; i++)
            gc[i] = gp.x[i];

        for (int w = wmin; w <= wmax; w++)
            for (int tr = 0; tr < trials; tr++)
            {
                memset(ev, 0, sizeof ev);
                int cnt = 0;
                while (cnt < w)
                {
                    int i = rand() % M;
                    if (ev[i])
                        continue;
                    ev[i] = (short)(1 + rand() % maxval);
                    cnt++;
                }
                short S[K];
                for (int j = 0; j < K; j++)
                {
                    short acc = 0;
                    for (int i = 0; i < M; i++)
                        if (ev[i])
                            acc = q_add(acc, q_mul(ev[i], mat[j][i]));
                    S[K - 1 - j] = acc;
                }
                int n = patterson_q(S, gc, 1, phimax, errs);
                if (n < 0)
                {
                    skipc[w]++;
                    continue;
                }
                int hit = 0;
                for (int c = 0; c < n; c++)
                    if (memcmp(errs[c], ev, sizeof ev) == 0)
                        hit = 1;
                if (hit && n == 1)
                    okc[w]++;
                else if (hit)
                    ambc[w]++;
                else
                    failc[w]++;
            }
    }

    printf("\n==== 一般化 Patterson 復号 (Pr=%d, E=%d, K=%d, 誤りの値 1..%d, φ=1..%d) ====\n",
           Pr, E, K, maxval, phimax);
    printf("  w : 成功 / 曖昧 / 失敗 / スキップ\n");
    for (int w = wmin; w <= wmax; w++)
        printf(" %2d : %d / %d / %d / %d\n", w, okc[w], ambc[w], failc[w], skipc[w]);
}

// synd() の出力と g の係数から Patterson 復号する。
// 戻り値: 候補数 (負なら復号不能)。
int patterson_check(const short *syndrome, const short *gcoef, int phimin, int phimax)
{
    static short errs[Q_MAXC][M];
    int n = patterson_q(syndrome, gcoef, phimin, phimax, errs);
    printf("Patterson: 候補数=%d", n);
    if (n == -1)
        printf(" (g が無平方でない / s が可逆でない)\n");
    else if (n == -2)
        printf(" (s=0: 誤りなし)\n");
    else
        printf("\n");
    if (n <= 0)
        return n;

    for (int c = 0; c < n; c++)
    {
        printf("--- 候補 %d\n", c);
        printf("pos  decoded\n");
        for (int i = 0; i < M; i++)
            if (errs[c][i])
                printf("i= %4d  %7d\n", i, errs[c][i]);
    }
    for (int i = 0; i < K; i++)
        printf("syn[%d]=%d\n", i, syndrome[i]);

    return n;
}
