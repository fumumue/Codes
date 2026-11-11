// 比較実験: 同じ既約 g から作る Goppa 符号を、一般化 Patterson 復号と alternant 復号 (BM + Forney) で復号して比べる
//   符号 1: Γ(L,g)            検査行 m t        Patterson (誤りの値が揃う: t 個まで / ばらばら: 約 2t/P 個まで) と alternant (R=t, t/2 個)
//   符号 2: Γ(L,g^(P-1)) = Γ(L,g^P)  (g 無平方)  検査行 m (P-1) t   alternant で R=(P-1)t なら (P-1)t/2 個、R=P t なら P t/2 個
//   台 L = F_q の非零元で g(α) != 0 のもの (既約 g なので n = q-1)。
// 元の説明:
// 一般化 Patterson 復号 (Barreto-Misoczki-Lindner, arXiv:1103.3296)
// 任意の素数 P と拡大次数 m (m=1 なら素体 F_P) に対応した版。
// 自己完結版: gcc -O2 -o gpatterson gpatterson.c
//
// 使い方: ./gpatterson P m t w trials mode [sign] [seed]
//
//   P      標数 (素数)。符号は F_P 上の符号。例: 2, 3, 5, 7, 31, 101
//   m      拡大次数。Goppa 多項式 g の係数体は F_q (q = P^m)。m=1 なら g も F_P 係数。
//   t      g の次数 (既約多項式)。符号の設計 (誤り訂正) パラメータ。
//   w      入れる誤りの個数 (重み)。
//   trials 試行回数 (試行ごとに g と誤りパターンを乱数で作り直す)。
//   mode   誤りの大きさ (値) の決め方
//            0: 全部 1。φ=1 だけ試す。
//            1: 全部同じ値 a (1..P-1 のどれか。復号側は a を知らない)。φ=1..P-1 を全部試す。
//            2: 位置ごとにばらばら (1..P-1 の一様乱数)。φ=1..P-1 を全部試す。
//   sign   (省略可) 鍵方程式の符号。既定 -1 (正しい)。+1 にすると失敗することの確認用。
//   seed   (省略可) 乱数の種。既定 12345。
//
// 符号: Γ(L,g) = { c in F_P^n | Σ c_i/(x-L_i) = 0 mod g }, L = {a in F_q | g(a) != 0}  (g 既約なので n = q)
// 鍵方程式: φ σ' = σ s_e (mod g),  σ_φ = Π (x-L_i)^(e_i/φ)  (指数は 0..P-1 に持ち上げる)
//
// 訂正できる誤り数の目安: mode 2 では約 (2/P) t 個まで。mode 0/1 (値が揃っている) では t 個まで。
// P=2 のときは通常の Patterson 復号 (t 個まで) に一致する。
//
// P 乗根 mod g は m=1 でも m>1 でも同じ方法: r = x^(1/P) mod g を事前計算し、
// z = Σ x^k z_k^P と分解して z^(1/P) = Σ r^k z_k とする (z_k の係数は F_q の P 乗根)。
// σ = Σ x^j a_j(x)^P を組み立てるとき、a_j の係数も P 乗する (m=1 では恒等写像)。

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXQ 16384 // q = P^m の上限
#define MAXD 1024 // 多項式の最大長。積の次数 2t-2 が入るよう 2t < MAXD
#define MAXN MAXQ
#define MAXC 14
#define KMAX 64 // 格子の行数 min(P-1,t)+1 の上限

static int P; // 標数
static int Mm, Q, Q1;
static int pw[16];
static int EXPT[2 * MAXQ], LOGT[MAXQ];
static int ZECH[MAXQ]; // ZECH[d] = log(1 + g^d), 和が 0 のとき -1
static int NEG[MAXQ], PPOW[MAXQ], PROOT[MAXQ];

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
static int dscale(int v, int c) // F_P の元 c 倍 (桁ごと)
{
    int r = 0;
    for (int i = 0; i < Mm; i++)
        r += ((((v / pw[i]) % P) * c) % P) * pw[i];
    return r;
}
static inline int mul(int a, int b) { return (a && b) ? EXPT[LOGT[a] + LOGT[b]] : 0; }
static inline int inv(int a) { return EXPT[Q1 - LOGT[a]]; }
// a + b = a (1 + b/a) = g^(log a + ZECH[log b - log a])  (表は大きさ q だけ)
static inline int fadd(int a, int b)
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
    // 原始多項式 x^m + Σ pc_i x^i を探す (x の位数が q-1)
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
    for (int a = 0; a < Q; a++)
        PPOW[a] = a ? EXPT[(int)((long)LOGT[a] * P % Q1)] : 0; // a^P
    for (int a = 0; a < Q; a++)
        PROOT[PPOW[a]] = a; // a^(1/P)
    return 1;
}

// ---------- F_q[x] の多項式 (d = 次数, 零多項式は d = -1) ----------
typedef struct
{
    int d;
    unsigned short c[MAXD];
} poly;

static void ptrim(poly *a)
{
    while (a->d >= 0 && a->c[a->d] == 0)
        a->d--;
}
static poly pzero(void)
{
    poly r;
    memset(&r, 0, sizeof r);
    r.d = -1;
    return r;
}
static poly pconst(int v)
{
    poly r = pzero();
    if (v)
    {
        r.c[0] = v;
        r.d = 0;
    }
    return r;
}
static poly padd(poly a, poly b)
{
    poly r = pzero();
    int m = a.d > b.d ? a.d : b.d;
    for (int i = 0; i <= m; i++)
        r.c[i] = fadd(i <= a.d ? a.c[i] : 0, i <= b.d ? b.c[i] : 0);
    r.d = m;
    ptrim(&r);
    return r;
}
static poly psub(poly a, poly b)
{
    poly r = pzero();
    int m = a.d > b.d ? a.d : b.d;
    for (int i = 0; i <= m; i++)
        r.c[i] = fadd(i <= a.d ? a.c[i] : 0, NEG[i <= b.d ? b.c[i] : 0]);
    r.d = m;
    ptrim(&r);
    return r;
}
static poly pscale(poly a, int s)
{
    poly r = pzero();
    if (s == 0 || a.d < 0)
        return r;
    for (int i = 0; i <= a.d; i++)
        r.c[i] = mul(a.c[i], s);
    r.d = a.d;
    ptrim(&r);
    return r;
}
static poly pshift(poly a, int n)
{
    if (a.d < 0)
        return a;
    if (a.d + n >= MAXD)
    {
        fprintf(stderr, "pshift overflow\n");
        exit(1);
    }
    poly r = pzero();
    for (int i = 0; i <= a.d; i++)
        r.c[i + n] = a.c[i];
    r.d = a.d + n;
    return r;
}
static poly pmul(poly a, poly b)
{
    poly r = pzero();
    if (a.d < 0 || b.d < 0)
        return r;
    if (a.d + b.d >= MAXD)
    {
        fprintf(stderr, "pmul overflow\n");
        exit(1);
    }
    for (int i = 0; i <= a.d; i++)
        for (int j = 0; j <= b.d; j++)
            r.c[i + j] = fadd(r.c[i + j], mul(a.c[i], b.c[j]));
    r.d = a.d + b.d;
    ptrim(&r);
    return r;
}
static void pdivmod(poly a, poly b, poly *q, poly *r)
{
    poly Qp = pzero(), R = a;
    int ib = inv(b.c[b.d]);
    while (R.d >= b.d)
    {
        int c = mul(R.c[R.d], ib);
        int sh = R.d - b.d;
        Qp.c[sh] = c;
        if (Qp.d < sh)
            Qp.d = sh;
        for (int i = 0; i <= b.d; i++)
            R.c[i + sh] = fadd(R.c[i + sh], NEG[mul(c, b.c[i])]);
        ptrim(&R);
    }
    ptrim(&Qp);
    *q = Qp;
    *r = R;
}
static poly pmod(poly a, poly g)
{
    poly q, r;
    pdivmod(a, g, &q, &r);
    return r;
}
static poly pmulmod(poly a, poly b, poly g) { return pmod(pmul(a, b), g); }
static int peval(poly a, int x)
{
    int r = 0;
    for (int i = a.d; i >= 0; i--)
        r = fadd(mul(r, x), a.c[i]);
    return r;
}
static poly pgcd(poly a, poly b)
{
    while (b.d >= 0)
    {
        poly t = pmod(a, b);
        a = b;
        b = t;
    }
    return a;
}
static poly ppowmod(poly base, long e, poly g)
{
    poly r = pconst(1);
    base = pmod(base, g);
    while (e)
    {
        if (e & 1)
            r = pmulmod(r, base, g);
        base = pmulmod(base, base, g);
        e >>= 1;
    }
    return r;
}
// a の g を法とする逆元。gcd != 1 なら 0 を返す
static int pinvmod(poly a, poly g, poly *out)
{
    poly r0 = g, r1 = pmod(a, g), t0 = pzero(), t1 = pconst(1);
    while (r1.d >= 0)
    {
        poly q, r;
        pdivmod(r0, r1, &q, &r);
        poly t2 = psub(t0, pmul(q, t1));
        r0 = r1;
        r1 = r;
        t0 = t1;
        t1 = t2;
    }
    if (r0.d != 0)
        return 0;
    *out = pmod(pscale(t0, inv(r0.c[0])), g);
    return 1;
}

// Ben-Or の既約性判定: 1 <= i <= t/2 で gcd(g, x^(q^i) - x) = 1
static int is_irreducible(poly g)
{
    int t = g.d;
    poly x = pmod(pshift(pconst(1), 1), g), h = x;
    for (int i = 1; i <= t / 2; i++)
    {
        h = ppowmod(h, Q, g);
        poly gg = pgcd(g, psub(h, x));
        if (gg.d > 0)
            return 0;
    }
    return 1;
}

// ---------- P 乗根 mod g (g は既約) ----------
// F_q[x]/g は F_{q^t} = F_{P^(m t)}。Frobenius f -> f^P の位数は m*t なので、
// r = x^(P^(m t - 1)) が x の P 乗根 (r^P = x)。
static poly rp[KMAX + 2]; // rp[k] = r^k mod g
static poly G;
static int tdeg;

static void setup_root(poly g)
{
    G = g;
    tdeg = g.d;
    poly r = pmod(pshift(pconst(1), 1), g);
    for (int i = 0; i < Mm * tdeg - 1; i++)
        r = ppowmod(r, P, g);
    int kmax = (P - 1 < tdeg - 1) ? P - 1 : tdeg - 1; // z の次数 < t なので x^k z_k^P の k は t-1 まで
    rp[0] = pconst(1);
    for (int k = 1; k <= kmax; k++)
        rp[k] = pmulmod(rp[k - 1], r, g);
}
// z (次数 < t) = Σ_k x^k z_k(x)^P と分解し、z^(1/P) = Σ r^k z_k mod g
static poly proot(poly z)
{
    poly res = pzero();
    int kmax = (P - 1 < tdeg - 1) ? P - 1 : tdeg - 1;
    for (int k = 0; k <= kmax; k++)
    {
        poly zk = pzero();
        for (int i = 0; P * i + k <= z.d; i++)
        {
            zk.c[i] = PROOT[z.c[P * i + k]];
            zk.d = i;
        }
        ptrim(&zk);
        if (zk.d >= 0)
            res = padd(res, pmulmod(rp[k], zk, G));
    }
    return res;
}

// ---------- 符号 ----------
static int n;
static int Ls[MAXN];
static int gL[MAXN];  // g(L_i)
static poly Gcode; // 符号の g (syndrome で (g(x) - g(L_i)) / (x - L_i) をその場で計算する)

static int SKIPZERO = 0; // 1: 台 L から 0 を除く (alternant 復号は α != 0 が必要)
static void setup_code(poly g)
{
    Gcode = g;
    n = 0;
    for (int a = 0; a < Q; a++)
    {
        int gv = peval(g, a);
        if (gv == 0 || (SKIPZERO && a == 0))
            continue;
        Ls[n] = a;
        gL[n] = gv;
        n++;
    }
}
// s_e(x) = Σ e_i/(x-L_i) mod g = - Σ e_i q_{L_i}(x) / g(L_i)   (e_i は F_P = {0..P-1} の元)
static poly syndrome(const int *e)
{
    poly s = pzero();
    for (int i = 0; i < n; i++)
        if (e[i])
        {
            int t = Gcode.d, a = Ls[i];
            poly q = pzero();
            q.c[t - 1] = Gcode.c[t];
            for (int j = t - 1; j >= 1; j--)
                q.c[j - 1] = fadd(Gcode.c[j], mul(a, q.c[j]));
            q.d = t - 1;
            ptrim(&q);
            s = padd(s, pscale(q, NEG[mul(e[i], inv(gL[i]))]));
        }
    return s;
}

// ---------- 弱 Popov 形式 (Mulders-Storjohann, 列ごとの重み付き次数) ----------
static poly M[KMAX][KMAX];

static void pivot(poly *row, int C, const int *w, int *I, int *D)
{
    *I = -1;
    *D = -1000000;
    for (int j = 0; j < C; j++)
        if (row[j].d >= 0)
        {
            int dd = row[j].d + w[j];
            if (dd >= *D)
            {
                *D = dd;
                *I = j;
            }
        }
}
static void weak_popov(int R, int C, const int *w, int *I, int *D)
{
    for (;;)
    {
        for (int r = 0; r < R; r++)
            pivot(M[r], C, w, &I[r], &D[r]);
        int found = 0;
        for (int k = 0; k < R && !found; k++)
        {
            if (I[k] < 0)
                continue;
            for (int l = 0; l < R; l++)
            {
                if (l == k || I[l] != I[k] || D[l] < D[k])
                    continue;
                int c = mul(M[l][I[l]].c[M[l][I[l]].d], inv(M[k][I[k]].c[M[k][I[k]].d]));
                int e = D[l] - D[k];
                for (int j = 0; j < C; j++)
                    M[l][j] = psub(M[l][j], pscale(pshift(M[k][j], e), c));
                found = 1;
                break;
            }
        }
        if (!found)
            break;
    }
}

// ---------- 復号 (Algorithm 1) ----------
// 格子の行 r から σ(x) = Σ_j x^j a_j(x)^P を作る。a_j の係数も P 乗する (F_q 上の Frobenius)
static poly build_sigma(int r, int C)
{
    poly sg = pzero();
    for (int j = 0; j < C; j++)
        for (int i = 0; i <= M[r][j].d; i++)
        {
            int idx = i * P + j;
            if (idx >= MAXD)
                return pzero();
            sg.c[idx] = fadd(sg.c[idx], PPOW[M[r][j].c[i]]);
        }
    sg.d = MAXD - 1;
    ptrim(&sg);
    return sg;
}

static long VALID_HITS = 0; // シンドロームまで一致した σ の個数 (重複含む)

// σ の根と重複度から e_i = φ μ_i を作り、シンドロームが一致すれば候補に加える
static void try_sigma(poly sg, int phi, poly s, int cand[][MAXN], int maxc, int *found)
{
    if (sg.d < 0)
        return;
    int e[MAXN], total = 0;
    for (int i = 0; i < n; i++)
    {
        e[i] = 0;
        if (peval(sg, Ls[i]) != 0)
            continue;
        poly cur = sg;
        int mu = 0;
        while (cur.d > 0 && peval(cur, Ls[i]) == 0)
        {
            poly q = pzero();
            q.c[cur.d - 1] = cur.c[cur.d];
            for (int j = cur.d - 1; j >= 1; j--)
                q.c[j - 1] = fadd(cur.c[j], mul(Ls[i], q.c[j]));
            q.d = cur.d - 1;
            cur = q;
            mu++;
        }
        if (mu >= P || mu == 0)
            return;
        e[i] = (phi * mu) % P;
        total += mu;
    }
    if (total != sg.d)
        return; // 根がすべて台の中にあるわけではない
    poly df = psub(syndrome(e), s);
    if (df.d >= 0)
        return;
    VALID_HITS++;
    for (int c = 0; c < *found; c++)
        if (memcmp(cand[c], e, sizeof(int) * n) == 0)
            return;
    if (*found < maxc)
    {
        memcpy(cand[*found], e, sizeof(int) * n);
        (*found)++;
    }
}

static int SIGN = -1; // v_k^P = x^k + SIGN*φ k x^(k-1)/s  (導出では -1)

// ---------- 線形代数版 ----------
// a_0 + Σ_{k=1}^{K} a_k v_k ≡ 0 (mod g) を、a_k の係数 (deg a_k <= (t-k)/P) を未知数とする
// F_q 上の連立一次方程式として解く。未知数は全部で t+1 個、式は t 本 (g で割った余りの係数)。
// 未知数 a_{k,i} (a_k の x^i の係数) には idx = i*P + k を割り当てる。σ の x^idx の係数が a_{k,i}^P なので、
// idx が小さい順にガウス・ジョルダン法をかけると、自由変数 f ごとの基底ベクトルは「最高次の項が x^f」の σ になる。
// 最小次数の解 (= 自由変数が最小の基底ベクトル) が正解の σ 候補。
static int ALGO = 0; // 0: 弱 Popov (格子簡約), 1: 線形代数 (ガウス・ジョルダン)
static poly VK[KMAX];
static unsigned short LA[MAXD / 2][MAXD / 2 + 1];
static int LA_pivcol[MAXD / 2], LA_rank, LA_free[MAXD / 2 + 1], LA_nfree;

static void la_solve(poly g, int K)
{
    int t = g.d, N = t + 1, cur[MAXD / 2];
    for (int k = 0; k <= K; k++)
    {
        for (int j = 0; j < t; j++)
            cur[j] = 0;
        if (k == 0)
            cur[0] = 1;
        else
            for (int j = 0; j <= VK[k].d; j++)
                cur[j] = VK[k].c[j];
        for (int i = 0; i * P + k <= t; i++)
        {
            int idx = i * P + k;
            for (int j = 0; j < t; j++)
                LA[j][idx] = (unsigned short)cur[j];
            int top = cur[t - 1]; // cur <- x * cur mod g (g は monic)
            for (int j = t - 1; j >= 1; j--)
                cur[j] = fadd(cur[j - 1], NEG[mul(top, g.c[j])]);
            cur[0] = NEG[mul(top, g.c[0])];
        }
    }
    int r = 0;
    LA_nfree = 0;
    for (int c = 0; c < N; c++)
    {
        int p = -1;
        for (int j = r; j < t; j++)
            if (LA[j][c])
            {
                p = j;
                break;
            }
        if (p < 0)
        {
            LA_free[LA_nfree++] = c;
            continue;
        }
        if (p != r)
            for (int cc = c; cc < N; cc++)
            {
                unsigned short tmp = LA[p][cc];
                LA[p][cc] = LA[r][cc];
                LA[r][cc] = tmp;
            }
        int iv = inv(LA[r][c]);
        for (int cc = c; cc < N; cc++)
            LA[r][cc] = (unsigned short)mul(LA[r][cc], iv);
        for (int j = 0; j < t; j++)
            if (j != r && LA[j][c])
            {
                int f = LA[j][c];
                for (int cc = c; cc < N; cc++)
                    LA[j][cc] = (unsigned short)fadd(LA[j][cc], NEG[mul(f, LA[r][cc])]);
            }
        LA_pivcol[r++] = c;
    }
    LA_rank = r;
}
// j 番目の自由変数に対応する基底ベクトルから σ を作る
static poly la_basis_sigma(int j)
{
    int f = LA_free[j];
    poly sg = pzero();
    sg.c[f] = 1;
    for (int r = 0; r < LA_rank; r++)
        if (LA_pivcol[r] < f && LA[r][f])
            sg.c[LA_pivcol[r]] = PPOW[NEG[LA[r][f]]]; // a^P (σ の係数は a_{k,i}^P)
    sg.d = f;
    ptrim(&sg);
    return sg;
}

// 戻り値: 候補数 (>=0)、s が可逆でなければ -1
static int decode(poly g, poly s, int phimin, int phimax, int cand[][MAXN], int maxc)
{
    int t = g.d, found = 0;
    poly sinv;
    if (!pinvmod(s, g, &sinv))
        return -1;

    int K = (P - 1 < t) ? P - 1 : t; // a_k は k > t なら 0 に決まるので行を落とせる
    int R = K + 1, C = K + 1;
    int b[KMAX] = {0}, w[KMAX], I[KMAX], D[KMAX];
    for (int j = 0; j < C; j++)
        b[j] = (t - j) / P;
    int B = b[0];
    for (int j = 0; j < C; j++)
        w[j] = B - b[j];

    for (int phi = phimin; phi <= phimax; phi++)
    {
        for (int k = 1; k <= K; k++)
        {
            poly xk = pmod(pshift(pconst(1), k), g);
            poly xk1 = pmod(pshift(pconst(1), k - 1), g);
            poly term = pscale(pmulmod(xk1, sinv, g), (int)((long)phi * k % P));
            poly u = (SIGN < 0) ? psub(xk, term) : padd(xk, term);
            VK[k] = proot(u);
        }

        if (ALGO == 1)
        {
            la_solve(g, K);
            if (LA_nfree == 0)
                continue;
            long before = VALID_HITS;
            try_sigma(la_basis_sigma(0), phi, s, cand, maxc, &found); // 最小次数の解
            if (VALID_HITS == before && LA_nfree > 1)                // 当たらなければ他の基底も試す
            {
                for (int j = 1; j < LA_nfree; j++)
                    try_sigma(la_basis_sigma(j), phi, s, cand, maxc, &found);
                if (LA_nfree == 2) // 2 次元なら F_q 倍の組み合わせも試す (弱 Popov 版の nr==2 と同じ)
                {
                    poly s0 = la_basis_sigma(0), s1 = la_basis_sigma(1);
                    for (int c = 1; c < Q; c++)
                        try_sigma(padd(s1, pscale(s0, c)), phi, s, cand, maxc, &found);
                }
            }
            continue;
        }

        for (int r = 0; r < R; r++)
            for (int j = 0; j < C; j++)
                M[r][j] = pzero();
        M[0][0] = g;
        for (int k = 1; k <= K; k++)
        {
            M[k][0] = psub(pzero(), VK[k]);
            M[k][k] = pconst(1);
        }
        weak_popov(R, C, w, I, D);

        int rows[KMAX], nr = 0;
        for (int r = 0; r < R; r++)
            if (I[r] >= 0 && D[r] <= B)
                rows[nr++] = r;
        for (int q = 0; q < nr; q++)
            try_sigma(build_sigma(rows[q], C), phi, s, cand, maxc, &found);
        if (nr == 2) // 正解が 2 行の一次結合のことがあるので、F_q 倍の組み合わせも試す
        {
            poly s0 = build_sigma(rows[0], C), s1 = build_sigma(rows[1], C);
            for (int c = 1; c < Q; c++)
                try_sigma(padd(s0, pscale(s1, c)), phi, s, cand, maxc, &found);
        }
    }
    return found;
}

// ---------- テスト ----------
static poly random_irreducible_g(int t)
{
    for (;;)
    {
        poly g = pzero();
        for (int i = 0; i < t; i++)
            g.c[i] = rand() % Q;
        g.c[t] = 1;
        g.d = t;
        if (g.c[0] == 0)
            continue;
        if (is_irreducible(g))
            return g;
    }
}

static int is_prime(int p)
{
    if (p < 2)
        return 0;
    for (int d = 2; d * d <= p; d++)
        if (p % d == 0)
            return 0;
    return 1;
}

// ---------- alternant 復号 (Berlekamp-Massey + Forney) ----------
#define MAXR 1024
static inline int fsub(int a, int b) { return fadd(a, NEG[b]); }
static int peval_arr(const int *c, int d, int x)
{
    int r = 0;
    for (int i = d; i >= 0; i--)
        r = fadd(mul(r, x), c[i]);
    return r;
}
// 検査行列 H = [ Yw_i α_i^j ] (j<R) に対するシンドローム
static void alt_syndromes(const int *v, int R, const int *Yw, int *s)
{
    for (int j = 0; j < R; j++)
        s[j] = 0;
    for (int i = 0; i < n; i++)
        if (v[i])
        {
            int Yk = mul(v[i], Yw[i]), pwr = 1;
            for (int j = 0; j < R; j++)
            {
                s[j] = fadd(s[j], mul(Yk, pwr));
                pwr = mul(pwr, Ls[i]);
            }
        }
}
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
static int alt_decode(const int *s, int R, const int *Yw, int *e)
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
        if (peval_arr(Lam, L, inv(Ls[i])) == 0)
        {
            if (cnt >= L)
                return 0;
            pos[cnt++] = i;
        }
    if (cnt != L)
        return 0;
    for (int k = 0; k < L; k++)
    {
        int v = 0;
        for (int i = 0; i <= k; i++)
            v = fadd(v, mul(Lam[i], s[k - i]));
        Om[k] = v;
    }
    for (int i = 1; i <= L; i++)
        dL[i - 1] = mul(i % P, Lam[i]);
    for (int k = 0; k < L; k++)
    {
        int i = pos[k], X = Ls[i], xi = inv(X);
        int dv = peval_arr(dL, L - 1, xi);
        if (dv == 0)
            return 0;
        int Y = NEG[mul(X, mul(peval_arr(Om, L - 1, xi), inv(dv)))];
        int ev = mul(Y, inv(Yw[i]));
        if (ev == 0 || ev >= P)
            return 0;
        e[i] = ev;
    }
    alt_syndromes(e, R, Yw, chk);
    for (int j = 0; j < R; j++)
        if (chk[j] != s[j])
            return 0;
    return 1;
}

// F_P 上の検査行列 (m R x n) の階数 (= n - 次元)。大きすぎるときは -1
static int rank_FP(int R, const int *Yw)
{
    size_t rows = (size_t)Mm * R;
    if (rows * (size_t)n > 6000000UL)
        return -1;
    unsigned short *H = calloc(rows * n, sizeof(unsigned short));
    for (int i = 0; i < n; i++)
    {
        int v = Yw[i];
        for (int j = 0; j < R; j++)
        {
            for (int d = 0; d < Mm; d++)
                H[((size_t)j * Mm + d) * n + i] = (unsigned short)((v / pw[d]) % P);
            v = mul(v, Ls[i]);
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
            for (int cc = col; cc < n; cc++)
            {
                unsigned short tmp = H[(size_t)p * n + cc];
                H[(size_t)p * n + cc] = H[(size_t)rank * n + cc];
                H[(size_t)rank * n + cc] = tmp;
            }
        long iv = 1, b = H[(size_t)rank * n + col];
        for (int ex = P - 2; ex; ex >>= 1)
        {
            if (ex & 1)
                iv = iv * b % P;
            b = b * b % P;
        }
        for (int j = rank + 1; j < (int)rows; j++)
            if (H[(size_t)j * n + col])
            {
                long f = (long)H[(size_t)j * n + col] * iv % P;
                for (int cc = col; cc < n; cc++)
                    H[(size_t)j * n + cc] = (unsigned short)((H[(size_t)j * n + cc] + (long)(P - f) * H[(size_t)rank * n + cc]) % P);
            }
        rank++;
    }
    free(H);
    return rank;
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "使い方: %s P m t trials [wmax] [seed]\n"
            "  P      標数 (素数)\n"
            "  m      拡大次数 (q=P^m <= %d)\n"
            "  t      Goppa 多項式 g の次数 (既約多項式を乱択, 2 <= t, 2t < %d)\n"
            "  trials 試行回数 (試行ごとに g を作り直す)\n"
            "  wmax   (省略可) 調べる誤り数の最大値。既定 P*t/2+2\n"
            "  seed   (省略可) 乱数の種。既定 12345\n"
            "例: %s 3 4 8 100      (F_81, deg g=8)\n"
            "    %s 2 5 6 200      (P=2: Patterson と alternant (R=2t) が同じ t 個を訂正する確認)\n",
            prog, MAXQ, MAXD, prog, prog);
}

static void rnd_error(int *e, int w, int mode) // mode 1: 値が全部同じ, 2: ばらばら
{
    memset(e, 0, sizeof(int) * n);
    int a = 1 + rand() % (P - 1), cnt = 0;
    while (cnt < w)
    {
        int i = rand() % n;
        if (e[i])
            continue;
        e[i] = (mode == 1) ? a : 1 + rand() % (P - 1);
        cnt++;
    }
}

int main(int argc, char **argv)
{
    if (argc < 5)
    {
        usage(argv[0]);
        return 1;
    }
    int p = atoi(argv[1]), m = atoi(argv[2]), t = atoi(argv[3]), trials = atoi(argv[4]);
    int wmax = argc > 5 ? atoi(argv[5]) : 0;
    unsigned seed = argc > 6 ? (unsigned)atoi(argv[6]) : 12345u;
    if (!is_prime(p) || m < 1 || t < 2 || 2 * t >= MAXD || trials < 1)
    {
        usage(argv[0]);
        return 1;
    }
    if (!init_field(p, m))
    {
        fprintf(stderr, "q=P^m が大きすぎます (上限 %d)\n", MAXQ);
        return 1;
    }
    if ((P - 1 < t ? P - 1 : t) + 1 > KMAX)
    {
        fprintf(stderr, "min(P-1,t)+1 が KMAX=%d を超えています\n", KMAX);
        return 1;
    }
    if (wmax <= 0)
        wmax = P * t / 2 + 2;
    srand(seed);
    ALGO = (P >= 5) ? 1 : 0; // P が小さいときは弱 Popov、大きいときは線形代数が速い
    SKIPZERO = 1;
    int Rs[3] = {t, (P - 1) * t, P * t};
    static int Yc[3][MAXN];
    int *ok = calloc((size_t)(wmax + 1) * 5, sizeof(int));
    int *e = malloc(sizeof(int) * MAXN), *er = malloc(sizeof(int) * MAXN), sa[MAXR + 2];
    static int cand[MAXC][MAXN];
    clock_t T0 = clock();
    int nn = 0, kk[3] = {-1, -1, -1};
    for (int tr = 0; tr < trials; tr++)
    {
        poly g = random_irreducible_g(t);
        setup_root(g);
        setup_code(g);
        nn = n;
        int ee[3] = {1, P - 1, P};
        for (int c = 0; c < 3; c++)
            for (int i = 0; i < n; i++)
            {
                int pe = 1;
                for (int q = 0; q < ee[c]; q++)
                    pe = mul(pe, gL[i]);
                Yc[c][i] = inv(pe);
            }
        if (tr == 0)
            for (int c = 0; c < 3; c++)
                if (Rs[c] <= MAXR)
                {
                    int rk = rank_FP(Rs[c], Yc[c]);
                    kk[c] = rk < 0 ? -1 : n - rk;
                }
        for (int w = 1; w <= wmax && w <= n; w++)
        {
            for (int mode = 1; mode <= 2; mode++) // Patterson: 1 = 値が揃う, 2 = ばらばら
            {
                rnd_error(e, w, mode);
                int f = decode(g, syndrome(e), 1, P - 1, cand, MAXC);
                int hit = 0;
                for (int c = 0; c < f; c++)
                    if (memcmp(cand[c], e, sizeof(int) * n) == 0)
                        hit = 1;
                if (hit && f == 1)
                    ok[w * 5 + (mode == 1 ? 1 : 0)]++;
                if (mode == 2) // 同じ誤りで alternant 復号 (3 通りの R)
                    for (int c = 0; c < 3; c++)
                        if (Rs[c] <= MAXR)
                        {
                            alt_syndromes(e, Rs[c], Yc[c], sa);
                            if (alt_decode(sa, Rs[c], Yc[c], er) && memcmp(e, er, sizeof(int) * n) == 0)
                                ok[w * 5 + 2 + c]++;
                        }
            }
        }
    }
    printf("P=%d m=%d t=%d  q=%d 台の大きさ n=%d  試行=%d  (Patterson の格子簡約: %s)\n", P, Mm, t, Q, nn, trials, ALGO ? "線形代数" : "弱Popov");
    printf("符号 (検査行 / 次元 k):\n");
    const char *nm[3] = {"Γ(L,g)", "Γ(L,g^(P-1))", "Γ(L,g^P)"};
    for (int c = 0; c < 3; c++)
        if (kk[c] >= 0)
            printf("  %-14s 検査行 R=%-4d (F_P 上 %d 行)  次元 k=%d\n", nm[c], Rs[c], Mm * Rs[c], kk[c]);
        else
            printf("  %-14s 検査行 R=%-4d (次元は未計算)\n", nm[c], Rs[c]);
    if (kk[1] >= 0 && kk[2] >= 0)
        printf("  Γ(L,g^P) ⊆ Γ(L,g^(P-1)) で次元が %s -> %s\n", kk[1] == kk[2] ? "等しい" : "異なる", kk[1] == kk[2] ? "同じ符号 (wild Goppa)" : "同じ符号ではない");
    printf("\n成功率 (%%)。Pat = 一般化 Patterson (Γ(L,g) のみ), Alt = alternant 復号 (R は検査行の数, ⌊R/2⌋ 個まで)\n");
    printf("  w | Pat(値が揃う) Pat(ばらばら) | Alt Γ(L,g) R=%d | Alt R=%d (=(P-1)t) | Alt R=%d (=Pt)\n", Rs[0], Rs[1], Rs[2]);
    for (int w = 1; w <= wmax && w <= nn; w++)
    {
        printf(" %3d |", w);
        printf(" %12.1f %13.1f |", 100.0 * ok[w * 5 + 1] / trials, 100.0 * ok[w * 5 + 0] / trials);
        for (int c = 0; c < 3; c++)
            if (Rs[c] <= MAXR)
                printf(" %*.1f |", c == 0 ? 15 : (c == 1 ? 19 : 17), 100.0 * ok[w * 5 + 2 + c] / trials);
            else
                printf(" %*s |", c == 0 ? 15 : (c == 1 ? 19 : 17), "-");
        printf("\n");
    }
    printf("(経過時間 %.1f 秒)\n", (double)(clock() - T0) / CLOCKS_PER_SEC);
    return 0;
}
