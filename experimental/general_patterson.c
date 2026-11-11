// 一般化 Patterson 復号 (Barreto-Misoczki-Lindner, arXiv:1103.3296)
// 標数 3 の拡大体 F_{3^m} 上の既約 Goppa 符号 Γ(L,g) (部分体 F_3 上の符号)
// 自己完結版: gcc -O2 -o gpatterson_3m gpatterson_3m.c
//
// 使い方: ./gpatterson_3m m t w trials mode [sign]
//   mode 0: 誤りの大きさは全部 1 (φ=1 だけ試す)
//   mode 1: 誤りの大きさは全部同じ値 a ∈ {1,2} (a は未知。φ=1,2 を試す)
//   mode 2: 誤りの大きさはばらばら (φ=1,2 を試す)
//
// 符号: Γ(L,g) = { c in F_3^n | Σ c_i/(x-L_i) = 0 mod g }, L = F_q (g は既約なので n = q)
// 鍵方程式: φ σ' = σ s_e (mod g),  σ_φ = Π (x-L_i)^(e_i/φ)
//
// m=1 版からの変更点:
//  (1) 3 乗根 mod g: F_q 係数では (Σ y_i x^i)^3 = Σ y_i^3 x^(3i) なので、Frobenius 行列は
//      F_q 線形でなく使えない。代わりに r = x^(1/3) mod g を事前計算し、
//      z = Σ x^k z_k^3 と分解して z^(1/3) = Σ r^k z_k とする (z_k の係数は F_q の 3 乗根)。
//  (2) σ = Σ x^j a_j(x)^3 を組み立てるとき、a_j の係数も 3 乗する。
//  (3) F_q の演算 (対数表) と、既約 g の生成 (Ben-Or の既約性判定)。

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PR 3
#define MAXQ 729 // m <= 6
#define MAXD 96
#define MAXN MAXQ
#define MAXC 8

// ---------- F_q (q = 3^m), 元は 3 進 m 桁の整数 (多項式基底) ----------
static int Mm, Q, Q1;
static int pw3[8];
static int EXPT[2 * MAXQ], LOGT[MAXQ];
static unsigned short ADD[MAXQ][MAXQ];
static int NEG[MAXQ], CUBE[MAXQ], CROOT[MAXQ];

static int dadd(int a, int b)
{
    int r = 0;
    for (int i = 0; i < Mm; i++)
        r += (((a / pw3[i]) % 3 + (b / pw3[i]) % 3) % 3) * pw3[i];
    return r;
}
static int dneg(int a)
{
    int r = 0;
    for (int i = 0; i < Mm; i++)
        r += ((3 - (a / pw3[i]) % 3) % 3) * pw3[i];
    return r;
}
static int dscale(int v, int c) { return c == 0 ? 0 : (c == 1 ? v : dadd(v, v)); }

static inline int mul(int a, int b) { return (a && b) ? EXPT[LOGT[a] + LOGT[b]] : 0; }
static inline int inv(int a) { return EXPT[Q1 - LOGT[a]]; }

static void init_field(int m)
{
    Mm = m;
    pw3[0] = 1;
    for (int i = 1; i < 8; i++)
        pw3[i] = pw3[i - 1] * 3;
    Q = pw3[m];
    Q1 = Q - 1;
    // 原始多項式 x^m + Σ pc_i x^i を探す (x の位数が q-1)
    for (int pc = 1; pc < Q; pc++)
    {
        int negpc = dneg(pc), a = 1, ord = 0;
        do
        {
            int top = a / pw3[m - 1], sh = (a % pw3[m - 1]) * 3;
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
                int top = a / pw3[m - 1], sh = (a % pw3[m - 1]) * 3;
                a = dadd(sh, dscale(negpc, top));
            }
            break;
        }
    }
    for (int a = 0; a < Q; a++)
    {
        for (int b = 0; b < Q; b++)
            ADD[a][b] = (unsigned short)dadd(a, b);
        NEG[a] = dneg(a);
    }
    for (int a = 0; a < Q; a++)
    {
        CUBE[a] = mul(a, mul(a, a));
        CROOT[CUBE[a]] = a;
    }
}

// ---------- F_q[x] の多項式 (d = 次数, 零多項式は d = -1) ----------
typedef struct
{
    int d;
    int c[MAXD];
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
        r.c[i] = ADD[i <= a.d ? a.c[i] : 0][i <= b.d ? b.c[i] : 0];
    r.d = m;
    ptrim(&r);
    return r;
}
static poly psub(poly a, poly b)
{
    poly r = pzero();
    int m = a.d > b.d ? a.d : b.d;
    for (int i = 0; i <= m; i++)
        r.c[i] = ADD[i <= a.d ? a.c[i] : 0][NEG[i <= b.d ? b.c[i] : 0]];
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
            r.c[i + j] = ADD[r.c[i + j]][mul(a.c[i], b.c[j])];
    r.d = a.d + b.d;
    ptrim(&r);
    return r;
}
static void pdivmod(poly a, poly b, poly *q, poly *r)
{
    poly Q_ = pzero(), R = a;
    int ib = inv(b.c[b.d]);
    while (R.d >= b.d)
    {
        int c = mul(R.c[R.d], ib);
        int sh = R.d - b.d;
        Q_.c[sh] = c;
        if (Q_.d < sh)
            Q_.d = sh;
        for (int i = 0; i <= b.d; i++)
            R.c[i + sh] = ADD[R.c[i + sh]][NEG[mul(c, b.c[i])]];
        ptrim(&R);
    }
    ptrim(&Q_);
    *q = Q_;
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
        r = ADD[mul(r, x)][a.c[i]];
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

// ---------- 3 乗根 mod g (g は既約) ----------
// F_q[x]/g は F_{q^t} なので、Frobenius f -> f^3 の位数は m*t。r = x^(3^(m t - 1)) が x の 3 乗根。
static poly rp[3];
static poly G;
static int tdeg;

static void setup_root(poly g)
{
    G = g;
    tdeg = g.d;
    poly r = pmod(pshift(pconst(1), 1), g);
    for (int i = 0; i < Mm * tdeg - 1; i++)
        r = pmulmod(pmulmod(r, r, g), r, g);
    rp[0] = pconst(1);
    rp[1] = r;
    rp[2] = pmulmod(r, r, g);
}
// z (次数 < t) = Σ_k x^k z_k(x)^3 と分解し、z^(1/3) = Σ r^k z_k mod g
static poly proot(poly z)
{
    poly res = pzero();
    for (int k = 0; k < PR; k++)
    {
        poly zk = pzero();
        for (int i = 0; 3 * i + k <= z.d; i++)
        {
            zk.c[i] = CROOT[z.c[3 * i + k]];
            zk.d = i;
        }
        ptrim(&zk);
        res = padd(res, pmulmod(rp[k], zk, G));
    }
    return res;
}

// ---------- 符号 ----------
static int n;
static int Ls[MAXN];
static int gL[MAXN];  // g(L_i)
static poly qL[MAXN]; // (g(x) - g(L_i)) / (x - L_i)

static void setup_code(poly g)
{
    int t = g.d;
    n = 0;
    for (int a = 0; a < Q; a++)
    {
        int gv = peval(g, a);
        if (gv == 0)
            continue;
        Ls[n] = a;
        gL[n] = gv;
        poly q = pzero();
        q.c[t - 1] = g.c[t];
        for (int j = t - 1; j >= 1; j--)
            q.c[j - 1] = ADD[g.c[j]][mul(a, q.c[j])];
        q.d = t - 1;
        ptrim(&q);
        qL[n] = q;
        n++;
    }
}
// s_e(x) = Σ e_i/(x-L_i) mod g = - Σ e_i q_{L_i}(x) / g(L_i)   (e_i は F_3 = {0,1,2} の元)
static poly syndrome(const int *e)
{
    poly s = pzero();
    for (int i = 0; i < n; i++)
        if (e[i])
            s = padd(s, pscale(qL[i], NEG[mul(e[i], inv(gL[i]))]));
    return s;
}

// ---------- 弱 Popov 形式 (Mulders-Storjohann, 列ごとの重み付き次数) ----------
#define KMAX 3
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
// 格子の行 r から σ(x) = Σ_j x^j a_j(x)^3 を作る。a_j の係数も 3 乗する (F_q 上の Frobenius)
static poly build_sigma(int r, int C, int t)
{
    poly sg = pzero();
    for (int j = 0; j < C; j++)
        for (int i = 0; i <= M[r][j].d; i++)
        {
            int idx = i * PR + j;
            if (idx >= MAXD)
                return pzero();
            sg.c[idx] = ADD[sg.c[idx]][CUBE[M[r][j].c[i]]];
        }
    sg.d = MAXD - 1;
    ptrim(&sg);
    (void)t;
    return sg;
}

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
                q.c[j - 1] = ADD[cur.c[j]][mul(Ls[i], q.c[j])];
            q.d = cur.d - 1;
            cur = q;
            mu++;
        }
        if (mu >= PR || mu == 0)
            return;
        e[i] = (phi * mu) % PR;
        total += mu;
    }
    if (total != sg.d)
        return; // 根がすべて台の中にあるわけではない
    poly df = psub(syndrome(e), s);
    if (df.d >= 0)
        return;
    for (int c = 0; c < *found; c++)
        if (memcmp(cand[c], e, sizeof(int) * n) == 0)
            return;
    if (*found < maxc)
    {
        memcpy(cand[*found], e, sizeof(int) * n);
        (*found)++;
    }
}

static int SIGN = -1; // v_k^3 = x^k + SIGN*φ k x^(k-1)/s  (導出では -1)

// 戻り値: 候補数 (>=0)、s が可逆でなければ -1
static int decode(poly g, poly s, int phimin, int phimax, int cand[][MAXN], int maxc)
{
    int t = g.d, found = 0;
    poly sinv;
    if (!pinvmod(s, g, &sinv))
        return -1;

    int K = (PR - 1 < t) ? PR - 1 : t;
    int R = K + 1, C = K + 1;
    int b[KMAX], w[KMAX], I[KMAX], D[KMAX];
    for (int j = 0; j < C; j++)
        b[j] = (t - j) / PR;
    int B = b[0];
    for (int j = 0; j < C; j++)
        w[j] = B - b[j];

    static poly sigs[MAXQ + 4];
    for (int phi = phimin; phi <= phimax; phi++)
    {
        for (int r = 0; r < R; r++)
            for (int j = 0; j < C; j++)
                M[r][j] = pzero();
        M[0][0] = g;
        for (int k = 1; k <= K; k++)
        {
            poly xk = pmod(pshift(pconst(1), k), g);
            poly xk1 = pmod(pshift(pconst(1), k - 1), g);
            poly term = pscale(pmulmod(xk1, sinv, g), (phi * k) % PR);
            poly u = (SIGN < 0) ? psub(xk, term) : padd(xk, term);
            M[k][0] = psub(pzero(), proot(u));
            M[k][k] = pconst(1);
        }
        weak_popov(R, C, w, I, D);

        int rows[KMAX], nr = 0;
        for (int r = 0; r < R; r++)
            if (I[r] >= 0 && D[r] <= B)
                rows[nr++] = r;
        int nsig = 0;
        for (int q = 0; q < nr; q++)
            sigs[nsig++] = build_sigma(rows[q], C, t);
        if (nr == 2) // 正解が 2 行の一次結合のことがあるので、F_q 倍の組み合わせも試す
            for (int c = 1; c < Q; c++)
                sigs[nsig++] = padd(sigs[0], pscale(sigs[1], c));
        for (int q = 0; q < nsig; q++)
            try_sigma(sigs[q], phi, s, cand, maxc, &found);
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

int main(int argc, char **argv)
{
    int m = argc > 1 ? atoi(argv[1]) : 3;
    int t = argc > 2 ? atoi(argv[2]) : 8;
    int w = argc > 3 ? atoi(argv[3]) : 6;
    int trials = argc > 4 ? atoi(argv[4]) : 100;
    int mode = argc > 5 ? atoi(argv[5]) : 1;
    if (argc > 6)
        SIGN = atoi(argv[6]);
    if (m < 1 || m > 6 || t < 3 || t >= (MAXD - 1) / 2)
    {
        fprintf(stderr, "m in 1..6, 3 <= t < %d\n", (MAXD - 1) / 2);
        return 1;
    }
    init_field(m);
    srand(12345);

    static int cand[MAXC][MAXN];
    int okc = 0, amb = 0, fail = 0, nonin = 0;
    for (int tr = 0; tr < trials; tr++)
    {
        poly g = random_irreducible_g(t);
        setup_root(g);
        setup_code(g);
        int e[MAXN] = {0}, cnt = 0, a = 1 + rand() % (PR - 1);
        if (mode == 0)
            a = 1;
        while (cnt < w)
        {
            int i = rand() % n;
            if (e[i])
                continue;
            e[i] = (mode == 2) ? 1 + rand() % (PR - 1) : a;
            cnt++;
        }
        poly s = syndrome(e);
        int lo = 1, hi = (mode == 0) ? 1 : PR - 1;
        int f = decode(g, s, lo, hi, cand, MAXC);
        if (f < 0)
        {
            nonin++;
            continue;
        }
        int hit = 0;
        for (int c = 0; c < f; c++)
            if (memcmp(cand[c], e, sizeof(int) * n) == 0)
                hit = 1;
        if (hit && f == 1)
            okc++;
        else if (hit)
            amb++;
        else
            fail++;
    }
    double pred = 1.0;
    for (int i = 0; i < t + 1 - w; i++)
        pred /= Q;
    printf("q=3^%d t=%d w=%d mode=%d sign=%d: 成功=%d 曖昧(候補複数)=%d 失敗=%d s不可逆=%d (試行 %d)  論文の予測成功率 1-1/q^(t+1-w)=%.6f\n",
           m, t, w, mode, SIGN, okc, amb, fail, nonin, trials, 1.0 - pred);
    return 0;
}