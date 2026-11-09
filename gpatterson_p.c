// 一般化 Patterson 復号 (Barreto-Misoczki-Lindner, arXiv:1103.3296)
// 素体 F_P 上の無平方 Goppa 符号 (m=1)。自己完結版: gcc -O2 -o gpatterson_p gpatterson_p.c
//
// 使い方: ./gpatterson_p P t w trials mode
//   mode 0: 誤りの大きさは全部 1 (φ=1 だけ試す)
//   mode 1: 誤りの大きさは全部同じ値 a (a は未知。φ=1..P-1 を試す)
//   mode 2: 誤りの大きさはばらばら (φ=1..P-1 を試す)
//
// 符号: Γ(L,g) = { c in F_P^n | Σ c_i/(x-L_i) = 0 mod g }, L = { a in F_P | g(a) != 0 }
// 鍵方程式: φ σ' = σ s_e (mod g),  σ_φ = Π (x-L_i)^(e_i/φ)  (指数は 0..P-1 に持ち上げる)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXD 128
#define MAXN 1400
#define MAXK 40
#define MAXC 8

static int P;

static inline int md(long a)
{
    a %= P;
    if (a < 0)
        a += P;
    return (int)a;
}
static int ipow(int a, long e)
{
    long r = 1, b = md(a);
    while (e)
    {
        if (e & 1)
            r = r * b % P;
        b = b * b % P;
        e >>= 1;
    }
    return (int)r;
}
static int inv(int a) { return ipow(a, P - 2); }

// ---------- F_P[x] の多項式 (d = 次数, 零多項式は d = -1) ----------
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
    v = md(v);
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
        r.c[i] = md((long)(i <= a.d ? a.c[i] : 0) + (i <= b.d ? b.c[i] : 0));
    r.d = m;
    ptrim(&r);
    return r;
}
static poly psub(poly a, poly b)
{
    poly r = pzero();
    int m = a.d > b.d ? a.d : b.d;
    for (int i = 0; i <= m; i++)
        r.c[i] = md((long)(i <= a.d ? a.c[i] : 0) - (i <= b.d ? b.c[i] : 0));
    r.d = m;
    ptrim(&r);
    return r;
}
static poly pscale(poly a, int s)
{
    poly r = pzero();
    s = md(s);
    if (s == 0 || a.d < 0)
        return r;
    for (int i = 0; i <= a.d; i++)
        r.c[i] = md((long)a.c[i] * s);
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
            r.c[i + j] = (int)((r.c[i + j] + (long)a.c[i] * b.c[j]) % P);
    r.d = a.d + b.d;
    ptrim(&r);
    return r;
}
static void pdivmod(poly a, poly b, poly *q, poly *r)
{
    poly Q = pzero(), R = a;
    int ib = inv(b.c[b.d]);
    while (R.d >= b.d)
    {
        int c = md((long)R.c[R.d] * ib);
        int sh = R.d - b.d;
        Q.c[sh] = c;
        if (Q.d < sh)
            Q.d = sh;
        for (int i = 0; i <= b.d; i++)
            R.c[i + sh] = md(R.c[i + sh] - (long)c * b.c[i]);
        ptrim(&R);
    }
    ptrim(&Q);
    *q = Q;
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
    long r = 0;
    for (int i = a.d; i >= 0; i--)
        r = (r * x + a.c[i]) % P;
    return (int)r;
}
static poly pderiv(poly a)
{
    poly r = pzero();
    for (int i = 1; i <= a.d; i++)
        r.c[i - 1] = md((long)a.c[i] * i);
    r.d = a.d - 1;
    ptrim(&r);
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

// ---------- p 乗根 mod g (g は無平方) ----------
// F_P 係数では (Σ y_i x^i)^P = Σ y_i x^(iP)。Frobenius 行列 F (列 i = x^(iP) mod g) の
// 逆行列を掛ければ P 乗根になる。g が無平方 <=> F が可逆。
static int Finv[MAXK + 2][MAXK + 2];
static int tdeg;

static int build_frobenius_inverse(poly g)
{
    int t = g.d;
    static int A[MAXK + 2][2 * (MAXK + 2)];
    poly xp = ppowmod(pshift(pconst(1), 1), P, g);
    poly col = pconst(1);
    for (int i = 0; i < t; i++)
    {
        for (int j = 0; j < t; j++)
        {
            A[j][i] = (j <= col.d) ? col.c[j] : 0;
            A[j][t + i] = (i == j);
        }
        col = pmulmod(col, xp, g);
    }
    for (int i = 0; i < t; i++)
    {
        int p = -1;
        for (int j = i; j < t; j++)
            if (A[j][i])
            {
                p = j;
                break;
            }
        if (p < 0)
            return 0; // 特異 = g が平方因子を持つ
        if (p != i)
            for (int k = 0; k < 2 * t; k++)
            {
                int tmp = A[i][k];
                A[i][k] = A[p][k];
                A[p][k] = tmp;
            }
        int iv = inv(A[i][i]);
        for (int k = 0; k < 2 * t; k++)
            A[i][k] = md((long)A[i][k] * iv);
        for (int j = 0; j < t; j++)
            if (j != i && A[j][i])
            {
                int f = A[j][i];
                for (int k = 0; k < 2 * t; k++)
                    A[j][k] = md(A[j][k] - (long)f * A[i][k]);
            }
    }
    for (int i = 0; i < t; i++)
        for (int j = 0; j < t; j++)
            Finv[i][j] = A[i][t + j];
    tdeg = t;
    return 1;
}
static poly proot(poly z)
{
    poly r = pzero();
    for (int i = 0; i < tdeg; i++)
    {
        long s = 0;
        for (int j = 0; j <= z.d && j < tdeg; j++)
            s += (long)Finv[i][j] * z.c[j];
        r.c[i] = (int)(s % P);
    }
    r.d = tdeg - 1;
    ptrim(&r);
    return r;
}

// ---------- 符号 ----------
static int n;
static int Ls[MAXN];
static int gL[MAXN]; // g(L_i)
static poly qL[MAXN]; // (g(x) - g(L_i)) / (x - L_i)

static void setup_code(poly g)
{
    int t = g.d;
    n = 0;
    for (int a = 0; a < P; a++)
    {
        int gv = peval(g, a);
        if (gv == 0)
            continue;
        Ls[n] = a;
        gL[n] = gv;
        poly q = pzero();
        q.c[t - 1] = g.c[t];
        for (int j = t - 1; j >= 1; j--)
            q.c[j - 1] = md(g.c[j] + (long)a * q.c[j]);
        q.d = t - 1;
        ptrim(&q);
        qL[n] = q;
        n++;
    }
}
// s_e(x) = Σ e_i/(x-L_i) mod g = - Σ e_i q_{L_i}(x) / g(L_i)
static poly syndrome(const int *e)
{
    poly s = pzero();
    for (int i = 0; i < n; i++)
        if (e[i])
        {
            int f = md(-(long)e[i] * inv(gL[i]));
            s = padd(s, pscale(qL[i], f));
        }
    return s;
}

// ---------- 弱 Popov 形式 (Mulders-Storjohann, 列ごとの重み付き次数) ----------
static poly M[MAXK][MAXK];

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
                int c = md((long)M[l][I[l]].c[M[l][I[l]].d] * inv(M[k][I[k]].c[M[k][I[k]].d]));
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
// 格子の行 r から σ(x) = Σ_j x^j a_j(x^P) を作る
static poly build_sigma(int r, int C, int t)
{
    poly sg = pzero();
    for (int j = 0; j < C; j++)
        for (int i = 0; i <= M[r][j].d; i++)
            sg.c[i * P + j] = md(sg.c[i * P + j] + M[r][j].c[i]);
    sg.d = t;
    ptrim(&sg);
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
                q.c[j - 1] = md(cur.c[j] + (long)Ls[i] * q.c[j]);
            q.d = cur.d - 1;
            cur = q;
            mu++;
        }
        if (cur.d == 0 && peval(cur, Ls[i]) == 0)
            mu++;
        if (mu >= P)
            return;
        e[i] = md((long)phi * mu);
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

// 戻り値: 候補数 (>=0)、s が可逆でなければ -1。cand[i][0..n-1] が誤りベクトル。
static int SIGN = -1; // v_k^p = x^k + SIGN*φ k x^(k-1)/s  (導出では -1)

static int decode(poly g, poly s, int phimin, int phimax, int cand[][MAXN], int maxc)
{
    int t = g.d, found = 0;
    poly sinv;
    if (!pinvmod(s, g, &sinv))
        return -1;

    int K = (P - 1 < t) ? P - 1 : t; // a_k は k > t なら 0 に決まるので行を落とせる
    int R = K + 1, C = K + 1;
    int b[MAXK], w[MAXK], I[MAXK], D[MAXK];
    for (int j = 0; j < C; j++)
        b[j] = (t - j) / P;
    int B = b[0];
    for (int j = 0; j < C; j++)
        w[j] = B - b[j];

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
            poly term = pscale(pmulmod(xk1, sinv, g), (long)phi * k % P);
            poly u = (SIGN < 0) ? psub(xk, term) : padd(xk, term);
            poly v = proot(u);
            M[k][0] = psub(pzero(), v);
            M[k][k] = pconst(1);
        }
        weak_popov(R, C, w, I, D);

        // D <= B の行 (短いベクトル) を集める。解の空間が 2 次元以上のときは
        // 正解が行の一次結合になっていることがあるので、2 本までは組み合わせも試す。
        int rows[MAXK], nr = 0;
        for (int r = 0; r < R; r++)
            if (I[r] >= 0 && D[r] <= B)
                rows[nr++] = r;
        poly sigs[MAXK + 2];
        int nsig = 0;
        for (int q = 0; q < nr; q++)
            sigs[nsig++] = build_sigma(rows[q], C, t);
        if (nr == 2)
            for (int c = 1; c < P; c++)
                sigs[nsig++] = padd(sigs[0], pscale(sigs[1], c));
        for (int q = 0; q < nsig; q++)
            try_sigma(sigs[q], phi, s, cand, maxc, &found);
    }
    return found;
}

// ---------- テスト ----------
static poly random_squarefree_g(int t)
{
    for (;;)
    {
        poly g = pzero();
        for (int i = 0; i < t; i++)
            g.c[i] = rand() % P;
        g.c[t] = 1;
        g.d = t;
        poly gd = pderiv(g);
        poly h = pgcd(g, gd);
        if (h.d == 0)
            return g;
    }
}

int main(int argc, char **argv)
{
    P = argc > 1 ? atoi(argv[1]) : 31;
    int t = argc > 2 ? atoi(argv[2]) : 10;
    int w = argc > 3 ? atoi(argv[3]) : 5;
    int trials = argc > 4 ? atoi(argv[4]) : 100;
    int mode = argc > 5 ? atoi(argv[5]) : 0;
    if (argc > 6)
        SIGN = atoi(argv[6]);
    srand(12345);

    static int cand[MAXC][MAXN];
    int okc = 0, amb = 0, fail = 0, nonin = 0, badg = 0;
    for (int tr = 0; tr < trials; tr++)
    {
        poly g = random_squarefree_g(t);
        if (!build_frobenius_inverse(g))
        {
            badg++;
            continue;
        }
        setup_code(g);
        int e[MAXN] = {0}, cnt = 0, a = 1 + rand() % (P - 1);
        if (mode == 0)
            a = 1;
        while (cnt < w)
        {
            int i = rand() % n;
            if (e[i])
                continue;
            e[i] = (mode == 2) ? 1 + rand() % (P - 1) : a;
            cnt++;
        }
        poly s = syndrome(e);
        int lo = (mode == 0) ? 1 : 1, hi = (mode == 0) ? 1 : P - 1;
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
    printf("P=%d t=%d w=%d mode=%d sign=%d: 成功=%d 曖昧(候補複数)=%d 失敗=%d s不可逆=%d g不適=%d (試行 %d)\n",
           P, t, w, mode, SIGN, okc, amb, fail, nonin, badg, trials);
    return 0;
}
