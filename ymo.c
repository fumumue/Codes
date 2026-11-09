#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include "global-p.h"
#include "struct.h"
#include "chash-p.c"

#define SEPARABLE 0

// extern int mlt(int x, int y);
// extern int mltn(int n, int x);

unsigned short g[K + 1] = {0};

// ランダム多項式の生成
static void
ginit(void)
{
    int j, count = 0, k = 0;
    unsigned short gg[K + 1] = {0};

    printf("in ginit\n");

    g[K] = 1;          // xor128();
    g[0] = rand() % N; // or N
    k = rand() % (K - 1);
    if (k > 0)
    {
        while (count < k)
        {
            printf("in whule\n");
            j = rand() % (K);
            if (j < K && j > 0 && g[j] == 0)
            {
                g[j] = rand() % N; // or N;
                count++;
            }
        }
    }

    for (j = 0; j < K + 1; j++)
        gg[j] = g[K - j];

    memcpy(g, gg, sizeof(g));
}

// OP型からベクトル型への変換
vec o2v(OP f)
{
    vec a = {0};
    int i;

    for (i = 0; i < K * E; i++)
    {
        if (f.t[i].a > 0 && f.t[i].n < K * E)
            a.x[f.t[i].n] = f.t[i].a;
    }

    return a;
}

// ベクトル型からOP型への変換
OP v2o(vec a)
{
    int i, j = 0;
    OP f = {0};

    // #pragma omp parallel for
    for (i = 0; i < K * E; i++)
    {
        if (a.x[i] > 0)
        {
            f.t[j].n = i;
            f.t[j++].a = a.x[i];
        }
    }

    return f;
}

unsigned short oinv(unsigned short a, unsigned short n)
{
    unsigned short i;

    if (a == 0)
        return 0;
    if (a == 1)
        return 1;
    for (i = 2; i < n; i++)
    {
        if ((i * a) % n == 1)
            return i;
    }
    printf("no return\n");
    exit(1);
}

// 多項式の次数(default)
int deg(vec a)
{
    int i, n = 0, flg = 0;

    // #pragma omp parallel for
    for (i = 0; i < DEG; i++)
    {
        if (a.x[i] > 0)
        {
            n = i;
            flg = 1;
        }
    }
    if (flg == 0)
        return 0;

    return n;
}

// 多項式を表示する(default)
void printpol(vec a)
{
    int i, n;

    n = deg(a);

    // printf ("baka\n");
    //  assert(("baka\n", n >= 0));

    for (i = n; i > -1; i--)
    {
        if (a.x[i] > 0)
        {
            printf("%u*", a.x[i]);
            // if (i > 0)
            printf("x^%d", i);
            // if (i > 0)
            printf("+");
        }
    }
    //  printf("\n");

    return;
}

vec kof2(unsigned short c, vec f)
{
    int i, k;
    vec b = {0}, h = {0};

    printf("c=%d\n", c);
    // exit(1);
    b = f; // o2v(f);
    k = deg(b);
    printpol(b);
    printf(" =b debugi\n");
    for (i = 0; i < k + 1; i++)
    {
        h.x[i] = (c * b.x[i]) % N;
    }
    // g = v2o(h);
    printpol(h);
    printf(" =h in kof2\n");
    return h;
}

vec vadd(vec a, vec b)
{
    int i;
    vec c = {0};

    // printf("deg=%d %d\n",deg(a),deg(b));

    for (i = 0; i < DEG; i++)
        c.x[i] = (a.x[i] + b.x[i]) % N;

    return c;
}

int mul = 0, mul2 = 0;
vec vmul(vec a, vec b)
{
    int i, j, k, l;
    vec c = {0};
    if (deg(a) > 128 && deg(b) > 128)
        mul++;
    mul2++;

    k = deg(a);
    l = deg(b);

    for (i = 0; i < k + 1; i++)
    {
        for (j = 0; j < l + 1; j++)
        // if (a.x[i] > 0)
        {
            c.x[i + j] += (a.x[i] * b.x[j]) % N;
            c.x[i + j] %= N;
        }
    }

    return c;
}

unsigned short vb[K * 2][N] = {0};
unsigned short gt[K * 2][K * 2] = {0};

void van(int kk)
{
    int i, j;

    printf("van der\n");
    
    for (i = 0; i < N; i++){
        mat[i][0] = vb[0][i] = 1;
        printf("%d,", vb[0][i]);
    }
    printf("\n");
    
    // #pragma omp parallel for private(i, j)
    for (i = 1; i < kk; i++)
    {
        for (j = 0; j < N; j++)
        {
            vb[i][j] = mltn(i , j);
            printf("g%d,", vb[i][j]);
            mat[j][i] = vb[i][j];
        }
        printf("\n");
    }
}

void ogt(unsigned short pp[], int kk)
{
    int i, j;

    // #pragma omp parallel for private(i, j)
    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < kk - i; j++)
        {
            gt[i][j + i] = g[j];
        }
    }
    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < kk; j++)
            printf("h%d,", gt[i][j]);
        printf("\n");
    }
    // exit(1);
}

// 配列の値を係数として多項式に設定する
OP setpol(unsigned short f[], int n)
{
    OP g;
    vec v = {0};
    int i;

    for (i = 0; i < n; i++)
    {
        v.x[n - 1 - i] = f[i];
    }

    g = v2o(v);

    return g;
}

OP mkpol()
{
    int i, j, k, flg, ii = 0;
    OP w = {0};

    do
    {
        // fail = 0;
        j = 0;
        k = 0;
        flg = 0;
        // l = 0;
        memset(g, 0, sizeof(g));
        // memset(ta, 0, sizeof(ta));
        memset(w.t, 0, sizeof(w));
        ginit();
        ii++;
        if (ii > 100)
        {
            printf("erro=%d\n", ii);
            exit(1);
        }

        for (i = 0; i < K; i++)
        {
            if (g[K - 1] > 0)
                flg = 1;
            if (i % 2 == 1 && g[i] > 0 && i < K)
                k++;
        }

        // 偶数項だけにならないようにする
        if ((k > 0 && flg == 0) || (k > 1 && flg == 1))
        // if(k>0)
        {
            w = setpol(g, K + 1);
            j = 1;
            // if(isquad(w)==-1)
            // exit(1);
        }
        // exit(1);

    } while (j == 0);

    printpol(o2v(w));
    printf(" ==g\n");
    // exit(1);

    return w;
}

unsigned short
v2a(oterm a)
{
    int j;

    if (a.a == 0)
        return 0;

    // printf("aa=%d\n",a.a);
    for (j = 0; j < M; j++)
    {
        if (j == a.a && a.a > 0)
        {
            // printf("j==%d\n",j);
            return j - 1;
        }
    }
    return 0;
}

void printsage(vec a)
{
    int i, j;
    oterm b;

    printf("poly=");
    for (i = 0; i < DEG; i++)
    {
        if (a.x[i] > 0)
        {
            b.a = a.x[i];
            b.n = i;
            j = v2a(b);
            printf("B('a^%d')*X**%d+", j, i); // for GF(2^m)
        }
    }
}

// 多項式の代入値
unsigned short
trace(OP f, unsigned short x)
{
    unsigned short u = 0;
    vec v = o2v(f);
    int d = deg((v)) + 1;

    for (int i = 0; i < d; i++)
    {
        if (v.x[i] > 0)
            u = (u + (v.x[i] * mltn(i, x))) % N;
    }

    return u;
}

// リーディグタームを抽出(default)
oterm vLT(vec f)
{
    int i;
    oterm t = {0};

    // k = deg (o2v (f));
    for (i = 0; i < DEG; i++)
    {
        // printf("a=%d %d\n",f.t[i].a,f.t[i].n);
        if (f.x[i] > 0)
        {
            t.n = i;
            t.a = f.x[i];
        }
    }

    return t;
}

// aに何をかけたらbになるか
unsigned short
equ(unsigned short a, unsigned short b)
{
    return (oinv(a, N) * b);
}

// 多項式を単行式で割る
oterm vLTdiv(vec f, oterm t)
{
    oterm tt = {0}, s = {
                        0};

    tt = vLT(f);
    if (tt.n < t.n)
    {
        s.n = 0;
        s.a = 0;
    }
    else if (tt.n == t.n)
    {
        s.n = 0;
        s.a = equ(t.a, tt.a);
    }
    else if (tt.n > t.n)
    {
        s.n = tt.n - t.n;
        s.a = equ(t.a, tt.a);
        // printf("%u\n",s.a);
    }
    else if (t.n == 0 && t.a > 0)
    {
        s.a = (tt.a * oinv(t.a, N)) % N;
        s.n = tt.n;
    }

    return s;
}

// 多項式を項ずつ掛ける
vec vterml(vec f, oterm t)
{
    // f = conv(f);
    // ssert(op_verify(f));
    int i;
    vec h = {0};

    // f=conv(f);
    // k = deg (o2v(f));

    for (i = 0; i < DEG; i++)
    {
        // h.t[i].n = f.t[i].n + t.n;
        if (f.x[i] > 0)
            h.x[i + t.n] = (f.x[i] * t.a) % N;
    }

    // h = conv(h);
    //  assert(op_verify(h));
    return h;
}

int vm = 0;
// 多項式の剰余を取る
vec vmod(vec f, vec g)
{
    vec h = {0};
    oterm b = {0}, c = {0};

    vm++;
    // printf("vmod-bl=%d k=%d\n",deg(f),deg(g));
    if (vLT(f).n < vLT(g).n)
    {
        //    exit(1);
        return f;
    }

    b = vLT(g);

    // printpol(f);
    // printf(" ==f\n");
    while (1)
    {

        c = vLTdiv(f, b);
        h = vterml(g, c);
        f = vadd(f, h);
        if (deg((f)) == 0 || deg((g)) == 0)
        {
            break;
        }

        if (c.n == 0)
            break;
    }
    // printf("vmod-baka== %d %d\n",deg(f),deg(g));
    return f;
}

// int mul = 0, mul2 = 0;
vec vmul_2(vec a, vec b)
{
    int i, j, k, l;
    vec c = {0};
    if (deg(a) > 128 && deg(b) > 128)
        mul++;
    mul2++;

    k = deg(a);
    l = deg(b);

    for (i = 0; i < k + 1; i++)
    {
        for (j = 0; j < l + 1; j++)
        // if (a.x[i] > 0)
        {
            c.x[i + j] += (a.x[i] * b.x[j]) % N;
            c.x[i + j] %= N;
        }
    }

    return c;
}

int cnty = 0;
vec vpp(vec f, vec mod)
{
    int i;
    vec s = {0};
    // t = f;
    s = f;

    // 繰り返し２乗法
    for (i = 1; i < E + 1; i++)
    {
        s = vmod(vmul_2(s, s), mod);
    }

    return s;
}

// gcd
vec vgcd(vec xx, vec yy)
{
    vec tt = {0}, tmp, h = {0};
    // ee.x[K] = 1;

    h.x[0] = 1;
    // h.x[0] = 0;
    if (deg((xx)) < deg((yy)))
    {
        tmp = xx;
        xx = yy;
        yy = tmp;
    }
    // tt = vmod(xx, yy);
    tt = vmod(xx, yy);
    while (deg(tt) > 0)
    {
        xx = yy;
        yy = tt;
        if (deg(yy) > 0)
            tt = vmod(xx, yy);
        if (vLT(tt).a == 0)
            return yy;
    }
    if (vLT(tt).a == 0)
    {
        return yy;
    }
    else
    {
        return h;
    }
    //  return yy;
}

// #define NN 16
vec renritu(MTX a)
{
    unsigned short p, d;
    int i, j, k;
    vec v = {0};

    for (i = 0; i < K; i++)
    {
        p = a.x[i][i];

        for (j = 0; j < (K + 1); j++)
        {
            a.x[i][j] = (a.x[i][j] * oinv(p, N)) % N;
        }

        for (j = 0; j < K; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (k = i; k < (K + 1); k++)
                {
                    a.x[j][k] = (a.x[j][k] + d * a.x[i][k]) % N;
                }
            }
        }
    }
    for (i = 0; i < K; i++)
    {
        if (a.x[i][i] != 1)
        {
            for (j = 0; j < K + 1; j++)
                printf("a%d,", a.x[i][j]);
            printf("\n");
            exit(1);
        }
    }
    printf("\n");
    vec x = {0};
    for (i = 0; i < K; i++)
    {
        v.x[i] = a.x[i][K];
        // v.x[128]=1;
        printf(" x%d = %d\n", i, v.x[i]);
        x.x[i + 1] = v.x[i];
    }
    x.x[0] = 1;

    OP pol = {0};
    pol = setpol(x.x, K + 1);
    printpol(o2v(pol));
    printf(" ==key\n");
    for (i = 0; i < N; i++)
    {
        if (trace(pol, i) % N == 0)
            printf("uz%d i=%d\n", i, i);
    }

    return v;
}

// #define NN 16
vec sol(MTX a)
{
    unsigned int p, d;
    int i, j, k;
    vec v = {0};

    for (i = 0; i < K / 2; i++)
    {
        p = a.x[i][i];

        for (j = 0; j < (K / 2 + 1); j++)
        {
            a.x[i][j] = (a.x[i][j] * oinv(p, N)) % N;
        }

        for (j = 0; j < K / 2; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (k = i; k < (K / 2 + 1); k++)
                {
                    if (a.x[j][k] > (d * a.x[i][k]) % N)
                    {
                        a.x[j][k] -= (d * a.x[i][k]) % N;
                    }
                    else
                    {
                        a.x[j][k] = (N + (a.x[j][k] - (d * a.x[i][k]) % N)) % N;
                    }
                }
            }
        }
    }
    for (i = 0; i < K / 2; i++)
    {
        for (j = 0; j < K / 2 + 1; j++)
            printf("%d,", a.x[i][j]);
        printf("\n");
        if (a.x[i][i] != 1)
        {
            for (j = 0; j < K / 2 + 1; j++)
                printf("a%d,", a.x[i][j]);
            printf("\n");
        }
    }
    printf("\n");
    vec x = {0};
    for (i = 0; i < K / 2; i++)
    {
        if (N > a.x[i][K / 2])
        {
            x.x[K / 2 - i] = (N - a.x[i][K / 2]) % N;
        }
        else
        {
            x.x[K / 2 - i] = a.x[i][K / 2] % N;
        }
        // x.x[i+1] = N-a.x[i][K / 2];
    }

    x.x[0] = 1;

    vec vv = {0};
    OP pol = {0};
    pol = setpol(x.x, K / 2 + 1);
    printpol(o2v(pol));
    printf(" ==key\n");
    for (i = 0; i < N; i++)
    {
        // v.x[i] = 0;
        if (trace(pol, i) % N == 0)
        {
            printf("i=%d\n", i);
            vv.x[i] = 1;
        }
    }

    return vv;
}

void vv(int kk)
{
    int i, j;
    OP r = mkpol();
    unsigned short tr[N];
    unsigned short ta[N] = {0};

    printf("van der\n");
    van(K);

    int l = -1;
    vec pp = {0}, tt = {0};

aa:
    // exit(1);
    r = mkpol();

    for (i = 0; i < N; i++)
    {
        ta[i] = trace(r, i);
        if (ta[i] == 0)
        {
            printf("trace 0 @ %d\n", i);
            // fail = 1;
            goto aa;
        }
    }

    for (i = 0; i < N; i++)
    {
        tr[i] = oinv(ta[i], N);
        // printf("%d,", tr[i]);
    }

    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < kk; j++)
        {
            mat[i][j] = (vb[j][i] * tr[i]) % N;
        }
    }
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < N; j++)
            printf("c%d,", mat[j][i]);
        printf("\n");
    }
}

void mkerr(unsigned short *z1, int num)
{
    int j, l;

    j = 0;

    memset(z1, 0, sizeof(2 * N));

    while (j < num)
    {
        l = rand() % (N - 1);
        // printf ("l=%d\n", l);
        if (0 == z1[l] && l > 0)
        {
            z1[l] = 1;
            // printf("l=%d\n", l);
            j++;
        }
    }
}

OP synd(unsigned short zz[], int kk)
{
    unsigned short syn[K] = {0}, s = 0;
    int i, j;
    OP f = {0};

    printf("in synd2\n");

    for (i = 0; i < kk; i++)
    {
        syn[i] = 0;
        s = 0;
        // #pragma omp parallel num_threads(16)
        for (j = 0; j < N; j++)
        {
            s = (s + (zz[j] * mat[j][i])) % N;
        }
        syn[i] = s;
        // printf ("syn%d,", syn[i]);
    }
    // printf ("\n");

    f = setpol(syn, kk);
    printpol(o2v(f));
    printf(" syn=============\n");
    //  exit(1);

    return f;
}

// chen探索
vec chen(vec f)
{
    vec e = {0};
    int i, n, x = 0, count = 0;
    unsigned short z;

    n = deg((f));
    for (x = 0; x < N; x++)
    {
        z = 0;
        for (i = 0; i < n + 1; i++)
        {
            if (f.x[i] > 0)
                z += (mltn(i, x) * f.x[i]) % N;
        }
        if (z % N == 0)
        {
            e.x[count] = x;
            count++;
            printf("change %d\n", (x));
        }
    }

    return e;
}

int main()
{
    int i;
    unsigned short s[K + 1] = {0}, z1[N] = {0};
    vec v = {0}, x = {0};
    OP f = {0};

    printf("%d %d %d\n", 3, oinv(3, N), 3 * oinv(3, N) % N);
    // exit(1);
    srand(clock());
    // mkg(K); // Goppa Code (EEA type)
    //van(K); // RS-Code generate
     vv(K);           // Goppa Code's Parity Check (Berlekamp type)
    for(i=0;i<T;i++)
    z1[i]=1;
    //mkerr(z1, T);    // generate error vector
    f = synd(z1, K); // calc syndrome
    x = o2v(f);      // transorm to vec
    // r = bma(x.x);    // Berlekamp-Massey Algorithm
    MTX b = {0};

    for (i = 0; i < K; i++)
        v.x[K - 1 - i] = x.x[i];

    for (i = 0; i < K / 2; i++)
    {
        for (int j = 0; j < K / 2 + 1; j++)
        {
            b.x[i][j] = v.x[i + j];
            // printf("%d,",b.x[i][i+j]);
        }
        // printf("\n");
    }
    printf("\n");
    for (i = 0; i < K / 2; i++)
    {
        for (int j = 0; j < K / 2 + 1; j++)
            printf("e%d,", b.x[i][j]);
        printf("\n");
    }
    x = sol(b);

    for (i = 0; i < N; i++)
    {
        if (z1[i] != x.x[i]){
            printf("baka=%d\n", i);
            exit(1);
        }
    }
    // printf("\n");

    return 0;
}

