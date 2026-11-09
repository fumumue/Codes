#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// #include "8192.h"
#include "global.h"
// #include "param.h"
#include "struct.h"
#include "chash.c"

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

// 有限体の元の逆数
unsigned short
oinv(unsigned short a)
{

    if (a == 0)
        return 0;
    if (a == 1)
        return a;
    return N - fg[a] + 1;

    printf("no return \n");

    //  exit(1);
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
            printf("%u", a.x[i]);
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

    c = fg[c];
    printf("c=%d\n", c);
    // exit(1);
    b = f; // o2v(f);
    k = deg(b);
    printpol(b);
    printf(" =b debugi\n");
    for (i = 0; i < k + 1; i++)
    {
        h.x[i] = gf[mlt(c, fg[b.x[i]])];
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
        c.x[i] = a.x[i] ^ b.x[i];

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
            c.x[i + j] ^= gf[mlt(fg[a.x[i]], fg[b.x[j]])];
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
    /*
    for (i = 0; i < N; i++){
        mat[i][0] = vb[0][i] = 1;
        printf("%d,", vb[0][i]);
    }
    printf("\n");
    */
    // #pragma omp parallel for private(i, j)
    for (i = 0; i < kk; i++)
    {
        for (j = 1; j < N; j++)
        {
            vb[i][j] = gf[mltn(i , fg[j])];
            printf("%d,", vb[i][j]);
            mat[j ][i] = vb[i][j];
        }
        printf("\n");
    }
    //exit(1);
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
            printf("%d,", gt[i][j]);
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
        if (gf[j] == a.a && a.a > 0)
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
    int i, d;
    unsigned short u = 0;

    d = deg(o2v(f));

    for (i = 0; i < d + 1; i++)
    {
        u ^= gf[mlt(fg[f.t[i].a], mltn(f.t[i].n, fg[x]))];
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
    return gf[mlt(oinv(a), fg[b])];
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
        s.a = gf[mlt(fg[tt.a], oinv(t.a))];
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
            h.x[i + t.n] = gf[mlt(fg[f.x[i]], fg[t.a])];
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
            c.x[i + j] ^= gf[mlt(fg[a.x[i]], fg[b.x[j]])];
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

// GF(2^m) then set m in this function.
int ben_or(vec f)
{
    int i, n; //, pid;
    vec s = {0}, u = {0}, r = {0};
    vec v = {0}; //, ff=o2v(f);
    // if GF(8192) is 2^m and m==13 or if GF(4096) and m==12 if GF(16384) is testing
    // int m = E;
    //  m=12 as a for GF(4096)=2^12 defined @ gloal.h or here,for example m=4 and GF(16)

    v.x[1] = 1;
    s = (v);
    // for (i = 0; i < K / 2; i++)
    r = s;
    n = deg((f));

    if (vLT(f).n == 0 && vLT(f).a == 1)
    {
        printf("f==0\n");
        exit(1);
    }
    if (n == 0)
        return -1;

    i = 0;

    // r(x)^{q^i} square pow mod
    for (i = 0; i < K / 2; i++)
    {
        printf(":i=%d", i);
        // irreducible over GH(8192) 2^13
        r = vpp(r, f);
        // if(r.x[0]==65535)
        // return -1;
        u = vadd(r, (s));
        u = vgcd(f, u);

        if (deg(u) > 0 || vLT(u).a == 0)
        {
            // flg[i]= -1;
            printf("ae\n");
            return -1;
        }
    }

    return 0;
}

unsigned short gf_mul(unsigned short in0, unsigned short in1)
{
    int i;

    uint32_t tmp;
    uint32_t t0;
    uint32_t t1;
    uint32_t t;

    t0 = in0;
    t1 = in1;

    tmp = t0 * (t1 & 1);

    for (i = 1; i < 12; i++)
        tmp ^= (t0 * (t1 & (1 << i)));

    t = tmp & 0x7FC000;
    tmp ^= t >> 9;
    tmp ^= t >> 12;

    t = tmp & 0x3000;
    tmp ^= t >> 9;
    tmp ^= t >> 12;

    return tmp & ((1 << 12) - 1);
}

/* input: in0, in1 in GF((2^m)^t)*/
/* output: out = in0*in1 */
void GF_mul(unsigned short *out, unsigned short *in0, unsigned short *in1)
{
    int i, j;

    unsigned short prod[K * 2 - 1] = {0};

    for (i = 0; i < K * 2 - 1; i++)
        prod[i] = 0;

    for (i = 0; i < K; i++)
    {
        for (j = 0; j < K; j++)
            prod[i + j] ^= gf[mlt(fg[in0[i]], fg[in1[j]])];
    }
    //

    for (i = (K - 1) * 2; i >= K; i--)
    {
        if (K == 512)
        {
            // GF(2^512) from sage
            prod[i - K + 8] ^= prod[i];
            prod[i - K + 5] ^= prod[i];
            prod[i - K + 2] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 256)
        {
            // GF(2^256) from sage
            prod[i - K + 10] ^= prod[i];
            prod[i - K + 5] ^= prod[i];
            prod[i - K + 2] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 128)
        {
            // 128
            prod[i - K + 7] ^= prod[i];
            prod[i - K + 2] ^= prod[i];
            prod[i - K + 1] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 32)
        {
            // 32
            prod[i - K + 15] ^= prod[i];
            prod[i - K + 9] ^= prod[i];
            prod[i - K + 7] ^= prod[i];
            prod[i - K + 4] ^= prod[i];
            prod[i - K + 3] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 16)
        {
            // 16
            prod[i - K + 5] ^= prod[i];
            prod[i - K + 3] ^= prod[i];
            prod[i - K + 2] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 8)
        {
            // 8
            prod[i - K + 4] ^= prod[i];
            prod[i - K + 3] ^= prod[i];
            prod[i - K + 2] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 6)
        {
            // 8
            prod[i - K + 4] ^= prod[i];
            prod[i - K + 3] ^= prod[i];
            prod[i - K + 1] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
        if (K == 4)
        {
            // 4
            prod[i - K + 1] ^= prod[i];
            prod[i - K + 0] ^= prod[i];
        }
    }

    for (i = 0; i < K; i++)
        out[i] = prod[i];
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
            a.x[i][j] = gf[mlt(fg[a.x[i][j]], oinv(p))];
        }

        for (j = 0; j < K; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (k = i; k < (K + 1); k++)
                {
                    a.x[j][k] = a.x[j][k] ^ gf[mlt(fg[d], fg[a.x[i][k]])];
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
        if (trace(pol, i) == 0)
            printf("%d i=%d\n", i, fg[i] - 1);
    }

    return v;
}

// #define NN 16
vec renritz(MTX a)
{
    unsigned short p, d;
    int i, j, k;
    vec v = {0};
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < K / 2 + 1; j++)
            printf("%d,", a.x[i][j]);
        printf("\n");
    }
    // exit(1);
    for (i = 0; i < K / 2; i++)
    {
        p = a.x[i][i];

        for (j = 0; j < (K / 2 + 1); j++)
        {
            a.x[i][j] = gf[mlt(fg[a.x[i][j]], oinv(p))];
        }

        for (j = 0; j < K / 2; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (k = i; k < (K / 2 + 1); k++)
                {
                    a.x[j][k] = a.x[j][k] ^ gf[mlt(fg[d], fg[a.x[i][k]])];
                }
            }
        }
    }
    for (i = 0; i < K / 2; i++)
    {
        if (a.x[i][i] != 1)
            exit(1);
        for (j = 0; j < K / 2 + 1; j++)
            printf("%d,", a.x[i][j]);
        printf("\n");
    }
    printf("\n");
    // exit(1);

    for (i = 0; i < K / 2; i++)
    {
        v.x[i] = a.x[i][K / 2];
        // v.x[128]=1;
        printf(" x%d = %d\n", i, fg[v.x[i]]);
    }

    return v;
}

// #define NN 16
vec sol(MTX a)
{
    unsigned short p, d;
    int i, j, k;
    vec v = {0};

    for (i = 0; i < K / 2; i++)
    {
        p = a.x[i][i];

        for (j = 0; j < (K / 2 + 1); j++)
        {
            a.x[i][j] = gf[mlt(fg[a.x[i][j]], oinv(p))];
        }

        for (j = 0; j < K / 2; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (k = i; k < (K / 2 + 1); k++)
                {
                    a.x[j][k] = a.x[j][k] ^ gf[mlt(fg[d], fg[a.x[i][k]])];
                }
            }
        }
    }
    for (i = 0; i < K / 2; i++)
    {
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
        x.x[K / 2 - i] = a.x[i][K / 2];
        // printf(" x%d = %d\n", i, v.x[i]);
    }

    x.x[0] = 1;
    int count = 0;
    OP pol = {0};
    pol = setpol(x.x, K / 2 + 1);
    printpol(o2v(pol));
    printf(" ==key\n");
    for (i = 0; i < N; i++)
    {
        if (trace(pol, i) == 0)
        {
            printf("%d i=%d\n", i, gf[i-1]);
            v.x[count++] = i;
        }
    }
    for(i=0;i<count;i++)
    printf("%d,",v.x[i]);
    printf("\n");
    return v;
}

/* input: f, element in GF((2^m)^t) */
/* output: out, minimal polynomial of f */
/* return: 0 for success and -1 for failure */
int mykey(unsigned short *out, vec x)
{
    unsigned short mat[K + 1][K] = {0};
    MTX a = {0};
    int i, j;

    // fill matrix

    mat[0][0] = 1;

    for (i = 1; i < K; i++)
        mat[0][i] = 0;

    for (i = 0; i < K; i++)
        mat[1][i] = x.x[i];

    for (j = 2; j <= K; j++)
    {
        GF_mul(mat[j], mat[j - 1], x.x);
    }
    // exit(1);
    //
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < K + 1; j++)
        {
            a.x[i][j] = mat[j][i];
            printf("%d,", mat[j][i]);
        }
        printf("\n");
    }
    printf("\n");
    // exit(1);

    vec v = {0};
    v = renritu(a);

    for (i = 0; i < K; i++)
    {
        out[i] = v.x[i];
        printf("%d,", out[i]);
    }
    printf("\n");
    return 0;
}

// 多項式の代入値
unsigned short eval(OP f, unsigned short x)
{
    OP g = {0};
    unsigned short u = 0, s = 0;
    vec v = o2v(f), h = {0};
    int d = deg((v)) + 1;

    for (int i = 0; i < d; i++)
    {
        if (v.x[i] > 0)
        {
            u ^= gf[mlt(fg[v.x[i]], mltn(i, fg[x]))];
            //
            // vtrace((oadd(dick[fg[u]], dick[mlt(fg[v.x[i]], mltn(i, x))])),x);
        }
    }

    return u;
}


OP vv(int kk)
{
    int i, j;
    OP r = mkpol();
    unsigned short tr[N];
    unsigned short ta[N] = {0};

    printf("van der\n");

    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < N; j++)
        {
            vb[i][j] = gf[mltn(i, j)];
        }
        // printf("\n");
    }

    int l = -1;
    vec pp = {0}, tt = {0};

aa:
    if (SEPARABLE == 0)
    {
        while (l < 0)
        {
            for (i = 0; i < K; i++)
                pp.x[i] = rand() % N;
            mykey(tt.x, pp);
            tt.x[K] = 1;
            l = ben_or(tt);
            if (l == 0)
            {
                printf("\n");
                printsage(tt);
                printf(" ==irr\n");
                // exit(1);
            }
        }
        r = v2o(tt);
    }
    // exit(1);
    if (SEPARABLE == 1)
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
        tr[i] = oinv(ta[i]);
        // printf("%d,", tr[i]);
    }

    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < kk; j++)
        {
            mat[i][j] = gf[mlt(fg[vb[j][i]], tr[i])];
        }
    }
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < N; j++)
            printf("%d,", mat[j][i]);
        printf("\n");
    }

return r;
}

// Patterson & EEA 用（ランダム多項式、次元指定）
OP mkg(int kk)
{
    int i, j, l, ii = 0;
    OP w = {0};
    unsigned short tr[N] = {0};
    unsigned short ta[N] = {0};
    unsigned short po[K + 1] = {0}; //{1, 0, 1, 0, 5};

aa:

    // printf("\n");

    // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
    // 既約多項式しか使わない。

    l = -1;
    ii = 0;

    // while (l == -1)
    {
        w = mkpol();
        // l = ben_or(o2v(w));
        printf("irr=%d\n", l);
        if (ii > 300)
        {
            printf("too many error\n");
            exit(1);
        }
        ii++;
        //
    }

    // w = mkpol();
    memset(ta, 0, sizeof(ta));
    // g[0]=1;

    // 多項式の値が0でないことを確認
    for (i = 0; i < N; i++)
    {
        ta[i] = trace(w, i);
        if (ta[i] == 0)
        {
            printf("trace 0 @ %d\n", i);
            // fail = 1;
            goto aa;
        }
    }

    // 多項式を固定したい場合コメントアウトする。
    /*
  memset(ta, 0, sizeof(ta));
  w = setpol(po, K + 1);
  printpol(o2v(w));
  printf(" =poly\n");
//  exit(1);
  for (i = 0; i < N; i++)
  {    ta[i] = trace(w, i);
    if (ta[i] == 0)
    {
      printf("trace 0 @ %d\n", i);
      exit(1);
    }
  }
  */
    printpol(o2v(w));
    printf("\n");
    printsage(o2v(w));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");

    for (i = 0; i < N; i++)
    {
        tr[i] = oinv(ta[i]);
        // printf("%d,", tr[i]);
    }

    printpol(o2v(w));
    printf(" =irreducible\n");
    printsage(o2v(w));
    printf("\n");
    // wait();

    memset(vb, 0, sizeof(vb));
    memset(gt, 0, sizeof(gt));
    van(kk);
    ogt(po, kk);
    memset(mat, 0, sizeof(mat));

    printf("\nすげ、オレもうイキそ・・・\n");

    for (j = 0; j < N; j++)
    {
        for (i = 0; i < K; i++)
        {
            ma[j][i] = gf[mlt(fg[vb[i][j]], tr[j])];
        }
    }

    unsigned short s;
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < N; j++)
        {
            s = 0;

            for (int k = 0; k < K; k++)
                s ^= gf[mlt(fg[gt[k][i]], fg[ma[j][k]])];
            // printf("%d,",s);
            mat[j][i] = s;
        }
        printf("\n");
    }
    printf("\n");
    // exit(1);
    /*
        for (j = 0; j < N; j++)
        {
            for (i = 0; i < K; i++)
                printf("%d,", mat[j][i]);
            printf("\n");
        }
        printf("\n");
        // wait();
    */
    return w;
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
            s ^= gf[mlt(fg[zz[j]], fg[mat[j][i]])];
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
                z ^= gf[mlt(mltn(i, fg[x]), fg[f.x[i]])];
        }
        if (z == 0)
        {
            e.x[count] = x;
            count++;
            printf("change %d\n", (fg[x]));
        }
    }

    return e;
}

// Input:符号の次元をKとすると、K個のシンドロームを要素として持つ配列ｓ
// Output:誤り位置決定多項式（この多項式が０になる値を探すことでエラーの位置を求める）
vec bma(unsigned short s[]) // sはシンドロームの値
{
    int L = 0, m = -1, d[K] = {0}, k = 0, i, e;
    vec f = {0}, g = {0}, h, v;

    f.x[0] = g.x[0] = 1;

    while (k <= (2 * T - 1))
    {
        e = 0;
        for (i = 0; i < L; i++)
            e ^= gf[mlt(fg[f.x[i]], fg[s[k - i]])];

        d[k] = gf[mlt(fg[f.x[i]], fg[s[k - i]])] ^ e;
        if (d[k] > 0)
        {
            h = f;
            memset(v.x, 0, sizeof(v.x));
            v.x[k - m] = 1;

            unsigned short a;
            a = (m < 0) ? 1 : oinv(d[m]);
            f = vadd(f, vmul(kof2(gf[mlt(fg[d[k]], a)], g), v));
            if (L <= k / 2)
            {
                L = k + 1 - L;
                m = k;
                g = h;
            }
        }
        k++;
    }

    return f;
}

vec bms(unsigned short s[], int kk)
{
    int i, j, k, ll = 0, l, d[2 * K + 1] = {0};
    vec lo[2 * K + 1] = {0}, b[2 * K + 1] = {0}, t[2 * K + 1] = {0}, h = {0}, g = {0};
    vec v = {0}, x = {0};

    x.x[1] = 1;
    h = (x);
    v.x[0] = 1;
    // f = v2o(x);
    lo[0] = (v);
    b[0] = lo[0];
    ll = 0;
    for (j = 1; j < T * 2 + 1; j++)
    {
        v = (lo[j - 1]);
        k = 0;

        l = deg((lo[j - 1]));
        for (i = 1; i < l + 1; i++)
        {
            k ^= gf[mlt(fg[v.x[i]], fg[s[j - i]])];
            printf("v[%d]=%d %d %d\n", i, v.x[i], s[j - i], k);
        }
        d[j] = s[j] ^ k;
        printf("d[%d]=%d\n", j, d[j]);
        if (d[j] == 0)
        {
            lo[j] = lo[j - 1];
            b[j] = (vmul((b[j - 1]), (h)));
            // ll=j-1;
        }
        else // if (d[j] > 0)
        {
            g = (vmul(kof2(d[j], (h)), (b[j - 1])));
            t[j] = (vadd((lo[j - 1]), (g)));
            lo[j] = t[j];
            if (ll * 2 > (j - 1))
            {
                // lo[j]=t[j];
                b[j] = (vmul((b[j - 1]), (h)));
            }
            else // if(2*ll <= j)
            {
                // printpol(o2v(t[j]));
                // printf("==t[%d]\n", j);
                b[j] = (kof2(gf[oinv(d[j])], (lo[j - 1])));
                // lo[j]=t[j];
                ll = j - ll;

                if (j == 2 * T)
                {
                    lo[j - 1] = (vmul((lo[j - 2]), (h)));
                    break;
                }
            }
        }
        printf("l=%d\n", ll);
        k = 0;
        // printpol(o2v(b[j]));
        // printf(" ==b[%d]\n", j);
    }

    k = 0;
    // printpol(o2v(lo[j - 1]));
    // printf(" ==coef\n");
    if (deg((lo[j - 1])) == T)
    {
        chen((lo[j - 1]));
    }
    else
    {
        printf("baka==\n");
        exit(1);
    }

    return (lo[j - 1]);
}

// inverse matrix
MTX matinv(MTX a, int n)
{

    // unsigned short a[F][F];     //={{1,2,0,1},{1,1,2,0},{2,0,1,1},{1,2,1,1}}; //入力用の配列
    unsigned short inv_a[N][N]; // ここに逆行列が入る
    unsigned short buf;         // 一時的なデータを蓄える
    unsigned short b[N][N] = {0}, dd[N][N] = {0};
    int i, j, k, count; // カウンタ
    // MTX a={0};
    unsigned short c[N][N] = {0};
    MTX z = {0};
    unsigned short cc[N][N] = {0};

lab:
    memset(b, 0, sizeof(b));
    // memset(a.x, 0, sizeof(a.x));
    /*
      for (i = 0; i < n; i++)
      {
        for (j = 0; j < n; j++)
        {
          a.x[i][j] = rand() % 256;
          printf("%d,", a.x[i][j]);
        }
        printf("\n");
      }
      */
    // exit(1);
    //  printf("\n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            c[i][j] = a.x[i][j];
    }
    // 単位行列を作る
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            inv_a[i][j] = (i == j) ? 1 : 0;
        }
    }
    // 掃き出し法
    // #pragma omp parallel for num_threads(omp_get_max_threads()) //private(i,j,k)
    for (i = 0; i < n; i++)
    {
        buf = gf[oinv(a.x[i][i])];
        for (j = 0; j < n; j++)
        {
            a.x[i][j] = gf[mlt(fg[buf], fg[a.x[i][j]])];
            inv_a[i][j] = gf[mlt(fg[buf], fg[inv_a[i][j]])];
        }
        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                buf = a.x[j][i];
                for (k = 0; k < n; k++)
                {
                    a.x[j][k] ^= gf[mlt(fg[a.x[i][k]], fg[buf])];
                    inv_a[j][k] ^= gf[mlt(fg[inv_a[i][k]], fg[buf])];
                }
            }
        }
    }

    // printf("\n\n逆行列を出力\n");
    for (i = 0; i < n; i++)
    {
        count = 0;
        for (j = 0; j < n; j++)
        {
            if (inv_a[i][j] == 0)
                count++;
            if (count == n)
            {
                printf("\nbaka\n\n");
                goto lab;
            }
            printf(" %d", inv_a[i][j]);
            z.x[i][j] = inv_a[i][j];
        }
        printf("\n");
    }
    // exit(1);

    printf("行列を出力\n ={\n");
    // #pragma omp parallel for num_threads(omp_get_max_threads()) //private(i,j,k)
    for (i = 0; i < n; i++)
    {
        printf("{");
        for (j = 0; j < n; j++)
        {
            // a[i][j]=rand()%N;
            printf("%3d,", c[i][j]);
        }
        printf("},\n");
    }
    printf("};");

    printf("\n逆行列を出力\n ={\n");
    for (i = 0; i < n; i++)
    {
        count = 0;
        printf("{");
        for (j = 0; j < n; j++)
        {
            if (inv_a[i][j] == 0)
                count++;
            if (count == n)
            {
                printf("\nbaka\n\n");
                goto lab;
            }
            printf("%3d,", inv_a[i][j]);
        }
        printf("},\n");
    }
    printf("};\n");
    // exit(1);

    memset(b, 0, sizeof(b));
    printf("検算\n");
    //   #pragma omp parallel for num_threads(omp_get_max_threads()) //private(i,j,k)
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
                b[i][j] ^= gf[mlt(fg[c[i][k]], fg[inv_a[k][j]])];

            printf("%d,", b[i][j]);
        }
        printf("\n");
    }
printf("\n");

    int flg = 0;
    for (i = 0; i < n; i++)
    {
        //   printf("%d",b[i][i]);
        // printf("==\n");
        if (b[i][i] == 1)
        {
            // printf("baka");
            //    exit(1);
            flg++;
        }
    }
    count = 0;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (b[i][j] == 0 && i != j)
                count++;
        }
    }
    if (flg == n && n * n - n == count)
        return z;

    goto lab;
}

int main()
{
    int i;
    unsigned short s[K + 1] = {0}, z1[N] = {0};
    vec v = {0}, x = {0} ,xx={0};
    OP f = {0},g={0};

    srand(clock());
     //mkg(K); // Goppa Code (EEA type)
    //van(K);          // RS-Code generate
    g=vv(K);           // Goppa Code's Parity Check (Berlekamp type)
    //mkerr(z1, T);    // generate error vector
    //z1[0]=gf[1];
    z1[5]=1;
    z1[6]=1;
    f = synd(z1, K); // calc syndrome
    x = o2v(f);      // transorm to vec
    // r = bma(x.x);    // Berlekamp-Massey Algorithm
    MTX b = {0};
    for(i=0;i<N;i++)
    printf("%d,",x.x[i]);
    printf("\n");
    //exit(1);
    /*
    x.x[0]=7;
    x.x[1]=3;
    x.x[2]=7;
    x.x[3]=7;
    */
    for (i = 0; i < K; i++){
        //v.x[i] = x.x[i];
        v.x[K - 1 - i] = x.x[i]; //gf[mlt(fg[x.x[i]],fg[trace(g,i)])];
    }
    for(i=0;i<N;i++){
        xx.x[i]=gf[mlt(fg[v.x[i]],fg[trace(g,gf[i])])];
    printf("%d,",v.x[i]);
    }
    printf("\n");
    //exit(1);

    for (i = 0; i < K / 2; i++)
    {
        for (int j = 0; j < K / 2 + 1; j++)
        {
            b.x[i][j] = v.x[i + j];
             printf("%d,",b.x[i][j]);
        }
         printf("\n");
    }
    printf("\n");

    MTX ee = {0}, ff = {0};
    
    x=sol(b);
    //exit(1);
    ee.x[0][0]=gf[x.x[0]+1];
    ee.x[0][1]=gf[x.x[1]+1];
    ee.x[1][0]=gf[mltn(2,fg[ee.x[0][0]])];//gf[2];
    ee.x[1][1]=gf[mltn(2,fg[ee.x[0][1]])]; ////gf[6];
    ff=matinv(ee,2);
    printpol(x);
    sol(ee);
    //exit(1);
    vec c={0};
    unsigned short or=0;
    c.x[0]=5;
    c.x[1]=4; //gf[4];
    for(int i=0;i<K/2;i++){
    or=0;
    for(int j=0;j<K/2;j++)
    or^=gf[mlt(fg[ff.x[i][j]],c.x[j])];
    printf("e=%d\n",or-1);
    }
    exit(1);
    
    for (i = 0; i < K / 2; i++)
    {
        for (int j = 0; j < K / 2 + 1; j++)
        {
            printf("%d,", fg[b.x[i][j]]);
            if (j < K / 2)
                ee.x[i][j] = fg[b.x[i][j]];
        }
        printf("\n");
    }
    
    v = sol(b);
    for(i=0;i<N;i++)
    if(v.x[i]>0)
    printf("zz=%d %d\n",fg[v.x[i]],v.x[i]);
    printf("\n");
    MTX pter={0};
    for(int i=0;i<T;i++){
    for(int j=0;j<T;j++){
    pter.x[i][j]=gf[mltn(i+1,v.x[j])];
    //
    }
    pter.x[i][T]=x.x[i];
    }
    for(int i=0;i<T;i++){
    for(int j=0;j<T+1;j++)
    printf("%d,",pter.x[i][j]);
    printf("\n");
    }
    renritz(pter);
//exit(1);
    //exit(1);
    v=renritz(pter);
    for (i = 0; i < N; i++)
        if (z1[i] > 0)
            printf("anaser=%d %d %d\n", i, z1[i],v.x[i]);
    printf("\n");
// exit(1);
    for(int j=0;j<T;j++){
        ff.x[j][0]=xx.x[j];
    for(i=0;i<T;i++){
        ee.x[0][i]=v.x[i];
        ee.x[j][i]=gf[mltn(j+1,fg[v.x[i]])];
    }
    }
    //exit(1);
    ee = matinv(ee, K/2);
    for(i=0;i<K/2;i++){
        //for(int j=0;j<2;j++)
        printf("%d,",ff.x[i][0]);
        printf("\n");
    }
    exit(1);
    //
    //
    //renritz(ee);
    // exit(1);
    for (int i = 0; i < K / 2; i++)
    {
        unsigned short k=0,l=0,n=0;
        //printf("xx=%d\n", ff.x[i][0]);
        ff.x[K/2][i]=0;
        for (int j = 0; j < K / 2; j++)
        {
            //printf("%d,",ee.x[i][j]);
            ff.x[K/2][i] ^= gf[mlt(fg[ee.x[i][j]], fg[ff.x[j][0]])];
        }
        k=eval(g,v.x[i]);
        l=ff.x[K/2][i];
        n=gf[mlt(fg[k],fg[l])];
        //printf("ff %d %d\n",l, k);
        printf("ff %d %d %d %d\n",gf[k+1], n, l,gf[oinv(l)]);
    }
    for (i = 0; i < N; i++)
        if (z1[i] > 0)
            printf("anaser=%d %d\n", i, z1[i]);
    printf("\n");
    exit(1);

    for (i = 0; i < K; i++)
        s[i + 1] = x.x[i];
    v = bms(s, K + 1);
    printf("the errors below\n");
    for (i = 0; i < N; i++)
    {
        if (z1[i] > 0)
            printf("i=%d\n", i); // print answer
    }

    printf("%d\n", deg(v));
    printpol(v); // print error locater polynomial
    printf("==bms\n");

    chen((v)); // searching error position

    /*
    printf("%d\n", deg(r));
    printpol(r); // print error locater polynomial
    printf("\n");
    chen((r)); // searching error position
    printpol(r);
    printf("==bma\n");
    */

    return 0;
}
