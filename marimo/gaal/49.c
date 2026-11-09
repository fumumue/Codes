#include <stdio.h>

#include "global-p.h"
#include "struct-p.h"

#define Pr 7
#define O 49
#define ORD 49
#define EXP 2

unsigned short pp[13][4] = {{0, 0, 9, 2}, {0, 0, 11, 2}, {0, 0, 16, 3}, {0, 0, 15, 2}, {0, 0, 1, 2}, {0, 1, 0, 2}, {0, 0, 1, 1}, {0, 0, 1, 2}, {1, 1, 2, 2}, {0, 0, 1, 2}, {0, 0, 21, 5}, {0, 0, 30, 3}, {0, 0, 1, 4}};
unsigned short gf[O], fg[O];

unsigned short plus(unsigned short a, unsigned short b);

typedef struct
{
    unsigned long long re;
    unsigned long long im;
} com;

typedef struct
{
    com x;
    com y;
} PO;

// 多項式の次数(default)
int deg(vec a)
{
    int n = 0, flg = 0;

    // #pragma omp parallel for
    for (int i = 0; i < DEG; i++)
    {
        if (a.x[i] > 0 && a.x[i] <= N)
        {
            n = i;
            // flg = 1;
        }
    }

    return n;
}

// 多項式を表示する（vec型）
void printpol(vec g)
{
    int n;
    vec f = {0};
    f = g;
    // f = conv(f);
    n = deg(g);
    // printf("n=%d\n", n);
    //  printf("terms=%d\n", terms(f));
    //  printf("deg=%d\n", deg(f));

    for (int i = n; i > -1; i--)
    {
        if (f.x[i] > 0)
            printf("%dx^%d+", f.x[i], i);
    }
    // exit(1);
}

// 多項式を項ずつ掛ける
vec oterml(vec f, oterm t)
{
    // f=conv(f);
    //  assert (vec_verify (f));
    // int i, k, j;
    vec h = {0};
    vec test;
    unsigned short n;

    printpol(f);
    printf(" ==f\n");
    printf("t.a=%d\n", t.a);
    printf("t.n=%d\n", t.n);
    // exit(1);
    //
    int k = deg(f) + 1;
    printf("k=%d\n", k);
    int j = 0;
    for (int i = 0; i < k; i++)
    {
        if (f.x[i] > 0)
        {
            h.x[i + t.n] = (f.x[i] * t.a) % Pr;
            printf("t+n=%d\n", i + t.n);
        }
    }

    // h=conv(h);
    // assert (vec_verify (h));
    return h;
}

// 多項式の足し算
vec vadd(vec a, vec b)
{
    vec c = {0};
    // int i, j, k, l = 0;
    vec h = {0}, f2 = {0}, g2 = {0};

    for (int i = 0; i < DEG; i++)
    {
        c.x[i] = (a.x[i] + b.x[i]) % Pr;
    }

    return c;
}

// 20200816:正規化したいところだがうまく行かない
// 多項式の足し算
vec vadd2(vec f, vec g)
{
    // assert(vec_verify(f));
    // assert(vec_verify(g));

    vec a = {0}, b = {0}, c = {0};
    // int i, j, k, l = 0;
    vec h = {0}, f2 = {0}, g2 = {0};

    a = (f);
    b = (g);

    // k=deg((f));
    // l=deg((g));

    for (int i = 0; i < DEG; i++)
    {
        c.x[i] = plus(a.x[i], b.x[i]);
        // h.x[i]=f.x[i]^g.x[i];
    }
    h = (c);
    // h=conv(h);
    // assert(vec_verify(h));
    return h;
}

// 20200816:正規化したいところだがうまく行かない

// 多項式の掛け算
vec vmul(vec f, vec g)
{

    // assert (vec_verify (f));
    // assert (vec_verify (g));
    int count = 0, k, l, m;
    oterm t = {0};
    vec h = {0}, e = {0}, r = {0};
    vec c = {0};

    l = deg((f));
    m = deg((g));

    for (int i = 0; i < 2; i++)
        printf("%d\n", g.x[i]);
    // exit(1);

    if (l >= m)
    {
        k = l;
    }
    else
    {
        k = m;
    }
    // printf("l=%d",l);
    // printf("m=%d",m);
    printpol((f));
    printf(" =f\n");
    printpol((g));
    printf(" =g\n");
    // exit(1);

    for (int i = 0; i < k + 1; i++)
    {
        t.a = (g.x[i]);
        t.n = i;
        // t.n=i;
        if (t.a > 0 && t.a < O)
        {
            printf("t[%d]=%d,%d\n", i, t.a, t.n);
            e = oterml(f, t);
            printpol((e));
            printf(" =e\n");
            // exit(1);
            printpol((h));
            printf(" =h\n");
            h = vadd(h, e);
        }
    }
    printpol((h));
    printf(" =h2\n");

    // printpol(o2v(g));
    // printf(" =g\n");
    //    exit(1);
    // assert (op_verify (h));
    return h;
}

vec confer(vec f, int a)
{
    vec r;
    int n;
    vec g;

    r = (f);
    n = deg(r);
    for (int i = 0; i < n + 1; i++)
        r.x[i] = (r.x[i] * a) % Pr;
    g = (r);

    return g;
}

unsigned short xtrace(vec f, unsigned short x)
{
    int d, z = 0;
    vec g = (f);
    unsigned short u = 0;

    d = deg(g) + 1;
    // if(g.x[0]>0)
    // z+=g.x[0];

    for (int i = 0; i < d; i++)
    {
        u = 1;
        if (g.x[i] > 0)
        {
            for (int j = 0; j < i; j++)
                u = (u * x) % O;
            u = u * g.x[i] % O;
            z += u % O;
        }

        // plus(u,gf[mlt(fg[g.x[i]],mltn(i,fg[Pr]))]);
    }
    // if(g.x[0]>0)
    // z+=g.x[0];

    return z % O;
}

void makefg()
{
    // unsigned short i, j, count = 0;

    // for( int i= 0; i < O; i++)
    {
        for (int j = 0; j < O; j++)
        {
            fg[gf[j]] = j;
        }
    }

    // exit(1);

    printf("unsigned short fg[%d]={", O);
    for (int i = 0; i < O; i++)
        printf("%d,", fg[i]);
    printf("};\n");
    // printf("count=%d\n", count);
    //  exit(1);

    return;
}

// 配列の値を係数として多項式に設定する
vec setpol(unsigned short f[], int n)
{
    vec g;
    vec a = {0};

    // memset (a, 0, sizeof (c));
    // memcpy (c, f, n);
    for (int i = 0; i < n; i++)
    {
        a.x[n - 1 - i] = f[i];
    }

    // exit(1);
    // a = Setvec (n);

    g = (a);

    return a;
}

// リーディグタームを抽出(defauvLT)
oterm vLT(vec f)
{
    int k;
    oterm t = {0};

    // k = deg ( (f));
    for (int i = 0; i < DEG; i++)
    {
        // printf("a=%d %d\n",f.x[i],f.x[i].n);
        if (f.x[i] > 0)
        {
            t.n = i;
            t.a = f.x[i];
        }
    }

    return t;
}

vec dick[O] = {0};
unsigned short val[N] = {0}; //, a = {0}, cc = {0};
void mkmf()
{
    // int i, j, k, count = 0;
    vec f = {0}, g = {0}, h = {0}, w = {0}, s = {0}, u = {0};
    vec b = {0}, a = {0}, d = {0}, t = {0}, v = {0};
    oterm o = {0};
    unsigned short ccp[4] = {0};

    if (O == 1331)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[0][i];
    }
    if (O == 2197)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[1][i];
    }
    if (O == 4913)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[2][i];
    }
    if (O == 6859)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[3][i];
    }
    if (O == 3125)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[4][i];
    }
    if (O == 2187)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[5][i];
    }
    if (O == 9)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[6][i];
    }
    if (O == 27)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[7][i];
    }
    if (O == 243)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[9][i];
    }
    if (O == 19683)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[8][i];
    }
    if (O == 12167)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[10][i];
    }
    if (O == 29791)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[11][i];
    }
    if (O == 49)
    {
        for (int i = 0; i < 4; i++)
            ccp[i] = pp[12][i];
    }

    g = (setpol(ccp, 4));
    for (int i = 0; i < 4; i++)
        printf("%d,", g.x[i]);
    printf("\n");
    // exit(1);
    printpol((g));
    printf("\n");
    // exit(1);
    //  b.x[0]=2;
    //  b.x[1]=9;
    b.x[0] = 0;
    a.x[0] = 1;
    vec vv = {0};
    v.x[1] = 1;
    v.x[0] = 0;
    // u = (v);
    //  d.x[EXP] = 1;
    printpol((g));
    printf(" =g\n");
    // exit(1);

    // g=v2o(b);
    s = (v);
    // s.t[1].a=1;
    // s.t[1].n=1;
    // gf[12]=P;
    // gf[13]=P*P;
    w = g;
    printpol((w));
    printf(" =w\n");
    printpol((g));
    printf(" =g\n");
    printpol((v));
    printf(" =s\n");
    // exit(1);
    vec c = {0};
    c.x[0] = 2;
    c.x[1] = 0;
    c.x[2] = 2;
    printf("err=%d\n", xtrace((c), Pr));
    // exit(1);
    printf("\n");
    // for(i=0;i<P;i++)
    gf[0] = 0;
    gf[1] = 1;
    //  gf[2] = Pr;
    dick[0] = (b);
    // gf[EXP]=Pr;
    dick[1] = (a);
    // dick[2]=v2o(v);
    for (int i = 2; i < EXP + 1; i++)
    {
        gf[i] = (gf[i - 1] * Pr);
        dick[i] = vmul(v, dick[i - 1]);
        printpol(dick[i]);
        printf(" %d dick\n", gf[i]);
    }
    // exit(1);

    // gf[2]=Pr;
    // gf[3]=Pr*Pr;
    gf[EXP + 1] = xtrace(w, Pr);
    printf("\naa=%d\n", gf[Pr]);
    // exit(1);
    // w=omul(w,s);
    // gf[12]=1111;
    dick[EXP + 1] = w;
    int count = EXP + 1;
    for (int i = 0; i < EXP + 1; i++)
    {
        printpol(dick[i]);
        printf(" ==beef\n");
    }
    // g=w;
    while (1)
    {
        printpol((g));
        printf(" ==beef %d\n", count);
        dick[count] = g;
        gf[count] = xtrace(g, Pr);
        printf("count2=%d %d ", count, gf[count]);
        printpol((g));
        printf(" =gg\n\n");
        if (gf[count] == 1)
        {
            printf("count!=%d\n", count);
            // break;
        }

        g = vmul(g, s);
        printpol((g));
        printf(" =g\n\n");
        printf(" papaya\n");

        // exit(1);

        o = vLT(g);
        memset(d.x, 0, sizeof(d));

        if (o.n == EXP)
        {
            vec xx = (g);
            vec ww = (w);
            xx.x[EXP] = 0;
            // printpol(v2o(xx));
            // exit(1);
            // d.x[o.n] = o.a;
            // xx=vmul(xx,v);
            f = w;
            g = (xx);
            // g = osub(g, h);
            if (o.a > 0)
                f = confer(f, o.a);
            g = vadd(g, f);
            // w=omod(w,u);
        }

        count++;
        if (count == O)
            break;
    }

    printf("unsigned short gf[%d]={", O);
    for (int i = 0; i < O; i++)
        printf("%d,", gf[i]);
    printf("};");
    printf("\n");

    // exit(1);
}

void de()
{

    int i, j;
    for (int i = 0; i < O; i++)
    {
        printpol((dick[i]));
        printf(", %d, %d lueda\n", xtrace(dick[i], Pr), i);
    }
}

unsigned short plus(unsigned short a, unsigned short b)
{
    unsigned short u;
    if (a == 0 && b > 0)
        return b;
    if (b == 0 && a > 0)
        return a;
    if (a == 0 && b == 0)
        return 0;
    u = xtrace(vadd(dick[fg[a]], dick[fg[b]]), Pr);

    return u;
}

// unsigned long long p = 1234567891;
unsigned short p = 49;
com cadd(com a, com b)
{
    com c;

    c.re = plus(a.re, b.re);
    c.im = plus(a.im, b.im);
    c.re = c.re;
    c.im = c.im;

    return c;
}

com inv_add(com a)
{ // -a
    com c;

    c.re = -1;
    c.im = -1;
    c.re = c.re * a.re % p;
    c.im = c.im * a.im % p;

    return c;
}

com csub(com a, com b)
{
    com c;

    c.re = (a.re - b.re);
    c.im = (a.im - b.im);

    return c;
}

int mlt(int x, int y)
{

    if (x == 0 || y == 0)
        return 0;

    return ((x + y - 2) % (N - 1)) + 1;
}

int mltn(int n, int x)
{
    int i, j;

    if (n == 0)
    {
        return 1;
    }
    if (x == 0)
    {
        return 0;
    }
    else
    {
        i = x;
        for (j = 0; j < n - 1; j++)
            i = mlt(i, x);

        return i;
    }
}

com cmul(com a, com b)
{
    com c;
    vec g = {0};
    unsigned long long d, e;
    // printf("a=%d %d ,b=%d b=%d\n",a.im,a.re,b.im,b.re);
    g = dick[mlt(fg[a.im], fg[b.im])];
    d = deg(g) + 1;
    for (int i = 0; i < d; i++)
        g.x[i] = (Pr - g.x[i]) % Pr;
    unsigned f = xtrace(g, Pr);
    // printf("min=%d\n",f);
    c.re = plus(gf[mlt(fg[a.re], fg[b.re])], f);
    d = gf[mlt(fg[a.re], fg[b.im])];
    e = gf[mlt(fg[b.re], fg[a.im])];

    c.im = plus(d, e);

    return c;
}

com cmod(com a, com b)
{
    a.im %= b.im;
    a.re %= b.re;

    return a;
}

int count = 0;
void hermite(com x, com y)
{
    int n;
    com c = {0};

    c = cmul(cmul(cmul(y, y), cmul(y, y)), cmul(cmul(y, y), y));
    c = cadd(c, y);
    com b = cmul(cmul(cmul(x, x), cmul(x, x)), cmul(cmul(x, x), cmul(x, x)));

    if (c.re % p == b.re % p && c.im % p == b.im % p)
    {
        printf("%d+%di,%d+%di\n", x.re, x.im, y.re, y.im);
        printf("count==%d c=%d %d b=%d %d\n", count, c.re, b.re, c.im, b.im);
        count++;
    }
}

int main()
{
    com x, y;
    mkmf();
    makefg();
    de();
    // exit(1);
    for (int i = 0; i < Pr; i++)
    {
        x.re = i;
        for (int j = 0; j < Pr; j++)
        {
            x.im = j;
            for (int k = 0; k < Pr; k++)
            {
                y.re = k;
                for (int l = 0; l < Pr; l++)
                {
                    y.im = l;
                    hermite(x, y);
                }
            }
        }
    }

    return 0;
}