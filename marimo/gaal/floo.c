// date 20200331:xgcdの終了条件が２つになってしまう。ogcdとxgcdで使い分ける。
// date : 20200326 鍵生成が det と deta で高い確率で一致する。detaは並列処理。
// date 20200229 : pattarson algorithm implementation ver 1.0
//  xgcd & osqrtを追加した
// date      :  20160310,20191218,20191220,20191221,20191223,20191224,20191225,20191229,20191230
// auther    : the queer who thinking about cryptographic future
// code name :  一変数多項式演算ライブラリのつもり
// status    : now in debugging (ver 0.8)
//  0ベクトルが出ないように生成多項式のトレースチェックを入れた。
// date      :  20160310,20210419
// auther    : the queer who thinking about cryptographic future
// code name : OVP - One Variable Polynomial library with vecenMP friendly
// status    : now in debugging (ver 0.9)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <assert.h>
#include <execinfo.h>

#include "global-p.h"
#include "struct-p.h"

#include "debug.c"
#include "chash-p.c"
#include "sha3.c"
#include "gaus.c"
#include "inv_mat.c"

#define TH omp_get_max_threads()
// #define O 27 // 1331,2197,4913,6859,3125,2187,19683,29791,12167
//  #define F 16
extern unsigned long xor128(void);
extern int mlt(int x, int y);
// extern int mltn(int n, int a);
extern void makeS();

typedef struct
{
    vec f;
    vec g;
    vec h;
} ymo;

vec vv(int kk);
unsigned short plus(unsigned short a, unsigned short b);
unsigned short bexp(unsigned short a);
// unsigned short minus(unsigned short a);
unsigned short eval(vec f, unsigned short x);
void printpol(vec a);
ymo bm_itr(unsigned short s[]);
// #pragma omp threadprivate(mat)
// シンドロームのコピー
unsigned short sy[K] = {0};
unsigned short pp[13][4] = {{0, 0, 9, 2}, {0, 0, 11, 2}, {0, 0, 16, 3}, {0, 0, 15, 2}, {0, 0, 1, 2}, {0, 1, 0, 2}, {0, 0, 1, 1}, {0, 0, 1, 2}, {1, 1, 2, 2}, {0, 0, 1, 2}, {0, 0, 21, 5}, {0, 0, 30, 3}, {0, 0, 1, 4}};
unsigned short pt[3] = {0, 1, 2};
static CTX S, inv_S;
static CTX big = {0};
static CTX fls = {0};
static CTX lng ={0};
// Gvecpa多項式
static unsigned short g[K + 1] = {1, 2, 3, 4, 6};
//{1,0,0,0,0};
//{1,0,1,0,5};
//{1};
//{1, 0, 0, 0, 1, 0, 1};
//

// 有限体の元の逆数
unsigned short oinv(unsigned short a)
{

    if (a == 0)
        return 0;
    if (a == 1)
        return a;

    return N - fg[a] + 1;
}

unsigned short oinv2(unsigned short a, unsigned short n)
{
    unsigned short i;

    if (a == 0)
        return 0;
    if (a == 1)
        return 1;
    for (int i = 2; i < N; i++)
    {
        if (gf[mlt(i, fg[a])] == 1)
            return i;
    }
    printf("no return\n");
    exit(1);
}

// 有限体の元の逆数
unsigned short oinv3(unsigned short a)
{

    if (a == 0)
        return 0;

    for (int i = 0; i < N; i++)
    {
        if (gf[mlt(fg[a], i)] == 1)
            return (unsigned short)i;
    }

    printf("no return3 \n");
    //  exit (1);
}

// aに何をかけたらbになるか
unsigned short equ(unsigned short a, unsigned short b)
{
    int i = 0;
    for (int i = 0; i < M; i++)
    {
        if (gf[mlt(fg[a], fg[i])] == b)
            return i;
    }
    return i;
}

// 停止コマンド
void wait(void)
{
    int a;                                     // 読み込む変数はローカルに取るべし
    printf(" (enter number and hit return) "); // 何か表示させたほうが良いだろう
    fflush(stdout);                            // just in case
    scanf("%d", &a);                           // fgets(line, LINESIZE, stdin); という手も
}

int deg2(vec v)
{
    int n, i;

    for (int i = 0; i < DEG; i++)
        if (v.x[i] > 0)
            n = i;

    return i;
}

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

void vec_print_raw(const vec f)
{
    puts("vec_print_raw:");
    for (int i = 0; i < DEG; i++)
    {
        if (f.x[i] > 0)
            printf("[%d] %ux^%u\n", i, f.x[i], i);
    }
}

bool vec_verify(const vec f)
{
    bool end = false;
    unsigned short n_max = 0;
    for (int i = 0; i < DEG; i++)
    {
        if (end && (i != 0 || f.x[i] != 0))
        {
            vec_print_raw(f);
            printf("found data after end: i=%d\n", i);
            print_trace();
            fflush(stdout);
            return false;
        }
        if (f.x[i] == 0)
        {
            end = true;
            continue;
        }
        if (i + 1 <= n_max)
        {
            vec_print_raw(f);
            printf("found invalid order: i=%d\n", i);
            print_trace();
            fflush(stdout);
            return false;
        }
        n_max = i + 1;
    }
    return true;
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

// 20200816:正規化したいところだがうまく行かない
// 多項式の足し算
vec vsub(vec a, vec b)
{
    vec c = {0};
    // int i, j, k, l = 0;
    vec h = {0}, f2 = {0}, g2 = {0};

    for (int i = 0; i < DEG; i++)
    {
        c.x[i] = (a.x[i] - b.x[i]) % Pr;
    }

    return c;
}

// 20200816:正規化したいところだがうまく行かない
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

// 多項式を項ずつ掛ける
vec vterml(vec f, oterm t)
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
            unsigned short tmp = equ(t.a, f.x[i]);
            h.x[i - t.n] = (f.x[i] - gf[mlt(fg[tmp], fg[t.a])]) % Pr;
            printf("t+n=%d\n", i - t.n);
        }
    }

    // h=conv(h);
    // assert (vec_verify (h));
    return h;
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

// 多項式を項ずつ掛ける
vec oterml2(vec f, oterm t)
{

    // assert (vec_verify (f));
    // int i, k, j;
    vec h = {0};
    vec test;
    unsigned short n;

    // f=conv(f);
    int k = deg(f);
    int j = 0;
    for (int i = 0; i < k + 1; i++)
    {
        h.x[i + t.n] = gf[mlt(fg[f.x[i]], fg[t.a])];
    }

    // h=conv(h);
    // assert (vec_verify (h));
    return h;
}

// 多項式を項ずつ掛ける
vec otermla(vec f, oterm t)
{
    // f = conv(f);
    // assert(vec_verify(f));
    // int i, k, j;
    vec h = {0};
    vec test;
    unsigned int n;

    // f=conv(f);
    // k = deg ((f));
    int j = 0;
    for (int i = 0; i < DEG; i++)
    {
        if (f.x[i] > 0)
            h.x[i + t.n] = gf[mlt(fg[f.x[i]], fg[t.a])];
    }

    // h = conv(h);
    // assert(vec_verify(h));
    return h;
}

vec padd(vec f, vec g)
{
    vec c = {0}, d = (f), e = (g);
    int k;

    if (deg(d) >= deg(e))
    {
        k = deg(d) + 1;
    }
    else
    {

        k = deg(e) + 1;
    }
    for (int i = 0; i < k; i++)
    {
        c.x[i] = plus(d.x[i], e.x[i]);
    }
    // c.x[i] = a.x[i] ^ b.x[i];
    //  h =  (c);

    return (c);
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

// 項の数
int terms(vec f)
{
    int count = 0;

    for (int i = 0; i < DEG; i++)
        if (f.x[i] > 0)
            count++;

    return count;
}

// 項の順序を降順に揃える
vec sort(vec f)
{
    oterm o = {0};
    // int i, j, k;

    int k = terms(f);
    for (int i = 0; i < k + 1; ++i)
    {
        for (int j = i + 1; j < k + 1; ++j)
        {
            if (f.x[i] > f.x[j])
            {
                o = vLT(f);
                f.x[i] = f.x[j];
                f.x[o.n] = o.a;
            }
        }
    }

    return f;
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

// モニック多項式にする
vec coeff(vec f, unsigned short d)
{
    int j, k;
    vec a, b;

    k = deg((f)) + 1;
    for (int i = 0; i < k; i++)
        f.x[i] = gf[mlt(fg[f.x[i]], oinv(d))];

    return f;
}

// 多項式の剰余を取る
vec vmod(vec f, vec g)
{
    // int i = 0, j, n, k;
    vec h = {0}, e = {0};
    vec x, y;
    oterm a, b = {0}, c = {0};

    int n = vLT(g).n;

    //  assert (("baka^\n", vLT (f).n != 0));

    //  assert (("baka(A)\n", vLT (g).n != 0));

    if (vLT(f).n < vLT(g).n)
    {
        //    exit(1);
        return f;
    }

    // printf ("in vmod\n");
    // exit(1);

    int k = vLT(g).n;
    b = vLT(g);

    x = (g);
    assert(("double baka\n", b.a != 0 && b.n != 0));
    while (vLT(f).n >= 0 && vLT(g).n > 0)
    {

        c = vLTdiv(f, b);
        h = oterml((x), c);
        f = vsub(f, h);
        if (deg((f)) == 0 || deg((g)) == 0)
        {
            //      printf("blake1\n");
            break;
        }

        if (c.n == 0 || b.n == 0)
            break;
    }

    return f;
}

// 多項式の商を取る
vec vdiv(vec f, vec g)
{

    assert(vec_verify(f));
    assert(vec_verify(g));
    // int i = 0, j, n, k;
    vec h = {0}, e = {0}, tt = {0};
    oterm a, b = {0}, c = {0};
    vec x, y;

    if (vLT(f).n == 0 && vLT(g).a == 0)
    {
        printf("baka^\n");
        // return f;
        exit(1);
    }
    if (vLT(g).a == 0)
    {
        print_trace();
        exit(1);
    }
    if (vLT(g).n == 0 && vLT(g).a > 1)
        return coeff(f, vLT(g).a);

    int k = deg(g);
    b = vLT(g);
    if (b.a == 1 && b.n == 0)
        return f;
    if (b.a == 0 && b.n == 0)
    {
        printf("baka in vdiv\n");
        exit(1);
    }
    if (deg((f)) < deg((g)))
    {
        return f;
        //  a=vLT(f);
    }

    x = (g);
    int i = 0;
    while (vLT(f).n > 0 && vLT(g).n > 0)
    {
        c = vLTdiv(f, b);
        assert(c.n < DEG);
        tt.x[c.n] = c.a;
        i++;

        h = (oterml((x), c));

        f = vadd(f, h);
        if (deg((f)) == 0 || deg((g)) == 0)
        {
            // printf ("blake2\n");
            break;
        }

        if (c.n == 0)
            break;
    }

    // tt は逆順に入ってるので入れ替える
    vec ret = {0};
    int tt_terms = terms(tt);
    for (int i = 0; i < tt_terms; i++)
    {
        ret.x[i] = tt.x[tt_terms - i - 1];
    }

    assert(vec_verify(ret));
    return ret;
}

// 多項式のべき乗
vec opow(vec f, int n)
{
    // int i;
    vec g = {0};

    g = f;

    for (int i = 1; i < n; i++)
        g = vmul(g, f);

    return g;
}

// 多項式のべき乗余
vec opowmod(vec f, vec mod, int n)
{
    // int i, j = 0;

    // 繰り返し２乗法
    for (int i = 1; i < n + 1; i++)
    {
        f = vmul(f, f);
        if (deg(f) > deg(mod))
            f = vmod(f, mod);
    }

    return f;
}

vec osub(vec f, vec g)
{
    vec a = {0}, b = {0}, d = {0};
    // int i, k, l, m;
    vec ans = {0};

    a = (f);
    b = (g);
    int l = deg(a);
    int m = deg(b);
    int k = 0;
    if (l >= m)
    {
        int k = l;
    }
    else
    {
        k = m;
    }
    for (int i = 0; i < k + 1; i++)
    {
        if (a.x[i] >= b.x[i])
        {
            d.x[i] = (a.x[i] - b.x[i]) % Pr;
        }
        else
        {
            d.x[i] = (Pr + (a.x[i] - b.x[i])) % Pr;
        }
    }
    ans = (d);

    return ans;
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

    // exit(1);

    // printpol(o2v(f));
    // printf(" =f\n");
    // printf("gf[%d]={\n",O);
    printf("unsigned short gf[%d]={", O);
    for (int i = 0; i < O; i++)
        printf("%d,", gf[i]);
    printf("};");
    printf("\n");

    // exit(1);
}

unsigned short v2c(vec v){
uni w={0};

    w.fugo.b0 = v.x[0];
    w.fugo.b1 = v.x[1];
    w.fugo.b2 = v.x[2];

return w.u;
}

vec c2v(uni in){
vec x={0};

x.x[0]=in.fugo.b0;
x.x[1]=in.fugo.b1;
x.x[2]=in.fugo.b2;

return x;
}

uni d2c(unsigned short a)
{
    vec v = dick[fg[a]];
    uni w = {0};

    w.fugo.b0 = v.x[0];
    w.fugo.b1 = v.x[1];
    w.fugo.b2 = v.x[2];
    // for (int i = 0; i < EXP; i++)
     printf("d%d,", w.u);
    // printf("\n");
    //  exit(1);
    return w;
}

// union to vec 2 short
unsigned short c2d(uni w)
{
    vec x = {0};
    x.x[0] = w.fugo.b0;
    x.x[1] = w.fugo.b1;
    x.x[2] = w.fugo.b2;
    return xtrace((x), Pr);
}

vec o2c(unsigned short a)
{
    vec v = dick[fg[a]];

    for (int i = 0; i < EXP; i++)
        printf("%d,", v.x[i]);
    // printf("\n");

    return v;
}

unsigned short c2o(vec x)
{

    return xtrace((x), Pr);
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

vec pmul(vec a, vec b)
{
    int i, j, k, l;
    vec c = {0};

    k = deg(a) + 1;
    l = deg(b) + 1;

    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < l; j++)
            if (a.x[i] > 0)
            {
                c.x[i + j] = plus(c.x[i + j], gf[mlt(fg[a.x[i]], fg[b.x[j]])]);
                // printf("%d=c ",c.x[i+j]);
            }
        // printf("\n");
    }
    /*
    printf("\n");
    printpol(v2o(c));
    printf(" ==c\n");
    printpol(v2o(a));
    printf(" ==a\n");
    printpol(v2o(b));
    printf(" ==b\n");
    // exit(1);
    */
    return c;
}

unsigned short v2u(vec u)
{
    unsigned short s = 0;

    for (int j = 0; j < N; j++)
    {
        if (u.x[j] > 0)
            s = (s + gf[mltn(j, fg[Pr])]) % O;
    }
    return s;
}

// 多項式の代入値
unsigned short eval(vec f, unsigned short x)
{
    vec g = {0};
    unsigned short u = 0, s = 0;
    vec v = (f), h = {0};
    int d = deg((v)) + 1;

    for (int i = 1; i < d; i++)
    {
        if (v.x[i] > 0)
        {
            u = plus(u, gf[mlt(fg[v.x[i]], mltn(i, fg[x]))]);
            //
            // vtrace((vadd(dick[fg[u]], dick[mlt(fg[v.x[i]], mltn(i, x))])),x);
        }
    }
    if (v.x[0] > 0)
        u = plus(u, v.x[0]);

    return u;
}

unsigned short trace(vec g, unsigned short u)
{
    int n = 1, d = deg(g) + 1;

    for (int i = 0; i < d; i++)
    {
        if (g.x[i] > 0)
        {
            for (int j = 0; j < i; j++)
                n *= u % Pr;
            n += (n * g.x[i]) % Pr;
        }
        // plus(u,gf[mlt(fg[g.x[i]],mltn(i,fg[Pr]))]);
    }
    return n;
}

void de()
{

    int i, j;
    for (int i = 0; i < O; i++)
    {   j=v2c(dick[i]);
        printpol((dick[i]));
        printf(", %d, %d lueda %d\n", xtrace(dick[i], Pr), i,j);
    }
}

unsigned short plus2(unsigned short a, unsigned short b)
{
    unsigned short u;
    if (a == 0)
        return b;
    if (b == 0)
        return a;
    u = eval(vadd(dick[fg[a]], dick[fg[b]]), Pr);

    return u;
}

unsigned short hiku(unsigned short a, unsigned short b)
{
    unsigned short u;

    u = eval(osub(dick[fg[a]], dick[fg[b]]), Pr);

    return u;
}

unsigned short hiku2(unsigned short a, unsigned short b)
{
    unsigned short u;

    u = eval(osub(dick[fg[a]], dick[fg[b]]), Pr);

    return u;
}

// invert of polynomial
vec inv(vec a, vec n)
{
    vec d = {0}, x = {0}, s = {0}, q = {0}, r = {0}, t = {0}, u = {0}, v = {0}, w = {0}, tt = {0}, gcd = {0}, tmp = {0};
    oterm b = {0};
    vec vv = {0}, xx = {
                      0};

    if (deg((a)) > deg((n)))
    {
        tmp = a;
        a = n;
        n = tmp;
        printf("baka_i\n");
        // exit (1);
    }
    if (vLT(a).a == 0)
    {
        printf(" a ga 0\n");
        exit(1);
    }

    tt = n;

    d = n;
    x.x[0] = 0;
    s.x[0] = 1;
    s.x[0] = 0;
    while (deg((a)) > 0)
    {
        if (deg((a)) > 0)
            r = vmod(d, a);
        if (vLT(a).a == 0)
            break;
        if (vLT(a).a > 0)
            q = vdiv(d, a);

        d = a;
        a = r;
        t = vadd(x, vmul(q, s));
        ////printpol ( (a));
        // printf ("\nin rovec a==================%d\n", deg ((a)));
        // printf ("\n");

        x = s;
        s = t;
    }
    // exit(1);
    //  if(vLT(a).a>0){
    d = a;
    a = r;
    ////printpol ( (a));
    // printf ("\nin rovec a|==================%d\n", deg ((a)));
    // printf ("\n");

    x = s;
    s = t;

    ////printpol ( (d));
    // printf ("\nout1================\n");
    gcd = d; // $\gcd(a, n)$
    printpol(gcd);
    printf(" =========gcd\n");
    // exit(1);
    // printf ("\n");
    ////printpol ( (n));
    // printf ("\n");
    // printf ("out2===============\n");

    printf("before vdiv\n");
    // w=tt;

    b = vLT(w);
    ////printpol ( (w));
    // printf ("\nw=======%d %d\n", b.a, b.n);
    // w=tt;
    v = vadd(x, n);
    ////printpol ( (v));
    // printf ("\n");
    /*
     if (vLT (v).a == 0)
     {
     printf ("v=============0\n");
     }
     printf ("d==============\n");
   */
    //  } //end of a>0
    w = tt;
    ////printpol ( (v));
    // printf ("\n");
    // printf ("ss==============\n");
    //        exit(1);
    //  if(deg((w))>0)
    if (vLT(v).n > 0 && vLT(w).n > 0)
    {
        u = vmod(v, w);
    }
    else
    {
        // printpol ( (v));
        printf(" v===========\n");
        // printpol ( (x));
        printf(" x==0?\n");
        // printpol ( (n));
        printf(" n==0?\n");

        exit(1);
    }
    // caution !!
    if (vLT(u).a > 0 && vLT(d).a > 0)
    {
        u = vdiv(u, d);
    }

    if (vLT(u).a == 0 || vLT(d).a == 0)
    {
        printf("inv div u or d==0\n");
        // exit(1);
    }
    // u=coeff(u,d.x[0]);
    ////printpol ( (u));
    // printf ("\nu==================\n");
    if (vLT(u).a == 0)
    {
        printf("no return at u==0\n");
        exit(1);
    }

    return u;
}

unsigned short v2a(oterm a)
{
    int i, j;

    if (a.a == 0)
        return 0;

    // printf("aa=%d\n",a.a);
    for (int j = 0; j < M; j++)
    {
        if (gf[j] == a.a && a.a > 0)
        {
            // printf("j==%d\n",j);
            return j - 1;
        }
    }
}

void printsage(vec a)
{
    int i, j, k;
    oterm b;

    printf("poly=");
    for (int i = 0; i < DEG; i++)
    {
        if (a.x[i] > 0)
        {
            b.a = a.x[i];
            b.n = i;
            j = v2a(b);
            // printf("%d,==ba\n",b.a);
            // printf ("X**%d+", i); //for GF2
            printf("B('a^%d')*x**%d+", j, i); // for GF(2^m)
        }
    }
}

vec gcd(vec a, vec b)
{
    vec r = {0}, h = {0}, tmp = {0};

    h.x[0] = 1;

    if (deg(a) < deg(b))
    {
        tmp = a;
        a = b;
        b = tmp;
    }
    /*
  printpol((a));
  printf(" ========f\n");
  printpol((b));
  printf(" ========g\n");
*/
    /* 自然数 a > b を確認・入替 */
    if (deg(a) < deg(b))
    {
        tmp = a;
        a = b;
        b = tmp;
    }

    r = vmod(a, b);
    while (deg(r) > 0)
    {
        a = b;
        b = r;
        r = vmod(a, b);
        if (vLT(r).a == 0)
            return b;
    }

    if (vLT(r).a == 0)
    {
        return b;
    }
    else
    {
        // if(vLT(r).a>0)
        return h;
    }
}

// gcd for pattarson
vec zgcd(vec a, vec n)
{
    vec d = {0}, x = {0}, s = {0}, q = {0}, r = {0}, t = {0}, u = {0}, v = {0}, w = {0}, tt = {0}, gcd = {0}, rt = {0};
    oterm b = {0};
    vec vv = {0}, xx = {
                      0};

    if (deg(a) > deg(n))
    {
        rt = a;
        a = n;
        n = rt;
        printf("big is good\n");
        // exit (1);
    }
    if (vLT(a).a == 0)
    {
        printf(" a ga 0\n");
        exit(1);
    }

    tt = n;

    d = n;
    x.x[0] = 0;
    s.x[0] = 1;
    while (vLT(a).n > T)
    {

        r = vmod(d, a);
        q = vdiv(d, a);

        d = a;
        a = r;
        t = vadd(x, vmul(q, s));

        x = s;
        s = t;
    }

    d = a;
    a = r;

    x = s;
    s = t;
    gcd = d; // $\gcd(a, n)$

    printpol((x));
    printf(" =======x\n");
    printpol((a));
    printf(" =======a\n");
    printpol((s));
    printf(" =======s\n");
    printpol((r));
    printf(" =======r\n");

    return x;
}

// GCD for decode
vec ogcd(vec xx, vec yy)
{
    vec tt;

    while (deg(yy) > T - 1)
    {
        tt = vmod(xx, yy);
        xx = yy;
        yy = tt;
    }

    printpol((yy));
    printf(" =========yy\n");
    printpol((tt));
    printf(" =========tt\n");

    return tt;
}

// test gcd
vec agcd(vec xx, vec yy)
{
    vec tt = {0}, tmp;

    if (deg((xx)) < deg((yy)))
    {
        tmp = xx;
        xx = yy;
        yy = tmp;
    }
    tt = vmod(xx, yy);
    while (vLT(tt).n > 0)
    {
        xx = yy;
        yy = tt;
        tt = vmod(xx, yy);
    }

    return yy;
}

// error locater for decode
vec vx(vec f, vec g)
{
    vec h = {0}, ww = {
                     0};
    vec v[K] = {0}, vv = {
                        0};
    oterm a, b;
    int i, j;

    v[0].x[0] = 0;
    v[0].x[1] = 0;
    v[1].x[0] = 1;
    v[1].x[1] = 0;

    i = 0;

    while (1)
    {
        if (deg((g)) == 0)
            break;
        h = vmod(f, g);
        if (vLT(g).a == 0)
            break;
        ww = vdiv(f, g);
        v[i + 2] = vadd(v[i], vmul(ww, v[i + 1]));
        f = g;
        g = h;

        vv = v[i + 2];

        if (deg((vv)) == T)
            break;
        i++;
    }

    return vv;
}

// error locater for decode
vec sabun(vec f, vec g)
{
    vec h = {0}, ww = {
                     0};
    vec v[K] = {0}, vv = {
                        0};
    oterm a, b;
    int i, j;

    v[0].x[0] = 0;
    v[0].x[1] = 0;
    v[1].x[0] = 1;
    v[1].x[1] = 0;

    i = 0;

    while (1)
    {
        if (deg((g)) == 0)
            break;
        h = vmod(f, g);
        if (vLT(g).a == 0)
            break;
        ww = vdiv(f, g);
        v[i + 2] = vadd(v[i], vmul(ww, v[i + 1]));
        f = g;
        g = h;

        vv = v[i + 2];

        if (deg((vv)) == T * 2)
            break;
        i++;
    }

    return vv;
}

// 最終の項までの距離
int distance(vec f)
{
    int i, j = 0;

    for (int i = 0; i < DEG; i++)
    {
        if (f.x[i] > 0)
            j = i;
    }

    return j;
}

// ０多項式かどうかのチェック
unsigned char
chk(vec f)
{
    int i, flg = 0;
    vec x = {0};

    x = (f);
    for (int i = 0; i < DEG; i++)
    {
        if (x.x[i] > 0)
        {
            flg = 1;
            return 1;
        }
    }
    if (flg == 0)
        return 0;

    exit(1);
}

vec kof(unsigned short c, vec f)
{
    int i, j, k;
    vec b = {0}, h = {0};
    vec g = {0};

    c = fg[c];
    b = (f);
    k = deg(b);
    for (int i = 0; i < k + 1; i++)
    {
        h.x[i] = gf[mlt(c, fg[b.x[i]])];
    }
    g = (h);

    return g;
}

/*
// 拡張ユークリッドアルゴリズム
EX xgcd(vec f, vec g)
{
    vec h = {0}, ww = {0}, *v, *u;
    oterm a, b;
    int i = 0, j, k, flg = 0, m = deg(f);
    EX e = {0}, ee = {0};

    v = (vec *)malloc(sizeof(vec) * (DEG));
    u = (vec *)malloc(sizeof(vec) * (DEG));
    memset(v, 0, sizeof(vec) * DEG);
    memset(u, 0, sizeof(vec) * DEG);

    u[0].x[0] = 1;
    u[0].x[0] = 0;
    u[1].x[0] = 0;
    u[1].x[0] = 0;
    u[2].x[0] = 1;
    u[2].x[0] = 0;

    v[0].x[0] = 0;
    v[0].x[0] = 0;
    v[1].x[0] = 1;
    v[1].x[0] = 0;

    printpol((f));
    printf(" f===============\n");
    printpol((g));
    printf(" s===============\n");
    //  exit(1);

    k = 0;
    i = 0;
    while (vLT(g).n > 0)
    // for( int i= 0; i < T * 2; i++)
    {

        if ((deg((g)) == 0 && vLT(g).a == 0) || deg((f)) == 0)
        {
            flg = 1;
            printf("v[%d]=%d skipped deg(g)==0!\n", i, deg((v[i])));
            printpol((g));
            printf(" g========\n");
            e.d = f;
            e.v = v[i];
            e.u = u[i];

            free(v);
            free(u);
            // wait();
            return e;

            // exit (1);
            // return e;

            // break;
        }

        if (vLT(g).n > 0)
            h = vmod(f, g);
        printpol((h));
        printf(" ===hh\n");
        if (vLT(g).a > 0)
            ww = vdiv(f, g);
        printpol((ww));
        printf(" ===ww\n");

        v[i + 2] = vadd(v[i], vmul(ww, v[i + 1]));
        u[i + 2] = vadd(u[i], vmul(ww, u[i + 1]));
        printf("i+1=%d %d %d g=%d\n", i + 1, deg((v[i + 2])), T - 1, deg((g)));
        printpol((v[i + 2]));
        printf(" ==vvv\n");
        f = g;
        g = h;

        if (deg(v[i + 2]) == m - 1)
        {
            // printf("vaka\n");
            // wait();
            e.d = f;
            e.v = kof(gf[oinv(vLT(h).a)], v[i + 2]);
            e.u = u[i + 2];

            free(v);
            free(u);
            // wait();
            return e;
        }
        i++;
    }

    // printf ("i=%d\n", i);
    // wait();
    // printpol ((v[i]));
    printf("deg(v)=%d\n", deg((v[i + 2])));
    printf(" v=============\n");
    printf("deg(u)=%d\n", deg((u[i])));
    // printpol ( (u[i]));
    printf(" u=============\n");
    printf("deg(f)=%d\n", deg((f)));
    printf(" f=============\n");
    // exit(1);
    //   if(deg(v[i])==T-1){
    e.d = f;
    e.v = kof(gf[oinv(vLT(h).a)], v[i + 2]);
    e.u = u[i + 2];

    free(v);
    free(u);
    printf("end of fnc\n");
    //  wait ();

    return e;
}
*/
vec init_pol(vec f)
{
    int i;

    for (int i = 0; i < DEG; i++)
    {
        f.x[i] = 0;
    }

    return f;
}

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

    for (int j = 0; j < K + 1; j++)
        gg[j] = g[K - j];

    memcpy(g, gg, sizeof(g));
}

unsigned short pack(vec b)
{
    int i;
    unsigned short a = 0;

    for (int i = EXP - 1; i >= 0; i--)
        a ^= (b.x[EXP - i - 1] << 2);

    return a;
}

short unpack(unsigned short u)
{
    vec c = {0};
    c = dick[fg[u]];

    return xtrace(c, Pr);
}

// 整数からベクトル型への変換
vec i2v(unsigned int n)
{
    vec v = {0};
    int i = 0;

    while (n > 0)
    {
        v.x[i++] = n % 2;
        n = (n >> 1);
    }

    return v;
}

// ベクトル型から整数への変換
unsigned int
v2i(vec v)
{
    unsigned int d = 0, i, e = 0;

    for (int i = 0; i < deg(v) + 1; i++)
    {
        e = v.x[i];
        d ^= (e << i);
    }

    return d;
}

void printvec(vec v)
{
    int i, j;

    for (int i = 0; i < deg(v) + 1; i++)
    {
        // if (v.x[i] > 0)
        printf("%d:%d\n", i, v.x[i]);
    }
}

// 整数のべき乗
unsigned int
ipow(unsigned int q, unsigned int u)
{
    unsigned int i, m = 1;

    for (int i = 0; i < u; i++)
        m *= q;

    printf("in ipow====%d\n", m);

    return m;
}

vec ww[T] = {0};

// chen探索
vec chen(vec f)
{
    vec e = {0};
    int i, count = 0, n, x = 0;
    unsigned short z;

    n = deg((f));
    // exit(1);
    // #pragma omp parallel for private(i)
    for (x = 0; x < M; x++)
    {
        z = 0;
        // #pragma omp parallel for reduction (^:z)
        for (int i = 0; i < n + 1; i++)
        {
            if (f.x[i] > 0)
                z = plus(z, gf[mlt(mltn(i, fg[x]), fg[f.x[i]])]) % O;
        }
        if (z % O == 0)
        {
            e.x[count] = x;
            count++;
            printf("chen = %d\n", x);
        }
    }
    if (count != T)
    {
        printpol((e));
        printf(" ==cheneee!\n");
        // exit(1);
    }
    return e;
}

int oequ(vec f, vec g)
{
    vec v, x;
    int i, flg = 0;

    v = (f);
    x = (g);
    for (int i = 0; i < DEG; i++)
    {
        if (v.x[i] != x.x[i])
            return -1;
    }

    return 0;
}

// GF(2^m) then set m in this function.
int ben_or(vec f)
{
    int i, n, flg = 0;
    vec s = {0}, u = {0}, r = {0};
    vec v = {0}, x = {0};
    // if GF(8192) is 2^m and m==13 or if GF(4096) and m==12 if GF(16384) is testing
    int m = E;
    // m=12 as a for GF(4096)=2^12 defined @ gloal.h or here,for example m=4 and GF(16)

    v.x[1] = 1;
    s = (v);
    r = s;
    n = deg((f));

    if (n == 0)
        return -1;

    i = 0;

    // r(x)^{q^i} square pow mod
    while (i < n / 2 + 1)
    {
        flg = 1;
        // irreducible over GH(8192) 2^13
        r = opowmod(r, f, m);

        // irreducible over GF2
        // r=vmod(vecow(r,2),f);

        u = vadd(r, s);
        if (deg((u)) == 0 && vLT(u).a == 0)
            return -1;
        if (deg((u)) == 0 && vLT(u).a == 1)
        {
            i++;
            flg = 0;
        }
        if (deg((u)) > 0)
            u = gcd(f, u);

        if (deg((u)) > 0)
            return -1;

        if (flg == 1)
            i++;
    }

    return 0;
}

// バイナリ型パリティチェック行列を生成する
void bdet()
{
    int i, j, k;
    vec l = {0};
    uni ll = {0};

    unsigned char dd[E * K] = {0};
    FILE *ff;
    // static CTX R = {0};

    // ff = fopen("Hb.key", "wb");

    for (int j = 0; j < K; j++)
    {
        for (int i = 0; i < M; i++)
        {
            big.d[j][i] = d2c(mat[j][i]);
            // big.x[j][i] = mat[j][i];
        }
    }

    printf("\n");
    /*
    // #pragma omp parallel for
    for (int j = 0; j < K; j++)
    {
        for (int i = 0; i < M; i++)
            printf("%d ", big.d[j][i].fugo.b0);
        printf("\n");
        for (int i = 0; i < M; i++)
            printf("%d ", big.d[j][i].fugo.b1);
        printf("\n");
        if(E==3){
        for (int i = 0; i < M; i++)
            printf("%d ", big.d[j][i].fugo.b2);
        printf("\n");
        }        // fwrite(dd, 1, E * K, ff);
    }
    */
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < M; j++)
        {
            printf("%d ", big.d[i][j].fugo.b0);
            lng.x[i*E+0][j] = big.d[i][j].fugo.b0;
        }
        printf("\n");
        for (j = 0; j < M; j++)
        {
            printf("%d ", big.d[i][j].fugo.b1);
            lng.x[i*E+1][j] = big.d[i][j].fugo.b1;
        }
        printf("\n");
        for (j = 0; j < M; j++)
        {
            lng.x[i*E+2][j] = big.d[i][j].fugo.b2;
            printf("%d ", big.d[i][j].fugo.b2);
        }
        printf("\n");
    }
    //exit(1);
    //  fclose(ff);
    return;
}

// バイナリ型パリティチェック行列を生成する
static CTX bdet2()
{
    int i, j, k;
    vec l = {0};
    uni ll = {0};
    unsigned char dd[E * K] = {0};
    FILE *ff;
    static CTX R = {0};

    // ff = fopen("Hb.key", "wb");

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < K; j++)
        {
            ll = d2c(mat[j][i]);
            // l = o2c(mat[i][j]);
            //  #pragma omp parallel for
            // for (int k= 0; k < E; k++)
            {
                R.d[j][i] = ll;
                // R.d[i][j*E+0].u=ll.fugo.b0;
                // R.d[i][j*E+1].u=ll.fugo.b1;
                // R.d[i][j*E+2].u=ll.fugo.b2;
                printf("V%d %d %d ", R.d[j][i].fugo.b0, R.d[j][i].fugo.b1, R.d[j][i].fugo.b2);
                //
            }
        }
    }
    printf("\n");
    for (int i = 0; i < M; i++)
    {
        // #pragma omp parallel for
        for (int j = 0; j < K; j++)
        {
            printf("Xx%d,", R.d[j][i].u);
            // printf("%d %d %d",R.d[j][i].,R.d[j][i].b1,R.d[j][i].b2);
            //  dd[j] = BH[j][i];
        }
        // fwrite(dd, 1, E * K, ff);
        printf("\n");
    }

    // fclose(ff);
    return R;
}

static MTX bd2()
{
    int i, j, k, l;
    unsigned char dd[E * K] = {0};
    FILE *ff;
    vec v = {0};
    static MTX R = {0};

    // ff = fopen("Hb.key", "wb");

    // memset(BB.z,0,sizeof(BB.z));
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < K / 2 + 1; j++)
        {
            l = bm[j][i];
            printf("bm==%d %d\n", l, j);

            v = i2v(l);
            // #pragma omp parallel for
            for (int k = 0; k < E; k++)
            {
                R.x[i][j * E + k] = v.x[k];
                // l = (l >> 1);
            }
        }
    }

    return R;
}

MT bin(unsigned short c[K])
{
    int i, j, count = 0;
    vec v = {0}, x = {0};
    MT n = {0};

    for (int i = 0; i < K; i++)
    {
        memset(v.x, 0, sizeof(v.x));
        v = i2v(c[i]);
        for (int j = 0; j < E; j++)
            n.v[j++] = v.x[i];
    }
    n.f = j;
    return n;
}

MT vin(unsigned short s[K * E])
{
    vec v = {0};
    int i, j, count = 0;
    MT n = {0};

    for (int j = 0; j < K; j++)
    {
        for (int i = 0; i < E; i++)
        {
            v.x[i] = s[count];
            count++;
        }
        n.v[i] = v2i(v);
    }
    n.f = K;

    return n;
}

static MTX BB = {0};
// バイナリ型パリティチェック行列を生成する
void toBit(MTX L)
{
    unsigned char dd[E * K] = {0};
    FILE *ff;

    // ff = fopen("Hb.key", "wb");

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < K; j++)
        {
            int l = L.x[i][j];
            printf("l=%d,", l);
            // #pragma omp parallel for
            for (int k = 0; k < E; k++)
            {
                BB.x[i][j * E + k] = l % 2;
                l = (l >> 1);
            }
        }
        printf("\n");
    }
    // exit(1);
    /*
    for( int i= 0; i < N; i++)
    {
        //#pragma omp parallel for
        for( int j = 0; j < E * K; j++)
        {
            printf("%d,", BB.z[i][j]);
            //dd[j] = BH[j][i];
        }
        //fwrite(dd, 1, E * K, ff);
        printf("\n");
    }
*/
    // fclose(ff);
}

unsigned short HH[N][K];
unsigned short TE[N][K / 2 + 1];

static MTX toByte(MTX SH, int kk)
{
    vec v = {0};
    int i, j, k, cnt;
    static MTX R = {0};

    memset(HH, 0, sizeof(HH));
    printf("HH=");
    // exit(1);
    for (int i = 0; i < N; i++)
    {
        printf("%d\n", i);
        // #pragma omp parallel for
        for (int j = 0; j < kk; j++)
        {
            cnt = 0;
            for (int k = j * E; k < j * E + E; k++)
                v.x[cnt++] = SH.x[i][k];

            HH[i][j] = v2i(v);
            R.x[i][j] = v2i(v);
            // printf("%d,", HH[i][j]);
            //= BH[j][i];
        }
        // fwrite(dd, 1, E * K, ff);
        // printf("\n");
    }
    printf("end of byte\n");
    // exit(1);
    // wait();

    return R;
}

// 秘密置換を生成する
void Pgen()
{
    unsigned int i, j;
    FILE *fp;

    // fp = fopen("P.key", "wb");
    for (int i = 0; i < M; i++)
        P[i] = i;
    random_shuffle(P, SIZE_OF_ARRAY(P));
    // random_permutation(P);

    printf("\n\n");
    for (int i = 0; i < N; i++)
        printf("%d,", P[i]);
    printf("\n");
    // exit(1);

    for (int i = 0; i < M; i++)
        inv_P[P[i]] = i;
    // fwrite(P, 2, N, fp);
    // fclose(fp);
    /*
    for( int i= 0; i < N; i++)
    printf ("%d,", P[i]);
    printf ("\n");
    exit(1);
    */
    // fp = fopen("inv_P.key", "wb");
    // fwrite(SP, 2, N, fp);
    // fclose(fp);
}

// ハッシュ１６進表示
static void
byte_to_hex(uint8_t b, char s[23])
{
    unsigned i = 1;
    s[0] = s[1] = '0';
    s[2] = '\0';
    while (b)
    {
        unsigned t = b & 0x0f;
        if (t < 10)
        {
            s[i] = '0' + t;
        }
        else
        {
            s[i] = 'a' + t - 10;
        }
        i--;
        b >>= 4;
    }
}

// 有限体の元の平方を計算する
int isqrt(unsigned short u)
{
    int i, j, k;

    for (int i = 0; i < N; i++)
    {
        if (gf[mlt(i, i)] == u)
            return i;
    }

    printf("来ちゃいけないところに来ました\n");
    exit(1);
}

// 512bitの秘密鍵を暗号化
void encrypt(char buf[], unsigned char sk[64])
{
    const uint8_t *hash = {0};
    sha3_context c = {0};
    int image_size = 512, i;
    FILE *fp;
    //  unsigned short dd=0;

    printf("plain text=");
    for (int i = 0; i < 64; i++)
        printf("%u,", sk[i]);
    printf("\n");
    //  puts(buf);
    // printf("\n");
    // exit(1);

    // scanf("%s",buf);
    sha3_Init256(&c);
    sha3_Update(&c, (char *)buf, strlen(buf));
    hash = sha3_Finalize(&c);

    // j=0;

    for (int i = 0; i < 64; i++)
    {
        printf("%d", hash[i]);
        // char s[3];
        // byte_to_hex(hash[i],s);

        sk[i] ^= hash[i];
    }
    printf("\nencrypt sk=");
    for (int i = 0; i < 64; i++)
        printf("%d,", sk[i]);
    printf("\n");

    fp = fopen("enc.sk", "wb");
    fwrite(sy, 2, K, fp);
    fwrite(sk, 1, 64, fp);
    fclose(fp);
}

void decrypt(vec w)
{
    FILE *fp;
    int i, j;
    unsigned char sk[64] = {0}, err[N] = {
                                    0};
    unsigned short buf[K] = {0}, tmp[K] = {
                                     0};
    vec f = {0}, r = {0};
    vec v = {0};
    const uint8_t *hash = {0};
    sha3_context c = {0};
    int image_size = 512;

    j = 0;
    fp = fopen("enc.sk", "rb");

    fread(tmp, 2, K, fp);
    fread(sk, 1, 64, fp);
    fclose(fp);

    for (int i = 0; i < K; i++)
        buf[i] = tmp[K - i - 1];

    printf("in decrypt\n");
    f = setpol(buf, K);
    // v = pattarson(w, f);

    // elo(r);
    // exit(1);
    // v=(r);

    j = 0;
    if (v.x[1] > 0 && v.x[0] == 0)
    {
        err[0] = 1;
        j++;
    }

    printf("j=%d\n", j);
    printf("after j\n");
    for (int i = j; i < 2 * T; i++)
    {
        if (v.x[i] > 0 && v.x[i] < M)
        {
            err[v.x[i]] = 1;
        }
    }

    char buf0[8192] = {0}, buf1[10] = {
                               0};

    // #pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        snprintf(buf1, 10, "%d", err[i]);
        strcat(buf0, buf1);
    }
    // puts (buf0);
    printf("vector=%d\n", strlen(buf0));
    // exit(1);
    printf("cipher sk2=");
    for (int i = 0; i < 64; i++)
        printf("%u,", sk[i]);
    printf("\n");

    sha3_Init256(&c);
    sha3_Update(&c, (char *)buf0, strlen(buf0));
    hash = sha3_Finalize(&c);

    j = 0;
    printf("hash=");
    for (int i = 0; i < 64; i++)
    {
        printf("%d", hash[i]);
        // char s[3];
        // byte_to_hex(hash[i],s);

        sk[i] ^= hash[i];
    }
    printf("\ndecript sk=");
    for (int i = 0; i < 64; i++)
        printf("%u,", sk[i]);
    printf("\n");
    //  exit(1);

    return;
}

unsigned short vb[K * 2][N] = {0};
unsigned short gt[K * 2][K * 2] = {0};
void van(int kk)
{
    int i, j;

    printf("van der\n");
    /*
    for( int i= 0; i < N; i++){
        mat[i][0] = vb[0][i] = 1;
        printf("%d,", vb[0][i]);
    }
    printf("\n");
    */
    // #pragma omp parallel for private(i, j)
    for (int i = 0; i < kk; i++)
    {
        for (int j = 1; j < M; j++)
        {
            vb[i][j] = gf[mltn(i, fg[j])];
            printf("g%d,", vb[i][j]);
            mat[i][j] = vb[i][j];
        }
        printf("\n");
    }
    // exit(1);
}

void vanshe(int kk)
{
    int i, j, k;

    printf("van der\n");

    for (int i = 0; i < N; i++)
        vb[0][i] = 1;
    // #pragma omp parallel for private(i, j)
    for (int i = 1; i < kk; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mat[i][j] = gf[mltn(i, j)];
        }
        // printf("\n");
    }
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < N; j++)
            printf("%d,", mat[i][j]);
        printf("\n");
    }
    // exit(1);
}

void ogt(unsigned short pp[], int kk)
{
    int i, j, k;
    vec w = {0};

#pragma omp parallel for private(i, j)
    for (int i = 0; i < kk; i++)
    {
        for (int j = 0; j < kk - i; j++)
        {
            gt[i][j + i] = pp[j];
        }
    }
    for (int i = 0; i < kk; i++)
    {
        for (int j = 0; j < kk; j++)
            printf("%d,", gt[i][j]);
        printf("\n");
    }
    // exit(1);
}

vec synd(unsigned short zz[], int kk)
{
    unsigned short syn[K] = {0}, s = 0;
    int i, j, t1;
    vec f = {0};

    printf("in synd2\n");

    for (int i = 0; i < K; i++)
    {
        // syn[i] = 0;
        s = 0;
        // #pragma omp parallel num_threads(16)
        for (int j = 0; j < M; j++)
        {
            if (zz[j] > 0)
            {
                s = plus(s, gf[mlt(fg[zz[j]], fg[mat[i][j]])]);
                // printf("%d,%d %d %d", s, mat[j][i], zz[j], j);
                // printf(" ==sind\n");
            }
        }
        syn[K - 1 - i] = s;
        printf("syn%d,\n", syn[K - 1 - i]);
    }
    // printf ("\n");
    for (int i = 0; i < K; i++)
        printf("%d,", bexp(syn[i]));
    printf("\n");
    // exit(1);

    f = setpol(syn, K);
    printpol((f));
    printf(" syn=============\n");
    printf("%d %d\n", bexp(2), bexp(26));
    // exit(1);

    return f;
}

unsigned short tas(uni a, uni x)
{
    unsigned short c = {0};

    c = plus(xtrace(dick[fg[c2d(a)]], Pr), xtrace(dick[fg[c2d(x)]], Pr));
    return c;
}

vec sind(unsigned short zz[], int kk)
{
    // uni syn[K] = {0};
    unsigned short s[K * E] = {0};
    unsigned short syn[K * E] = {0};
    int i, j, t1;
    vec ww = {0};
    uni om = {0}, u[K] = {0};

    for (int i = 0; i < kk; i++)
    {
        t1 = 0;
        for (int j = 0; j < M; j++)
        {
            if (zz[j] > 0)
            {
                t1 += fls.x[i][j];
                printf("s=%d,", fls.x[i][j]);
            }
        }
        s[i] = t1 % Pr;
        // ww.x[i]=c2d(om);
        /*
        ww.x[i*E+0]=om.fugo.b0;
        ww.x[i*E+1]=om.fugo.b1;
        if(E==3)
        ww.x[i*E+2]=om.fugo.b2;
        */
        printf("w2 %d %d %d %d \n", s[i], om.fugo.b0, om.fugo.b1, om.fugo.b2);
        // s[i*E+0]=om.fugo.b0;
        // s[i*E+1]=om.fugo.b1;
        om.u = 0;
    }
    // exit(1);

    unsigned a = 0;
    unsigned s2[K * E] = {0};

    for (int i = 0; i < K * E; i++)
    {
        a = 0;
        for (int j = 0; j < K * E; j++)
        {
            a += (inv_S.x[i][j] * s[j]) % Pr;
            printf("%d,", inv_S.x[i][j]);
        }
        printf("\n");
        s2[i] = a % Pr;
    }

    for (int i = 0; i < K; i++)
    {
        u[i].fugo.b0 = s2[i * E + 0];
        u[i].fugo.b1 = s2[i * E + 1];
        u[i].fugo.b2 = s2[i * E + 2];
    }

    // memset(ww.x,0,sizeof(ww.x));
    for (i = 0; i < K; i++)
    {
        ww.x[i] = c2d(u[i]);
        printf("ww2 %d %d %d %d\n", ww.x[i], u[i].fugo.b0, u[i].fugo.b1, u[i].fugo.b2);
    }
    // exit(1);
    /*
    unsigned a=0;
    printf("\n");
    printf(" syn=============\n");
    //
*/
    return ww;
}

vec d2v(unsigned short d[K])
{
    vec v = {0};
    uni e = {0};

    for (int i = 0; i < K; i++)
    {
        e.u = d[i];
        v.x[i * E] = e.fugo.b0;
        v.x[i * E + 1] = e.fugo.b1;
        v.x[i * E + 2] = e.fugo.b2;
    }

    return v;
}

// 64バイト秘密鍵の暗号化と復号のテスト
void test(vec w, unsigned short zz[])
{
    int i;
    vec v = {0};
    const uint8_t *hash;
    sha3_context c;
    // int image_size=512;
    vec f = {0};
    FILE *fp;

    fp = fopen("aes.key", "rb");
    /*
     static char base64[] = {
     'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
     'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
     'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
     'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
     'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
     'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
     'w', 'x', 'y', 'z', '0', '1', '2', '3',
     '4', '5', '6', '7', '8', '9', '+', '/',
     };
   */
    char buf[8192] = {0}, buf1[10] = {
                              0};
    unsigned char sk[64] = {0};
    // unsigned short s[K]={0};
    // fread(sk,1,32,fp);
    for (int i = 0; i < 64; i++)
        sk[i] = i + 1;

    for (int i = 0; i < N; i++)
    {
        snprintf(buf1, 10, "%d", zz[i]);
        strcat(buf, buf1);
    }
    // puts (buf);
    printf("vector=%u\n", strlen(buf));
    // exit(1);

    printf("sk0=");
    for (int i = 0; i < 64; i++)
        printf("%u,", sk[i]);
    printf("\n");
    // exit(1);

    f = synd(zz, K);
    v = (f);
    // printf("v=");
    for (int i = 0; i < K; i++)
    {
        sy[i] = v.x[i];
        printf("%d,", sy[i]);
    }
    printf("\n");

    encrypt(buf, sk);
    decrypt(w);

    sha3_Init256(&c);
    sha3_Update(&c, (char *)buf, strlen(buf));
    hash = sha3_Finalize(&c);
}

void readkey()
{
    FILE *fp, *fq;
    unsigned short dd[K * N] = {0};
    int i, j;

    // 鍵をファイルに書き込むためにはkey2を有効にしてください。
    //   key2(g);
    fp = fopen("sk.key", "rb");
    fread(g, 2, K + 1, fp);
    fclose(fp);
    // 固定した鍵を使いたい場合はファイルから読み込むようにしてください。
    fq = fopen("H.key", "rb");
    fread(dd, 2, K * N, fq);
    // #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < K; j++)
            mat[j][i] = dd[K * i + j];
    }
    fclose(fq);
}

// vec sx={0},ty={0};
int isquad(vec w)
{
    int i, j, flg = 0;
    vec b = {0};

    b = (w);
    for (int i = 0; i < DEG; i++)
    {
        if (b.x[i] > 0 && i % 2 == 1)
            return 0;
    }

    return -1;
}

vec mkpol()
{
    int i, j, k, fail, flg, l, ii = 0;
    vec w = {0};
    if (K == 2)
    {
        ginit();
        w = setpol(g, K + 1);
        return w;
    }

    do
    {
        fail = 0;
        j = 0;
        k = 0;
        flg = 0;
        l = 0;
        memset(g, 0, sizeof(g));
        // memset(ta, 0, sizeof(ta));
        memset(w.x, 0, sizeof(w));
        ginit();
        ii++;
        if (ii > 100)
        {
            printf("erro=%d\n", ii);
            exit(1);
        }

        for (int i = 0; i < K; i++)
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

    printpol((w));
    printf(" ==g\n");
    // exit(1);

    return w;
}

unsigned short dd[N][N] = {0};

// Patterson & EEA 用（ランダム多項式、次元指定）
vec mkg(int kk)
{
    int i, j, k, l, ii = 0;
    vec w = {0};
    unsigned short tr[N] = {0};
    unsigned short ta[N] = {0};
    unsigned short po[K + 1] = {1, 0, 1, 0, 5};

aa:

    // printf("\n");

    // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
    // 既約多項式しか使わない。

    l = -1;
    ii = 0;
    /*
    while (l == -1)
    {
        w = mkpol();
        l = ben_or(w);
        printf("irr=%d\n", l);
        if (ii > 300)
        {
            printf("too many error\n");
            exit(1);
        }
        ii++;
        //
    }
    */
    w = mkpol();
    memset(ta, 0, sizeof(ta));
    // g[0]=1;

    // 多項式の値が0でないことを確認
    for (int i = 0; i < N; i++)
    {
        ta[i] = eval(w, i);
        if (ta[i] == 0)
        {
            printf("eval 0 @ %d\n", i);
            // fail = 1;
            goto aa;
        }
    }

    // 多項式を固定したい場合コメントアウトする。
    /*
  memset(ta, 0, sizeof(ta));
  w = setpol(po, K + 1);
  printpol((w));
  printf(" =poly\n");
//  exit(1);
  for( int i= 0; i < N; i++)
  {    ta[i] = trace(w, i);
    if (ta[i] == 0)
    {
      printf("trace 0 @ %d\n", i);
      exit(1);
    }
  }
  */
    printpol(w);
    printf("\n");
    printsage((w));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");

    for (int i = 0; i < N; i++)
    {
        tr[i] = oinv(ta[i]);
        // printf("%d,", tr[i]);
    }

    printpol((w));
    printf(" =irreducible\n");
    printsage((w));
    printf("\n");
    // wait();

    memset(vb, 0, sizeof(vb));
    memset(gt, 0, sizeof(gt));
    van(kk);
    ogt(po, kk);
    memset(mat, 0, sizeof(mat));

    //  wait();

    // #pragma omp parallel for

    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < K; i++)
        {
            ma[i][j] = gf[mlt(fg[vb[i][j]], tr[j])];
        }
        // printf("tr[%d]=%d\n",j,tr[j]);
    }

    unsigned short s;
#pragma omp parallel for
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < N; j++)
        {
            s = 0;

            for (int k = 0; k < K; k++)
                s ^= gf[mlt(fg[gt[k][i]], fg[ma[k][j]])];
            // printf("%d,",s);
            mat[i][j] = s;
        }
        printf("\n");
    }
    printf("\n");
    // exit(1);

    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < K; i++)
            printf("%d,", mat[i][j]);
        printf("\n");
    }
    printf("\n");
    // wait();

    return w;
}

// BMA 専用（多項式と次元指定）
vec mkc()
{
    int i, j, k, l, ii = 0;

    unsigned short tr[N] = {0};
    unsigned short ta[N] = {0};
    vec v = {0},w;
    unsigned short po[K + 1] = {1, 0, 1, 0, 5};
    // vec w={0};
    vec r = {0};

aa:

    // printf("\n");
    memset(mat, 0, sizeof(mat));
    // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
    // 既約多項式しか使わない。

    l = -1;
    ii = 0;
    // irreducible gvecpa code (既役多項式が必要なら、ここのコメントを外すこと。)

    while (l == -1)
    {
        w = mkpol();
        l = ben_or(w);
        printf("irr=%d\n", l);
        if (ii > 300)
        {
            printf("too many error\n");
            exit(1);
        }
        ii++;
        //
    }

    // separable gvecpa code
    // w = mkpol();
    // r = vmul(w, w);

    memset(ta, 0, sizeof(ta));
    // w = setpol(g, K + 1);
    printpol((r));
    // printf(" =poly\n");

    // 多項式の値が0でないことを確認
    for (int i = 0; i < N; i++)
    {
        ta[i] = eval(r, i);
        if (ta[i] == 0)
        {
            printf("eval 0 @ %d\n", i);
            // fail = 1;
            goto aa;
        }
    }
    for (int i = 0; i < N; i++)
    {
        tr[i] = oinv(ta[i]);
        // printf("%d,", tr[i]);
    }
    memset(g, 0, sizeof(g));
    // g[0] = 1;

    // 多項式を固定したい場合コメントアウトする。
    printpol(r);
    printf("\n");
    printsage((r));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");
    memset(v.x, 0, sizeof(v.x));
    //  v=rev(w);
    van(K);
    //  v=(w);
    ogt(r.x, K);

    // wait();

    // #pragma omp parallel for

    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < K; i++)
        {
            mat[i][j] = gf[mlt(fg[vb[i][j]], tr[j])];
        }
        // printf("tr[%d]=%d\n",j,tr[j]);
    }

    // printf("\n");
    // exit(1);
    /*
    for( int j = 0; j < N; j++)
    {
        for( int i= 0; i < kk; i++)
            printf("%d,", mat[j][i]);
        printf("\n");
    }
    //exit(1);
    //wait();
*/

    return w;
}

//
vec mkd(vec w, int kk)
{
    int i, j, k, l, ii = 0;

    unsigned short tr[N] = {0};
    unsigned short ta[N] = {0};
    vec v = {0};
    unsigned short po[K + 1] = {1, 0, 1, 0, 5};
    // vec w={0};
    vec r = {0};

aa:

    // printf("\n");
    memset(mat, 0, sizeof(mat));
    // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
    // 既約多項式しか使わない。

    l = -1;
    ii = 0;
    // irreducible gvecpa code (既役多項式が必要なら、ここのコメントを外すこと。)

    while (l == -1)
    {
        w = mkpol();
        l = ben_or(w);
        printf("irr=%d\n", l);
        if (ii > 300)
        {
            printf("too many error\n");
            exit(1);
        }
        ii++;
        //
    }

    // separable gvecpa code
    // w = mkpol();
    r = w;
    //  r=vmul(w,w);
    memset(ta, 0, sizeof(ta));
    // w = setpol(g, K + 1);
    printpol((r));
    // printf(" =poly\n");

    // 多項式の値が0でないことを確認
    for (int i = 0; i < N; i++)
    {
        ta[i] = eval(r, i);
        if (ta[i] == 0)
        {
            printf("eval 0 @ %d\n", i);
            // fail = 1;
            goto aa;
        }
    }
    for (int i = 0; i < N; i++)
    {
        tr[i] = oinv(ta[i]);
        // printf("%d,", tr[i]);
    }
    memset(g, 0, sizeof(g));
    g[0] = 1;

    // 多項式を固定したい場合コメントアウトする。
    printpol(r);
    printf("\n");
    printsage((r));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");
    memset(v.x, 0, sizeof(v.x));
    //  v=rev(w);
    van(kk);
    //  v=(w);
    ogt(g, kk);

    // wait();

    // #pragma omp parallel for

    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < kk; j++)
        {
            mat[j][i] = vb[j][i];
        }
    }

    // printf("\n");
    // exit(1);
    /*
    for( int j = 0; j < N; j++)
    {
        for( int i= 0; i < kk; i++)
            printf("%d,", mat[j][i]);
        printf("\n");
    }
    //exit(1);
    //wait();
*/

    return w;
}

void half(int kk)
{
    int i, j, k;

    for (int i = 0; i < N; i++)
        bm[0][i] = mat[0][i];
    for (int i = 1; i < kk; i++)
    {
        for (int j = 0; j < N; j++)
            bm[i][j] = mat[i * 2 - 1][j];
    }
    /*
    for( int i= 0; i < N; i++)
    {
        for( int j = 0; j < kk+1; j++)
            printf("%d,", bm[i][j]);
        printf("  ==bm\n");
    }
*/
    // exit(1);
}

static CTX transi(CTX H)
{
    CTX X = {0};

    int count = 0;
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < M; j++)
        {
            big.d[i][j].fugo.b0 = H.x[i * E + 0][j];
            big.d[i][j].fugo.b1 = H.x[i * E + 1][j];
            big.d[i][j].fugo.b2 = H.x[i * E + 2][j];
            // printf("%d ",big.d[i][j].fugo.b2);
        }
        printf("\n");
    }
    // exit(1);

    printf("\n");
    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < M; j++)
            printf("%d ", big.d[i][j].fugo.b0);
        printf("\n");
        for (int j = 0; j < M; j++)
            printf("%d ", big.d[i][j].fugo.b1);
        printf("\n");
        for (int j = 0; j < M; j++)
        {
            // big.d[i][j].fugo.b2 = 0;
            printf("%d ", big.d[i][j].fugo.b2);
        }

        printf("\n");
    }
    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < M; j++)
        {
            X.x[i][j] = c2d(big.d[i][j]);
            printf("%d ", X.x[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    // exit(1);
    /*
    //R = toByte(O_bin, K);
    for( int i= 0; i < N; i++)
    {
        for( int j = 0; j < K; j++)
        {
            R_bin.x[i][j] = O_bin.x[i][j];
            // printf("%d,",R_bin.x[i][j]);
        }
    }
    */
    // exit(1);

    return X;
}

vec zind(unsigned short zz[], int kk)
{
    // uni syn[K] = {0};
    unsigned short s[K * E] = {0};
    int i, j, t1;
    vec ww = {0};
    uni om = {0};

    for (int i = 0; i < kk; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (zz[j] > 0)
            {
                om = d2c(tas(om, big.d[i][j]));
            }
        }
        ww.x[i] = c2d(om);
        printf("w2%d,%d %d %d\n", ww.x[i], om.fugo.b0, om.fugo.b1, om.fugo.b2);
        om.u = 0;
    }
    printf("\n");
    printf(" syn=============\n");
    // exit(1);

    return ww;
}

void printuni(uni a)
{
    printf("%d %d %d\n", a.fugo.b0, a.fugo.b1, a.fugo.b2);
}

// Niederreiter暗号の公開鍵を作る(Gvecpa)
void pk_gen()
{
    int i, j, k, l;
    FILE *fp;
    unsigned char dd[E * K] = {0};
    vec w = {0};
    static CTX Q = {0}, O_bin = {0};
    // R = {0}, R_bin = {0},
    unsigned short zz[M] = {0};
    vv(K);
    //vv(K);
    // zz[1] = 1;
    // zz[2] = 1;
   
    mkerr(zz, T);
    vec x = synd(zz, K);
    ymo m = bm_itr(x.x);
    chen(m.f);
    for(i=0;i<N;i++)
    if(zz[i]>0)
    printf("i=%d %d\n",i,zz[i]);
    //exit(1);
    
    //  w = mkd(w, K);
    //   w = mkg(K);
    //   half(K / 2 + 1);

    printpol(w);
    printf("\n");
    printsage((w));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");

    bdet();
    CTX D = {0};
    uni test[M] = {0};

    // exit(1);
    unsigned short ww[K] = {0}, nn[K] = {0};
    uni om = {0};
    static CTX H = {0};

    memset(zz, 0, sizeof(zz));
    // zz[1] = 1;
    // zz[3] = 1;
    /*
    mkerr(zz, T);
    // exit(1);
    x = zind(zz, K);
    ymo tt = bm_itr(x.x);
    chen(tt.f);
    // exit(1);
    */
    vec v = {0};
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < M; j++)
        {
            D.x[i * E + 0][j] = big.d[i][j].fugo.b0;
            D.x[i * E + 1][j] = big.d[i][j].fugo.b1;
            if (E == 3)
                D.x[i * E + 2][j] = big.d[i][j].fugo.b2;
        }
    }
    Pgen();
    printf("pre\n");
inv_S.b=0;
    
    //do
    {
        memset(inv_S.x, 0, sizeof(inv_S.x));
        memset(S.x, 0, sizeof(S.x));
        
        for (int i = 0; i < K * E; i++)
        {
            for (int j = 0; j < K * E; j++)
                S.x[i][j] = random() % Pr;
        }
        inv_S = genS(S);
    } //while (inv_S.b == 0);

    Pgen();
    for (int i = 0; i < M; i++)
        printf("i%d,", P[i]);
    printf("\n");
    CTX temp = {0};
    temp = kenzan(S, D, 2);

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < K * E; j++)
            fls.x[j][i] = temp.x[j][P[i]];
    }
    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < M; j++)
        {
            printf("%d,", fls.x[i][j]);
            // printuni((H.d[i][j]));
        }
        printf("\n");
    }
    printf("\n\n");

    H = kenzan(inv_S, fls, 2);
    memset(zz, 0, sizeof(zz));
    int y;
    vec z = {0};
    uni in[K] = {0};
    int a = 0;
    unsigned short ss[M] = {0};

     mkerr(zz,T);
    z = sind(zz, K * E);
    ymo o = bm_itr(z.x);
    v = chen(o.f);
    unsigned short rr[M] = {0};
    for (int i = 0; i < T; i++)
        ss[v.x[i]] = 1;
    for (int i = 0; i < M; i++)
        if (ss[i] > 0)
            printf("%d,", i);
    printf("\n");
    for (int i = 0; i < M; i++)
        rr[i] = ss[P[i]];
    for (int i = 0; i < M; i++){
        if (rr[i] > 0)
            printf("%d,", i);
    }
    printf("\n");
    for (int i = 0; i < M; i++){
        if (zz[i] > 0)
            printf("%d,", i);
    }
    printf("\n");
    exit(1);

    return;
}

vec dec(unsigned short ss[])
{
    int i, j, k;
    vec v = {0};
    vec s = {0};
    uni d[K] = {0};
    unsigned ch[K * E] = {0};
    uni h2o[K] = {0};

    printf("!1\n");
    for (int i = 0; i < K; i++)
    {
        // v = i2v(ss[i]);
        // for( int j = 0; j < E; j++)
        d[i] = d2c(ss[i]);
    }
    for (int i = 0; i < K * E; i++)
        printf("%d", ch[i]);
    printf("\n");

    unsigned short uk[K] = {0};

    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < K; j++)
        {
            h2o[i].fugo.b0 += (d[j].fugo.b0 * inv_S.x[i + E][j]);
            h2o[i].fugo.b1 += (d[j].fugo.b1 * inv_S.x[i * E + 1][j]);
            h2o[i].fugo.b2 += (d[j].fugo.b2 * inv_S.x[i * E + 2][j]);
        }
    }

    // for( int i= 0; i < K * E; i++)
    // printf("%d,", h2o[i]);
    // printf("\n");
    memset(v.x, 0, sizeof(v.x));
    for (int i = 0; i < K; i++)
    {
        v.x[j] = h2o[i].u;
    }
    for (int i = 0; i < K; i++)
        printf("%d,", v.x[i]);
    printf("\n");
    //    exit(1);

    return v;
}

// 鍵生成
void key2(unsigned short g[])
{
    FILE *fp;
    unsigned short dd[K] = {0};
    int i, j, k;

    printf("鍵を生成中です。４分程かかります。\n");
    fp = fopen("H.key", "wb");
    i = 0;

    mkg(K);

    // exit(1);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < K; j++)
            dd[j] = mat[j][i];
        fwrite(dd, 2, K, fp);
    }
    fclose(fp);
    fp = fopen("sk.key", "wb");
    fwrite(g, 2, K + 1, fp);
    fclose(fp);
}

int printerr(vec r)
{
    int count, i, j, k;

    count = 0;

    for (int i = 0; i < T; i++)
    {

        if (i > 0 && i == 0)
        {
            printf("err baka-z\n");
            // return -1;
            exit(1);
        }

        if (r.x[i] > 0 && i > 0) // == r.x[i].n)
        {
            printf("err=%d %d %s\n", r.x[i], i, "お");
            count++;
        }
        if (i == 0)
        {
            printf("\nerr=%d %d %s\n", r.x[i], i, "う");
            count++;
        }
        // zz[r.x[i].n]=r.x[i];
    }
    if (count != T)
    {
        printf("error pattarn too few %d\n", count);
        return -1;
    }
    // exit(1);

    return count;
}

int elo2(vec r)
{
    int count, i, j, k;

    count = 0;

    unsigned char x[N] = {0}, yy[N] = {0};
    for (int i = 0; i < T; i++)
    {

        if (i > 0 && i == 0)
        {
            printf("err baka-z\n");
            count++;
            // return -1;
            // exit(1);
        }

        if (r.x[i] > 0 && i > 0) // == r.x[i].n)
        {
            x[i] = r.x[i];
            // printf("err=%d %d %s\n", r.x[i], r.x[i].n, "お");
            count++;
        }
        if (i == 0)
        {
            x[i] = r.x[i];
            // printf("\nerr=%d %d %s\n", r.x[i], r.x[i].n, "う");
            count++;
        }
        // zz[r.x[i].n]=r.x[i];
    }
    if (count < T)
    {
        printf("err is too few\n");
        exit(1);
    }
    for (int i = 0; i < N; i++)
        yy[i] = x[P[i]];
    for (int i = 0; i < N; i++)
    {
        if (yy[i] > 0)
            printf("err= %d\n", i);
    }

    return count;
}

int ero(vec v)
{
    int i, j, count = 0;

    for (int i = 0; i < T * 2; i++)
    {
        if (i == 0)
        {
            printf("error position=%d %d う\n", i, v.x[i]);
            count++;
        }
        if (i > 0 && v.x[i] > 0)
        {
            printf("error position=%d %d お\n", i, v.x[i]);
            count++;
        }
        if (i > 0 && v.x[i] == 0)
        {
            printf("baka %d %d\n", i, v.x[i]);
            printf("v.x[K-1]=%d\n", v.x[K - 1]);
            /*
        printpol((w));
        printf(" ============gvecpa\n");
        printsage((w));
        printf(" ============sage\n");
        printsage((f));
        printf(" ============syn\n");
        printpol((f));
        printf(" ==========synd\n");
        printf("{");
        for (int k= 0; k < N; k++)
        {
          if (z1[k] > 0)
            printf("%d,", z1[k]);
        }
        printf("};\n");
        //AA++;
        //wait();
        */
            break;
            //
            // exit (1);
        }
        int cnt = 0;
        /*
      for (int k= 0; k < N; k++)
      {
        if (z1[k] > 0)
        {
          if (k != v.x[cnt])
          {
            printf("%d,%d\n", k, v.x[cnt]);
            printsage((w));
            printf(" ========w\n");
            AA++;
            break;
            //exit(1);
          }
          cnt++;
        }
      }
*/
    }

    if (count == T * 2)
    {
        printf("err=%dっ!! \n", count);
    }
    if (count < T * 2)
    {
        printf("error is too few\n");

        // memcpy (zz, z1, sizeof (zz));
        /*
      printf("{");
      for( int i= 0; i < D; i++)
        printf("%d,", z1[i]);
      printf("};\n");
      printpol((w));
      printf(" =========gvecpa\n");
      printsage((w));
      printf(" =========sage\n");
      printsage((f));
      printf(" =========syn\n");
      printpol((f));
      printf(" ==========synd\n");
      */
        printf("へげえええーっ\n");
        // exit(1);
        exit(1);
    }

    return count;
}

int ero2(vec v)
{
    int i, j, count = 0;
    unsigned short ya[N] = {0}, xa[N] = {0};

    for (int i = 0; i < T; i++)
    {
        if (i == 0)
        {
            xa[v.x[i]] = 1;
            // printf("error position=%d %d う\n", i, v.x[i]);
            count++;
        }
        if (i > 0 && v.x[i] > 0)
        {
            xa[v.x[i]] = 1;
            // printf("error position=%d %d お\n", i, v.x[i]);
            count++;
        }
        if (i > 0 && v.x[i] == 0)
        {
            printf("baka %d %d\n", i, v.x[i]);
            printf("v.x[K-1]=%d\n", v.x[K - 1]);
            break;
            //
            // exit (1);
        }
    }

    int cnt = 0;
    for (int i = 0; i < M; i++)
        ya[i] = xa[P[i]];
    for (int i = 0; i < M; i++)
    {
        if (ya[i] > 0 && i == 0)
        {
            printf("error position=%d う\n", i);
            cnt = 1;
        }
        else if (ya[i] > 0)
        {
            if (cnt == 0)
            {
                printf("error position=%d う\n", i);
                cnt = 1;
            }
            else
            {
                printf("error position=%d お\n", i);
            }
        }
    }
    // exit(1);

    if (count == T)
    {
        printf("err=%dっ!! \n", count);
    }
    if (count < T)
    {
        printf("error is too few\n");

        // memcpy (zz, z1, sizeof (zz));
        /*
      printf("{");
      for( int i= 0; i < D; i++)
        printf("%d,", z1[i]);
      printf("};\n");
      printpol((w));
      printf(" =========gvecpa\n");
      printsage((w));
      printf(" =========sage\n");
      printsage((f));
      printf(" =========syn\n");
      printpol((f));
      printf(" ==========synd\n");
      */
        printf("へげえええーっ\n");
        // exit(1);
        exit(1);
    }

    return count;
}

int ero3(vec v)
{
    int i, j, count = 0;
    unsigned short ya[N] = {0}, xa[N] = {0};

    for (int i = 0; i < T; i++)
    {
        if (i == 0)
        {
            xa[v.x[i]] = 1;
            // printf("error position=%d %d う\n", i, v.x[i]);
            count++;
        }
        if (i > 0 && v.x[i] > 0)
        {
            xa[v.x[i]] = 1;
            // printf("error position=%d %d お\n", i, v.x[i]);
            count++;
        }
        if (i > 0 && v.x[i] == 0)
        {
            printf("baka %d %d\n", i, v.x[i]);
            printf("v.x[K-1]=%d\n", v.x[K - 1]);
            break;
            //
            // exit (1);
        }
    }

    int cnt = 0;
    for (int i = 0; i < N; i++)
        ya[i] = xa[i];
    for (int i = 0; i < N; i++)
    {
        if (ya[i] > 0 && i == 0)
        {
            printf("error position=%d う\n", i);
            cnt = 1;
        }
        else if (ya[i] > 0)
        {
            if (cnt == 0)
            {
                printf("error position=%d う\n", i);
                cnt = 1;
            }
            else
            {
                printf("error position=%d お\n", i);
            }
        }
    }
    // exit(1);

    if (count == T)
    {
        printf("err=%dっ!! \n", count);
    }
    if (count < T)
    {
        printf("error is too few\n");

        // memcpy (zz, z1, sizeof (zz));
        /*
      printf("{");
      for( int i= 0; i < D; i++)
        printf("%d,", z1[i]);
      printf("};\n");
      printpol((w));
      printf(" =========gvecpa\n");
      printsage((w));
      printf(" =========sage\n");
      printsage((f));
      printf(" =========syn\n");
      printpol((f));
      printf(" ==========synd\n");
      */
        printf("へげえええーっ\n");
        // exit(1);
        exit(1);
    }

    return count;
}

void mkerr(unsigned short *z1, int num)
{
    int j, l;

    j = 0;

    memset(z1, 0, sizeof(z1[N]));

    while (j < num)
    {
        l = rand() % (M);
        // printf ("l=%d\n", l);
        if (0 == z1[l])
        {
            z1[l] = 2;
            // printf("l=%d\n", l);
            j++;
        }
    }
}

void fun()
{
    unsigned short i, k;

    vec qq = {0};
    for (int i = 0b10001; i < 0b11111 + 1; i++)
    {
        qq = (i2v(i));
        k = ben_or(qq);
        if (k == 0)
        {
            printpol((qq));
            printf(" =irreducible\n");
        }
    }
}

vec newhalf(unsigned short e[K])
{
    int i, j, k;
    vec v = {0};
    unsigned short t[K] = {0};

    for (int i = 0; i < K / 2 + 1; i++)
    {
        t[i] = e[i];
        printf("e=%d\n", e[i]);
    }
    for (int i = 0; i < K / 2 + 1; i++)
    {
        printf("t=%d\n", t[i]);
    }
    // exit(1);

    v.x[0] = t[0];
    v.x[1] = t[1];
    k = 2;
    for (int i = 2; i < K; i++)
    {
        if (i % 2 == 1)
        {
            v.x[i] = t[k];
            k++;
        }
        if (i % 2 == 0)
            v.x[i] = gf[mlt(fg[v.x[i / 2]], fg[v.x[i / 2]])];
    }

    return v;
}
/*
vec bfd(unsigned short ss[])
{
    int i, j, k, count = 0;
    vec v = {0};
    vec s = {0};
    unsigned ch[K * E * 2] = {0};
    unsigned char h2o[K * E * 2] = {0};

    // count=(K/2+1)*E-1;
    for( int i= 0; i < (K / 2) + 1; i++)
    {
        v = i2v(ss[i]);
        for( int j = 0; j < E; j++)
        {
            ch[count] = v.x[j];
            count++;
        }
    }
    printf("bfd_bin=\n");
    for( int i= 0; i < (K / 2 + 1) * E; i++)
        printf("%d", ch[i]);
    printf("\n");
    // exit(1);

    unsigned short uk[K] = {0};

    for( int i= 0; i < (K / 2 + 1) * E; i++)
    {
        for( int j = 0; j < (K / 2 + 1) * E; j++)
            h2o[i] ^= (ch[j] & inv_S.x[i][j]);
    }

    printf("deco_bin=\n");
    for( int i= 0; i < (K / 2 + 1) * E; i++)
        printf("%d,", h2o[i]);
    printf("\n");

    // count=(K/2+1)*E-1;

    count = 0;
    for( int i= 0; i < (K / 2) + 1; i++)
    {
        memset(v.x, 0, sizeof(v.x));
        for( int j = 0; j < E; j++)
        {
            v.x[j] = h2o[count];
            count++;
        }
        uk[i] = v2i(v);
    }

    printf("bm_int=\n");
    for( int i= 0; i < K / 2 + 1; i++)
        printf("%d,", bm[1][i] ^ bm[2][i] ^ bm[3][i]);
    printf("\n");

    printf("u_int=\n");
    for( int i= 0; i < K / 2 + 1; i++)
        printf("%d,", uk[i]);
    printf("\n");
    // exit(1);

    s = setpol(uk, K / 2 + 1);
    v = (s);

    return v;
}
*/

unsigned short logx(unsigned short u)
{
    unsigned short i;
    if (u == 0)
        return 1;
    for (int i = 1; i < O; i++)
        if (gf[i] == u)
            return i - 1;

    printf("baka-von\n");
}

vec kof2(unsigned short c, vec f)
{
    int i, j, k;
    vec b = {0}, h = {0};
    vec g = {0};

    c = fg[c];
    printf("c=%d\n", c);
    // exit(1);
    b = f; // (f);
    k = deg(b);
    printpol((b));
    printf(" =b debugi\n");
    for (int i = 0; i < k + 1; i++)
    {
        h.x[i] = gf[mlt(c, fg[b.x[i]])];
    }
    // g = (h);
    printpol((h));
    printf(" =h in kof2\n");
    return h;
}

unsigned short minus(unsigned short a)
{
    int i;
    vec c = dick[fg[a]];
    vec b = {0}, d = (c);
    int k = deg((c)) + 1;
    unsigned short u = 0;
    vec f;
    printf("a=%d\n", a);
    for (int i = 0; i < k; i++)
    {
        b.x[i] = (Pr - d.x[i]) % Pr;
        printf("chen_b=%d\n", b.x[i]);
    }
    u = eval((b), Pr);

    return u;
}

ymo bm_itr(unsigned short s[])
{
    vec U1[2][2] = {0}, U2[2][2][2] = {0}, null = {0};
    int i, j, k;
    ymo t = {0};

    U2[0][0][0].x[0] = 1;        // f[0];
    U2[0][0][1].x[0] = 0;        // fai[0];
    U2[0][1][0].x[0] = 0;        // g[0];
    U2[0][1][1].x[0] = minus(1); // thi[0];
    int m = 0, d = 0, p = 2 * d - m - 1, myu = 0;
    printf("m=%d d=%d myu=%d p=%d\n", m, d, myu, p);
    for (m = 0; m < K; m++)
    {
        d = deg(U2[0][0][0]);
        p = 2 * d - m - 1;
        myu = 0;
        for (int i = 0; i <= d; i++)
            myu = plus(myu, gf[mlt(fg[U2[0][0][0].x[i]], fg[s[i + (m - d)]])]);

        printf("m=%d ad=%d myu=%d p=%d\n", m, d, myu, p);
        memset(U1, 0, sizeof(U1));
        if (myu == 0 || p >= 0)
        {
            U1[0][0].x[0] = 1;
            U1[0][1].x[p] = minus(myu);
            U1[1][0].x[0] = 0;
            U1[1][1].x[0] = 1;
            // exit(1);
        }
        else if (myu > 0 && p < 0)
        {
            if (p < 0)
            {
                p = -1 * (p);
            }
            U1[0][0].x[p] = 1;
            U1[0][1].x[0] = minus(myu);
            U1[1][0].x[0] = gf[oinv(myu)];
            U1[1][1].x[0] = 0;
        }
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                    U2[1][i][j] = (vadd2((U2[1][i][j]), (pmul(U1[i][k], U2[0][k][j]))));
            }
        }
        memcpy(U2[0], U2[1], sizeof(U2[0]));
        memset(U2[1], 0, sizeof(U2[1]));
    }
    t.f = U2[0][0][0];
    t.g = U2[0][1][0];
    t.h = U2[0][0][1];
    if (deg(t.f) == T)
    {
        printsage((t.f));
        printf(" ==chen00\n");
        return t;
    }
    else
    {
        t.f = U2[1][0][0];
        printsage((t.f));
        printf("baka\n");
        exit(1);
    }
}

vec lag(vec v)
{
    vec w[T] = {0}, t, l, tmp;

    // #pragma omp parallel for private(i,j)
    for (int i = 0; i < T; i++)
    {
        t = (tmp);
        //
        for (int j = 0; j < T; j++)
        {
            // #pragma omp critical
            if (i != j)
            {
                t = (vmul((w[i]), (w[j])));
                printf("wa=%d %d\n", i, j);
                printpol(t);
                printf("  ==t\n");
            }
            l = vadd(l, t);
            printpol(l);
            printf(" ===l\n");
        }

        printpol((l));
        printf(" ==snach\n");
        // exit(1);
        if (deg((l)) == 0)
        {
            printf("baka9\n");
            // exit(1);
        }
    }
}

vec bib(vec x)
{
    vec f = {0};

    for (int i = 1; i < deg(x) + 1; i++)
        f.x[i - 1] = gf[mlt(fg[i], fg[x.x[i]])];

    return f;
}

// 多項式の形式的微分
vec bibun_old(vec a)
{
    vec w[T * 2] = {0};
    vec l = {0}, t = {
                     0};
    int i, j, n, id;
    vec tmp = {0};

    n = deg(a);
    printf("n=%d\n", n);
    if (n == 0)
    {
        printf("baka8\n");
        //  exit(1);
    }

    //
    for (int i = 0; i < T; i++)
    {
        w[i].x[0] = minus(a.x[i]);
        w[i].x[1] = 1;
        printpol((w[i]));
        // printf("d=%d\n",a.x[i]);
        printf("\n");
    }
    tmp.x[0] = 1;
    //
    // exit(1);
    // #pragma omp parallel for private(i,j)
    for (int i = 0; i < T; i++)
    {
        t = (tmp);
        //
        for (int j = 0; j < T; j++)
        {
            // #pragma omp critical
            if (i != j)
            {
                t = (pmul((t), (w[j])));
                printf("wa=%d %d\n", i, j);
                printpol(t);
                printf("  ==t\n");
            }
            printpol(l);
            printf(" ===l\n");
        }
        l = vadd2(l, t);
        printpol((l));
        printf(" ==jaak\n");
        // exit(1);
        if (deg((l)) == 0)
        {
            printf("baka9\n");
            // exit(1);
        }
    }
    // exit(1);
    return (l);
}

vec yes(ymo t)
{
    vec e = vmul(t.g, t.h);
    vec b2 = bib((t.f));
    vec b = bibun_old((t.f));
    printsage((t.f));
    printf("\n");
    printsage((b));
    printf("\n");
    vec ixi = {0};
    int count = 0;
    for (int i = 0; i < N; i++)
    {
        if (eval((e), i) == 1)
        // if(a[i]>0)
        {
            printf("ixi=%d %d\n", gf[oinv(gf[mlt(fg[eval((b), i)], fg[eval((t.g), i)])])], i); // gf[mlt(oinv(eval((b),i)),fg[eval((t.h),i)])],i);
            ixi.x[i] = gf[mlt(oinv(eval((b), i)), fg[eval((t.h), i)])];
        }
    }

    return ixi;
}

vec bms(unsigned short s[])
{
    int L = 0, m = -1, d[K] = {0}, k = 0, i, e;
    vec f = {0}, g = {0}, h, v, c = {0}, cc = {0};

    // c.x[1] = 1;
    f.x[0] = g.x[0] = 1;
    d[m] = 1;
    printf("\nUh! %d\n", s[0]);
    // exit(1);

    while (k <= (2 * T - 1))
    {
        e = s[k];
        for (int i = 1; i < L; i++)
            e = plus(e, gf[mlt(fg[f.x[i]], fg[s[k - i]])]);
        printf("L=%d i=%d\n", L, i);
        d[k] = e; // plus(s[k],e); // plus(gf[mlt(fg[f.x[L]], fg[s[k - i]])], e); // s[k] ^ e;
        printf("d[%d]=%d s[k]=%d\n", k, d[k], s[k]);

        // printpol((f));
        printf("%d,%d====f\n", k, d[k]);
        if (d[k] > 0)
        {
            h = f;
            memset(v.x, 0, sizeof(v.x));
            v.x[k - m] = 1;
            unsigned short a;
            a = (m < 0) ? 1 : oinv(d[m]);
            cc = (kof2(gf[mlt(fg[d[k]], a)], (g)));
            f = (vadd((f), (vmul((cc), (v)))));
            if (L <= k / 2)
            {
                L = k + 1 - L;
                m = k;
                g = h;
            }
        }
        k++;
        // if(deg(f)==T)
        // break;
    }
    // exit(1);
    int j;
    vec ff = {0};
    if (deg(f) == T + 1)
    {
        for (int i = 0; i < T; i++)
            ff.x[i + 1] = f.x[i];
        f = ff;
        printpol((f));
        printf("  ===f\n");
        return f;
    }
    // for(i=0;i<O;i++){
    //     printf("ooo %d,",eval((f),i));
    // }

    if (deg(f) != T)
    {
        printpol((f));
        printf("\nI am a baka\n");
        exit(1);
    }

    return f;
}

vec bma(unsigned short s[], int kk)
{
    int i, j, k, ll = 0, l, d[2 * K + 1] = {0};
    vec lo[2 * K + 1] = {0}, b[2 * K + 1] = {0}, t[2 * K + 1] = {0}, a = {0}, f = {0}, h = {0}, g = {0}, hh = {0};
    vec v = {0}, x = {0}, w = {0};

    x.x[1] = 1;
    h = (x);
    v.x[0] = 1;
    f = (x);
    lo[0] = (v);
    b[0] = lo[0];
    ll = 0;
    for (int j = 1; j < T * 2 + 1; j++)
    {
        v = (lo[j - 1]);
        k = 0;
        // printpol(v);
        // printf(" ==lo\n");

        l = deg((lo[j - 1]));
        for (int i = 1; i < l + 1; i++)
        {
            k = plus(k, gf[mlt(fg[v.x[i]], fg[s[j - i]])]);
            // printf("v[%d]=%d\n", i, v.x[i]);
        }
        d[j] = plus(s[j], k);
        // printf("d[%d]=%d\n", j, d[j]);
        if (d[j] == 0)
        {
            lo[j] = lo[j - 1];
            b[j] = vmul(b[j - 1], h);
            // ll=j-1;
        }
        else // if (d[j] > 0)
        {
            g = vmul((kof(d[j], (h))), b[j - 1]);
            t[j] = (vadd((lo[j - 1]), (g)));
            lo[j] = t[j];
            if (ll * 2 > (j - 1))
            {
                // lo[j]=t[j];
                b[j] = vmul(b[j - 1], h);
            }
            else // if(2*ll <= j)
            {
                // printpol((t[j]));
                // printf("==t[%d]\n", j);
                b[j] = (kof(gf[oinv(d[j])], (lo[j - 1])));
                // lo[j]=t[j];
                ll = j - ll;

                if (j == 2 * T)
                {
                    if (!(d[T * 2 - 1] == 0 && d[T * 2 - 3] == 0 && deg(lo[j - 1]) == T) || !(deg(lo[j - 1]) == T))
                    {
                        if ((d[T * 2 - 1] == 0 && deg(lo[j - 2]) == T - 1))
                        {
                            lo[j - 1] = vmul(lo[j - 2], h);
                            // printpol((lo[j - 1]));
                            // printf("\n");
                        }
                    }
                    break;
                }
            }
        }
        printf("l=%d\n", ll);
        k = 0;
        // printpol((b[j]));
        // printf(" ==b[%d]\n", j);
    }

    k = 0;
    int count = 0;
    // printpol((lo[j - 1]));
    // printf(" ==coef\n");
    if (deg(lo[j - 1]) == T)
    {
        x = chen((lo[j - 1]));
    }
    else
    {
        printf("baka\n");
        exit(1);
        // return -1;
    }
    // exit(1);
    for (int i = 0; i < deg(x) + 1; i++)
    {
        if (x.x[i] >= 0)
        {
            printf("xx[%d]=1\n", (x.x[i]));
            count++;
        }
        //

        if (x.x[i] == 0)
            k++;
        if (k > 1)
        {
            printf("baka0\n");
            // printvec((x));
            // for( int i= 0; i < N; i++)
            // printf("%d,", zz[i]);
            exit(1);
            // return f;
        }
    }
    if (count < T)
    {
        printf("vaka in bms %d\n", count);
        exit(1);
    }

    // return count;
    return (lo[j - 1]);
}

vec rev(vec f)
{
    unsigned short i, tmp, j = 0, c[512] = {0}, d[512] = {0}, count = 0;
    vec v = {0}, x = {0};
    vec w = {0};

    v = (f);

    j = deg(f) + 1;
    printf("d=");
    for (int i = 0; i < j; i++)
    {
        // d[f.x[i].n]=f.x[i];
        // c[count]=f.x[i];
        d[j - 1 - i] = v.x[i];
        printf("%d,", v.x[i]);
        count++;
    }
    printf("\n");
    printf("c=");
    // memset(v.x,0,sizeof(v.x));
    for (int i = 0; i < count; i++)
    {
        x.x[i] = d[i];
        printf("%d,", d[i]);
    }
    printf("\n");
    // for(i=0;i<count;i++)
    // v.x[d[count-i-1]]=d[i];
    w = setpol(v.x, K + 1);
    // v=(w);
    printpol((w));
    printf(" ==rev?\n");
    // exit(1);
    // f=(v);

    return x;
}

vec sendrier(unsigned short zz[N], int kk)
{
    unsigned short syn[K / 2 + 1] = {0}, s = 0, rt[K * 3] = {0};
    int i, j, k;
    vec f = {0};
    vec v = {0}, x[K * 2] = {0};

    for (int j = 0; j < N; j++)
    {
        if (zz[j] > 0)
        {
            // memcpy(syn, bm[j], sizeof(syn));
            printf("bm_in_sen= %d || ", j);
            for (int i = 0; i < K / 2 + 1; i++)
            {
                syn[i] ^= bm[i][j];
                printf("%d,", syn[i]);
            }
        }
    }
    printf("\n");
    v = newhalf(syn);
    // printf("%d\n",j);
    for (int k = 0; k < kk; k++)
        rt[k] = v.x[k];
    //}
    // exit(1);
    // printf ("%d\n", j);
    // printf ("\n");
    //}

    f = setpol(rt, kk);
    // printpol((f));
    // printf(" syn=============\n");

    return f;
}

unsigned short bexp(unsigned short a)
{
    int i;
    /*
    if (a == 0)
        return 0;
    if (a == 1)
        return 1;
    */
    for (int i = 0; i < N; i++)
    {
        if (gf[i] == a)
            return i;
    }
}

// #define NN 16
vec sol(MTX a)
{
    unsigned int p, d;
    // int i, j, k;
    vec v = {0};

    for (int i = 0; i < K / 2; i++)
    {
        p = a.x[i][i];

        for (int j = 0; j < (K / 2 + 1); j++)
        {
            a.x[i][j] = gf[mlt(fg[a.x[i][j]], oinv(p))];
        }

        for (int j = 0; j < K / 2; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (int k = i; k < (K / 2 + 1); k++)
                {
                    if (a.x[j][k] > gf[mlt(fg[d], fg[a.x[i][k]])])
                    {
                        a.x[j][k] = plus(a.x[j][k], minus(gf[mlt(fg[d], fg[a.x[i][k]])]));
                    }
                    else
                    {
                        a.x[j][k] = (plus(a.x[j][k], minus(gf[mlt(fg[d], fg[a.x[i][k]])])));
                    }
                }
            }
        }
    }
    for (int i = 0; i < K / 2; i++)
    {
        for (int j = 0; j < K / 2 + 1; j++)
            printf("@%d,", a.x[i][j]);
        printf("\n");
        if (a.x[i][i] != 1)
        {
            for (int j = 0; j < K / 2 + 1; j++)
                printf("a%d,", a.x[i][j]);
            printf("\n");
        }
    }
    printf("\n");
    vec x = {0};
    for (int i = 0; i < K / 2; i++)
    {
        if (N > a.x[i][K / 2])
        {
            x.x[K / 2 - i] = minus(a.x[i][K / 2]);
        }
        else
        {
            x.x[K / 2 - i] = a.x[i][K / 2] % N;
        }
        // x.x[i+1] = N-a.x[i][K / 2];
    }

    x.x[0] = 1;

    vec vv = {0};
    vec pol = {0};
    pol = setpol(x.x, K / 2 + 1);
    printpol((pol));
    printf(" ==key\n");
    for (int i = 0; i < M; i++)
    {
        // v.x[i] = 0;
        if (xtrace(pol, i) == 0)
        {
            printf("i=%d\n", i);
            vv.x[i] = 1;
        }
    }
    // exit(1);

    return vv;
}

// #define NN 16
vec renri2(MTX a)
{
    unsigned short p, d;
    // int i, j, k;
    vec v = {0};

    for (int i = 0; i < K; i++)
    {
        p = a.x[i][i];

        for (int j = 0; j < (K + 1); j++)
        {
            a.x[i][j] = gf[mlt(fg[a.x[i][j]], oinv(p))];
        }

        for (int j = 0; j < K; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (int k = i; k < (K + 1); k++)
                {
                    a.x[j][k] = a.x[j][k] ^ gf[mlt(fg[d], fg[a.x[i][k]])];
                }
            }
        }
    }
    for (int i = 0; i < K; i++)
    {
        if (a.x[i][i] != 1)
            // exit(1);
            for (int j = 0; j < K + 1; j++)
                printf("%d,", a.x[i][j]);
        printf("\n");
    }
    printf("\n");

    for (int i = 0; i < K; i++)
    {
        v.x[i] = a.x[i][K];
        // v.x[128]=1;
        printf(" x%d = %d\n", i, v.x[i]);
    }

    return v;
}

vec vv(int kk)
{
    int i, j;
    vec r; //= mkpol();
    unsigned short tr[N];
    unsigned short ta[N] = {0};

    printf("van der\n");

    for (int i = 0; i < kk; i++)
    {
        for (int j = 0; j < M; j++)
        {
            vb[i][j] = gf[mltn(i, fg[j])];
            // printf("%d,", mltn(i, j));
        }
        printf("\n");
    }

    int l = -1;
    vec pp = {0}, tt = {0};

aa:
    // exit(1);
    r = mkpol();
    // for( int i= 0; i < K; i++)
    // pp.x[i] = rand() % N;
    // mykey(r.x, pp);
    // r.x[K] = 1;
    // printsage(r);
    // r=(tt);
    for (int i = 0; i < M; i++)
    {
        ta[i] = eval(r, i);
        if (ta[i] == 0)
        {
            printf("eval 0 @ %d\n", i);
            // fail = 1;
            goto aa;
        }
    }

    for (int i = 0; i < M; i++)
    {
        tr[i] = oinv(ta[i]);
        // printf("%d,", tr[i]);
    }
    printf("%d %d\n", r.x[0], r.x[K]);
    // exit(1);
    //memset(pp.x, 0, sizeof(pp.x));
    // for(int i=1;i<K+1;i++)
    //     pp.x[i]=r.x[i];
    ogt(r.x, K);
    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < K; j++)
        {
            mat[j][i] = gf[mlt(fg[vb[j][i]], tr[i])];
        }
    }
    /*
        for (int i = 0; i < K; i++)
        {
            for (int j = 0; j < M; j++)
            {
                for (int k = 0; k < K; k++)
                {
                    mat[i][j] = plus(mat[i][j], gf[mlt(fg[gt[k][i]],fg[ma[k][j]])]);
                }
                printf("c%d,", mat[i][j]);
            }
            printf("\n");
        }
    */
    // exit(1);
    //  bdet();
    return r;
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
            prod[i + j] = plus(prod[i + j], gf[mlt(fg[in0[i]], fg[in1[j]])]);
    }
    //

    for (i = (K - 1) * 2; i >= K; i--)
    {
        if (K == 512)
        {
            // GF(2^512) from sage
            prod[i - K + 8] = plus(prod[i - K + 8], prod[i]);
            prod[i - K + 5] = plus(prod[i - K + 5], prod[i]);
            prod[i - K + 2] = plus(prod[i - K + 2], prod[i]);
            prod[i - K + 0] = plus(prod[i - K + 0], prod[i]);
        }
        if (K == 256)
        {
            // GF(2^256) from sage
            prod[i - K + 10] = plus(prod[i - K + 10], prod[i]);
            prod[i - K + 5] = plus(prod[i - K + 5], prod[i]);
            prod[i - K + 2] = plus(prod[i - K + 2], prod[i]);
            prod[i - K + 0] = plus(prod[i - K + 0], prod[i]);
        }
        if (K == 238)
        {
            // GF(2^256) from sage
            prod[i - K + 73] = plus(prod[i - K + 10], prod[i]);
            prod[i - K + 0] = plus(prod[i - K + 0], prod[i]);
        }
        if (K == 128)
        {
            // 128
            prod[i - K + 7] = plus(prod[i - K + 7], prod[i]);
            prod[i - K + 2] = plus(prod[i - K + 2], prod[i]);
            prod[i - K + 1] = plus(prod[i - K + 1], prod[i]);
            prod[i - K + 0] = plus(prod[i - K + 0], prod[i]);
        }
        if (K == 96)
        {
            // 96
            prod[i - K + 63] = plus(prod[i - K + 63], prod[i]);
            prod[i - K + 61] = plus(prod[i - K + 61], prod[i]);
            prod[i - K + 60] = plus(prod[i - K + 60], prod[i]);
            prod[i - K + 59] = plus(prod[i - K + 59], prod[i]);
            prod[i - K + 54] = plus(prod[i - K + 54], prod[i]);
            prod[i - K + 50] = plus(prod[i - K + 50], prod[i]);
            prod[i - K + 46] = plus(prod[i - K + 46], prod[i]);
            prod[i - K + 44] = plus(prod[i - K + 44], prod[i]);
            prod[i - K + 40] = plus(prod[i - K + 40], prod[i]);
            prod[i - K + 37] = plus(prod[i - K + 37], prod[i]);
            prod[i - K + 33] = plus(prod[i - K + 33], prod[i]);
            prod[i - K + 32] = plus(prod[i - K + 32], prod[i]);
            prod[i - K + 28] = plus(prod[i - K + 28], prod[i]);
            prod[i - K + 24] = plus(prod[i - K + 24], prod[i]);
            prod[i - K + 23] = plus(prod[i - K + 23], prod[i]);
            prod[i - K + 21] = plus(prod[i - K + 21], prod[i]);
            prod[i - K + 20] = plus(prod[i - K + 20], prod[i]);
            prod[i - K + 19] = plus(prod[i - K + 19], prod[i]);
            prod[i - K + 18] = plus(prod[i - K + 18], prod[i]);
            prod[i - K + 16] = plus(prod[i - K + 16], prod[i]);
            prod[i - K + 15] = plus(prod[i - K + 15], prod[i]);
            prod[i - K + 14] = plus(prod[i - K + 14], prod[i]);
            prod[i - K + 13] = plus(prod[i - K + 13], prod[i]);
            prod[i - K + 11] = plus(prod[i - K + 11], prod[i]);
            prod[i - K + 9] = plus(prod[i - K + 9], prod[i]);
            prod[i - K + 8] = plus(prod[i - K + 8], prod[i]);
            prod[i - K + 7] = plus(prod[i - K + 7], prod[i]);
            prod[i - K + 6] = plus(prod[i - K + 6], prod[i]);
            prod[i - K + 4] = plus(prod[i - K + 4], prod[i]);
            prod[i - K + 3] = plus(prod[i - K + 3], prod[i]);
            prod[i - K + 2] = plus(prod[i - K + 2], prod[i]);
            prod[i - K + 0] = plus(prod[i - K + 0], prod[i]);
        }
        if (K == 64)
        {
            // 64
            prod[i - K + 33] = plus(prod[i - K + 33], prod[i]);
            prod[i - K + 30] = plus(prod[i - K + 30], prod[i]);
            prod[i - K + 26] = plus(prod[i - K + 26], prod[i]);
            prod[i - K + 25] = plus(prod[i - K + 25], prod[i]);
            prod[i - K + 24] = plus(prod[i - K + 24], prod[i]);
            prod[i - K + 23] = plus(prod[i - K + 23], prod[i]);
            prod[i - K + 22] = plus(prod[i - K + 22], prod[i]);
            prod[i - K + 21] = plus(prod[i - K + 21], prod[i]);
            prod[i - K + 20] = plus(prod[i - K + 20], prod[i]);
            prod[i - K + 18] = plus(prod[i - K + 18], prod[i]);
            prod[i - K + 13] = plus(prod[i - K + 13], prod[i]);
            prod[i - K + 12] = plus(prod[i - K + 12], prod[i]);
            prod[i - K + 11] = plus(prod[i - K + 11], prod[i]);
            prod[i - K + 10] = plus(prod[i - K + 10], prod[i]);
            prod[i - K + 7] = plus(prod[i - K + 7], prod[i]);
            prod[i - K + 5] = plus(prod[i - K + 5], prod[i]);
            prod[i - K + 4] = plus(prod[i - K + 4], prod[i]);
            prod[i - K + 2] = plus(prod[i - K + 2], prod[i]);
            prod[i - K + 1] = plus(prod[i - K + 1], prod[i]);
            prod[i - K + 0] = plus(prod[i - K + 0], prod[i]);
        }
        
        if (K == 32)
        {
            // 32
            prod[i - K + 15] = plus(prod[i - K + 15], prod[i]);
            prod[i - K + 9] = plus(prod[i - K + 9], prod[i]);
            prod[i - K + 7] = plus(prod[i - K + 7], prod[i]);
            prod[i - K + 4] = plus(prod[i - K + 4], prod[i]);
            prod[i - K + 3] = plus(prod[i - K + 3], prod[i]);
            prod[i - K + 0] = plus(prod[i - K + 0], prod[i]);
        }
        if (K == 16)
        {
            // 16
            prod[i - K + 5] = plus(prod[i - K + 5], prod[i]);
            prod[i - K + 3] = plus(prod[i - K + 3], prod[i]);
            prod[i - K + 2] = plus(prod[i - K + 2], prod[i]);
            prod[i - K + 0] = plus(prod[i - K + 0], prod[i]);
        }
        if (K == 8)
        {
            // 8
            prod[i - K + 4] = plus(prod[i - K + 4], prod[i]);
            prod[i - K + 3] = plus(prod[i - K + 3], prod[i]);
            prod[i - K + 2] = plus(prod[i - K + 2], prod[i]);
            prod[i - K + 0] = plus(prod[i - K + 0], prod[i]);
        }
        if (K == 6)
        {
            // 8
            prod[i - K + 4] = plus(prod[i - K + 4], prod[i]);
            prod[i - K + 3] = plus(prod[i - K + 3], prod[i]);
            prod[i - K + 1] = plus(prod[i - K + 1], prod[i]);
            prod[i - K + 0] = plus(prod[i - K + 0], prod[i]);
        }
        if (K == 4)
        {
            // 4
            prod[i - K + 1] = plus(prod[i - K + 1], prod[i]);
            prod[i - K + 0] = plus(prod[i - K + 0], prod[i]);
        }
    }

    for (i = 0; i < K; i++)
        out[i] = prod[i];
}

/* input: f, element in GF((2^m)^t) */
/* output: out, minimal polynomial of f */
/* return: 0 for success and -1 for failure */
int mykey(unsigned short *out, vec x)
{
    unsigned short m[K + 1][K] = {0};
    MTX a = {0};
    int i, j;

    // fill matrix

    m[0][0] = 1;

    for (i = 1; i < K; i++)
        m[0][i] = 0;

    for (i = 0; i < K; i++)
        m[1][i] = x.x[i];

    for (j = 2; j <= K; j++)
    {
        GF_mul(m[j], m[j - 1], x.x);
    }
    // exit(1);
    //
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < K + 1; j++)
        {
            a.x[i][j] = m[j][i];
            printf("%d,", m[j][i]);
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

void v2(int kk)
{
    int i, j, l = -1;
    vec f, r, tt, pp;
    unsigned short tr[N], te[M] = {0};
    ;
    unsigned short ta[N] = {0};

    printf("van der\n");

    // while (l < 0)
    {
        for (i = 0; i < K; i++)
            pp.x[i] = rand() % N;
        mykey(tt.x, pp);
        tt.x[K] = 1;
        // l = ben_or(tt);
        if (l == 0)
        {
            printf("\n");
            printsage(tt);
            printf(" ==irr\n");
            // exit(1);
        }
    }
    // f = (tt);
    // exit(1);

    l = 0;
aa:
    f = mkpol();
    for (int j = 0; j < M; j++)
    {
        te[j] = eval(f, fg[j]);
        if (te[j] == 0)
            goto aa;
    }

    vec h = {0}, g = {0};
    unsigned short tb[M] = {0}, tc[M] = {0}, td[M] = {0};

bb:
    g = mkpol();
    for (int j = 0; j < M; j++)
    {
        td[i] = eval(g, fg[i]);

        if (td[i] == 0)
            goto bb;
    }

cc:
    h = mkpol();

    for (int j = 0; j < M; j++)
    {
        ta[i] = eval(h, fg[i]);
        if (ta[i] == 0)
            goto cc;
    }

    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < kk; i++)
        {
            vb[i][j] = gf[mlt(mltn(i, fg[j]), mltn(j, fg[te[j]]))];
        }
    }
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < K; j++)
        { // この辺の添字は適当にやっても大丈夫
            tc[i] = plus(gf[mltn(M - j - 1, fg[plus(random() % O, minus(td[i]))])], mltn(i, fg[ta[j]]));
        }
        tr[i] = gf[oinv(tc[i])];
    }

    for (int i = 0; i < M; i++)
    {
        int zero = 0;
        for (int j = 0; j < K; j++)
        {
            mat[j][i] = gf[mlt(fg[vb[j][i]], fg[tr[i]])];
            if (mat[j][i] == 0)
                zero++;
        }
        if (zero == K)
            goto aa;
    }
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < M; j++)
            printf("c%d,", mat[i][j]);
        printf("\n");
    }
    printf("\n");
    unsigned short zz[M] = {0};
    mkerr(zz, T);
    vec v = synd(zz, K);
    ymo y = bm_itr(v.x);
    chen(y.f);
    for(i=0;i<N;i++)
    if(zz[i]>0)
    printf("i=%d %d\n",i,zz[i]);

    exit(1);
    //  bdet();
}

void v3(int kk)
{
    int i, j;
    vec f;
    unsigned short tr[N], te[M] = {0};
    ;
    unsigned short ta[N] = {0};

    printf("van der\n");

aa:
    f = mkpol();
    for (int j = 0; j < M; j++)
    {
        te[j] = eval(f, fg[j]);
        for (int i = 0; i < kk; i++)
        {
            vb[i][j] = gf[mlt(mltn(i, fg[j]), mltn(j, fg[te[j]]))];
        }
    }

    int l = -1;
    vec pp = {0}, tt = {0}, h = {0}, g = {0};
    unsigned short tb[M] = {0}, tc[M] = {0}, td[M] = {0};

    h = mkpol();
    g = mkpol();

    for (int i = 0; i < M; i++)
    {
        td[i] = eval(g, fg[i]);
        ta[i] = eval(h, fg[i]);
        if (td[i] == 0 || ta[i] == 0)
            goto aa;
        for (int j = 0; j < K; j++)
        { // この辺の添字は適当にやっても大丈夫
            tc[i] = plus(gf[mltn(M - j - 1, fg[plus(random() % O, minus(td[i]))])], mltn(i, fg[ta[j]]));
        }
        tr[i] = gf[oinv(tc[i])];
    }

    for (int i = 0; i < M; i++)
    {
        int zero = 0;
        for (int j = 0; j < K; j++)
        {
            mat[j][i] = gf[mlt(fg[vb[j][i]], fg[tr[i]])];
            if (mat[j][i] == 0)
                zero++;
        }
        if (zero == K)
            goto aa;
    }
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < M; j++)
            printf("c%d,", mat[i][j]);
        printf("\n");
    }
    printf("\n");
    unsigned short zz[M] = {0};
    mkerr(zz, T);
    vec v = synd(zz, K);
    ymo y = bm_itr(v.x);
    chen(y.f);
    // exit(1);
    //  bdet();
}

int badd(int i, int j)
{
    int z = (i % 2 | (i >> 1)) & (j % 2 | (j >> 1));
    int c = z ^ (i % 2 | j % 2);
    int d = z ^ (i >> 1 | j >> 1);
    c ^= (d << 1);

    return c;
}

int bsub(int i, int j)
{
    int z = (i >> 1 | i % 2) & (j & 2 | j >> 1);
    int c = z ^ (i % 2 | j >> 1);
    int d = z ^ (i >> 1 | j % 2);
    c ^= (d << 1);

    return c;
}

vec to_vec(unsigned short a)
{
    vec f = (dick[fg[a]]);

    return f;
}

unsigned short to_val(vec a)
{
    unsigned short c = 0;

    c = eval(a, Pr);

    return c;
}

unsigned int gf_div(unsigned int a, unsigned int b)
{
    int i = fg[a] - fg[b];

    if (a == 0 || b == 0)
        return 0;
    if (i == 0)
        return 1;
    if (i > 0)
        return (N - i);
    if (i < 0)
        return 1 - i;
}

uni u2d(unsigned short s)
{
    vec v = dick[fg[s]];
    uni d = {s};

    d.fugo.flag = __builtin_parity(d.u) ^ 1;

    return d;
}

unsigned short d2u(uni d)
{
    vec v = {d.fugo.b0, d.fugo.b1, d.fugo.b2, d.fugo.flag};
    bool c = __builtin_parity(d.u) ^ 1;

    if ((d.fugo.flag == 1 && c == 0) || (d.fugo.flag == 0 && c == 1))
    {
        printf("error\n");
        exit(1);
    }

    return xtrace(v, Pr);
}

vec m2v(short a)
{
    short u = a ^ (1 << 15);
    vec v = {0};

    v.x[0] = u & 0x1f;
    u = (u >> 5);
    v.x[1] = u & 0x1f;
    u = (u >> 5);
    v.x[2] = u & 0x1f;

    return v;
}

void forney(ymo t)
{
    vec v = {0}, tmp[T] = {0}, w = {0}, kai = {0}, y = {0}, x;
    unsigned short u = 0, z, a, n;
    x = t.f;
    if (x.x[0] > 0)
    {
        u = gf[oinv(x.x[0])];
        y = kof(u, x);
        printpol(y);
        printf("\n");
    }
    else
    {
        y = x;
    }
    v = chen(y);
    for (int i = 0; i < T; i++)
        w.x[i] = gf[O - v.x[i] - 1];
    printpol(w);
    printf("\n");

    kai.x[0] = 1;
    kai.x[1] = minus(w.x[0]);
    u = eval(kai, v.x[1]);
    z = gf[mlt(fg[u], fg[w.x[1]])];
    n = eval(t.g, v.x[1]);
    a = gf[mlt(fg[n], oinv(z))];
    printpol(x);
    printf("%d\n", (a));
    exit(1);

    return;
}

vec deli(vec a, vec b)
{
    vec c = {0};

    for (int i = 0; i < deg(b); i++)
        c.x[i] = a.x[i];

    return c;
}

void test1()
{
    unsigned short zz[M] = {0};

    printf("aa %d \n", plus(6, 5));
    // exit(1);

    memset(zz, 0, sizeof(zz));
    mkerr(zz, T);
    // zz[1] = 1;
    // zz[3] = 1;
    //  zz[aa]=i;
    vec cc = {0};
    cc.x[T * 2] = 1;
    // zz[1]=4;
    vec x = synd(zz, K);
    ymo oo = bm_itr(x.x);
    vec v = chen(oo.f);
    cc = deli(vmul(oo.f, x), cc);
    printsage(x);
    printf(" ====syn\n");
    printsage(oo.f);
    printf(" ====ffx\n");
    exit(1);
    for (int i = 0; i < O; i++)
    {
        int qt = eval(cc, gf[O - 1 - v.x[0]]);
        printf("qt=%d\n", qt);
        oo.f = bib(oo.f);
        printpol(oo.f);
        printf("\n");
        int qt2 = eval(oo.f, gf[O - 1 - v.x[0]]);
        printf("qt2=%d %d\n", qt2, gf[O - 1 - v.x[0]]);
        printf("%d\n", gf[mlt(fg[qt], oinv(qt2))]);
        exit(1);
        // gf[mlt(fg[(xtrace(oo.h, fg[i]))], oinv(xtrace(oo.f,fg[i])))];
        if (qt == zz[1] || qt == zz[3])
            printf("ee %d %d\n", i, qt);
    }
    exit(1);
    int uu = gf[mlt(fg[eval(oo.g, v.x[0])], fg[eval(oo.h, v.x[0])])];
    printf("dd %d\n", uu);
    printpol(oo.g);
    printf("\n");
    if (uu == 1)
    {
        uu = gf[oinv(gf[mlt(fg[eval(oo.g, fg[v.x[1]])], fg[eval(oo.f, fg[v.x[1]])])])];
        printf("%d %d \n", uu, gf[oinv(gf[mlt(fg[eval(oo.f, v.x[0])], fg[eval(oo.h, v.x[0])])])]);
        exit(1);
        oo.h = kof(equ(vLT(oo.h).a, 2), oo.h);
        printpol(oo.h);
        printf("\n");
        printpol(oo.f);
        printf("\n");
        printpol(oo.g);
        printf("\n");
    }
    // printpol(oo.g);
    // printf("\n");
    exit(1);
    for (int j = 1; j < O; j++)
    {
        // for (int i = 1; i < O; i++)
        {
            printpol(oo.f);
            printf(" =====oo.f\n");
            // exit(1);
            vec c = {0};
            c.x[T * 2] = 1;
            vec xv = (vmul(x, oo.f));
            printpol(xv);
            printf(" ==Uh!\n");
            for (int i1 = 0; i1 < T; i1++)
                printf("%d,", v.x[i1]);
            printf("\n");
            // exit(1);
            oo.h = kof((Pr), oo.h);
            printpol(oo.h);
            printf(" =====Ub2oo.h\n");
            printpol(oo.g);
            printf(" =====Ub2oo.g\n");
            printpol(bib(oo.f));
            printf(" =====Ub2bib.f\n");
            printpol(bibun_old(oo.f));
            printf(" =====Ub2ibun.f\n");
            vec ff = {0};
            ff.x[0] = 7;
            ff.x[1] = 3;
            printpol(xv);
            printf(" ==U\n");
            printpol(bibun_old(oo.f));
            printf(" ==U2\n");
            printpol(bib(oo.f));
            printf(" ==U3\n");
            printf("U%d\n", mlt(fg[eval(ff, 3)], oinv(eval(bib(oo.f), 3))));
            printf("U%d\n", mlt(fg[eval(ff, 5)], oinv(eval(bib(oo.f), 5))));
            printf("U%d\n", mlt(fg[eval(oo.h, 6)], oinv(eval(bib(oo.f), 6))));
            printf("U%d\n", mlt(fg[eval(oo.h, 7)], oinv(eval(bib(oo.f), 7))));

            exit(1);
            for (int k = 0; k < T; k++)
            {
                if (logx(zz[v.x[k]]) == ((logx(eval(oo.h, v.x[k]))) + j) % (O - 1))
                    printf("jaa===%d  %d  %d   %d  %d  %d   %d j=%d\n", zz[v.x[k]], logx(zz[v.x[k]]), gf[mlt(fg[eval(oo.h, v.x[k])], oinv(eval(bib(oo.f), v.x[k])))], logx(gf[oinv(xtrace(bibun_old(oo.f), v.x[k]))]) % (O - 1), logx(gf[oinv(xtrace(bib(oo.f), v.x[k]))]), (logx(gf[oinv(xtrace((oo.f), v.x[k]))])), logx(gf[oinv(xtrace(bib(oo.f), v.x[k]))]) % (O - 1), j);
            }
            // printf("%d  %d\n",gf
            // printf("%d  %d\n",gf[mlt(fg[eval(oo.h,gf[2])],oinv(eval((oo.f),gf[2])))],gf[oinv(eval((oo.f),gf[2]))]);
        }
    }
}

int ind(unsigned short a)
{
    return logx(a) + 1;
}

void stamena()
{

    do
    {
        memset(inv_S.x, 0, sizeof(inv_S.x));
        memset(S.x, 0, sizeof(S.x));
        for (int i = 0; i < K * E; i++)
        {
            for (int j = 0; j < K * E; j++)
                S.x[i][j] = random() % Pr;
        }
        inv_S = genS(S);
    } while (inv_S.b != 1);
    CTX Pub = {0};
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < M; j++)
        {
            Pub.x[i][j] = mat[i][j];
            printf("%d,", mat[i][j]);
        }
        printf("\n");
    }

    CTX R = kenzan(S, Pub, 2);
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < M; j++)
            printf("%d ", R.x[i][j]);
        printf("\n");
    }
    exit(1);

    vec v = {0};
    for (int i = 0; i < K * E; i++)
    {
        v.x[i] = (i + 1) % Pr;
        printf("%d,", v.x[i]);
    }
    printf("\n");

    vec cc = {0};
    // exit(1);
    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < K * E; j++)
            cc.x[i] += (S.x[i][j] * v.x[j]) % Pr;
        cc.x[i] %= Pr;
        printf("%d,", cc.x[i]);
    }
    printf("\n");
    // exit(1);
    memset(v.x, 0, sizeof(v.x));
    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < K * E; j++)
            v.x[i] += (inv_S.x[i][j] * cc.x[j]) % Pr;
        v.x[i] %= Pr;
        printf("%d ", v.x[i]);
    }
    printf("\n");
    // exit(1);
    memset(R.x, 0, sizeof(R.x));
    R = kenzan(S, inv_S, 1);
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            printf("%d,", R.x[i][j] % Pr);
        }
        printf("\n");
    }
    printf("\n");
    // exit(1);

    int u = 0, u1 = 0, u2 = 0, u3 = 0;
    // #pragma omp parallel for
    for (int i = 0; i < K; i++)
    {
        int k = 0;
        for (int j = 0; j < M; j++)
        {
            big.d[i][j].fugo.b0 = dick[fg[mat[i][j]]].x[2];
            big.d[i][j].fugo.b1 = dick[fg[mat[i][j]]].x[1];
            big.d[i][j].fugo.b2 = dick[fg[mat[i][j]]].x[0];
            vec xvi = {0};
            xvi = m2v(big.x[i][j]);
            printf("s1=%d %d %d %d %d %d %d \n", big.d[i][j].fugo.b0, big.d[i][j].fugo.b1, big.d[i][j].fugo.b2, big.x[i][j], xvi.x[0], xvi.x[1], xvi.x[2]);
            v.x[2] = big.d[i][j].fugo.b0;
            v.x[1] = big.d[i][j].fugo.b1;
            v.x[0] = big.d[i][j].fugo.b2;
            // printpol(v);
            // printf(" ====v\n");
        }
        // printf("\n");
    }

    // exit(1);
    printf("\n");
    printf("zもとの行列\n");
    static CTX C = {0};
    for (int j = 0; j < K; j++)
    {
        for (int i = 0; i < M; i++)
        {
            C.x[j * E][i] = (Pr + big.d[j][i].fugo.b0) % Pr;
            printf("%d ", big.d[j][i].fugo.b0);
        }
        printf("\n");
        for (int i = 0; i < M; i++)
        {
            C.x[j * E + 1][i] = (Pr + big.d[j][i].fugo.b1) % Pr;
            printf("%d ", big.d[j][i].fugo.b1);
        }
        printf("\n");
        for (int i = 0; i < M; i++)
        {
            C.x[j * E + 2][i] = (Pr + big.d[j][i].fugo.b2) % Pr;
            printf("%d ", big.d[j][i].fugo.b2);
        }
        printf("\n");
    }
    // exit(1);
    // inv_S.d[0][0].b = 0;
    //  OO=pk_gen();
    printf("over\n");
    printf("zかけた結果\n");
    R = kenzan(S, C, 2);
    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < M; j++)
        {
            // mat[i][j] = R.x[i][j];
            printf("%d,", (R.x[i][j]) % Pr);
        }
        printf("\n");
    }
    printf("\n");
    // exit(1);

    printf("over2\n");
    printf("zもとに戻す\n");
    CTX OO = kenzan(inv_S, R, 2);
    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < M; j++)
        {
            // mat[i][j] = R.x[i][j];
            printf("%d ", (OO.x[i][j]) % Pr);
        }
        printf("\n");
    }
    printf("\n");
    // exit(1);
}

// 言わずもがな
int main(void)
{
    int i, a, b, c;
    unsigned short zz[M] = {0};
    vec r = {0}, g = {0};
    vec v = {0}, x = {0}, f = {0};
    static CTX R = {0}, OO = {0}; //, S = {0}, inv_S = {0};
    unsigned short s[K + 1] = {0, 3, 0, 6, 0};
    static data d = {17, 23, 31}; //, dd[K * 3][N] = {0};
    uni e = {d};

    e.fugo.flag = __builtin_parity(e.u);

    if (K > N)
        printf("configuration error! K is bigger than N\n");
    d.b2 = (d.b0 + d.b1);
    printf("%d %d %d %d %u %d\n", d.b0, d.b1, d.b2, (17 + 23) % 23, e.u, e.fugo.flag);

    printf("%d %d %d %d\n", sizeof(uni), sizeof(MTX), sizeof(CTX), sizeof(data));
    // exit(1);
    // （謎）
    memset(mat, 0, sizeof(mat));
     srand(clock());
    /*
    R=lu();
    OO=genS(R);
    OO=kenzan(R,OO,1);
    for(int i=0;i<F;i++){
    for(int j=0;j<F;j++)
    printf("%d ",OO.x[i][j]);
    printf("\n");
    }
    printf("\n");
    exit(1);
    */
    mkmf();
    makefg();
    de();
    printf("%d %d %d\n", minus(1), minus(2), minus(3));
    //exit(1);

    unsigned short z1[M] = {0};
    z1[1] = 1;
    z1[2] = 1;
    pk_gen();
    exit(1);

    // big = fls;
    int count = 0;
    uni on[4] = {0};
    unsigned short es = 0;
    unsigned short ee[K] = {0};

    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < M; j++)
        {
            printf("%d ", OO.x[i][j]);
        }
        on[count].fugo.b0 = (OO.x[i * E + 0][3] + OO.x[i][6]) % Pr;
        on[count].fugo.b1 = (OO.x[i * E + 1][3] + OO.x[i + 1][6]) % Pr;
        on[count++].fugo.b2 = 0;
        printf("\n");
    }
    printf("aaaaaaaaaaaaaaaaaaaaaaaaaa\n");
    for (i = 0; i < K; i++)
        v.x[i] = c2d(on[i]);
    ymo u = bm_itr(v.x);
    chen(u.f);
    printf("\n");
    exit(1);

    transi(OO);
    // exit(1);
    unsigned short ss[K * E] = {0};
    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (z1[j] > 0)
                ss[i] = (ss[i] + big.x[i][j]) % Pr;
        }
        printf("%d,", ss[i]);
    }
    printf("\n");
    // exit(1);

    // big =
    pk_gen();
    uni op[K] = {0};
    memset(zz, 0, sizeof(zz));
    memset(x.x, 0, sizeof(x.x));
    memset(v.x, 0, sizeof(v.x));
    mkerr(zz, T);
    // zz[1] = 1;
    // zz[2] = 1;
    //  big=OO;
    v = sind(zz, K);
    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < K * E; j++)
        {
            x.x[i] = (x.x[i] + inv_S.x[i][j] * v.x[j]) % Pr;
        }
    }
    unsigned short oo[K] = {0};
    for (int i = 0; i < K; i++)
    {
        op[i].fugo.b0 = x.x[i * E];
        op[i].fugo.b1 = x.x[i * E + 1];
        op[i].fugo.b2 = 0;
        oo[i] = c2d(op[i]);
    }
    ymo yy = bm_itr(oo);
    chen(yy.f);
    exit(1);

    for (int i = 0; i < K * E; i++)
        printf("!%d,", op[i].fugo.b1);
    printf("\n");
    for (int i = 0; i < K * E; i++)
        printf("!%d,", op[i].fugo.b2);
    printf("\n");
    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < M; j++)
            printf("%d ", inv_S.x[i][j]);
        printf("\n");
    }
    printf("\n");
    exit(1);

    for (int i = 0; i < K * E; i++)
    {
        for (int j = 0; j < M; j++)
        {
            printf("%d ", (OO.x[i][j]));
        }
        printf("\n");
        // for(int j=0;j<M;j++)
        // printf("%d ",OO.d[i][j].fugo.b2);
        // printf("\n");
    }
    // exit(1);

    //---------------------------------------------------20230325
    // exit(1);
    for (int i = 0; i < O; i++)
        if (gf[mlt(fg[5], fg[i])] == 2)
            printf("%d\n", i);
    // printf("1331=%d\n",c2o(o2c(dick[fg[1331]])));
    // exit(1);
    //   次元KのGvecpa符号の検査行列を生成する
    // van(K);
    // v2(K);
    f = vv(K);
    // exit(1);
    int j;
    static vec cc = {0};
    printf("%d\n", minus(3));
    // exit(1);
    // test1();
    // exit(1);

    vec v1 = {0}, v2 = {0}, x1 = {0}, x2 = {0};
    int it = 0;
    while (1)
    {
        //    big=R;
        //   エラーベクトルの初期化
        memset(zz, 0, sizeof(zz));
        //   重みTのエラーベクトルを生成する
        // mkerr(zz, T);
        // for (int i = 0; i < K / 2; i++)
        //    zz[i] = i + 1;
        // zz[0] = 1;
        mkerr(zz, T);
        // zz[1] = 4;
        // zz[2] = 1;
        //  r = sind(zz, K);
        printpol(r);
        printf(" ===r\n");
        // x = synd(zz, K);
        printpol(x);
        printf(" ===x\n");
        // exit(1);

        int u;
        uni let = {0};
        unsigned short h[K * E][M] = {0};
        vec me = {0}, z;
        memset(v.x, 0, sizeof(v.x));
        memset(x.x, 0, sizeof(x.x));
        // r = synd(zz, K);

        // r=d2v(r.x);
        for (int i = 0; i < K; i++)
        {
            for (int j = 0; j < M; j++)
            {
                h[i * E][j] = inv_S.x[i * E][j];
                h[i * E + 1][j] = inv_S.x[i * E + 1][j];
                h[i * E + 2][j] = inv_S.x[i * E + 2][j];
            }
            let.u = r.x[i];
            me.x[i * E] = let.fugo.b0;
            me.x[i * E + 1] = let.fugo.b1;
            me.x[i * E + 2] = let.fugo.b2;
        }

        for (int i = 0; i < K * E; i++)
        {
            u = 0;
            for (int k = 0; k < K * E; k++)
                u = (u + me.x[k] * inv_S.x[k][i]) % Pr;
            z.x[i] = u % Pr;
        }
        // exit(1);
        uni dd[K] = {0};
        vec hi = {0};
        for (int i = 0; i < K; i++)
        {
            dd[i].fugo.b0 = z.x[i * E];
            dd[i].fugo.b1 = z.x[i * E + 1];
            dd[i].fugo.b2 = z.x[i * E + 2];
            hi.x[i] = dd[i].u;
        }
        x = synd(zz, K);
        printpol(x);
        printf(" ===x\n");
        printpol(hi);
        printf(" ===hi\n");
        // exit(1);

        // exit(1);
        /*
        printpol(bibun_old(oo.f));
        printf("   bibun\n");
        printpol(bib(oo.f));
        printf("  bib\n");

        exit(1);
        printf("   %d\n", eval(bib(oo.f), gf[v.x[0]]));
        printf("   %d\n", gf[mlt(fg[eval(bib(oo.f), v.x[0])],fg[eval(oo.g, v.x[0])])]);
        //exit(1);
        //forney(v);
        x=bibun_old(oo.f);
        printf("ixi=%d\n",gf[oinv(gf[mlt(fg[eval(x,v.x[0])],fg[eval(oo.g,v.x[0])])])]);
        exit(1);
        */
        /*
            uni vv[K]={0};
            for (int i = 0; i < K*3; i += 3)
            {
                vv[i].fugo.b0 = v.x[i + 0];
                vv[i].fugo.b1 = v.x[i + 1];
                vv[i].fugo.b2 = v.x[i + 2];
                r.x[i] = vv[i].u;
                printf("r%d", r.x[i]);
            }
            printf("\n");
        */
        // r=sin(zz,OO);
        //  exit(1);
        // r=dec(r.x);
        // printpol(v);
        // exit(1);

        memset(zz, 0, sizeof(zz));
        mkerr(zz, T);
        // for(i=1;i<T+1;i++)
        // zz[i]=i;
        /*
        zz[1] = 1;
        zz[2] = 1;
        zz[3] = 1;
        zz[4]=1;
        zz[5]=1;
        */
        vec c = {0};
        c.x[T * 2] = 1;
        r = synd(zz, K);
        for (int i = 0; i < M; i++)
        {
            if (zz[i] > 0)
                printf("aieu=%d %d\n", i, zz[i]);
        }
        printpol(r);
        printf("\n");
        // exit(1);
        ymo t = bm_itr(r.x);
        vec ll = r;
        // unsigned short u;
        printpol((t.h));
        printf("  trf\n");
        r = chen(t.f);
        for (int i = 0; i < M; i++)
            if (zz[i] > 0)
                printf("err=%d\n", i);
        exit(1);
        t.f = bibun_old(r);
        printpol(t.f);
        printf("\n");
        // exit(1);
        for (int i = 0; i < T; i++)
        {
            // u = gf[oinv(mlt(fg[eval(t.g,r.x[i])],fg[eval(t.f,r.x[i])])))];
            u = gf[mlt(fg[eval(t.h, r.x[i])], oinv(eval(t.f, r.x[i])))];
            printf("u=%d %d %d\n", u, gf[(logx(u) + 3) % (O - 1) + 1], r.x[i]);
        }
        exit(1);
        /*
        r=vmul(r,t.f);
        r=deli(r,c);
        printsage(r);
        printf("  --r\n");
        //exit(1);
        printpol(t.g);
        printf(" t.g\n");
        printf("%d ans1\n",gf[mlt(fg[eval(r,8)],oinv(7))]);
        printpol(t.h);
        printf(" t.h\n");
        printf("%d ans2\n",gf[mlt(fg[xtrace(r,5)],oinv(2))]);
        //forney(t);
        exit(1);
        */

        x = chen((t.f));
        // exit(1);
        vec ans = {0};
        MTX ww = {0};
        for (int i = 0; i < T + 1; i++)
        {
            printf("i=%d %d \n", x.x[i], x.x[i] + 1);
            ww.x[0][i] = gf[x.x[i] + 1];
            ww.x[1][i] = gf[mlt(x.x[i] + 1, x.x[i] + 1)];
            ww.x[i][2] = ll.x[i];
            printf("b= %d %d\n", ww.x[0][i], ww.x[1][i]);
        }
        printpol(ll);
        printf(" ===synd\n");
        sol(ww);
        exit(1);
        // yes(t);
        //  exit(1);

        vec rr = {0};
        for (int i = 0; i < T; i++)
        {
            rr.x[x.x[i]] = 1; // x.x[i];
            printf("x=%d\n", x.x[i]);
        }
        for (int i = 0; i < M; i++)
        {
            if ((zz[i] > 0) && zz[i] != rr.x[i])
            {
                printf("cheki=%d %d %d\n", i, zz[i], rr.x[i]);
                printf("baka\n");
                exit(1);
            }
        }
        exit(1);

        if (it == 1000)
        {
            printf("いいっ！\n");
            break;
        }
    }

    return 0;
}
