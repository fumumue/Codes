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

// #include "8192.h"
// #include "4096.h"
#include "global-p.h"
#include "struct.h"
#include "debug.c"
#include "chash-p.c"
#include "lu.c"
#include "sha3.c"
#include "inv_mat.c"
// #include "golay.c"

#define TH omp_get_max_threads()
#define O 27 // 1331,2197,4913,6859,3125,2187,19683
#define EXP 3
#define Pr 3
//#define F 16
extern unsigned long xor128(void);
extern int mlt(int x, int y);
// extern int mltn(int n, int a);
extern void makeS();

void vv(int kk);
unsigned short plus(unsigned short a, unsigned short b);
unsigned short bexp(unsigned short a);
unsigned short minus(unsigned short a);
unsigned short eval(vec f, unsigned short x);
void printpol(vec a);

// #pragma omp threadprivate(mat)
// シンドロームのコピー
unsigned short sy[K] = {0};
unsigned short pp[10][4] = {{0, 0, 9, 2}, {0, 0, 11, 2}, {0, 0, 16, 3}, {0, 0, 15, 2}, {0, 0, 1, 2}, {0, 1, 0, 2}, {0, 0, 1, 1}, {0, 0, 1, 2}, {1, 1, 2, 2}, {0, 0, 1, 2}};
unsigned short pt[3] = {0, 1, 2};
// Gvecpa多項式
static unsigned short g[K + 1] = {1, 2, 3, 4, 6};
//{1,0,0,0,0};
//{1,0,1,0,5};
//{1};
//{1, 0, 0, 0, 1, 0, 1};
//

MTX BB = {0};
MTX H = {0};

unsigned int AA = 0, B = 0; //, C = 0, A2 = 0;

// 有限体の元の逆数
unsigned short
oinv(unsigned short a)
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
    for (i = 2; i < N; i++)
    {
        if (gf[mlt(i, fg[a])] == 1)
            return i;
    }
    printf("no return\n");
    exit(1);
}

// 有限体の元の逆数
unsigned short
oinv3(unsigned short a)
{
    int i;

    if (a == 0)
        return 0;

    for (i = 0; i < N; i++)
    {
        if (gf[mlt(fg[a], i)] == 1)
            return (unsigned short)i;
    }

    printf("no return3 \n");
    //  exit (1);
}

// aに何をかけたらbになるか
unsigned short
equ(unsigned short a, unsigned short b)
{
    int i;

    for (i = 0; i < N; i++)
    {
        if (gf[mlt(fg[a], fg[i])] == b)
            break;
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

// vec型を正規化する
vec conv(vec f)
{
    vec v = {0};
    vec g = {0};

    v = (f);
    g = (v);

    return g;
}

int deg2(vec v){
    int i,n=0;

    for(i=0;i<DEG;i++)
    if(v.x[i]>0)
    n=i;

return i;
}

// 多項式の次数(default)
int deg(vec a)
{
    int i, n = 0, flg = 0;

    // #pragma omp parallel for
    for (i = 0; i < DEG; i++)
    {
        if (a.x[i] > 0 && a.x[i] <= O)
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
    int i, n;
    vec f = {0};
    f=g;
    // f = conv(f);
    n = deg(g);
    //printf("n=%d\n", n);
    // printf("terms=%d\n", terms(f));
    // printf("deg=%d\n", deg(f));

    for (i = n; i > -1; i--)
    {
        if (f.x[i] > 0)
            printf("%dx^%d+", f.x[i], i);
    }
//exit(1);
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
    int i;

    // memset (a, 0, sizeof (c));
    // memcpy (c, f, n);
    for (i = 0; i < n; i++)
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
vec vadd(vec f, vec g)
{
    vec a = {0}, b = {0}, c = {0};
    int i, j, k, l = 0;
    vec h = {0}, f2 = {0}, g2 = {0};

    // for(i=0;i<257;i++)
    //  printf("%d %d %d %d %d\n",i,f.x[i],f.x[i].n,g.x[i],g.x[i].n);

    //  exit(1);
    // printpol(f);
    // printf("\n");

    a = (f);
    // exit(1);
    b = (g);

    j = deg((f));
    l = deg((g));
    // printpol((f));
    // printf(" =f\n");
    // printpol((g));
    // printf(" =g\n");
    if (j >= l)
    {
        k = j + 1;
    }
    else
    {

        k = l + 1;
    }
    // for(i=0;i<k;i++)
    // printf("%d %d\n",i,b.x[i]);
    //   exit(1);

    for (i = 0; i < k; i++)
    {
        // if(f.x[i]>0 || g.x[i]>0)
        // h.x[i]=f.x[i]^g.x[i];
        c.x[i] = (a.x[i] + b.x[i]) % Pr;
        // if(a.x[i]>0 || b.x[i]>0)
        // printf("i=%d %d\n",a.x[i],b.x[i]);
    }
    //
    h = (c);
    // printpol((h));
    // printf(" ==cclemon\n");

    return h;
}

// 20200816:正規化したいところだがうまく行かない
// 多項式の足し算
vec vadd2(vec f, vec g)
{
    //assert(vec_verify(f));
    //assert(vec_verify(g));

    vec a = {0}, b = {0}, c = {0};
    int i, j, k, l = 0;
    vec h = {0}, f2 = {0}, g2 = {0};

    a = (f);
    b = (g);

    // k=deg((f));
    // l=deg((g));

    for (i = 0; i < DEG; i++)
    {
        c.x[i] = plus(a.x[i], b.x[i]);
        // h.x[i]=f.x[i]^g.x[i];
    }
    h = (c);
    // h=conv(h);
    //assert(vec_verify(h));
    return h;
}

// 多項式を項ずつ掛ける
vec oterml(vec f, oterm t)
{   
    //f=conv(f);
    // assert (vec_verify (f));
    int i, k, j;
    vec h = {0};
    vec test;
    unsigned short n;

printpol(f);
printf(" ==f\n");
printf("t.a=%d\n",t.a);
printf("t.n=%d\n",t.n);
//exit(1);
    // 
    k = deg(f)+1;
    j = 0;
    for (i = 0; i < k ; i++)
    {   if(f.x[i]>0){
        h.x[i+t.n] = (f.x[i] * t.a) % Pr;
        printf("t+n=%d\n",i+t.n);
    }
    }

     //h=conv(h);
    // assert (vec_verify (h));
    return h;
}

// 多項式を項ずつ掛ける
vec oterml2(vec f, oterm t)
{

    // assert (vec_verify (f));
    int i, k, j;
    vec h = {0};
    vec test;
    unsigned short n;

    // f=conv(f);
    k = deg(f);
    j = 0;
    for (i = 0; i < k + 1; i++)
    {
        h.x[i+t.n] = gf[mlt(fg[f.x[i]], fg[t.a])];
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
    int i, k, j;
    vec h = {0};
    vec test;
    unsigned int n;

    // f=conv(f);
    // k = deg ((f));
    j = 0;
    for (i = 0; i < DEG; i++)
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
    int i, k;

    if (deg(d) >= deg(e))
    {
        k = deg(d) + 1;
    }
    else
    {

        k = deg(e) + 1;
    }
    for (i = 0; i < k; i++)
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
    int i, k;
    oterm t = {0};

    // k = deg ( (f));
    for (i = 0; i < DEG; i++)
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
    // f=conv(f);
    // g=conv(g);

    // assert (op_verify (f));
    // assert (op_verify (g));
    int i, count = 0, k, l, m;
    oterm t = {0};
    vec h = {0}, e = {0}, r = {0};
    vec c = {0};

    l = deg((f));
    m = deg((g));

    //g = conv(g);
    for (i = 0; i < 2; i++)
        printf("%d\n", g.x[i]);
      //exit(1);

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

    for (i = 0; i < k+1; i++)
    {
        t.a = (g.x[i]);
        t.n=i;
        //t.n=i;
        if (t.a > 0)
        {
            printf("t[%d]=%d,%d\n", i, t.a, t.n);
            e = oterml(f, t);
            printpol((e));
            printf(" =e\n");
            //exit(1);
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


// 項の順序を降順に揃える
vec sort(vec f)
{
    oterm o = {0};
    int i, j, k;

    k = terms(f);
    for (i = 0; i < k + 1; ++i)
    {
        for (j = i + 1; j < k + 1; ++j)
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
    int i, j, k;
    vec a, b;

    f = conv(f);
    k = deg((f)) + 1;
    for (i = 0; i < k; i++)
        f.x[i] = gf[mlt(fg[f.x[i]], oinv(d))];

    return f;
}
/*
// 多項式を表示する(defauvLT)
void printpol(vec a)
{
    int i, n;
    vec v = {0};
    v = (a);
    n = deg((v));

    // printf ("baka\n");
    assert(("baka\n", n >= 0));

    for (i = n; i > -1; i--)
    {
        if (a.x[i] > 0)
        {
            printf("%u*", a.x[i]);
            printf("x^%d", a.x[i].n);
            printf("+");
        }
    }
    //  printf("\n");

    return;
}
*/
// 多項式の剰余を取る
vec vmod(vec f, vec g)
{
    int i = 0, j, n, k;
    vec h = {0}, e = {0};
    vec x, y;
    oterm a, b = {0}, c = {0};

    n = vLT(g).n;

    //  assert (("baka^\n", vLT (f).n != 0));

    //  assert (("baka(A)\n", vLT (g).n != 0));

    if (vLT(f).n < vLT(g).n)
    {
        //    exit(1);
        return f;
    }

    // printf ("in vmod\n");
    // exit(1);

    k = vLT(g).n;
    b = vLT(g);

    x = (g);
    assert(("double baka\n", b.a != 0 && b.n != 0));
    while (vLT(f).n > 0 && vLT(g).n > 0)
    {

        c = vLTdiv(f, b);
        h = oterml((x), c);
        f = vadd2(f, h);
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

    // f = conv(f);
    // g = conv(g);
    assert(vec_verify(f));
    assert(vec_verify(g));
    int i = 0, j, n, k;
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

    k = deg(g);
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
    i = 0;
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
    for (i = 0; i < tt_terms; i++)
    {
        ret.x[i] = tt.x[tt_terms - i - 1];
    }
    ret = conv(ret);
    assert(vec_verify(ret));
    return ret;
}

// 多項式のべき乗
vec opow(vec f, int n)
{
    int i;
    vec g = {0};

    g = f;

    for (i = 1; i < n; i++)
        g = vmul(g, f);

    return g;
}

// 多項式のべき乗余
vec opowmod(vec f, vec mod, int n)
{
    int i, j = 0;

    // 繰り返し２乗法
    for (i = 1; i < n + 1; i++)
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
    int i, k, l, m;
    vec ans = {0};

    a = (f);
    b = (g);
    l = deg(a);
    m = deg(b);
    if (l >= m)
    {
        k = l;
    }
    else
    {
        k = m;
    }
    for (i = 0; i < k + 1; i++)
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
    int n, i;
    vec g;

    r = (f);
    n = deg(r);
    for (i = 0; i < n + 1; i++)
        r.x[i] = (r.x[i] * a) % Pr;
    g = (r);

    return g;
}

unsigned short xtrace(vec f, unsigned short x)
{
    int i, d, j, z = 0;
    vec g = (f);
    unsigned short u = 0;

    d = deg(g) + 1;
    // if(g.x[0]>0)
    // z+=g.x[0];

    for (i = 0; i < d; i++)
    {
        u = 1;
        if (g.x[i] > 0)
        {
            for (j = 0; j < i; j++)
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
    unsigned short i, j, count = 0;

    // for (i = 0; i < O; i++)
    {
        for (j = 0; j < O; j++)
        {
            fg[gf[j]] = j;
        }
    }

    // exit(1);

    printf("unsigned short fg[%d]={", O);
    for (i = 0; i < O; i++)
        printf("%d,", fg[i]);
    printf("};\n");
    printf("count=%d\n", count);
    // exit(1);

    return;
}

vec o2c(vec a){
    vec v=(a);
    
    for(int i=0;i<EXP;i++)
    printf("%d,",v.x[i]);
    //printf("\n");

return v;
}

//vec dick[O] = {0};
//unsigned short val[N] = {0}; //, a = {0}, cc = {0};
vec dick[O] = {0};
unsigned short val[N] = {0}; //, a = {0}, cc = {0};
void mkmf()
{
    int i, j, k, count = 0;
    vec f = {0}, g = {0}, h = {0}, w = {0}, s = {0}, u = {0};
    vec b = {0}, a = {0}, d = {0}, t = {0}, v = {0};
    oterm o = {0};
    unsigned short ccp[4] = {0};

    if (O == 1331)
    {
        for (i = 0; i < 4; i++)
            ccp[i] = pp[0][i];
    }
    if (O == 2197)
    {
        for (i = 0; i < 4; i++)
            ccp[i] = pp[1][i];
    }
    if (O == 4913)
    {
        for (i = 0; i < 4; i++)
            ccp[i] = pp[2][i];
    }
    if (O == 6859)
    {
        for (i = 0; i < 4; i++)
            ccp[i] = pp[3][i];
    }
    if (O == 3125)
    {
        for (i = 0; i < 4; i++)
            ccp[i] = pp[4][i];
    }
    if (O == 2187)
    {
        for (i = 0; i < 4; i++)
            ccp[i] = pp[5][i];
    }
    if (O == 9)
    {
        for (i = 0; i < 4; i++)
            ccp[i] = pp[6][i];
    }
    if (O == 27)
    {
        for (i = 0; i < 4; i++)
            ccp[i] = pp[7][i];
    }
    if (O == 243)
    {
        for (i = 0; i < 4; i++)
            ccp[i] = pp[9][i];
    }
    if (O == 19683)
    {
        for (i = 0; i < 4; i++)
            ccp[i] = pp[8][i];
    }

    g = (setpol(ccp, 4));
    for(int i=0;i< 4;i++)
    printf("%d,",g.x[i]);
    printf("\n");
    //exit(1);
    printpol((g));
    printf("\n");
    //exit(1);
    // b.x[0]=2;
    // b.x[1]=9;
    b.x[0] = 0;
    a.x[0] = 1;
    vec vv ={0};
    v.x[1]=1;
    v.x[0]=0;
    //u = (v);
    // d.x[EXP] = 1;
    printpol((g));
    printf(" =g\n");
    //exit(1);

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
    //exit(1);
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
    for (i = 2; i < EXP + 1; i++)
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
    count = EXP + 1;
    for (i = 0; i < EXP + 1; i++)
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
    for (i = 0; i < O; i++)
        printf("%d,", gf[i]);
    printf("};");
    printf("\n");

    // exit(1);
}


unsigned short c2o(vec x){
    
    return xtrace((x),Pr);
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

    for (i = 0; i < k; i++)
    {
        for (j = 0; j < l; j++)
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

void de()
{

    int i, j;
    for (i = 0; i < O; i++)
    {
        printpol((dick[i]));
        printf(", %d, %d lueda\n", eval(dick[i], Pr), i);
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

unsigned short
v2a(oterm a)
{
    int i, j;

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
}

void printsage(vec a)
{
    int i, j, k;
    oterm b;

    printf("poly=");
    for (i = 0; i < DEG; i++)
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

    for (i = 0; i < DEG; i++)
    {
        if (f.x[i] > 0)
            j = i;
    }

    return j;
}

// 項の数
int terms(vec f)
{
    int i, count = 0;

    for (i = 0; i < DEG; i++)
        if (f.x[i] > 0)
            count++;

    return count;
}

// ０多項式かどうかのチェック
unsigned char
chk(vec f)
{
    int i, flg = 0;
    vec x = {0};

    x = (f);
    for (i = 0; i < DEG; i++)
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
    for (i = 0; i < k + 1; i++)
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
    // for (i = 0; i < T * 2; i++)
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

    for (i = 0; i < DEG; i++)
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

    for (j = 0; j < K + 1; j++)
        gg[j] = g[K - j];

    memcpy(g, gg, sizeof(g));
}

unsigned short
b2B(unsigned short b[E])
{
    int i;
    unsigned short a = 0;

    for (i = E - 1; i > -1; i--)
        a ^= (b[E - i - 1] << i);

    return a;
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

    for (i = 0; i < deg(v) + 1; i++)
    {
        e = v.x[i];
        d ^= (e << i);
    }

    return d;
}

void printvec(vec v)
{
    int i, j;

    for (i = 0; i < deg(v) + 1; i++)
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

    for (i = 0; i < u; i++)
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
    for (x = 0; x < N; x++)
    {
        z = 0;
        // #pragma omp parallel for reduction (^:z)
        for (i = 0; i < n + 1; i++)
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
    for (i = 0; i < DEG; i++)
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
MTX bdet()
{
    int i, j, k; 
    vec l={0};
    unsigned char dd[E * K] = {0};
    FILE *ff;
    MTX R = {0};

    // ff = fopen("Hb.key", "wb");

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K; j++)
        {
            l = o2c(dick[fg[mat[i][j]]]);
            // #pragma omp parallel for
            for (k = 0; k < E; k++)
            {
                R.x[i][j * E + k] = l.x[k];
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        // #pragma omp parallel for
        for (j = 0; j < E * K; j++)
        {
            printf("%d,", R.x[i][j]);
            // dd[j] = BH[j][i];
        }
        // fwrite(dd, 1, E * K, ff);
        printf("\n");
    }

    // fclose(ff);
    return R;
}

MTX bd2()
{
    int i, j, k, l;
    unsigned char dd[E * K] = {0};
    FILE *ff;
    vec v = {0};
    MTX R = {0};

    // ff = fopen("Hb.key", "wb");

    // memset(BB.z,0,sizeof(BB.z));
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K / 2 + 1; j++)
        {
            l = bm[i][j];
            printf("bm==%d %d\n", l, j);

            v = i2v(l);
            // #pragma omp parallel for
            for (k = 0; k < E; k++)
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

    for (i = 0; i < K; i++)
    {
        memset(v.x, 0, sizeof(v.x));
        v = i2v(c[i]);
        for (j = 0; j < E; j++)
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

    for (j = 0; j < K; j++)
    {
        for (i = 0; i < E; i++)
        {
            v.x[i] = s[count];
            count++;
        }
        n.v[i] = v2i(v);
    }
    n.f = K;

    return n;
}

// バイナリ型パリティチェック行列を生成する
void toBit(MTX L)
{
    int i, j, k, l;
    unsigned char dd[E * K] = {0};
    FILE *ff;

    // ff = fopen("Hb.key", "wb");

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K; j++)
        {
            l = L.x[i][j];
            printf("l=%d,", l);
            // #pragma omp parallel for
            for (k = 0; k < E; k++)
            {
                BB.x[i][j * E + k] = l % 2;
                l = (l >> 1);
            }
        }
        printf("\n");
    }
    // exit(1);
    /*
    for (i = 0; i < N; i++)
    {
        //#pragma omp parallel for
        for (j = 0; j < E * K; j++)
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

MTX toByte(MTX SH, int kk)
{
    vec v = {0};
    int i, j, k, cnt;
    MTX R = {0};

    memset(HH, 0, sizeof(HH));
    printf("HH=");
    // exit(1);
    for (i = 0; i < N; i++)
    {
        printf("%d\n", i);
        // #pragma omp parallel for
        for (j = 0; j < kk; j++)
        {
            cnt = 0;
            for (k = j * E; k < j * E + E; k++)
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
    for (i = 0; i < N; i++)
        P[i] = i;
    random_shuffle(P, SIZE_OF_ARRAY(P));
    // random_permutation(P);

    // for (i = 0; i < N; i++)
    //     P[i] = i;
    for (i = 0; i < N; i++)
        inv_P[P[i]] = i;
    // fwrite(P, 2, N, fp);
    // fclose(fp);
    /*
    for (i = 0; i < N; i++)
    printf ("%d,", P[i]);
    printf ("\n");
    exit(1);
    */
    // fp = fopen("inv_P.key", "wb");
    // fwrite(inv_P, 2, N, fp);
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

    for (i = 0; i < N; i++)
    {
        if (gf[mlt(i, i)] == u)
            return i;
    }

    printf("来ちゃいけないところに来ました\n");
    exit(1);
}

/*
EX extgcd(vec a, vec b)
{

    vec s = {0}, sx = {0}, sy = {0}, t = {0}, tx = {0}, ty = {0}, tmp = {0};
    EX c = {0};

    if (deg(b) > deg(a))
    {
        tmp = a;
        a = b;
        b = tmp;
    }
    s = a;
    t = b;
    sx.x[0] = 1;
    sx.x[0] = 0;
    ty.x[0] = 1;
    ty.x[0] = 0;

    //  vec temp={0};
    tmp = vmod(s, t);
    if (deg(tmp) == 0)
    {
        c.d = t;
        c.v = tx;
        c.u = ty;
        printf("ppp\n");
        return c;
    }
    while (deg(tmp) > 0)
    {
        printpol(((tmp)));
        printf(" ========vmod\n");
        vec temp = vdiv(s, t);
        vec u = vadd(s, vmul(t, temp));
        vec ux = vadd(sx, vmul(tx, temp));
        vec uy = vadd(sy, vmul(ty, temp));
        s = t;
        sx = tx;
        sy = ty;
        t = u;
        tx = ux;
        ty = uy;
        tmp = vmod(s, t);
    }
    printpol((tmp));
    printf(" ========vmod!\n");

    if (vLT(tmp).a == 1)
    {
        c.d.x[0] = 1;
        // c.d=t;
        c.v = tx;
        c.u = ty;
        printf("bbb\n");
        return c;
    }
    if (vLT(tmp).a == 0)
    {

        c.d = t;
        c.v = tx;
        c.u = ty;
        printf("ccc\n");

        return c;
    }
}
*/

// 512bitの秘密鍵を暗号化
void encrypt(char buf[], unsigned char sk[64])
{
    const uint8_t *hash = {0};
    sha3_context c = {0};
    int image_size = 512, i;
    FILE *fp;
    //  unsigned short dd=0;

    printf("plain text=");
    for (i = 0; i < 64; i++)
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

    for (i = 0; i < 64; i++)
    {
        printf("%d", hash[i]);
        // char s[3];
        // byte_to_hex(hash[i],s);

        sk[i] ^= hash[i];
    }
    printf("\nencrypt sk=");
    for (i = 0; i < 64; i++)
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

    for (i = 0; i < K; i++)
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
    for (i = j; i < 2 * T; i++)
    {
        if (v.x[i] > 0 && v.x[i] < N)
        {
            err[v.x[i]] = 1;
        }
    }

    char buf0[8192] = {0}, buf1[10] = {
                               0};

    // #pragma omp parallel for
    for (i = 0; i < N; i++)
    {
        snprintf(buf1, 10, "%d", err[i]);
        strcat(buf0, buf1);
    }
    // puts (buf0);
    printf("vector=%d\n", strlen(buf0));
    // exit(1);
    printf("cipher sk2=");
    for (i = 0; i < 64; i++)
        printf("%u,", sk[i]);
    printf("\n");

    sha3_Init256(&c);
    sha3_Update(&c, (char *)buf0, strlen(buf0));
    hash = sha3_Finalize(&c);

    j = 0;
    printf("hash=");
    for (i = 0; i < 64; i++)
    {
        printf("%d", hash[i]);
        // char s[3];
        // byte_to_hex(hash[i],s);

        sk[i] ^= hash[i];
    }
    printf("\ndecript sk=");
    for (i = 0; i < 64; i++)
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
            vb[i][j] = gf[mltn(i, fg[j])];
            printf("g%d,", vb[i][j]);
            mat[j - 1][i] = vb[i][j];
        }
        printf("\n");
    }
    // exit(1);
}

void vanshe(int kk)
{
    int i, j, k;

    printf("van der\n");

    for (i = 0; i < N; i++)
        vb[0][i] = 1;
    // #pragma omp parallel for private(i, j)
    for (i = 1; i < kk; i++)
    {
        for (j = 0; j < N; j++)
        {
            mat[j][i] = gf[mltn(i, j)];
        }
        // printf("\n");
    }
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < N; j++)
            printf("%d,", mat[j][i]);
        printf("\n");
    }
    // exit(1);
}

void ogt(unsigned short pp[], int kk)
{
    int i, j, k;
    vec w = {0};

#pragma omp parallel for private(i, j)
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

vec synd(unsigned short zz[], int kk)
{
    unsigned short syn[K] = {0}, s = 0;
    int i, j, t1;
    vec f = {0};

    printf("in synd2\n");

    for (i = 0; i < K; i++)
    {
        // syn[i] = 0;
        s = 0;
        // #pragma omp parallel num_threads(16)
        for (j = 0; j < N; j++)
        {
            if (zz[j] > 0)
            {
                s = plus(s, gf[mlt(fg[zz[j]], fg[mat[j][i]])]);
                // printf("%d,%d %d %d", s, mat[j][i], zz[j], j);
                // printf(" ==sind\n");
            }
        }
        syn[K - 1 - i] = s;
        printf("syn%d,\n", syn[K - 1 - i]);
    }
    // printf ("\n");
    for (i = 0; i < K; i++)
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

vec synd2(unsigned short zz[], int kk)
{
    unsigned short syn[K] = {0}, s = 0;
    int i, j, t1;
    vec f = {0};

    printf("in synd2\n");

    for (i = 0; i < kk; i++)
    {
        syn[i] = 0;
        s = 0;
        // #pragma omp parallel num_threads(16)
        for (j = 0; j < M; j++)
        {
            s ^= gf[mlt(fg[zz[j]], fg[mat[j][i]])];
        }
        syn[i] = s;
        // printf ("syn%d,", syn[i]);
    }
    // printf ("\n");

    f = setpol(syn, kk);
    printpol((f));
    printf(" syn=============\n");
    //  exit(1);

    return f;
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
    for (i = 0; i < 64; i++)
        sk[i] = i + 1;

    for (i = 0; i < N; i++)
    {
        snprintf(buf1, 10, "%d", zz[i]);
        strcat(buf, buf1);
    }
    // puts (buf);
    printf("vector=%u\n", strlen(buf));
    // exit(1);

    printf("sk0=");
    for (i = 0; i < 64; i++)
        printf("%u,", sk[i]);
    printf("\n");
    // exit(1);

    f = synd(zz, K);
    v = (f);
    // printf("v=");
    for (i = 0; i < K; i++)
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
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K; j++)
            mat[i][j] = dd[K * i + j];
    }
    fclose(fq);
}

// vec sx={0},ty={0};
int isquad(vec w)
{
    int i, j, flg = 0;
    vec b = {0};

    b = (w);
    for (i = 0; i < DEG; i++)
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
    for (i = 0; i < N; i++)
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
  for (i = 0; i < N; i++)
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

    for (i = 0; i < N; i++)
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

    for (j = 0; j < N; j++)
    {
        for (i = 0; i < K; i++)
        {
            ma[j][i] = gf[mlt(fg[vb[i][j]], tr[j])];
        }
        // printf("tr[%d]=%d\n",j,tr[j]);
    }

    unsigned short s;
#pragma omp parallel for defauvLT(none) private(i, j, k, s) shared(mat, gt, ma, gf, fg)
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < N; j++)
        {
            s = 0;

            for (k = 0; k < K; k++)
                s ^= gf[mlt(fg[gt[k][i]], fg[ma[j][k]])];
            // printf("%d,",s);
            mat[j][i] = s;
        }
        printf("\n");
    }
    printf("\n");
    // exit(1);

    for (j = 0; j < N; j++)
    {
        for (i = 0; i < K; i++)
            printf("%d,", mat[j][i]);
        printf("\n");
    }
    printf("\n");
    // wait();

    return w;
}

// BMA 専用（多項式と次元指定）
vec mkc(vec w, int kk)
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
    // r = vmul(w, w);

    memset(ta, 0, sizeof(ta));
    // w = setpol(g, K + 1);
    printpol((r));
    // printf(" =poly\n");

    // 多項式の値が0でないことを確認
    for (i = 0; i < N; i++)
    {
        ta[i] = eval(r, i);
        if (ta[i] == 0)
        {
            printf("eval 0 @ %d\n", i);
            // fail = 1;
            goto aa;
        }
    }
    for (i = 0; i < N; i++)
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
    van(kk);
    //  v=(w);
    ogt(g, kk);

    // wait();

    // #pragma omp parallel for

    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (j = 0; j < N; j++)
    {
        for (i = 0; i < kk; i++)
        {
            mat[j][i] = gf[mlt(fg[vb[i][j]], tr[j])];
        }
        // printf("tr[%d]=%d\n",j,tr[j]);
    }

    // printf("\n");
    // exit(1);
    /*
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < kk; i++)
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
    for (i = 0; i < N; i++)
    {
        ta[i] = eval(r, i);
        if (ta[i] == 0)
        {
            printf("eval 0 @ %d\n", i);
            // fail = 1;
            goto aa;
        }
    }
    for (i = 0; i < N; i++)
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

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < kk; j++)
        {
            mat[i][j] = vb[j][i];
        }
    }

    // printf("\n");
    // exit(1);
    /*
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < kk; i++)
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

    for (i = 0; i < N; i++)
        bm[i][0] = mat[i][0];
    for (i = 1; i < kk; i++)
    {
        for (j = 0; j < N; j++)
            bm[j][i] = mat[j][i * 2 - 1];
    }
    /*
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < kk+1; j++)
            printf("%d,", bm[i][j]);
        printf("  ==bm\n");
    }
*/
    // exit(1);
}

// Niederreiter暗号の公開鍵を作る(RS)
MTX mk_pub()
{
    int i, j, k, l;
    FILE *fp;
    unsigned char dd[E * K] = {0};
    vec w = {0};
    MTX Z = {0}, FX = {0}, G_bin = {0}, G_int = {0}, R = {0}; //, O = {0};

    w = mkd(w, K * 2);
    // w = mkg(K);
    half(K / 2 + 1);

    printpol(w);
    printf("\n");
    printsage((w));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");

    R = bd2();
    printf("deco_rev= ");
    for (i = 0; i < (K / 2 + 1) * E; i++)
        printf("%d,", R.x[1][i] ^ R.x[2][i] ^ R.x[3][i]);
    printf("\n");
    Pgen();
    do
    {
        memset(inv_S.x, 0, sizeof(inv_S.x));
        memset(S.x, 0, sizeof(S.x));
        for (i = 0; i < (K / 2 + 1) * E; i++)
        {
            for (j = 0; j < (K / 2 + 1) * E; j++)
                S.x[i][j] = xor128() % 2;
        }
    } while (mkRS(S, inv_S.x) == -1);
    // mkS();
    // O = toByte(R, K / 2 + 1);
    //   exit(1);
    // Z = mulmat(S, R, 2);
    printf("Zz=\n");
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < (K / 2 + 1) * E; i++)
        {
            // G.z[j][i] = Z.w[P[j]][i];
            G_bin.x[j][i] = Z.x[P[j]][i];
            // Z.x[j][i]=Z.x[j][i];
            // printf("%d.",inv_S.w[j][i]);
        }
        // printf("\n");
    }
    // printf("\n");

    FX = toByte(G_bin, K / 2 + 1);
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K / 2 + 1; j++)
            G_int.x[i][j] = FX.x[i][j];
    }

    return G_int;
}

// Niederreiter暗号の公開鍵を作る(Gvecpa)
MTX pk_gen()
{
    int i, j, k, l;
    FILE *fp;
    unsigned char dd[E * K] = {0};
    vec w = {0};
    MTX R = {0}, R_bin = {0}, Q = {0}, O_bin = {0};

    vv(K);
    //w = mkd(w, K);
    // w = mkg(K);
    // half(K / 2 + 1);

    printpol(w);
    printf("\n");
    printsage((w));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");

    //Q = bdet();
    //exit(1);
    //Pgen();
    printf("pre\n");
    /*
    do
    {
        memset(inv_S.x, 0, sizeof(inv_S.x));
        memset(S.x, 0, sizeof(S.x));
        for (i = 0; i < K * E; i++)
        {
            for (j = 0; j < K * E; j++)
                S.x[i][j] = xor128() % 3;
        }
    } while (is_reg2(S, inv_S.x) == -1);
    exit(1);
    */
   //do
   do
   {
   for(int i=0; i<F;i++){
    for (j = 0; j < F; j++){
        S.x[i][j] = xor128() % Pr;
        printf("%d ",S.x[i][j]);
    }
    printf("\n");
   }
   printf("\n");
   
   }while(mainga(S).col==0);
   inv_S=mainga(S);
   for(i=0;i<F;i++){
    for(j=0;j<F;j++)
    printf("%d,",inv_S.x[i][j]%Pr);
    printf("\n");
   }
   printf("\n");
   exit(1);
   
   MTX XA={0};
   unsigned short u;
   for(i=0;i<F;i++){
    for(j=0;j<F;j++){
        u=0;
        for(k=0;k<F;k++)
        u+=(S.x[i][k]*inv_S.x[k][j])%Pr;
        XA.x[i][j]=u%Pr;
    }
   }
   for(i=0;i<F;i++){
    for(j=0;j<F;j++)
    printf("%d ",XA.x[i][j]);
    printf("\n");
   }
   printf("\n");
   //while(is_reg(S,inv_S.x)==0);
    
    exit(1);
    
    // makeS();
    //   exit(1);
    // H = mulmat(S, Q, 1);
    for (i = 0; i < K * E; i++)
    {
        for (j = 0; j < N; j++)
        {
            // O.z[j][i] = H.z[P[j]][i];
            O_bin.x[j][i] = H.x[P[j]][i];
        }
    }
    R = toByte(O_bin, K);
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K * E; j++)
        {
            R_bin.x[i][j] = O_bin.x[i][j];
            // printf("%d,",R_bin.x[i][j]);
        }
    }

    // exit(1);

    return R;
}

vec dec(unsigned short ss[])
{
    int i, j, k;
    vec v = {0};
    vec s = {0};
    unsigned ch[K * E] = {0};
    unsigned char h2o[K * E] = {0};

    printf("!1\n");
    for (i = 0; i < K; i++)
    {
        v = i2v(ss[i]);
        for (j = 0; j < E; j++)
            ch[i * E + j] = v.x[j];
    }
    for (i = 0; i < K * E; i++)
        printf("%d", ch[i]);
    printf("\n");

    unsigned short uk[K] = {0};

    for (i = 0; i < K * E; i++)
    {
        for (j = 0; j < K * E; j++)
            h2o[i] ^= (ch[j] & inv_S.x[i][j]);
    }
    // for (i = 0; i < K * E; i++)
    // printf("%d,", h2o[i]);
    // printf("\n");

    for (i = 0; i < K; i++)
    {
        memset(v.x, 0, sizeof(v.x));
        for (j = 0; j < E; j++)
            v.x[j] = h2o[i * E + j];
        uk[i] = v2i(v);
    }
    for (i = 0; i < K; i++)
        printf("%d,", uk[i]);
    printf("\n");
    //    exit(1);
    s = setpol(uk, K);

    return s;
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
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K; j++)
            dd[j] = mat[i][j];
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

    for (i = 0; i < T; i++)
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
    for (i = 0; i < T; i++)
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
    for (i = 0; i < N; i++)
        yy[i] = x[P[i]];
    for (i = 0; i < N; i++)
    {
        if (yy[i] > 0)
            printf("err= %d\n", i);
    }

    return count;
}

int ero(vec v)
{
    int i, j, count = 0;

    for (i = 0; i < T * 2; i++)
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
        for (k = 0; k < N; k++)
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
      for (k = 0; k < N; k++)
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
        B++;
    }
    if (count < T * 2)
    {
        printf("error is too few\n");

        AA++;
        // memcpy (zz, z1, sizeof (zz));
        /*
      printf("{");
      for (i = 0; i < D; i++)
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

    for (i = 0; i < T; i++)
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
    for (i = 0; i < N; i++)
        ya[i] = xa[P[i]];
    for (i = 0; i < N; i++)
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
        B++;
    }
    if (count < T)
    {
        printf("error is too few\n");

        AA++;
        // memcpy (zz, z1, sizeof (zz));
        /*
      printf("{");
      for (i = 0; i < D; i++)
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

    for (i = 0; i < T; i++)
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
    for (i = 0; i < N; i++)
        ya[i] = xa[i];
    for (i = 0; i < N; i++)
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
        B++;
    }
    if (count < T)
    {
        printf("error is too few\n");

        AA++;
        // memcpy (zz, z1, sizeof (zz));
        /*
      printf("{");
      for (i = 0; i < D; i++)
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
        l = rand() % (N);
        // printf ("l=%d\n", l);
        if (0 == z1[l] && l>0)
        {
            z1[l] = l;
            // printf("l=%d\n", l);
            j++;
        }
    }
}

void fun()
{
    unsigned short i, k;

    vec qq = {0};
    for (i = 0b10001; i < 0b11111 + 1; i++)
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

    for (i = 0; i < K / 2 + 1; i++)
    {
        t[i] = e[i];
        printf("e=%d\n", e[i]);
    }
    for (i = 0; i < K / 2 + 1; i++)
    {
        printf("t=%d\n", t[i]);
    }
    // exit(1);

    v.x[0] = t[0];
    v.x[1] = t[1];
    k = 2;
    for (i = 2; i < K; i++)
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

vec bfd(unsigned short ss[])
{
    int i, j, k, count = 0;
    vec v = {0};
    vec s = {0};
    unsigned ch[K * E * 2] = {0};
    unsigned char h2o[K * E * 2] = {0};

    // count=(K/2+1)*E-1;
    for (i = 0; i < (K / 2) + 1; i++)
    {
        v = i2v(ss[i]);
        for (j = 0; j < E; j++)
        {
            ch[count] = v.x[j];
            count++;
        }
    }
    printf("bfd_bin=\n");
    for (i = 0; i < (K / 2 + 1) * E; i++)
        printf("%d", ch[i]);
    printf("\n");
    // exit(1);

    unsigned short uk[K] = {0};

    for (i = 0; i < (K / 2 + 1) * E; i++)
    {
        for (j = 0; j < (K / 2 + 1) * E; j++)
            h2o[i] ^= (ch[j] & inv_S.x[i][j]);
    }

    printf("deco_bin=\n");
    for (i = 0; i < (K / 2 + 1) * E; i++)
        printf("%d,", h2o[i]);
    printf("\n");

    // count=(K/2+1)*E-1;

    count = 0;
    for (i = 0; i < (K / 2) + 1; i++)
    {
        memset(v.x, 0, sizeof(v.x));
        for (j = 0; j < E; j++)
        {
            v.x[j] = h2o[count];
            count++;
        }
        uk[i] = v2i(v);
    }

    printf("bm_int=\n");
    for (i = 0; i < K / 2 + 1; i++)
        printf("%d,", bm[1][i] ^ bm[2][i] ^ bm[3][i]);
    printf("\n");

    printf("u_int=\n");
    for (i = 0; i < K / 2 + 1; i++)
        printf("%d,", uk[i]);
    printf("\n");
    // exit(1);

    s = setpol(uk, K / 2 + 1);
    v = (s);

    return v;
}

vec sin2(unsigned short zz[], MTX R)
{
    int i, j;
    vec s = {0};
    vec v = {0};
    unsigned short ss[K] = {0};

    for (i = 0; i < N; i++)
    {
        if (zz[i] > 0)
        {
            for (j = 0; j < K; j++)
            {
                // ss[j]
                v.x[j] ^= R.x[i][j];
                printf("%d,", R.x[i][j]);
            }
        }
        printf("\n");
    }

    return v;
}

unsigned short logx(unsigned short u)
{
    unsigned short i;

    return oinv(u);

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
    for (i = 0; i < k + 1; i++)
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

    for (i = 0; i < k; i++)
    {
        b.x[i] = (Pr - d.x[i]) % Pr;
        // printf("chen_b=%d\n",b.x[i]);
    }
    u = xtrace((b), Pr);

    return u;
}

typedef struct
{
    vec f;
    vec g;
    vec h;
} ymo;

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
        myu=0;
        for (i = 0; i <= d; i++)
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
        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < 2; j++)
            {
                for (k = 0; k < 2; k++)
                    U2[1][i][j] = (vadd2((U2[1][i][j]), (pmul(U1[i][k], U2[0][k][j]))));
                
            }
        }
        memcpy(U2[0], U2[1], sizeof(U2[0]));
        memset(U2[1], 0, sizeof(U2[1]));
    }
    t.f = U2[0][0][0];
    t.g = U2[0][1][0];
    t.h =U2[0][0][1];
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

vec bib(vec x){
vec f={0};

for(int i=1;i<deg(x)+1;i++)
f.x[i-1]=gf[mlt(fg[i],fg[x.x[i]])];

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
    for (i = 0; i < T; i++)
    {
        w[i].x[0] = 1;
        w[i].x[1] = a.x[i];
        printpol((w[i]));
        printf("\n");
    }
    tmp.x[0] = 1;
    //

    // #pragma omp parallel for private(i,j)
    for (i = 0; i < T; i++)
    {
        t = (tmp);
        //
        for (j = i; j < T; j++)
        {
            // #pragma omp critical
            if (i != j)
            {
                t = (vmul((w[i]), (w[j])));
                printf("wa=%d %d\n",i,j);
                printpol(t);
                printf("  ==t\n");
            }
        l = vadd(l, t);
        printpol(l);
        printf(" ===l\n");
        }

        printpol((l));
        printf(" ==snach\n");
        //exit(1);
        if (deg((l)) == 0)
        {
            printf("baka9\n");
            // exit(1);
        }
    }
    //exit(1);
    return (l);
}

/*
vec yes(ymo t){
    vec e=vmul(t.g,t.h);
    vec b=bib((t.f));
    vec b2 = bibun_old((t.f));
    printsage((t.f));
    printf("\n");
    printsage((b));
    printf("\n");
    vec ixi={0};
    int count=0;
    if(deg(e)==0){
        if(deg(t.g)==0)
        printf("muno\n");
        if(deg(t.h)==0)
        printf("devil\n");
        return t.g;
    }
    for(int i=0;i<N;i++){
    if(eval((e),i)==1)
    {
        printf("ixi=%d %d\n",gf[oinv(gf[mlt(fg[eval((b),i)],fg[eval((t.g),i)])])],i);  //gf[mlt(oinv(eval((b),i)),fg[eval((t.h),i)])],i);
        ixi.x[i]=gf[mlt(oinv(eval((b),i)),fg[eval((t.h),i)])];
    }
    }

return ixi;
}
*/

vec yes(ymo t){
    vec e=vmul(t.g,t.h);
    vec b2=bib((t.f));
    vec b=bibun_old((t.f));
    printsage((t.f));
    printf("\n");
    printsage((b));
    printf("\n");
    vec ixi={0};
    int count=0;
    for(int i=0;i<N;i++){
    if(eval((e),i)==1)
    //if(a[i]>0)
    {
        printf("ixi=%d %d\n",gf[oinv(gf[mlt(fg[eval((b),i)],fg[eval((t.g),i)])])],i);  //gf[mlt(oinv(eval((b),i)),fg[eval((t.h),i)])],i);
        ixi.x[i]=gf[mlt(oinv(eval((b),i)),fg[eval((t.h),i)])];
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
        for (i = 1; i < L; i++)
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
        for (i = 0; i < T; i++)
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
    for (j = 1; j < T * 2 + 1; j++)
    {
        v = (lo[j - 1]);
        k = 0;
        // printpol(v);
        // printf(" ==lo\n");

        l = deg((lo[j - 1]));
        for (i = 1; i < l + 1; i++)
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
    for (i = 0; i < deg(x) + 1; i++)
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
            // for (i = 0; i < N; i++)
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
    for (i = 0; i < j; i++)
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
    for (i = 0; i < count; i++)
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

    for (j = 0; j < N; j++)
    {
        if (zz[j] > 0)
        {
            // memcpy(syn, bm[j], sizeof(syn));
            printf("bm_in_sen= %d || ", j);
            for (i = 0; i < K / 2 + 1; i++)
            {
                syn[i] ^= bm[j][i];
                printf("%d,", syn[i]);
            }
        }
    }
    printf("\n");
    v = newhalf(syn);
    // printf("%d\n",j);
    for (k = 0; k < kk; k++)
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

vec sendrier2(unsigned short zz[N], MTX L)
{
    unsigned short syn[K + 1] = {0}, s[K + 1] = {0}, rt[K / 2 + 1] = {0}, uu[(K / 2 + 1) * E] = {0}, es[(K / 2 + 1) * E] = {0};
    int i, j, k, count = 0;
    vec f = {0}, w = {0};
    vec v = {0}, x = {0}, u = {0}, t = {0};
    unsigned short tmp[(K / 2 + 1) * E] = {0}, m[K + 1] = {0};

    for (j = 0; j < N; j++)
    {
        if (zz[j] > 0)
        {
            // for(i=0;i<(K/2)+1;i++){

            // memcpy(syn, L.w[j], sizeof(syn));

            for (k = 0; k < K / 2 + 1; k++)
            {
                rt[k] ^= L.x[j][k];
                // rt[k] = bm[j][k];
            }
        }
    }
    for (i = 0; i < K / 2 + 1; i++)
        printf("%d,", rt[i]);
    printf("\n");
    // exit(1);

    x = bfd(rt);
    for (i = 0; i < K / 2 + 1; i++)
    {
        u.x[K / 2 - i] = x.x[i];
    }
    for (i = 0; i < K / 2 + 1; i++)
        printf("%d,%d\n", u.x[i], x.x[i]);
    printf("\n");
    // exit(1);

    printf("rt=\n");
    for (i = 0; i < K / 2 + 1; i++)
        printf("%d,", u.x[i]);
    printf("\n");
    printf("bm_in se2 == %d || ", j);
    // exit(1);
    v = newhalf(u.x);

    printf("P= ");
    for (i = 0; i < N; i++)
        printf("%d,", P[i]);
    printf("\n");

    memset(s, 0, sizeof(s));
    for (i = 0; i < K; i++)
        s[i + 1] = v.x[i];

    printf("rt_deco= ");
    f = bma(s, K);
    x = chen(f);
    ero3(x);
    exit(1);
    wait();

    f = setpol(syn, K);

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
    for (i = 0; i < N; i++)
    {
        if (gf[i] == a)
            return i;
    }
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
            a.x[i][j] = gf[mlt(fg[a.x[i][j]], oinv(p))];
        }

        for (j = 0; j < K / 2; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (k = i; k < (K / 2 + 1); k++)
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
    for (i = 0; i < K / 2; i++)
    {
        for (j = 0; j < K / 2 + 1; j++)
            printf("@%d,", a.x[i][j]);
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
    for (i = 0; i < N; i++)
    {
        // v.x[i] = 0;
        if (xtrace(pol, i) == 0)
        {
            printf("i=%d\n", i);
            vv.x[i] = 1;
        }
    }
    exit(1);

    return vv;
}

void vv(int kk)
{
    int i, j;
    vec r; //= mkpol();
    unsigned short tr[N];
    unsigned short ta[N] = {0};

    printf("van der\n");

    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < N; j++)
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
    for (i = 0; i < N; i++)
    {
        ta[i] = eval(r, i);
        if (ta[i] == 0)
        {
            printf("eval 0 @ %d\n", i);
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
        for (j = 0; j < K; j++)
        {
            mat[i][j] = gf[mlt(fg[vb[j][i]], tr[i])];
        }
    }
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < N; j++)
            printf("c%d,", mat[j][i]);
        printf("\n");
    }
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

MTX A2M(unsigned short A[N][K])
{
    int i, j;
    MTX J = {0};

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K; j++)
            J.x[i][j] = A[i][j];
    }

    return J;
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

// 言わずもがな
int main(void)
{
    int i, a, b, c;
    unsigned short zz[N] = {0};
    vec f = {0}, r = {0}, w = {0}, r1 = {0};
    vec v = {0}, x = {0};
    MTX R = {0}, OO = {0};
    unsigned short s[K + 1] = {0, 3, 0, 6, 0};

    if (K > N)
        printf("configuration error! K is bigger than N\n");

    // （謎）
    memset(mat, 0, sizeof(mat));
    srand(clock());

    mkmf();
    makefg();
      de();
    //exit(1);
    
    //printf("1331=%d\n",c2o(o2c(dick[fg[1331]])));
    //exit(1);
    //  次元KのGvecpa符号の検査行列を生成する
    //van(K);
    //vv(K);
    //OO=pk_gen();
    vec v1 = {0}, v2 = {0}, x1 = {0}, x2 = {0};
        int it=0;
        while (1)
        {
        //   エラーベクトルの初期化
            memset(zz, 0, sizeof(zz));
        //   重みTのエラーベクトルを生成する
            mkerr(zz, T);
            //zz[1]=1;
            //zz[2]=1;

            r = synd(zz, K);
            // exit(1);
            printpol(r);
            //exit(1);
            x = (r);
            for (i = 0; i < N; i++)
            {
                if (zz[i] > 0)
                    printf("aie=%d\n", i);
            }
            ymo t = bm_itr(x.x);
            printpol((t.f));
            printf("\n");
 
            x = chen((t.f));
            //exit(1);
            yes(t);
            //exit(1);
            
            vec rr = {0};
            for (i = 0; i < T; i++)
            {
                rr.x[x.x[i]] = 1; //x.x[i];
                printf("x=%d\n", x.x[i]);
            }
            for (i = 0; i < N; i++)
            {
                if ((zz[i] > 0) && zz[i] != rr.x[i])
                {
                    printf("cheki=%d %d %d\n", i, zz[i], rr.x[i]);
                        printf("baka\n");
                        exit(1);
                }
            }
            exit(1);
            it++;
            if(it==1000){
                printf("いいっ！\n");
            break;
            }
        }
    return 0;
}
