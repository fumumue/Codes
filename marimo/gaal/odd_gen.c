#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "global-p.h"
#include "struct.h"
#include "debug.c"
#include "chash-p.c"

#define O 9 // 1331,2197,4913,6859,3125,2187
#define EXP 2
#define Pr 3

// sagemath上での原始多項式
unsigned short pp[8][4] = {{0, 0, 9, 2}, {0, 0, 11, 2}, {0, 0, 16, 3}, {0, 0, 15, 2}, {0, 0, 1, 2}, {0, 1, 0, 2}, {0, 0, 1, 1}, {0, 0, 1, 2}};
// {0,0,9,2}, {1,0,11,2}, {1,0,16,3}, {1,0,15,2};
// GF(11^3,13^3,17^3,19^3)
unsigned short ff[2][4] = {0, 2, 0, 2}; // GF(3^7,5^5)
static unsigned short g[K + 1] = {0};
// unsigned short gf[O] = {0}, fg[O] = {0};
//  int N =0,M=0;
//  unsigned short c[K+1]={0};

// OP型からベクトル型への変換
vec o2v(OP f)
{
  vec a = {0};
  int i;

  // #pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0 && f.t[i].n < DEG)
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
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      f.t[j].n = i;
      f.t[j++].a = a.x[i];
    }
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

void op_print_raw(const OP f)
{
  puts("op_print_raw:");
  for (int i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0)
      printf("[%d] %ux^%u\n", i, f.t[i].a, f.t[i].n);
  }
}

bool op_verify(const OP f)
{
  bool end = false;
  unsigned short n_max = 0;
  for (int i = 0; i < DEG; i++)
  {
    if (end && (f.t[i].n != 0 || f.t[i].a != 0))
    {
      op_print_raw(f);
      printf("found data after end: i=%d\n", i);
      print_trace();
      fflush(stdout);
      return false;
    }
    if (f.t[i].a == 0)
    {
      end = true;
      continue;
    }
    if (f.t[i].n + 1 <= n_max)
    {
      op_print_raw(f);
      printf("found invalid order: i=%d\n", i);
      print_trace();
      fflush(stdout);
      return false;
    }
    n_max = f.t[i].n + 1;
  }
  return true;
}

// 多項式の次数(default)
int deg(vec a)
{
  int i, n = 0;

  // #pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
      n = i;
  }

  return n;
}

// 配列からベクトル表現の多項式へ変換する

// 配列の値を係数として多項式に設定する
OP setpol(unsigned short f[], int n)
{
  OP g;
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

  g = v2o(a);

  return g;
}

// 多項式を表示する(default)
void printpol(vec a)
{
  int i, n;

  n = deg(a);

  // printf ("baka\n");
  assert(("baka\n", n >= 0));

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

// 多項式の代入値
unsigned short
xtrace(OP f, unsigned short x)
{
  int i, d;
  unsigned short u = 0, v = 1;

  d = deg(o2v(f));
  // printpol(o2v(f));
  // printf(" =ff\n");

  for (i = 0; i < d + 1; i++)
  {
    v = 0;
    if (f.t[i].a > 0)
    {
      v = 1;
      for (int j = 0; j < f.t[i].n; j++)
      {
        v = (v * x);
      }
      v = (v * f.t[i].a) % O;

      // printf("\nv=%d",v);
    }
    u = (u + v);
  }
  // printf("u=%d\n",u%O);

  return u % O;
}

void makefg(int n)
{
  unsigned short i, j, count = 0;

  for (i = 0; i < O; i++)
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

// 多項式の次数(degのOP型)
int odeg(OP f)
{
  int i, j = 0, k;

  // k=terms(f);
  for (i = 0; i < 512; i++)
  {
    if (j < f.t[i].n && f.t[i].a > 0)
      j = f.t[i].n;
  }

  return j;
}

// 20200816:正規化したいところだがうまく行かない
// 多項式の足し算
OP oadd(OP f, OP g)
{
  vec a = {0}, b = {0}, c = {0};
  int i, j, k, l = 0;
  OP h = {0}, f2 = {0}, g2 = {0};

  // for(i=0;i<257;i++)
  //  printf("%d %d %d %d %d\n",i,f.t[i].a,f.t[i].n,g.t[i].a,g.t[i].n);

  //  exit(1);
  // f=conv(f);
  // g=conv(g);

  a = o2v(f);
  // exit(1);
  b = o2v(g);

  j = deg(o2v(f));
  l = deg(o2v(g));
  printpol(o2v(f));
  printf(" =f\n");
  printpol(o2v(g));
  printf(" =g\n");
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
    // if(f.t[i].a>0 || g.t[i].a>0)
    // h.t[i].a=f.t[i].a^g.t[i].a;
    c.x[i] = (a.x[i] + b.x[i]) % Pr;
  }
  //
  h = v2o(c);
  printpol(o2v(h));
  printf("\n");

  return h;
}

// 多項式を項ずつ掛ける
OP oterml(OP f, oterm t)
{

  // assert (op_verify (f));
  int i, k, j;
  OP h = {0};
  vec test;
  unsigned short n;

  // f=conv(f);
  k = odeg(f);
  j = 0;
  for (i = 0; i < k + 1; i++)
  {
    h.t[i].n = f.t[i].n + t.n;
    h.t[i].a = (f.t[i].a * t.a) % Pr;
  }

  // h=conv(h);
  // assert (op_verify (h));
  return h;
}

// OP型を正規化する
OP conv(OP f)
{
  vec v = {0};
  OP g = {0};

  v = o2v(f);
  g = v2o(v);

  return g;
}

// 多項式の掛け算
OP omul(OP f, OP g)
{
  // f=conv(f);
  // g=conv(g);

  // assert (op_verify (f));
  // assert (op_verify (g));
  int i, count = 0, k, l, m;
  oterm t = {0};
  OP h = {0}, e = {0}, r = {0};
  vec c = {0};

  l = deg(o2v(f));
  m = deg(o2v(g));

  g = conv(g);
  for (i = 0; i < 2; i++)
    printf("%d\n", o2v(g).x[i]);
  //  exit(1);

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
  printpol(o2v(f));
  printf(" =f\n");
  printpol(o2v(g));
  printf(" =g\n");
  // exit(1);

  for (i = 0; i < k; i++)
  {
    t = g.t[i];
    if (t.a > 0)
    {
      printf("t[%d]=%d,%d\n", i, t.a, t.n);
      e = oterml(f, t);
      printpol(o2v(e));
      printf(" =e\n");
      printpol(o2v(h));
      printf(" =h\n");
      h = oadd(h, e);
    }
  }
  printpol(o2v(h));
  printf(" =h2\n");

  // printpol(o2v(g));
  // printf(" =g\n");
  //    exit(1);
  // assert (op_verify (h));
  return h;
}

/*
//多項式の掛け算
OP
omul (OP f, OP g)
{
  f = conv (f);
  g = conv (g);
  assert (op_verify (f));
  assert (op_verify (g));
  int i, count = 0, k;
  oterm t = { 0 };
  OP h = { 0 }, e = {
    0
  }, r = {
    0
  };
  vec c = { 0 };
  if (odeg ((f)) > odeg ((g)))
    {
      k = odeg ((f));
    }
  else
    {
      k = odeg ((g));
    }
  for (i = 0; i < k + 1; i++)
    {
      t = g.t[i];
      e = oterml (f, t);
      h = oadd (h, e);
    }
  assert (op_verify (h));
  return h;
}
*/

// リーディグタームを抽出(default)
oterm LT(OP f)
{
  int i, k;
  oterm t = {0};

  k = deg(o2v(f));
  for (i = 0; i < DEG; i++)
  {
    // printf("a=%d %d\n",f.t[i].a,f.t[i].n);
    if (f.t[i].a > 0)
    {
      t.n = f.t[i].n;
      t.a = f.t[i].a;
    }
  }

  return t;
}

OP osub(OP f, OP g)
{
  vec a = {0}, b = {0}, d = {0};
  int i, k, l, m;
  OP ans = {0};

  a = o2v(f);
  b = o2v(g);
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
    d.x[i] = (a.x[i] - b.x[i]) % Pr;
    // if(d.x[i]<0)
    // d.x[i]+=P;
  }
  ans = v2o(d);

  return ans;
}

OP confer(OP f, int a)
{
  vec r;
  int n, i;
  OP g;

  r = o2v(f);
  n = deg(r);
  for (i = 0; i < n + 1; i++)
    r.x[i] = (r.x[i] * a) % Pr;
  g = v2o(r);

  return g;
}

int oequ(OP f, OP g)
{
  vec v, x;
  int i, flg = 0;

  v = o2v(f);
  x = o2v(g);
  for (i = 0; i < 512; i++)
  {
    if (v.x[i] != x.x[i])
      return -1;
  }

  return 0;
}

OP mkpol()
{
  int i, j, k, fail, flg, l, ii = 0;
  OP w = {0};

  do
  {
    fail = 0;
    j = 0;
    k = 0;
    flg = 0;
    l = 0;
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

//
unsigned short vb[K * 2][N] = {0};
unsigned short gt[K * 2][K * 2] = {0};

void van(int kk)
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
      vb[i][j] = gf[mltn(i, fg[j])];
      printf("%d,", vb[i][j]);
    }
    printf("\n");
  }
}

void ogt(unsigned short pp[], int kk)
{
  int i, j, k;
  OP w = {0};

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

// aに何をかけたらbになるか
unsigned short
equ(unsigned short a, unsigned short b)
{
  int i;

  return gf[mlt(oinv(a), fg[b])];
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
      printf("B('a^%d')*X**%d+", j, i); // for GF(2^m)
    }
  }
}

OP mkd(OP w, int kk)
{
  int i, j, k, l, ii = 0;

  unsigned short tr[N] = {0};
  unsigned short ta[N] = {0};
  vec v = {0};
  unsigned short po[K + 1] = {1, 0, 1, 0, 5};
  // OP w={0};
  OP r = {0};

aa:

  // printf("\n");
  memset(mat, 0, sizeof(mat));
  // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
  // 既約多項式しか使わない。

  l = -1;
  ii = 0;
  // irreducible goppa code (既役多項式が必要なら、ここのコメントを外すこと。)
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
  // separable goppa code
  w = mkpol();
  r = w;
  //  r=omul(w,w);
  memset(ta, 0, sizeof(ta));
  // w = setpol(g, K + 1);
  printpol(o2v(r));
  // printf(" =poly\n");

  // 多項式の値が0でないことを確認
  for (i = 0; i < N; i++)
  {
    ta[i] = xtrace(r, i);
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
  memset(g, 0, sizeof(g));
  g[0] = 1;

  // 多項式を固定したい場合コメントアウトする。
  printpol(o2v(r));
  printf("\n");
  printsage(o2v(r));
  printf("\n");
  printf("sagemath で既約性を検査してください！\n");
  memset(v.x, 0, sizeof(v.x));
  //  v=rev(w);
  van(kk);
  //  v=o2v(w);
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

OP dick[N] = {0}, func[N], a = {0}, cc = {0};
void mkmf()
{
  int i, j, k, count = 0;
  OP f = {0}, g = {0}, h = {0}, w = {0}, s = {0}, u = {0};
  vec b = {0}, a = {0}, d = {0}, t = {0}, v = {0};
  oterm o;
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

  g = setpol(ccp, 4);
  printpol(o2v(g));
  printf("\n");
  // exit(1);
  // b.x[0]=2;
  // b.x[1]=9;
  b.x[0] = 0;
  a.x[0] = 1;
  v.x[1] = 1;
  u = v2o(v);
  d.x[EXP] = 1;
  printpol(o2v(g));
  printf(" =g\n");
  // exit(1);

  // g=v2o(b);
  s = v2o(v);
  // s.t[1].a=1;
  // s.t[1].n=1;
  // gf[12]=P;
  // gf[13]=P*P;
  OP first;
  first.t[0].a = 1;
  first.t[0].n = 0;
  w = g;
  printpol(o2v(w));
  printf(" =w\n");
  printpol(o2v(g));
  printf(" =g\n");
  printpol(o2v(s));
  printf(" =s\n");

  printf("\n");
  // for(i=0;i<P;i++)
  gf[0] = 0;
  gf[1] = 1;
  dick[0] = v2o(b);
  // gf[EXP]=Pr;
  dick[1] = v2o(a);
  for (i = 2; i < EXP + 1; i++)
  {
    gf[i] = gf[i - 1] * Pr;
    dick[i] = omul(s, dick[i - 1]);
    printf("gf %d,",gf[i]);
  }

  // gf[2]=Pr;
  // gf[3]=Pr*Pr;
  gf[EXP + 1] = xtrace(g, Pr);
  printf("\naa=%d\n", gf[Pr]);
  // exit(1);
  // w=omul(w,s);
  // gf[12]=1111;
  dick[EXP + 1] = w;
  count = EXP + 2;
  while (1)
  {
    g = omul(g, s);
    printpol(o2v(g));
    printf(" =g\n\n");
    printf(" papaya\n");

    // exit(1);

    o = LT(g);
    memset(d.x, 0, sizeof(d));

    if (o.n == EXP)
    {
      d.x[o.n] = o.a;
      h = v2o(d);
      g = osub(g, h);
      f = confer(w, o.a);
      g = oadd(g, f);
      // w=omod(w,u);
      printpol(o2v(f));
      printf("\n");
    }
    dick[count] = g;
    gf[count] = xtrace(g, Pr);
    printf("count2=%d %d ", count, gf[count]);
    printpol(o2v(g));
    printf(" =gg\n\n");
    if (gf[count] == 1)
    {
      printf("count!=%d\n", count);
      break;
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

void de()
{

  int i, j;
  for (i = 0; i < O; i++)
  {
    printpol(o2v(dick[i]));
    printf(" %d\n", i);
  }
}

OP kei(unsigned short u, OP g)
{
  vec v = {0};
  int j, i = 0;
  OP f = {0};

  v = o2v(g);
  i = deg(v);
  printpol(v);
  printf("\n");
  for (j = 0; j < i + 1; j++)
    v.x[j] = (v.x[j] * u) % Pr;
  f = v2o(v);
  printpol(v);
  printf("\n");
  // exit(1);

  return f;
}

/*
//aに何をかけたらbになるか
unsigned short
equ (unsigned short a, unsigned short b)
{
  int i;


  for (i = 0; i < N; i++)
    {
      if ((a*i)%Pr == b)
  break;
    }
  return i;
}
*/

// invert of integer
unsigned short inv(unsigned short a, unsigned short n)
{
  unsigned short d, x, s, q, r, t, gcd;
  d = n;
  x = 0;
  s = 1;
  while (a != 0)
  {
    q = d / a;
    r = d % a;
    d = a;
    a = r;
    t = (x - q * s) % Pr;
    x = s;
    s = t;
  }
  gcd = d % Pr;

  return ((x + n) % (n / d)) % Pr;
}

// 多項式を単行式で割る
oterm LTdiv(OP f, oterm t)
{
  oterm tt = {0}, s = {
                      0};

  tt = LT(f);
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
    s.a = (tt.a * inv(t.a, Pr)) % Pr;
    s.n = tt.n;
  }
  else
  {
    printf("debug in LTdiv\n");
    exit(1);
  }

  return s;
}

// 多項式の剰余を取る
OP omod(OP f, OP g)
{
  int i = 0, j, n, k;
  OP h = {0}, e = {
                  0};
  oterm a, b = {0}, c = {0};

  n = LT(g).n;

  //  assert (("baka^\n", LT (f).n != 0));

  //  assert (("baka(A)\n", LT (g).n != 0));

  if (LT(f).n < LT(g).n)
  {
    //    exit(1);
    return f;
  }

  // printf ("in omod\n");
  // exit(1);

  k = LT(g).n;
  b = LT(g);
  OP ll;

  assert(("double baka\n", b.a != 0 && b.n != 0));
  while (LT(f).n > 0 && LT(g).n > 0)
  {

    c = LTdiv(f, b);
    h = oterml(g, c);
    printpol(o2v(f));
    printf("======f_before_omod\n");
    printpol(o2v(h));
    printf("======h_before_omod\n");
    f = osub(f, (h));
    printpol(o2v((h)));
    printf(" =====h_minus_omod\n");
    printpol(o2v(f));
    printf(" =====f_after_omod\n");
    // exit(1);
    if (odeg((f)) == 0 || odeg((g)) == 0)
    {
      //      printf("blake1\n");
      break;
    }

    if (c.n == 0 || b.n == 0)
      break;
    if (LT(f).a == 4 && deg(o2v(f)) == 0)
      exit(1);
  }
  printpol(o2v(f));
  printf("\n");
  // exit(1);

  return f;
}

OP hyoe[1331] = {0};
unsigned short sue[1331] = {0};
unsigned short sie[1331] = {0};
void tas()
{
  int i, j, k;
  unsigned short ccp[4] = {0, 0, 9, 2};
  OP g = {0}, h = {0}, f = {0}, m = {0};
  vec v = {0}, x = {0}, w = {0};

  v.x[1] = 9;
  v.x[0] = 2;
  m = v2o(v);
  printpol(v);
  printf("\n");
  g = m;
  w.x[1] = 1;
  printpol(w);
  printf("\n");

  h = v2o(w);
  // g=omul(g,h);
  printpol(o2v(g));
  printf("\n");
  //  exit(1);

  memset(v.x, 0, sizeof(v.x));
  for (i = 0; i < 2; i++)
  {
    v.x[0] = i;
    hyoe[i] = v2o(v);
  }
  memset(v.x, 0, sizeof(v.x));
  v.x[1] = 1;
  hyoe[2] = v2o(v);
  memset(v.x, 0, sizeof(v));
  v.x[2] = 1;
  hyoe[3] = v2o(v);
  x.x[3] = 1;
  f = v2o(x);
  hyoe[4] = g;
  g = omul(g, h);
  printpol(o2v(g));
  hyoe[5] = g;
  // exit(1);

  for (i = 6; i < 1331; i++)
  {
    g = omul(h, g);
    if (deg(o2v(g)) == 3)
    {
      g = oadd(g, kei(LT(g).a, m));
      g = omod(g, f);
      printpol(o2v(g));
      printf("====3\n");
      // exit(1);
      printpol(o2v(m));
      printf("\n");
    }
    //  exit(1);
    //  g=oadd(g,

    printpol(o2v(g));
    printf("\n");
    hyoe[i] = g;
    printpol(o2v(g));
    printf("\n");
  }
  // exit(1);
  for (i = 0; i < 2; i++)
    sue[i] = i;
  for (i = 2; i < 1331; i++)
    sue[i] = xtrace(hyoe[i], Pr);
  for (i = 0; i < 1331; i++)
    sie[sue[i]] = i;
}

int main()
{
  unsigned short i, count = 0, c2 = 0;

  // tas();
  // exit(1);

  mkmf();
  // exit(1);
  makefg(O);
  de();
  exit(1);

  count = 0;
  c2 = 0;
  for (i = 0; i < O; i++)
  {
    if (gf[i] >= 0)
      count++;
    if (fg[i] >= 0)
    {
      c2++;
    }
  }

  printf("%d %d\n", count, c2);
  // exit(1);

  return 0;
}