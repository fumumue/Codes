// date : 20210325（ver.1）
// date      : 20160310,20191218,20191220,20191221,20191223,20191224,20191225,20191229,20191230
// auther    : the queer who thinking about cryptographic future
// code name : Ben_Or's irreducibly test for polynomial over GF(2^m)
// code name : OVP - One Variable Polynomial library
// status    : release (ver 1.2) 20211201

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
// #include <assert.h>
#include <execinfo.h>
// #include <omp.h>
#include <unistd.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sys/types.h>
// #include "1024.h"
// #include "8192.h"
// #include "2048.h"
// #include "4096.h"
// #include "512.h"
#include "global-p.h"
#include "struct.h"
#include "debug.c"
#include "chash-p.c"

#include <pthread.h>
#include <err.h>
#include <errno.h>

extern int mltu();
extern int mlt();
// extern int mltn();
extern void print_trace();

#define DAT 5
int num = 0;

// nomal bases
// unsigned short gf[M]={0,1,2,4,8,9,11,15,7,14,5,10,13,3,6,12};
// unsigned short fg[M]={0,1,2,13,3,10,14,8,4,5,11,6,15,12,9,7};

// sage比較用
// static const unsigned short gf[16]={0,1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};
// static const unsigned short fg[16]={0,1,2,5,3,9,6,11,4,15,10,8,7,14,12,13};

// GF(2^8)
// static const unsigned short gf[N] = {0, 1, 2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76, 152, 45, 90, 180, 117, 234, 201, 143, 3, 6, 12, 24, 48, 96, 192, 157, 39, 78, 156, 37, 74, 148, 53, 106, 212, 181, 119, 238, 193, 159, 35, 70, 140, 5, 10, 20, 40, 80, 160, 93, 186, 105, 210, 185, 111, 222, 161, 95, 190, 97, 194, 153, 47, 94, 188, 101, 202, 137, 15, 30, 60, 120, 240, 253, 231, 211, 187, 107, 214, 177, 127, 254, 225, 223, 163, 91, 182, 113, 226, 217, 175, 67, 134, 17, 34, 68, 136, 13, 26, 52, 104, 208, 189, 103, 206, 129, 31, 62, 124, 248, 237, 199, 147, 59, 118, 236, 197, 151, 51, 102, 204, 133, 23, 46, 92, 184, 109, 218, 169, 79, 158, 33, 66, 132, 21, 42, 84, 168, 77, 154, 41, 82, 164, 85, 170, 73, 146, 57, 114, 228, 213, 183, 115, 230, 209, 191, 99, 198, 145, 63, 126, 252, 229, 215, 179, 123, 246, 241, 255, 227, 219, 171, 75, 150, 49, 98, 196, 149, 55, 110, 220, 165, 87, 174, 65, 130, 25, 50, 100, 200, 141, 7, 14, 28, 56, 112, 224, 221, 167, 83, 166, 81, 162, 89, 178, 121, 242, 249, 239, 195, 155, 43, 86, 172, 69, 138, 9, 18, 36, 72, 144, 61, 122, 244, 245, 247, 243, 251, 235, 203, 139, 11, 22, 44, 88, 176, 125, 250, 233, 207, 131, 27, 54, 108, 216, 173, 71, 142};
// static const unsigned short fg[N] = {0, 1, 2, 26, 3, 51, 27, 199, 4, 224, 52, 239, 28, 105, 200, 76, 5, 101, 225, 15, 53, 142, 240, 130, 29, 194, 106, 249, 201, 9, 77, 114, 6, 139, 102, 48, 226, 37, 16, 34, 54, 148, 143, 219, 241, 19, 131, 70, 30, 182, 195, 126, 107, 40, 250, 186, 202, 155, 10, 121, 78, 229, 115, 167, 7, 192, 140, 99, 103, 222, 49, 254, 227, 153, 38, 180, 17, 146, 35, 137, 55, 209, 149, 207, 144, 151, 220, 190, 242, 211, 20, 93, 132, 57, 71, 65, 31, 67, 183, 164, 196, 73, 127, 111, 108, 59, 41, 85, 251, 134, 187, 62, 203, 95, 156, 160, 11, 22, 122, 44, 79, 213, 230, 173, 116, 244, 168, 88, 8, 113, 193, 248, 141, 129, 100, 14, 104, 75, 223, 238, 50, 198, 255, 25, 228, 166, 154, 120, 39, 185, 181, 125, 18, 69, 147, 218, 36, 33, 138, 47, 56, 64, 210, 92, 150, 189, 208, 206, 145, 136, 152, 179, 221, 253, 191, 98, 243, 87, 212, 172, 21, 43, 94, 159, 133, 61, 58, 84, 72, 110, 66, 163, 32, 46, 68, 217, 184, 124, 165, 119, 197, 24, 74, 237, 128, 13, 112, 247, 109, 162, 60, 83, 42, 158, 86, 171, 252, 97, 135, 178, 188, 205, 63, 91, 204, 90, 96, 177, 157, 170, 161, 82, 12, 246, 23, 236, 123, 118, 45, 216, 80, 175, 214, 234, 231, 232, 174, 233, 117, 215, 245, 235, 169, 81, 89, 176};

// static unsigned short c[E * K + 1] = {0};

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
    tmp = plus(tmp, (t0 * (t1 & (1 << i))));

  t = tmp & 0x7FC000;
  tmp ^= t >> 9;
  tmp ^= t >> 12;

  t = tmp & 0x3000;
  tmp ^= t >> 9;
  tmp ^= t >> 12;

  return tmp & ((1 << 12) - 1);
}

// 整数のべき乗
unsigned int
ipow(unsigned int q, unsigned int u)
{
  unsigned int i, m = 1;

  for (i = 0; i < u; i++)
    m = gf[mlt(fg[m], fg[q])];

  printf("in ipow====%d\n", m);

  return m;
}

// a<b : gf[b]%gf[a]
unsigned short gf_mod(unsigned short a, unsigned short b)
{
  int i = fg[b] % fg[a];
  int k = 0;

  printf("a=%d b=%d\n", fg[a], fg[b]);
  if (i > 0)
  {
    return i + 1;
  }
  if (i == 0)
    return 1;
  if (fg[b] < fg[a])
    return fg[b];
}

unsigned short gcd(unsigned short a, unsigned short b)
{
  int r, tmp;

  /* 自然数 a > b を確認・入替 */
  if (fg[a] < fg[b])
  {
    tmp = a;
    a = b;
    b = tmp;
  }

  /* ユークリッドの互除法 */
  r = gf[gf_mod(a, b)];
  while (r != 0)
  {
    a = b;
    b = r;
    r = gf[gf_mod(a, b)];
  }

  /* 最大公約数を出力 */
  printf("最大公約数 = %d\n", b);

  return b;
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
    if (K == 64)
    {
      prod[i - K + 33] ^= prod[i];
      prod[i - K + 30] ^= prod[i];
      prod[i - K + 26] ^= prod[i];
      prod[i - K + 25] ^= prod[i];
      prod[i - K + 24] ^= prod[i];
      prod[i - K + 23] ^= prod[i];
      prod[i - K + 22] ^= prod[i];
      prod[i - K + 21] ^= prod[i];
      prod[i - K + 20] ^= prod[i];
      prod[i - K + 18] ^= prod[i];
      prod[i - K + 13] ^= prod[i];
      prod[i - K + 12] ^= prod[i];
      prod[i - K + 11] ^= prod[i];
      prod[i - K + 10] ^= prod[i];
      prod[i - K + 7] ^= prod[i];
      prod[i - K + 5] ^= prod[i];
      prod[i - K + 4] ^= prod[i];
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
    if (K == 4)
    {
      // 16
      prod[i - K + 3] ^= prod[i];
      prod[i - K + 0] ^= prod[i];
    }
  }

  for (i = 0; i < K; i++)
    out[i] = prod[i];
}

// 有限体の元の逆数
unsigned short
oinv(unsigned short a)
{
  int i;

  if (a == 0)
    return 0;

  return N - fg[a] + 1;

  printf("no return \n");

  exit(1);
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
      // exit(1);
      for (j = 0; j < K + 1; j++)
        printf("%d,", a.x[i][j]);
    printf("\n");
  }
  printf("\n");

  for (i = 0; i < K; i++)
  {
    v.x[i] = a.x[i][K];
    // v.x[128]=1;
    printf(" x%d = %d\n", i, v.x[i]);
  }

  return v;
}

/* input: f, element in GF((2^m)^t) */
/* output: out, minimal polynomial of f */
/* return: 0 for success and -1 for failure */
int mykey(unsigned short *out, vec x)
{
  unsigned short matt[K + 1][K] = {0};
  MTX a = {0};
  int i, j, k;

  // fill matrix

  matt[0][0] = 1;

  for (i = 1; i < K; i++)
    matt[0][i] = 0;

  for (i = 0; i < K; i++)
    matt[1][i] = x.x[i];

  for (j = 2; j <= K; j++)
  {
    // for(i=0;i<128;i++)
    // mat[j][i]=gf[mlt(fg[mat[j-1][i]],fg[x.x[i]])];
    GF_mul(matt[j], matt[j - 1], x.x);

    // for(i=0;i<K;i++)
    // printf("%d,",mat[j][i]);
    // printf("\n");
  }
  // exit(1);
  //
  for (i = 0; i < K; i++)
  {
    for (j = 0; j < K + 1; j++)
    {
      a.x[i][j] = matt[j][i];
      printf("%d,", matt[j][i]);
    }
    printf("\n");
  }
  printf("\n");
  // exit(1);

  vec v = {0};
  v = renritu(a);
  // printsage(x);
  // printf("\n");
  // for(i=0;i<T;i++)
  // v.x[T-i-1]=x.x[i];
  // printsage(v);
  // printf("\n");

  for (i = 0; i < K; i++)
  {
    out[i] = v.x[i];
    printf("%d,", out[i]);
  }
  printf("\n");
  // exit(1);
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

static uint8_t xtime(uint8_t x)
{
  return ((x << 1) ^ (((x >> 7) & 1) * 0x1b));
}

// Multiply is used to multiply numbers in the field GF(2^8)
// Note: The last call to xtime() is unneeded, but often ends up generating a smaller binary
//       The compiler seems to be able to vectorize the operation better this way.
//       See https://github.com/kokke/tiny-AES-c/pull/34
static uint8_t Multiply(uint8_t x, uint8_t y)
{
  return (((y & 1) * x) ^
          ((y >> 1 & 1) * xtime(x)) ^
          ((y >> 2 & 1) * xtime(xtime(x))) ^
          ((y >> 3 & 1) * xtime(xtime(xtime(x)))) ^
          ((y >> 4 & 1) * xtime(xtime(xtime(xtime(x)))))); /* this last call to xtime() can be omitted */
}

// aに何をかけたらbになるか
unsigned short
equ(unsigned short a, unsigned short b)
{
  int i;

  return gf[mlt(oinv(a), fg[b])];
}

// 停止コマンド
void wait2(void)
{
  printf(" (enter number and hit return) "); // 何か表示させたほうが良いだろう
  fflush(stdout);                            // just in case
  getchar();
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

int primitive(vec x)
{
  int i, j, k = deg(x), l;

  for (i = 0; i < k; i++)
  {
    for (j = i; j < k; j++)
    {
      if (i != j)
      {
        l = gcd(x.x[i], x.x[j]);
        if (l > 1)
          return -1;
      }
    }
  }
  return 0;
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

vec vadd(vec a, vec b)
{
  int i;
  vec c = {0};

  // printf("deg=%d %d\n",deg(a),deg(b));

  for (i = 0; i < DEG; i++)
    c.x[i] = a.x[i] ^ b.x[i];

  return c;
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

// 多項式を表示する(default)
void printpol(vec a)
{
  int i, n;

  n = deg(a);

  for (i = n; i > -1; i--)
  {
    if (a.x[i] > 0)
    {
      printf("%u", a.x[i]);
      // if (i > 0)
      printf("x^%d", i);
      if (i > 0)
        printf("+");
    }
  }
  // printf("\n");

  return;
}

int mul = 0, mul2 = 0;
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

typedef struct
{
  vec a;
  vec b;
  vec c;
} tri;

tri wake(vec x, int D)
{
  int i;
  tri v = {0};
  for (i = 0; i < D; i++)
  {
    v.a.x[i] ^= x.x[i];
    v.b.x[i] ^= x.x[i + D];
    v.c.x[i] ^= x.x[i] ^ x.x[i + D];
  }

  return v;
}

tri kake(tri x, tri y)
{
  tri t = {0};

  t.a = vmul_2(x.a, y.a);
  t.b = vmul_2(x.b, y.b);
  t.c = vadd(vadd(vmul_2(x.c, y.c), t.a), t.b);
  // printf("%d^\n",deg(t.a));
  // exit(1);

  return t;
}

vec tas(tri n, int D)
{
  int i;
  vec e = {0};

  for (i = 0; i < D; i++)
  {
    e.x[i] ^= n.a.x[i];
    e.x[i + D] ^= n.b.x[i];
    e.x[i + D / 2] ^= n.c.x[i];
  }
  // printpol(e);
  // printf("==eeeee\n");
  // exit(1);

  return e;
}

// カラツバ法：3/4だけ早くなるｗ
vec karatuba(vec f, vec g)
{
  int i;
  vec f1 = {0}, f2 = {0}, g1 = {0}, g2 = {0}, f3 = {0}, g3 = {0}, a = {0}, b = {0}, c = {0}, d = {0}, e;
  vec o1, o2, o3, o4, o5, o6, x4, y1, y2, y3, y4, y5, y6, y7 = {0}, z1, z2, z3, z4, z5, z6, z7 = {0}, w1, w2, w3, w4, w5, w6, w7, v1, v2, v3, v4, v5, v6, v7, x1, x2, x3, x5, x6, x7 = {0};
  if (deg(f) > 255 || deg(g) > 255)
    exit(1);

  int j = deg(f) + deg(g);

  for (i = 0; i < 128; i++)
  {
    f1.x[i] = f.x[i];
    g1.x[i] = g.x[i];
  }
  for (i = 128; i < 256; i++)
  {
    f2.x[i - 128] = f.x[i];
    g2.x[i - 128] = g.x[i];
  }
  for (i = 0; i < 128; i++)
  {
    f3.x[i] = f.x[i] ^ f.x[i + 128];
    g3.x[i] = g.x[i] ^ g.x[i + 128];
  }

  for (i = 0; i < 64; i++)
  {
    x1.x[i] = f1.x[i];
    x2.x[i] = f1.x[i + 64];
    x3.x[i] = f1.x[i] ^ f1.x[i + 64];
    y1.x[i] = f2.x[i];
    y2.x[i] = f2.x[i + 64];
    y3.x[i] = f2.x[i] ^ f2.x[i + 64];
    v1.x[i] = f3.x[i];
    v2.x[i] = f3.x[i + 64];
    v3.x[i] = f3.x[i] ^ f3.x[i + 64];

    z1.x[i] = g1.x[i];
    z2.x[i] = g1.x[i + 64];
    z3.x[i] = g1.x[i] ^ g1.x[i + 64];
    w1.x[i] = g2.x[i];
    w2.x[i] = g2.x[i + 64];
    w3.x[i] = g2.x[i] ^ g2.x[i + 64];
    o1.x[i] = g3.x[i];
    o2.x[i] = g3.x[i + 64];
    o3.x[i] = g3.x[i] ^ g3.x[i + 64];
  }
  x4 = vmul_2(x1, z1);
  y4 = vmul_2(x2, z2);
  z4 = vadd(vadd(vmul_2(x3, z3), x4), y4);
  x5 = vmul_2(y1, w1);
  y5 = vmul_2(y2, w2);
  z5 = vadd(vadd(vmul_2(y3, w3), x5), y5);
  x6 = vmul_2(v1, o1);
  y6 = vmul_2(v2, o2);
  z6 = vadd(vadd(vmul_2(v3, o3), x6), y6);

  for (i = 0; i < 128; i++)
  {
    x7.x[i] ^= x4.x[i];
    x7.x[i + 128] ^= y4.x[i];
    x7.x[i + 64] ^= z4.x[i];
    y7.x[i] ^= x5.x[i];
    y7.x[i + 128] ^= y5.x[i];
    y7.x[i + 64] ^= z5.x[i];
    z7.x[i] ^= x6.x[i];
    z7.x[i + 128] ^= y6.x[i];
    z7.x[i + 64] ^= z6.x[i];
  }

  for (i = 0; i < 256; i++)
  {
    a.x[i] ^= x7.x[i];
    a.x[i + 256] ^= y7.x[i];
    a.x[i + 128] ^= z7.x[i] ^ x7.x[i] ^ y7.x[i];
  }

  return a;
}

// カラツバ法：3/4だけ早くなるｗ
vec karatubaka(vec f, vec g)
{
  int i, k;
  vec f1 = {0}, f2 = {0}, g1 = {0}, g2 = {0}, f3 = {0}, g3 = {0}, a = {0}, b = {0}, c = {0}, d = {0}, e;
  vec o1 = {0}, o2 = {0}, o3 = {0}, o4 = {0}, o5 = {0}, o6 = {0}, y1 = {0}, y2 = {0}, y3 = {0}, y4, y5, y6, y7 = {0}, z1 = {0}, z2 = {0}, z3 = {0}, z4, z5, z6, z7 = {0}, w1 = {0}, w2 = {0}, w3 = {0}, w4, w5, w6, w7, v1 = {0}, v2 = {0}, v3 = {0}, v4, v5, v6, v7, x1 = {0}, x2 = {0}, x3 = {0}, x5, x6, x7 = {0};
  if (deg(f) > 255 || deg(g) > 255)
  {
    printf("too big! %d %d\n", deg(f), deg(g));
    printpol(f);
    printf("\n");
    printpol(g);
    printf("\n");
    // exit(1);
  }
  int j = deg(f) + deg(g);

  for (i = 0; i < 128; i++)
  {
    f1.x[i] = f.x[i];
    g1.x[i] = g.x[i];
  }
  for (i = 128; i < 256; i++)
  {
    f2.x[i - 128] = f.x[i];
    g2.x[i - 128] = g.x[i];
  }
  for (i = 0; i < 128; i++)
  {
    f3.x[i] = f.x[i] ^ f.x[i + 128];
    g3.x[i] = g.x[i] ^ g.x[i + 128];
  }
  vec xx[9], yy[9], vv[9], zz[9], ww[9], oo[9], pp[9], qq[9], rr[9], ss[9], xa[9], xb[9], xc[9], ya[9], yb[9], yc[9], za[9], zb[9], zc[9];

  for (i = 0; i < 64; i++)
  {
    x1.x[i] = f1.x[i];
    x2.x[i] = f1.x[i + 64];
    x3.x[i] = f1.x[i] ^ f1.x[i + 64];
    y1.x[i] = f2.x[i];
    y2.x[i] = f2.x[i + 64];
    y3.x[i] = f2.x[i] ^ f2.x[i + 64];
    v1.x[i] = f3.x[i];
    v2.x[i] = f3.x[i + 64];
    v3.x[i] = f3.x[i] ^ f3.x[i + 64];

    z1.x[i] = g1.x[i];
    z2.x[i] = g1.x[i + 64];
    z3.x[i] = g1.x[i] ^ g1.x[i + 64];
    w1.x[i] = g2.x[i];
    w2.x[i] = g2.x[i + 64];
    w3.x[i] = g2.x[i] ^ g2.x[i + 64];
    o1.x[i] = g3.x[i];
    o2.x[i] = g3.x[i + 64];
    o3.x[i] = g3.x[i] ^ g3.x[i + 64];
  }

  for (k = 0; k < 2; k++)
  {
    for (i = 0; i < 32; i++)
    {
      xx[k].x[i] = x1.x[i + 32 * k];
      xx[2].x[i] = x1.x[i] ^ x1.x[i + 32];
      yy[k].x[i] = x2.x[i + 32 * k];
      yy[2].x[i] = x2.x[i] ^ x2.x[i + 32];
      vv[k].x[i] = x3.x[i + 32 * k];
      vv[2].x[i] = x3.x[i] ^ x3.x[i + 32];
      zz[k].x[i] = y1.x[i + 32 * k];
      zz[2].x[i] = y1.x[i] ^ y1.x[32 + i];
      ww[k].x[i] = y2.x[i + 32 * k];
      ww[2].x[i] = y2.x[i] ^ y2.x[32 + i];
      oo[k].x[i] = y3.x[i + 32 * k];
      oo[2].x[i] = y3.x[i] ^ y3.x[32 + i];
      pp[k].x[i] = v1.x[i + 32 * k];
      pp[2].x[i] = v1.x[i] ^ v1.x[32 + i];
      qq[k].x[i] = v2.x[i + 32 * k];
      qq[2].x[i] = v2.x[i] ^ v2.x[32 + i];
      rr[k].x[i] = v3.x[i + 32 * k];
      rr[2].x[i] = v3.x[i] ^ v3.x[32 + i];

      xa[k].x[i] = z1.x[i + 32 * k];
      xa[2].x[i] = z1.x[i] ^ z1.x[i + 32];
      xb[k].x[i] = z2.x[i + 32 * k];
      xb[2].x[i] = z2.x[i] ^ z2.x[i + 32];
      xc[k].x[i] = z3.x[i + 32 * k];
      xc[2].x[i] = z3.x[i] ^ z3.x[i + 32];
      ya[k].x[i] = w1.x[i + 32 * k];
      ya[2].x[i] = w1.x[i] ^ w1.x[i + 32];
      yb[k].x[i] = w2.x[i + 32 * k];
      yb[2].x[i] = w2.x[i] ^ w2.x[i + 32];
      yc[k].x[i] = w3.x[i + 32 * k];
      yc[2].x[i] = w3.x[i] ^ w3.x[i + 32];
      za[k].x[i] = o1.x[i + 32 * k];
      za[2].x[i] = o1.x[i] ^ o1.x[i + 32];
      zb[k].x[i] = o2.x[i + 32 * k];
      zb[2].x[i] = o2.x[i] ^ o2.x[i + 32];
      zc[k].x[i] = o3.x[i + 32 * k];
      zc[2].x[i] = o3.x[i] ^ o3.x[i + 32];
    }
  }

  vec aa[3], bb[3], cc[3], dd[3], ee[3], ff[3], gg[3], hh[3], nn[3], mm[3];
  for (k = 0; k < 2; k++)
  {
    aa[k] = vmul_2(xx[k], xa[k]);
    bb[k] = vmul_2(yy[k], xb[k]);
    cc[k] = vmul_2(vv[k], xc[k]);
    dd[k] = vmul_2(zz[k], ya[k]);
    ee[k] = vmul_2(ww[k], yb[k]);
    ff[k] = vmul_2(oo[k], yc[k]);
    // gg[k]=vmul_2(pp[k],za[k]);
    // hh[k]=vmul_2(qq[k],zb[k]);
    // nn[k]=vmul_2(rr[k],zc[k]);
  }
  aa[2] = vadd(vadd(vmul_2(xx[2], xa[2]), aa[0]), aa[1]);
  bb[2] = vadd(vadd(vmul_2(yy[2], xb[2]), bb[0]), bb[1]);
  cc[2] = vadd(vadd(vmul_2(vv[2], xc[2]), cc[0]), cc[1]);
  dd[2] = vadd(vadd(vmul_2(zz[2], ya[2]), dd[0]), dd[1]);
  ee[2] = vadd(vadd(vmul_2(ww[2], yb[2]), ee[0]), ee[1]);
  ff[2] = vadd(vadd(vmul_2(oo[2], yc[2]), ff[0]), ff[1]);
  // gg[2]=vadd(vadd(vmul_2(pp[2],za[2]),gg[0]),gg[1]);
  // hh[2]=vadd(vadd(vmul_2(qq[2],zb[2]),hh[0]),hh[1]);
  // nn[2]=vadd(vadd(vmul_2(rr[2],zc[2]),nn[0]),nn[1]);

  vec oo1 = {0}, oo2 = {0}, oo3 = {0}, oo4 = {0}, oo5 = {0}, oo6 = {0}, oo7 = {0}, oo8 = {0}, oo9 = {0};
  // exit(1);
  for (i = 0; i < 64; i++)
  {
    x7.x[i] = aa[0].x[i];
    x7.x[i + 64] = aa[1].x[i];
    x7.x[i + 32] ^= aa[2].x[i];

    y7.x[i] = bb[0].x[i];
    y7.x[i + 64] = bb[1].x[i];
    y7.x[i + 32] ^= bb[2].x[i];

    z7.x[i] = cc[0].x[i];
    z7.x[i + 64] = cc[1].x[i];
    z7.x[i + 32] ^= cc[2].x[i];

    oo1.x[i] = dd[0].x[i];
    oo1.x[i + 64] = dd[1].x[i];
    oo1.x[i + 32] ^= dd[2].x[i];

    oo2.x[i] = ee[0].x[i];
    oo2.x[i + 64] = ee[1].x[i];
    oo2.x[i + 32] ^= ee[2].x[i];

    oo3.x[i] = ff[0].x[i];
    oo3.x[i + 64] = ff[1].x[i];
    oo3.x[i + 32] ^= ff[2].x[i];

    oo4.x[i] = gg[0].x[i];
    oo4.x[i + 64] = gg[1].x[i];
    // oo4.x[i+32]^=gg[2].x[i];
    oo5.x[i] = hh[0].x[i];
    oo5.x[i + 64] = hh[1].x[i];
    // oo5.x[i+32]^=hh[2].x[i];
    oo6.x[i] = nn[0].x[i];
    oo6.x[i + 64] = nn[1].x[i];
    // oo6.x[i+32]^=nn[2].x[i];
  }
  vec oe[3] = {0}, ooo = {0};
  for (i = 0; i < 128; i++)
  {
    oe[0].x[i] = x7.x[i];
    oe[0].x[i + 128] = y7.x[i];
    // oe[0].x[i+64]^=z7.x[i];
    oe[1].x[i] = oo1.x[i];
    oe[1].x[i + 128] = oo2.x[i];
    // oe[1].x[i+64]^=oo3.x[i];
    oe[2].x[i] = oo4.x[i];
    oe[2].x[i + 128] = oo5.x[i];
    // oe[2].x[i+64]^=oo6.x[i];
  }
  for (i = 0; i < 256; i++)
  {
    ooo.x[i] = oe[0].x[i];
    ooo.x[i + 256] = oe[1].x[i];
    // ooo.x[i+128]^=oe[2].x[i];
  }
  // printpol(oe[2]);
  // exit(1);

  return ooo;
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

vec vcoef(vec v)
{
  unsigned short n, k = deg(v);
  int i;

  // if(v.x[0]==0)
  // return v;

  if (v.x[0] > 1)
    n = oinv(v.x[0]);
  for (i = 0; i < k + 1; i++)
    v.x[i] = gf[mlt(n, fg[v.x[i]])];

  return v;
}

vec rev(vec a, int n)
{
  vec vx = (a), tmp = {0};
  int k = deg((a)), i;
  OP t = {0};

  // printf("k=%d\n", k);
  //  exit(1);
  for (i = 0; i < k + 1; i++)
    tmp.x[n - i] = a.x[i];
  // printpol(tmp);
  // printf("\n");
  //  exit(1);

  return tmp;
}

vec deli(vec a, vec b)
{
  int i = 0;
  OP t;
  vec v = {0};

  // printpol(a);
  // printf(" ==a\n");
  // printpol(b);
  /*
  if(deg(b)>deg(a)){
  printf(" ==b\n");
  printf("deg(b)=%d %d\n",deg(b),deg(a));
  return a;
  //exit(1);
  }
  */
  for (i = 0; i < deg(b); i++)
    v.x[i] = a.x[i];

  return v;
}

// モニック多項式にする
vec coeff(vec f)
{
  unsigned short d = vLT(f).a;
  int i, k;

  k = deg((f)) + 1;

  for (i = 0; i < k; i++)
    f.x[i] = gf[mlt(fg[f.x[i]], oinv(d))];

  return f;
}

// fg=1
vec kof(vec a, vec b)
{
  unsigned short u = oinv(b.x[0]), k = deg(a);
  int i;

  for (i = 0; i < k + 1; i++)
    a.x[i] = gf[mlt(u, fg[a.x[i]])];

  return a;
}

vec vinv(vec a)
{
  vec v = {0}, x = {0};
  int i;

  x.x[2] = 1;
  v.x[0] = 1;

  if (a.x[0] > 1)
    a = vcoef(a);

  i = 1;
  while (i < K + 1)
  {
    v = vmul_2(vmul_2(v, v), a);
    if (i > 2)
      x = vmul_2(x, x);
    v = deli(v, x);
    i *= 2;
  }

  return v;
}

int chkinv(vec b, vec e, vec d)
{
  e = deli(vmul_2(b, e), d);
  if (vLT(e).a != 1)
  {
    printpol(e);
    printf(" fail-1!\n");
    exit(1);
  }
  return 0;
}

vec kara(vec a, vec b)
{
  tri z, w, c, za, zb, zc, a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3;
  tri za1, za2, za3, zb1, zb2, zb3, zc1, zc2, zc3;
  tri aa1, aa2, aa3, bb1, bb2, bb3, cc1, cc2, cc3;
  tri x1, x2, x3, v1, v2, v3, w1, w2, w3, s1, s2, s3;
  tri x4, x5, x6;
  vec d;
  int i;

  z = wake(a, 128);

  za = wake(z.a, 64);
  zb = wake(z.b, 64);
  zc = wake(z.c, 64);

  w = wake(b, 128);

  aa1 = wake(w.a, 64);
  aa2 = wake(w.b, 64);
  aa3 = wake(w.c, 64);

  x1 = kake(aa1, za);
  x2 = kake(aa2, zb);
  x3 = kake(aa3, zc);

  c.a = tas(x1, 128);
  c.b = tas(x2, 128);
  c.c = vadd(vadd(tas(x3, 128), c.a), c.b);

  // c=kake(z,w);
  d = tas(c, 256);

  /*
    a1=wake(za.a,32);
    a2=wake(za.b,32);
    a3=wake(za.c,32);
    b1=wake(zb.a,32);
    b2=wake(zb.b,32);
    b3=wake(zb.c,32);
    c1=wake(zc.a,32);
    c2=wake(zc.b,32);
    c3=wake(zc.c,32);


    za1=wake(aa1.a,32);
    za2=wake(aa1.b,32);
    za3=wake(aa1.c,32);
    zb1=wake(aa2.a,32);
    zb2=wake(aa2.b,32);
    zb3=wake(aa2.c,32);
    zc1=wake(aa3.a,32);
    zc2=wake(aa3.b,32);
    zc3=wake(aa3.c,32);

    x1=kake(a1,za1);
    x2=kake(a2,za2);
    x3=kake(a3,za3);
    v1=kake(b1,zb1);
    v2=kake(b2,zb2);
    v3=kake(b3,zb3);
    w1=kake(c1,zc1);
    w2=kake(c2,zc2);
    w3=kake(c3,zc3);
    tri o1,o2,o3,p1,p2,p3,q1,q2,q3;
    vec ans={0},cf={0};

    o1.a=tas(x1,64);
    o1.b=tas(x2,64);
    o1.c=vadd(vadd(tas(x3,64),o1.a),o1.b);
    o2.a=tas(v1,64);
    o2.b=tas(v2,64);
    o2.c=vadd(vadd(tas(v3,64),o2.a),o2.b);
    o3.a=tas(w1,64);
    o3.b=tas(w2,64);
    o3.c=vadd(vadd(tas(w3,64),o3.a),o3.b);

    p1.a=tas(o1,128);
    p1.b=tas(o2,128);
    p1.c=vadd(vadd(tas(o3,128),p1.a),p1.b);

    d=tas(p1,256);
    */
  /*
  cf=vmul_2(a,b);
  for(i=0;i<512;i++){
  if(ans.x[i]!=cf.x[i]){
    printf("i=%d %d %d\n",i,ans.x[i],cf.x[i]);
  }
  }
  exit(1);
  */

  return d;
}

static inline unsigned int gf_pow(unsigned int n, unsigned int u)
{

  if (n % N == 0)
    return 1;
  return (u * n - n) % (N - 1) + 1;
}

int atom(unsigned short a)
{
  if (gf_pow(N - 1, a) == 1)
  {
    return 0;
  }
  else
  {
    return -1;
  }
}

vec vpowmod(vec f, vec mod)
{
  vec v = {0};
  vec ret = {0}, s = {0}, t = {0};

  int n = N;
  v.x[0] = 1;
  ret = v;
  while (n > 0)
  {
    // s=inv()
    if (n & 1)
    {
      ret = (karatuba(ret, f)); // n の最下位bitが 1 ならば x^(2^i) をかける
      if (deg(ret) > deg(mod))
        ret = karatuba(ret, mod);
    }
    f = (vmul_2(f, f)); // deg(f)>=256
    if (deg(f) > deg(mod))
      f = vmod(f, mod);
    n >>= 1; // n を1bit 左にずらす
  }
  // printpol((ret));
  // printf(" it be\n");
  // exit(1);

  return ret;
}

int fequ(vec a, vec b)
{
  int k = deg(a), l = deg(b), i;
  if (k != l)
    return 1;
  for (i = 0; i < k; i++)
    if (a.x[i] != b.x[i])
      return 1;

  return 0;
}

vec jorju(vec ww, vec xx)
{
  // unsigned short f[K + 1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0}; //big indian
  OP g, w;
  int i, count = 0;
  vec e[10] = {0}, v = {0}, x = {0}, z = {0}, ee = {0}, y = {0}, tt = {0};
  int l = -1, m, n;
  int ii = 0;
  // irreducible goppa code (既役多項式が必要なら、ここのコメントを外すこと。)
  vec q = {0}, r = {0}, c = {0}, dd = {0};

  if (deg(ww) < deg(xx))
    return ww;

  vec ff = ww, gg = xx, rr = {0}, vv = {0};

  // vv.x[0]=9;
  m = deg(ww);
  n = deg(xx);

  ww = rev(ww, m);
  xx = rev(xx, n);

  int o = 0;
  if (xx.x[0] > 1)
  {
    o = oinv(xx.x[0]);
    o = gf[o];
  }
  x = vinv(xx);
  v = vmul_2(ww, x);
  //  v=rev(v,m-n);
  //  y=vadd(ww,vmul_2(v,xx));

  q = rev(v, m - n);
  if (o > 0)
  {
    z.x[0] = o;
    r = vmul_2(vmul_2(gg, q), z);
  }
  else
  {
    r = vmul_2(gg, q);
  }
  vec www = {0};
  y = vadd(ff, r);

  return y;
}

int cnty = 0;
vec vpp(vec f, vec mod)
{
  int i;
  vec t = {0}, s = {0};
  t = f;
  s = f;

  // 繰り返し２乗法
  for (i = 1; i < E + 1; i++)
  {
    s = vmod(vmul_2(s, s), mod);
  }

  return s;
}

unsigned short v2a(oterm a)
{
  int j;

  if (a.a == 0)
    return 0;

  // printf("aa=%d\n",a.a);
  for (j = 0; j < N; j++)
  {
    if (gf[j] == a.a && a.a > 0)
    {
      // printf("j==%d\n",j);
      return j - 1;
    }
  }

  printf("baka-v2a\n");
  exit(1);
}

void printsage(vec a)
{
  int i, j, k = deg(a);
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
      printf("B('a^%d')*X**%d", j, i); // for GF(2^m)
      if (i != k)
        printf("+");
    }
  }
}

void fprintsage(vec a, FILE *fp)
{
  int i, j;
  oterm b;

  //  fp=fopen("dat","w");
  fprintf(fp, "poly=");
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      b.a = a.x[i];
      b.n = i;
      j = v2a(b);
      // printf("%d,==ba\n",b.a);
      // printf ("X**%d+", i); //for GF2
      if (i == K)
      {
        fprintf(fp, "B('a^%d')*X**%d;", j, i); // for GF(2^m)
      }
      else
      {
        fprintf(fp, "B('a^%d')*X**%d+", j, i); // for GF(2^m)
      }
    }
  }
}

// gcd
vec vgcd(vec xx, vec yy)
{
  vec tt = {0}, tmp, h = {0}, ss = {0}, ee = {0};
  ee.x[K] = 1;

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

// ０多項式かどうかのチェック
unsigned short
chk(vec f)
{
  int i, flg = 0;
  vec x = {0};

  // x = o2v(f);
  for (i = 0; i < DEG; i++)
  {
    if (f.x[i] > 0)
    {
      flg = 1;
      return 1;
    }
  }
  if (flg == 0)
    return 0;

  exit(1);
}

OP init_pol(OP f)
{
  int i;

  for (i = 0; i < DEG; i++)
  {
    f.t[i].a = 0;
    f.t[i].n = 0;
  }

  return f;
}

// ランダム多項式の生成
static void
ginit(unsigned short *g)
{
  int j, count = 0, k = 0;
  unsigned short gg[K + 1] = {0};

  // printf("in ginit\n");

  g[K] = 1;          // xor128();
  g[0] = rand() % 2; // or N
  k = rand() % 2;    //(K - 1);
  if (k > 0)
  {
    while (count < k)
    {
      // printf("in whule\n");
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

  memcpy(g, gg, sizeof(K + 1));
}

OP ww[T] = {0};

// 配列の値を係数として多項式に設定する
vec setpol(unsigned short f[], int n)
{
  OP g;
  vec a, v = {0};
  int i;

  for (i = 0; i < n; i++)
  {
    v.x[n - 1 - i] = f[i];
  }

  //  g = v2o(a);

  return v;
}

vec mkpol()
{
  int i, j, k, flg, ii = 0;
  vec w = {0};
  static unsigned short g[K + 1] = {0};

  do
  {
    j = 0;
    k = 0;
    flg = 0;

    memset(g, 0, sizeof(g));
    // memset(ta, 0, sizeof(ta));
    memset(w.x, 0, sizeof(w));
    ginit(g);
    ii++;
    if (ii > 100) // K / 2)
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

  // printpol(o2v(w));
  // printf(" ==g\n");
  // exit(1);

  return w;
}

// GF(2^m) then set m in this function.
int ben_or(vec f)
{
  int i, n; //, pid;
  OP uu;
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

// 多項式の代入値
unsigned short
trace(vec f, unsigned short x)
{
  int i, d;
  unsigned short u = 0;

  d = deg(f);

  for (i = 0; i < d + 1; i++)
  {
    u ^= gf[mlt(fg[f.x[i]], gf_pow(f.x[i], fg[x]))];
  }

  return u;
}

unsigned short dd[DEG][DEG] = {0};

void get_irrpoly(void)
{
  int i, j, l, ii = 0;
  vec w = {0};
  FILE *fp;
  unsigned short ta[DEG] = {0};

  j = 0;
  fp = fopen("dat.sage", "w");

aa:

  // printf("\n");

  // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
  l = -1;
  ii = 0;

  while (l == -1)
  {
    w = mkpol();
    l = ben_or(w);
    printf("irr=%d\n", l);
    if (ii > 300)
    {
      printf("too many tryal\n");
      goto aa;
      // exit(1);
    }
    printf("ben=%d\n", ii);
    ii++;
    //
  }

  // 多項式の値が0でないことを確認
  for (i = 0; i < N; i++)
  {
    ta[i] = trace(w, i);
    if (ta[i] == 0)
    {
      printf("trace 0 @ %d\n", i);
      // fail = 1;
      exit(1);
    }
  }

  printpol(w);
  printf(" =irreducible\n");
  fprintsage(w, fp);
  fprintf(fp, " print(poly.is_irreducible());\n");
  j++;

  if (j > 10)
  {
    fclose(fp);
    exit(1);
  }
  goto aa;
}

// sagemath仕様の既多項式をファイルに書き出す（マルチプロセス）
int irr_poly_to_file()
{
  vec w, x, y, z;
  int b = -1, l = -1, k1 = 0, k2 = 0, k3 = 0, k4 = 0, a = -1, c = -1, i = 0;
  FILE *f1, *f2, *f3, *f4;
  int pid[5] = {0};

  for (i = 0; i < 3 && (pid[i] = fork()) > 0; i++)
    printf("%d\n", pid[i]);
  // exit(1);
  f1 = fopen("dat.sage", "w");
  f2 = fopen("dat0.sage", "w");
  f3 = fopen("dat1.sage", "w");
  f4 = fopen("dat2.sage", "w");
  // f5=fopen("dat3.sage","w");
  // this is child process
  fprintf(f1, "B=GF(2^%d,'a')\n", E);
  fprintf(f1, "F.<X>=B[]\n");
  fprintf(f2, "B=GF(2^%d,'a')\n", E);
  fprintf(f2, "F.<X>=B[]\n");
  fprintf(f3, "B=GF(2^%d,'a')\n", E);
  fprintf(f3, "F.<X>=B[]\n");
  fprintf(f4, "B=GF(2^%d,'a')\n", E);
  fprintf(f4, "F.<X>=B[]\n");

  i = 0;
  while (1)
  {

    l = -1;
    w = mkpol();
    b = -1;
    x = mkpol();
    a = -1;
    y = mkpol();
    c = -1;
    z = mkpol();

    if (pid[2] > 0)
    {
      l = ben_or(w);
      if (l == 0 && k1 < DAT)
      {
        fprintsage(w, f1);
        fprintf(f1, " print(poly.is_irreducible());\n");
        l = -1;
        k1++;
        printf("k1=%d\n", k1);
        if (k1 == DAT)
        {
          fclose(f1);
          wait(&pid[0]);
          wait(&pid[1]);
          wait(&pid[2]);
          exit(1);
        }
        // goto aa;
      }
    }
    else if (pid[1] > 0)
    {
      // }//this is mother process

      b = ben_or(x);
      if (b == 0 && k2 < DAT)
      {
        fprintsage(x, f2);
        fprintf(f2, " print(poly.is_irreducible());\n");
        b = -1;
        k2++;
        printf("k2=%d\n", k2);

        if (k2 == DAT)
        {
          fclose(f2);
          _exit(0);
          // wait(&pid[1]);
        }
        // goto aa;
      }
    }
    //}//this is mother process

    else if (pid[0] > 0)
    {

      // if(pid[1] == 0){
      a = ben_or(y);
      if (a == 0 && k3 < DAT)
      {
        fprintsage(y, f3);
        fprintf(f3, " print(poly.is_irreducible());\n");
        a = -1;
        k3++;
        printf("k3=%d\n", k3);

        if (k3 == DAT)
        {
          fclose(f3);
          _exit(0);
          // wait(&pid[0]);
        }
        // goto aa;
      }
    }
    else if (pid[0] == 0)
    {
      // this is mother process

      c = ben_or(z);
      if (c == 0 && k4 < DAT)
      {
        fprintsage(z, f4);
        fprintf(f4, " print(poly.is_irreducible());\n");
        c = -1;
        k4++;
        printf("k4=%d\n", k4);
        if (k4 == DAT)
        {
          fclose(f4);
          _exit(0);
        }
        // goto aa;
      }
    }
  }

  return 0;
}

void compare()
{
  int i, j;
  MTX a = {0};
  vec pp = {0}, x = {0}, v = {0};

  do
  {
    for (i = 0; i < K; i++)
      pp.x[i] = rand() % N;
    MTX opu = {0};
    vec cc = {0};

    opu.x[0][0] = 1;
    for (i = 1; i < K; i++)
      opu.x[0][i] = 0; // rand() % N;

    for (i = 0; i < K; i++)
      opu.x[1][i] = pp.x[i];

    for (i = 0; i < K; i++)
    {
      for (j = 2; j <= K; j++)
      {
        opu.x[j][i] = gf[mltn(j, fg[pp.x[i]])];
        // GF_mul2(opu.x[j], opu.x[j - 1], pp.x);
      }
    }

    for (i = 0; i < K; i++)
    {
      for (j = 0; j < K + 1; j++)
      {
        a.x[i][j] = opu.x[j][i];
        printf("%d,", opu.x[j][i]);
      }
      printf(" ===\n");
    }
    printf("\n");

    v = renritu(a);

    for (i = 0; i < K; i++)
      x.x[K - i + 1] = v.x[i];
    x.x[K] = 1;

    v.x[K] = 1;

    printsage(v);
    printf("\n");
  } while (ben_or(v) == -1);
  exit(1);
}

void speed()
{
  int i;
  vec f, g, q, r;

  for (i = 0; i < 100000; i++)
  {
    // vmul_2(f,g);
    karatuba(f, g);
    // karatubaka(f,f);
    // kara(f,g);
  }
  exit(1);

  // q=karatuba(f,g);
  q = karatuba(f, g);
  r = vmul_2(f, g);
  printpol(q);
  printf("===q\n\n");
  printpol(r);
  printf("\naaa\n");
  for (i = 0; i < deg(r); i++)
    if (q.x[i] != r.x[i])
      printf("i=%d, %d %d\n", i, q.x[i], r.x[i]);
  exit(1);
}

/*
// 言わずもがな
int main(void)
{
  vec g, w;
  unsigned int i, j, count = 0;
  vec e[10] = {0}, v = {0}, x = {0}, z = {0}, ee = {0}, y = {0}, tt = {0}, ww, xx, f;
  int l = -1, m, n;
  int ii = 0;
  // irreducible goppa code (既役多項式が必要なら、ここのコメントを外すこと。)
  vec q = {0}, r = {0};
  unsigned short ff[256] = {1, 1, 1, 1};
  unsigned short gg[256] = {0, 0, 1, 1};
  MTX a = {0};
  char rr[16] = {0};

  srand(clock());
  i = 4;
  j = 8;
  printf("%d\n", gf[gf_mod(gf[5], gf[8])]);
  // exit(1);

  vec pp = {0};

if(K==4){
  l = -1;
  while (l < 0)
  {
    for (i = 0; i < K; i++)
      pp.x[i] = rand() % N;
    mykey(tt.x, pp);
    tt.x[K] = 1;
    if (ben_or(tt) == 0)
    {
      printf("\n");
      printsage(tt);
      printf(" ==irr\n");
      exit(1);
    }
  }
}else{
  l= -1;
  while (l<0)
  {
    l= -1;
    v = mkpol();
    l = ben_or(v);
    printf("irr=%d\n", l);
    if (l == 0)
    {
      printsage(v);
      printf(" ==irr\n");
      exit(1);
      ii++;
    }
    if (ii > 300)
    {
      printf("too many error\n");
      exit(1);
    }
    // ii++;
    //
  }
  printf("aa\n");
}
//    exit(1);

  // multi_process();
  // irr_poly_to_file();
  // get_irrpoly();

  return 0;
}
*/