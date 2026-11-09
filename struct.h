
/* -*- mode: C; coding:utf-8 -*- */

//monomial
typedef struct
{
   unsigned short n; //単項式の次数
   unsigned short a; //単項式の係数
} oterm;

//polynomial
typedef struct
{
  oterm t[DEG]; //単項式の配列として多項式を表現する
} OP;

typedef struct 
{
   unsigned short x[DEG]; //配列の添字を次数に、配列の値を係数に持つ多項式の表現
} vec;

typedef struct {
  OP q;
  OP r;
} rem;

typedef struct
{
   short v[N];
  int f;
} MT;

//extra gcd
typedef struct
{
  OP u; //inverse of polynomial?
  OP v; //error locater
 OP d; //gcd
} EX;

typedef union
{ //test(SIMN)
   long long int u[K/4];
   short s[K];
} SU;

typedef union
{
  __uint128_t z[K/16];
   long long int u[K/8];
   int t[K/4];
   short x[K/2];
   char d[K];
} arrayul;

typedef union
{
   long long int u[N];
   char d[N * 8];
} arrayull;

typedef union a4
{
   char ar[4];
   int n;
} array;

typedef struct a8
{
   char ar[8];
} array8;

typedef union
{
   int h[16];
   char c[64];
} array16;

typedef union aN
{
   int d[64];
   long long int u[32];
   char ar[256];
  //
} arrayn;

typedef struct pub
{
   char a[K];
   char b[K];
} set;


#define I8T char
#define U8C(v) (v##U)

#define U8V(v) (( char)(v)&U8C(0xFF))
#define ROTL8(v, n) \
  (U8V((v) << (n)) | ((v) >> (8 - (n))))

#define R(x, n) (((x) << (n)) | ((x) >> (32 - (n))))

/*
 int rotate_left( int x, int n)
{
  assert(0 < n && n < 32);
  return (x << n) | (x >> (32 - n));
}
*/
/*
typedef union
{
   long long int u[8];
   char d[64];
} arrayul;

typedef struct a4
{
   char ar[4];
} array;

typedef struct a8
{
   char ar[8];
} array8;

typedef struct
{
   int h[16];
} array16;

typedef struct aN
{
   char ar[N];
} arrayn;

typedef struct pub
{
   char a[N];
   char b[N];
} set;
*/

typedef struct {
   short x[N][N];
  OP f;
  int row; //行
  int col; //列
} MTX;

/*
typedef struct {
   short x[N][K];
   char z[N][K*E];
   char w[K*E][K*E];
  int i;
  int rank;
} MAT;
*/