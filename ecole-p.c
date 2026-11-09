#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "global-p.h"
#include "struct.h"
//#include "chash-p.c"
//#include "debug.c"
//#include "1331.h"
//#define K 5
#define K2 8
#define P 13
#define O 2197 //1331

//sagemath上での原始多項式
unsigned short pp[6][4]= {{0,0,9,2}, {0,0,11,2}, {0,0,16,3}, {0,0,15,2},{0,0,1,2},{0,1,0,3}};
 // {0,0,9,2}, {1,0,11,2}, {1,0,16,3}, {1,0,15,2};
 //GF(11^3,13^3,17^3,19^3)
//unsigned short ff[2][7]={{1,0,0,0,0,2,0,2},{0,0,1,0,0,0,1,2}}; //GF(3^7,5^5)

unsigned short gf[O]={0},fg[O]={0};
//int N =0,M=0;
unsigned short c[K2+1]={0};

//OP型からベクトル型への変換
vec
o2v (OP f)
{
  vec a = { 0 };
  int i;

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
    {
      if (f.t[i].a > 0 && f.t[i].n < DEG)
	      a.x[f.t[i].n] = f.t[i].a;
    }


  return a;
}


//ベクトル型からOP型への変換
OP
v2o (vec a)
{
  int i, j = 0;
  OP f = { 0 };

  //#pragma omp parallel for
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

void
op_print_raw (const OP f)
{
  puts ("op_print_raw:");
  for (int i = 0; i < DEG; i++)
    {
      if (f.t[i].a > 0)
	printf ("[%d] %ux^%u\n", i, f.t[i].a, f.t[i].n);
    }
}

//多項式の次数(default)
int
deg (vec a)
{
  int i, n = 0;

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
    {
      if (a.x[i] > 0)
	n = i;
    }


  return n;
}



int mlt(int x, int y)
{

  if (x == 0 || y == 0)
    return 0;

  return ((x + y - 2) % (O - 1)) + 1;
}


int mltn2(int n,int x){
  int i,j;

if(n==0)
return 1;
  i=x;
    for(j=0;j<n-1;j++)
      i=mlt(i,x);

  return i;
}



int mltn(int n, int x) {
    int ret = 1;
    while (n > 0) {
        if (n & 1) ret = gf[mlt(fg[ret] , fg[x])] ;  // n の最下位bitが 1 ならば x^(2^i) をかける
        x = (x * x)%P;
        n >>= 1;  // n を1bit 左にずらす
    }
    return ret
    ;
}

int nmul(int n,int b){
  int i=0,c=b;

  for(i=0;i<n;i++)
  c=c*b%P;
  
  return c;
}

int mlt2(int n, int x)
{
  int i, j;

  if (n == 0)
    return 1;
  i = x;
  for (j = 0; j < n - 1; j++)
    i = mlt(i, x);

  return i;
}


int qinv(unsigned short b)
{
  int i;

  if (b == 0)
    return 0;


  return (O-fg[b])%(O-1)+1;
}



//配列からベクトル表現の多項式へ変換する
vec
Setvec (int n)
{
  int i;
  vec v = { 0 };


  for (i = 0; i < n; i++)
    {
      v.x[n - 1 - i] = c[i];
    }


  return v;
}



//配列の値を係数として多項式に設定する
OP
setpol (unsigned short f[], int n)
{
  OP g;
  vec a;
  int i;
  
  memset (c, 0, sizeof (c));
  //memcpy (c, f, n);
  for(i=0;i<n;i++){
    c[i]=f[i];
    //printf("%d,",f[i]);
  }
//exit(1);  
  a = Setvec (n);

  g = v2o (a);


  return g;
}


//多項式を表示する(default)
void
printpol (vec a)
{
  int i, n;

  n = deg (a);

  //printf ("baka\n");
  assert (("baka\n", n >= 0));



  for (i = n; i > -1; i--)
    {
      if (a.x[i] > 0)
	{
	  printf ("%u", a.x[i]);
	  //if (i > 0)
	  printf ("x^%d", i);
	  //if (i > 0)
	  printf ("+");
	}
    }
  //  printf("\n");

  return;
}

//多項式の代入値
unsigned short
xtrace (OP f, unsigned short x)
{
  int i, d;
  unsigned short u = 0,v=1;
  

  d = deg (o2v (f));
//printpol(o2v(f));
//printf(" =ff\n");

  for (i = 0; i < d + 1; i++)
    {
      v=0;
      if(f.t[i].a>0){
        v=1;
        for(int j=0;j<f.t[i].n;j++){
          v=(v*x);
        }
          v=(v*f.t[i].a)%O;
        
      //printf("\nv=%d",v);
    }
    u=(u+v);
    }
  //printf("u=%d\n",u%O);

  return u%O;
}


void makefg(int n){
unsigned short i,j,count=0;

for(i=0;i<O;i++){

  for(j=0;j<O;j++){
    if(gf[i]==j){
      fg[j]=i;
      count++;
    }
  }
  
}
  printf("unsigned short fg[%d]={",O);
  for(i=0;i<O;i++)
  printf("%d,",fg[i]);
printf("};\n");
printf("count=%d\n",count);
//exit(1);
 
return;
}




//多項式の次数(degのOP型)
int
odeg (OP f)
{
  int i, j = 0, k;


  //k=terms(f);
  for (i = 0; i < 512; i++)
    {
      if (j < f.t[i].n && f.t[i].a > 0)
	j = f.t[i].n;
    }

  return j;
}

OP minus(OP f){
unsigned int i,j;

j=deg(o2v(f));
for(i=0;i<j+1;i++){
  if(f.t[i].a>0)
f.t[i].a=P-f.t[i].a;
}

return f;
}

//リーディグタームを抽出(default)
oterm
LT (OP f)
{
  int i, k;
  oterm t = { 0 };


  //k = deg (o2v (f));
  for (i = 0; i < DEG; i++)
    {
      //printf("a=%d %d\n",f.t[i].a,f.t[i].n);
      if (f.t[i].a > 0)
	{
	  t.n = f.t[i].n;
	  t.a = f.t[i].a%P;
	}
    }

  return t;
}

//OP型を正規化する
OP
conv (OP f)
{
  vec v = { 0 };
  OP g = { 0 };

  v = o2v (f);
  g = v2o (v);

  return g;
}


//20200816:正規化したいところだがうまく行かない
//多項式の足し算
OP
oadd (OP f, OP g)
{
  vec a = { 0 }
  , b =
  {
  0}
  , c =
  {
  0};
  int i, j, k, l = 0;
  OP h = { 0 },f2={0},g2={0};

  //for(i=0;i<257;i++)
  // printf("%d %d %d %d %d\n",i,f.t[i].a,f.t[i].n,g.t[i].a,g.t[i].n);

   //  exit(1);
  f=conv(f);
  g=conv(g);

  a = o2v (f);
  //exit(1);
  b = o2v (g);

  j=deg(o2v(f));
  l=deg(o2v(g));
    printpol(o2v(f));
    printf(" =f_& %d\n",j);
    printpol(o2v(g));
    printf(" =g_& %d\n",l);
    //exit(1);

  if (j >= l)
    {
      k = j + 1;
    }
  else
    {

      k = l + 1;

    }
  //for(i=0;i<k;i++)
  //printf("%d %d\n",i,b.x[i]);
  //  exit(1);
  
  for (i = 0; i < k; i++)
    {
      //if(a.x[i]>b.x[i])
      c.x[i]=(a.x[i]+b.x[i])%P;
      if(a.x[i]==b.x[i])
      c.x[i]=a.x[i]*2%P;

    }
  // 
  h = v2o (c);
  printpol(o2v(h));
  printf(" =====in oadd\n");

  return h;
}



//多項式を項ずつ掛ける
OP
oterml (OP f, oterm t)
{

  //assert (op_verify (f));
  int i, k,j;
  OP h = { 0 };
  vec test;
  unsigned short n;

  //f=conv(f);
  k = odeg (f);
  j=0;
  for (i = 0; i < k + 1; i++)
    {
      h.t[i].n = f.t[i].n + t.n;
      h.t[i].a = (f.t[i].a * t.a)%P;
    }
  
  //h=conv(h);
  //assert (op_verify (h));
  return h;
}


//多項式の掛け算
OP
omul (OP f, OP g)
{
  f=conv(f);
  g=conv(g);
  //assert (op_verify (f));
  //assert (op_verify (g));
  int i, count = 0, k, l;
  oterm t = { 0 };
  OP h = { 0 }, e = {
    0
  }, r = {
    0
  };
  vec c = { 0 };

  k = odeg (f);
  l = odeg (g);
  if (l > k)
    {
      k = l;
    }

  for (i = 0; i < k + 1; i++)
    {
      t = g.t[i];
      e = oterml (f, t);
      h = oadd (h, e);
    }
  //assert (op_verify (h));
  return h;
}




OP osub(OP f,OP g){
vec a={0},b={0},d={0};
int i,k,l,m;
OP ans={0};

  a=o2v(f);
  b=o2v(g);
  l=deg(a);
  m=deg(b);
  if(l>=m){
    k=l;
  } else{
    k=m;
  }
  for(i=0;i<k+1;i++){
    if(a.x[i]>=b.x[i]){
      d.x[i]=a.x[i]-b.x[i];
    }
    
    else{
    d.x[i]=(P-(a.x[i]+b.x[i]));
    }
    /*
    if(d.x[i]<0){
      printf("%d\n",d.x[i]);
      d.x[i]+=P;
    }
    */
  }
  
  ans=v2o(d);

return ans;
}

OP confer(OP f,int a){
  vec r;
  int n,i;
  OP g;

  r=o2v(f);
  n=deg(r);
  for(i=0;i<n+1;i++)
    r.x[i]=(r.x[i]*a)%P;
  g=v2o(r);

return g;
}

int oequ(OP f,OP g){
vec v,x;
int i,flg=0;

v=o2v(f);
x=o2v(g);
for(i=0;i<512;i++){
    if(v.x[i]!=x.x[i])
    return -1;
}

return 0;
}


//aに何をかけたらbになるか
unsigned short
equ (unsigned short a, unsigned short b)
{
  int i;


  for (i = 0; i < N; i++)
    {
      if ((a*i)%P == b)
	break;
    }
  return i;
}






void mkmf(){
  int i,j,k,count=0;
  OP f={0},g={0},h={0},w={0},s={0},u={0};
  vec b={0},a={0},d={0},t={0},v={0};
  oterm o;
  unsigned short ccp[4]={0};

  if(O==1331){
  for(i=0;i<K+1;i++)
  ccp[i]=pp[0][i];
  }
  if(O==2197){
  for(i=0;i<K+1;i++)
  ccp[i]=pp[1][i];
  }
  if(O==4913){
  for(i=0;i<K+1;i++)
  ccp[i]=pp[2][i];
  }
  if(O==6859){
  for(i=0;i<K+1;i++)
  ccp[i]=pp[3][i];
  }
  if(O==27){
  for(i=0;i<K+1;i++)
  ccp[i]=pp[4][i];
  }
  if(O==343){
  for(i=0;i<K+1;i++)
  ccp[i]=pp[5][i];
  }

g=setpol(ccp,4);
//b.x[0]=2;
//b.x[1]=9;
a.x[1]=1;
v.x[3]=1;
u=v2o(v);
d.x[3]=1;
printpol(o2v(g));
printf(" =g\n");
//exit(1);

//g=v2o(b);
s=v2o(a);
//s.t[1].a=1;
//s.t[1].n=1;
//gf[12]=P;
//gf[13]=P*P;

w=g;
printpol(o2v(w));
printf(" =w\n");
printpol(o2v(g));
printf(" =g\n");
printpol(o2v(s));
printf(" =s\n");

printf("\n");
//for(i=0;i<P;i++)
gf[0]=0;
gf[1]=1;
gf[2]=P;
gf[3]=P*P;
gf[4]=xtrace(g,P);
printf("\naa=%d\n",gf[4]);
//exit(1);
//w=omul(w,s);
//gf[12]=1111;
count=4+1;
while(1){
g=omul(g,s);
printpol(o2v(g));
printf(" =g\n\n");
printf(" papaya\n");

//exit(1);

o=LT(g);
memset(d.x,0,sizeof(d));
if(o.n==K){
  d.x[o.n]=o.a;
  h=v2o(d);
  g=osub(g,h);
  f=confer(w,o.a);
  g=oadd(g,f);
  //w=omod(w,u);
  printpol(o2v(f));
  printf("\n");
}

gf[count]=xtrace(g,P);
printf("count=%d %d ",count,gf[count]);
printpol(o2v(g));
printf(" =gg\n\n");
if(gf[count]==1){
printf("count=%d\n",count);
break;
}

count++;
}

//exit(1);

//printpol(o2v(f));
//printf(" =f\n");
//printf("gf[%d]={\n",O);
printf("unsigned short gf[%d]={",O);
for(i=0;i<O;i++)
printf("%d,",gf[i]);
printf("};");
printf("\n");

//exit(1);
}

// invert of integer
unsigned short inv(unsigned short a,unsigned short n)
{
  unsigned short d,x,s,q,r,t,gcd;
  d = n;
  x = 0;
  s = 1;
  while (a != 0){
    q = d / a;
    r = d % a;
    d = a;
    a = r;
    t = (x - q * s)%P;
    x = s;
    s = t;
  }
  gcd = d%P;  

  return ((x + n) % (n / d))%P;
}


//有限体の元の逆数
unsigned short
oinv(unsigned short a)
{
    int i;

    if (a == 0)
        return -1;

    for (i = 0; i < N; i++)
    {
        if ((a*i)%P == 1)
            return (unsigned short)i;
    }

    printf("no return \n");
    //  exit (1);
}


//多項式を単行式で割る
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
        //printf("%u\n",s.a);
    }
    else if (t.n == 0 && t.a > 0)
    {
        s.a = (tt.a*inv(t.a,P))%P;
        s.n = tt.n;
    }else{
      printf("debug in LTdiv\n");
      exit(1);
    }

    return s;
}


//多項式の剰余を取る
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

    //printf ("in omod\n");
    //exit(1);

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
        //exit(1);
        if (odeg((f)) == 0 || odeg((g)) == 0)
        {
            //      printf("blake1\n");
            break;
        }

        if (c.n == 0 || b.n == 0)
            break;
            if(LT(f).a==4 && deg(o2v(f))==0)
            exit(1);
    }
    printpol(o2v(f));
    printf("\n");
    //exit(1);

    return f;
}


//項の数
int
oterms (OP f)
{
  int i, count = 0;

  for (i = 0; i < DEG; i++)
    if (f.t[i].a > 0)
      count++;

  return count;
}


//モニック多項式にする
OP coeff(OP f)
{
    int i, j, k;
    vec a, b;
    oterm t;

    t=LT(f);
    //f = conv(f);
    k = odeg((f)) + 1;
    for (i = 0; i < k; i++)
        f.t[i].a = (f.t[i].a*inv(t.a,P))%P;

    return f;
}

OP kei(unsigned short u,OP g){
  vec v={0};
  int j,i=0;
  OP f={0};

  v=o2v(g);
  i=deg(v);
  printpol(v);
  printf("\n");
  for(j=0;j<i+1;j++)
  v.x[j]=(v.x[j]*u)%P;
  f=v2o(v);
  printpol(v);
  printf("vom\n");
  //exit(1);

  return f;
}

//多項式の商を取る
OP
odiv (OP f, OP g)
{

  f = conv (f);
  g = conv (g);
  //assert (op_verify (f));
  //assert (op_verify (g));
  int i = 0, j, n, k;
  OP h = { 0 }, e = { 0 }, tt = { 0 };
  oterm a, b = { 0 }, c = { 0 };

printpol(o2v(f));
printf("\n");
printpol(o2v(g));
printf("\n");
//exit(1);

  if (LT (f).n == 0 && LT (g).a == 0)
    {
      printf ("baka^\n");
      //return f;
      exit (1);
    }
  if (LT (g).a == 0)
    {
      printf("a==0\n");
      //print_trace ();
      exit (1);
    }
  if (LT (g).n == 0 && LT (g).a > 1)
    return g; //coeff (f);

  k = odeg (g);
  b = LT (g);
  if (b.a == 1 && b.n == 0)
    return f;
  if (b.a == 0 && b.n == 0)
    {
      printf ("baka in odiv\n");
      exit (1);
    }
  if (odeg ((f)) < odeg ((g)))
    {
      return f;
      //  a=LT(f);
    }
  OP null={0};
  i = 0;
  while (LT (f).n > 0 && LT (g).n > 0)
    {
      c = LTdiv (f, b);
      c.a=c.a%P;
      assert (c.n < DEG);
      tt.t[i] = c;
      i++;

      h = oterml (g, c);
      f = oadd (f, minus(h));
      printpol(o2v(h));
      printf(" ===h in_odiv\n");
      printpol(o2v(f));
      printf(" ===f in_odiv\n");
      if(deg(o2v(f))==0 && LT(f).a>0){
      return f;
      //exit(1);
      }else{
        //return g;
      }
      if (odeg ((f)) == 0 || odeg ((g)) == 0)
	  {
	    printf ("blake2\n");
	    break;
	  }
    if(oequ(f,g)==0){
      printpol(o2v(tt));
      printf("\n");
      break;
      //exit(1);
    }
      if (c.n == 0)
	break;
    }

  // tt は逆順に入ってるので入れ替える
  OP ret = { 0 };
  int tt_terms = oterms (tt);
  for (i = 0; i < tt_terms; i++)
    {
      ret.t[i] = tt.t[tt_terms - i - 1];
    }
  ret = conv (ret);
  printpol(o2v(ret));
  printf("\n");
//  exit(1);

  //assert (op_verify (ret));
  return ret;
}


// invert of polynomial
OP pinv(OP a, OP n)
{
    OP d = {0}, x = {0}, s = {0}, q = {0}, r = {0}, t = {0}, u = {0}, v = {0}, w = {0}, tt = {0}, gcd = {0}, tmp = {0};
    oterm b = {0};
    vec vv = {0}, xx = {
                      0};

    if (odeg((a)) > odeg((n)))
    {
        tmp = a;
        a = n;
        n = tmp;
        printf("baka_i\n");
        //exit (1);
    }
    if (LT(a).a == 0)
    {
        printf(" a ga 0\n");
        exit(1);
    }

    tt = n;

    d = n;
    x.t[0].a = 0;
    x.t[0].n = 0;
    s.t[0].a = 1;
    s.t[0].n = 0;
    while (odeg((a)) > 0)
    {
      printpol(o2v(a));
      printf("==========a\n");
        if (odeg((a)) > 0)
            r = omod(d, a);
        if (LT(a).a == 0)
            break;
        if (LT(a).n > 0)
            q = odiv(d, a);
        printpol(o2v(r));
        printf("=====r\n");
        printpol(o2v(q));
        printf("=====q\n");
        printpol(o2v(d));
        printf("=====d\n");
        printpol(o2v(a));
        printf("=====a\n");
        printpol(o2v(t));
        printf("=====t\n");
        printpol(o2v(x));
        printf("=====x\n");
        printpol(o2v(s));
        printf("=====s\n");
        //exit(1);

        d = a;
        a = r;
        t = osub(x, (omul(q, s)));
        ////printpol (o2v (a));
        //printf ("\nin roop a==================%d\n", odeg ((a)));
        //printf ("\n");

        x = s;
        s = t;
    }
    printpol(o2v(a));
    printf("=====final\n");

    // exit(1);
    //  if(LT(a).a>0){
    d = a;
    a = r;
    ////printpol (o2v (a));
    //printf ("\nin roop a|==================%d\n", odeg ((a)));
    //printf ("\n");

    x = s;
    s = t;

    ////printpol (o2v (d));
    //printf ("\nout1================\n");
    gcd = d; // $\gcd(a, n)$
    printpol(o2v(gcd));
    printf(" =========gcd\n");
    //exit(1);
    //printf ("\n");
    ////printpol (o2v (n));
    //printf ("\n");
    //printf ("out2===============\n");

    printf("before odiv\n");
    //w=tt;

    b = LT(w);
    v = oadd(x, n);
    printpol(o2v(v));
    printf("===========v\n");
    w = tt;

    if (LT(v).n > 0 && LT(w).n > 0)
    {
        u = omod(v, w);
        printpol(o2v(u));
        printf("========u\n");
        printpol(o2v(v));
        printf("========v2\n");
        printpol(o2v(w));
        printf("========w\n");
        
    }
    else
    {
        //printpol (o2v (v));
        printf(" v===========\n");
        //printpol (o2v (x));
        printf(" x==0?\n");
        //printpol (o2v (n));
        printf(" n==0?\n");

        exit(1);
    }
    //caution !!
    if (LT(u).a > 0 && LT(d).a > 0)
    {
        u = odiv(u, d);
        printpol(o2v(u));
        printf("========u2\n");
    }

    if (LT(u).a == 0 || LT(d).a == 0)
    {
        printf("inv div u or d==0\n");
        // exit(1);
    }

    if (LT(u).a == 0)
    {
        printf("no return at u==0\n");
        exit(1);
    }
    printpol(o2v(u));
    printf(" ==========u\n");
    exit(1);
    return u;
}

OP muri(OP x,OP m){
unsigned short s[K+1]={0};
OP t={0},u;

for(int i=0;i<P;i++){
  memset(s,0,sizeof(s));
  for(int j=0;j<P;j++){
    s[4]=j;
    s[5]=i;
  t=setpol(s,K+1);
  u=omod(omul(t,x),m);
  if(LT(u).n==0 && LT(u).a==1)
  break;
}
  if(LT(u).n==0 && LT(u).a==1)
  break;
}

return t;
}

OP hyoe[8192]={0};
unsigned short sue[8192]={0};
unsigned short sie[8192]={0};
void tas(){
  int i,j,k;
  //unsigned short ccp[4]={0,0,11,2};
  OP g={0},h={0},f={0},m={0};
  vec v={0},x={0},w={0};

  v.x[1]=11;
  v.x[0]=2;
  m=v2o(v);
  printpol(v);
  printf("\n");
  g=m;
  w.x[1]=1;
  printpol(w);
  printf("\n");

  h=v2o(w);
  //g=omul(g,h);
  printpol(o2v(g));
  printf("\n");
// exit(1);

  memset(v.x,0,sizeof(v.x));
for(i=0;i<2;i++){
  v.x[0]=i;
  hyoe[i]=v2o(v);
}
memset(v.x,0,sizeof(v.x));
v.x[1]=1;
hyoe[2]=v2o(v);
memset(v.x,0,sizeof(v));
v.x[2]=1;
hyoe[3]=v2o(v);
x.x[3]=1;
f=v2o(x);
  hyoe[4]=g;
  g=omul(g,h);
  if(deg(o2v(g))>=3){
  g=oadd(g,kei(LT(g).a,m));
  g=omod(g,f);
  printpol(o2v(g));
  printf("====3\n");
  }
//  exit(1);

  printpol(o2v(g));
  hyoe[5]=g;
  //exit(1);

for(i=6;i<O;i++){
  g=omul(h,g);
  printpol(o2v(g));
  //exit(1);
  if(deg(o2v(g))>=3){
  g=oadd(g,kei(LT(g).a,m));
  g=omod(g,f);
  printpol(o2v(g));
  printf("====3\n");
  //exit(1);
  printpol(o2v(m));
  printf("\n");
  }
  //exit(1);
//  g=oadd(g,

  printpol(o2v(g));
  printf("\n");
  hyoe[i]=g;
  printpol(o2v(g));
  printf("\n");
}
//exit(1);
for(i=0;i<2;i++)
sue[i]=i;
for(i=2;i<O;i++)
sue[i]=xtrace(hyoe[i],P);
for(i=0;i<O;i++)
sie[sue[i]]=i;

}

unsigned short plus(unsigned short a,unsigned short b){
  unsigned short u;
  
  u=xtrace(oadd(hyoe[fg[a]],hyoe[fg[b]]),P);

return u;
}

unsigned short hiku(unsigned short a,unsigned short b){
  unsigned short u;
  
  u=xtrace(osub(hyoe[fg[a]],hyoe[fg[b]]),P);

return u;
}


OP synd(unsigned short zz[], int kk)
{
    unsigned short syn[K2] = {0}, s = 0;
    int i, j, t1;
    OP f = {0};

    printf("in synd2\n");

    for (i = 0; i < kk; i++)
    {
        syn[i] = 0;
        s = 0;
        //#pragma omp parallel num_threads(16)
        for (j = 0; j < N; j++)
        {   if(zz[j]>0)
            s = plus(s, mat[j][i]);
            printf("%d,",mat[j][i]);
        }
        printf("\n");
        syn[i] = s;
         printf ("syn%d,", syn[i]);
    }
    // printf ("\n");

    f = setpol(syn, kk);
    printpol(o2v(f));
    printf(" syn=============\n");
      //exit(1);

    return f;
}


unsigned short vb[K * 2][N] = {0};
unsigned short gt[K * 2][K * 2] = {0};

void van(int kk)
{
    int i, j, k,t[K2+1]={0};

    printf("van der\n");

    for (i = 0; i < N; i++)
        vb[0][i] = 1;
    //#pragma omp parallel for private(i, j)
      for (j = 1; j < N; j++)
    {
            for (i = 0; i < kk; i++)
        {
            vb[i][j] = gf[mltn(i, j)];
            //vb[i][j] = gf[mltn2(i, fg[j])];
            mat[j][i]=vb[i][j];
            printf("%d,", vb[i][j]);
            //t[i+1]=plus2(t[i+1],vb[i][j]);
        }
        printf("\n");
    }
    /*
    for(i=0;i<K+1;i++)
    printf("%d,",t[i]);
    printf("\n");
    */
//exit(1);
}


// chen探索
vec chen(OP f)
{
    vec e = {0};
    int i, count = 0, n, x = 0;
    unsigned short z;

    n = deg(o2v(f));
    // exit(1);
    //#pragma omp parallel for private(i)
    for (x = 0; x < N; x++)
    {
        z = 0;
        //#pragma omp parallel for reduction (^:z)
        for (i = 0; i < n + 1; i++)
        {
            if (f.t[i].a > 0)
                z = (z+mlt(mltn(f.t[i].n, fg[x]), fg[f.t[i].a]))%P;
        }
        if (z == 0)
        {
            e.x[count] = x;
            count++;
            printf("%d\n", x);
        }
    }
    // printpol(e);
    // printf(" ==eee!\n");
    // exit(1);

    return e;
}



OP kof(unsigned short c, OP f)
{
    int i, j, k;
    vec b = {0}, h = {0};
    OP g = {0};

    c = fg[c];
    b = o2v(f);
    k = deg(b);
    for (i = 0; i < k + 1; i++)
    {
        h.x[i] = gf[mlt(c, fg[b.x[i]])];
    }
    g = v2o(h);

    return g;
}



OP bma(unsigned short s[], int kk)
{
    int i, j, k, ll = 0, l, d[2 * K + 1] = {0};
    OP lo[2 * K + 1] = {0}, b[2 * K + 1] = {0}, t[2 * K + 1] = {0}, a = {0}, f = {0}, h = {0}, g = {0}, hh = {0};
    vec v = {0}, x = {0}, w = {0};

    x.x[1] = 1;
    h = v2o(x);
    v.x[0] = 1;
    f = v2o(x);
    lo[0] = v2o(v);
    b[0] = lo[0];
    ll = 0;
    for (j = 1; j < T * 2 + 1; j++)
    {
        v = o2v(lo[j - 1]);
        k = 0;
        // printpol(v);
        // printf(" ==lo\n");

        l = deg(o2v(lo[j - 1]));
        for (i = 1; i < l + 1; i++)
        {
            k = plus(k,gf[mlt(fg[v.x[i]], fg[s[j - i]])]);
            // printf("v[%d]=%d\n", i, v.x[i]);
        }
        d[j] = plus(s[j] , k);
        // printf("d[%d]=%d\n", j, d[j]);
        if (d[j] == 0)
        {
            lo[j] = lo[j - 1];
            b[j] = omul(b[j - 1], h);
            // ll=j-1;
        }
        else // if (d[j] > 0)
        {
            g = omul(kof(d[j], h), b[j - 1]);
            t[j] = oadd(lo[j - 1], g);
            lo[j] = t[j];
            if (ll * 2 > (j - 1))
            {
                // lo[j]=t[j];
                b[j] = omul(b[j - 1], h);
            }
            else // if(2*ll <= j)
            {
                // printpol(o2v(t[j]));
                // printf("==t[%d]\n", j);
                b[j] = kof(gf[oinv(d[j])], lo[j - 1]);
                // lo[j]=t[j];
                ll = j - ll;

                if (j == 2 * T)
                {
                    if (!(d[T * 2 - 1] == 0 && d[T * 2 - 3] == 0 && deg(o2v(lo[j - 1])) == T) || !(deg(o2v(lo[j - 1])) == T))
                    {
                        if ((d[T * 2 - 1] == 0 && odeg(lo[j - 2]) == T - 1))
                        {
                            lo[j - 1] = omul(lo[j - 2], h);
                            // printpol(o2v(lo[j - 1]));
                            // printf("\n");
                        }
                    }
                    break;
                }
            }
        }
        printf("@l=%d\n", ll);
        k = 0;
        // printpol(o2v(b[j]));
        // printf(" ==b[%d]\n", j);
    }

    k = 0;
    int count = 0;
    // printpol(o2v(lo[j - 1]));
    // printf(" ==coef\n");
    if (odeg(lo[j - 1]) == T)
    {
        x = chen(lo[j - 1]);
    }
    else
    {
        printf("else baka\n");
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
        printf("something buggy in sendrier\n");
        exit(1);
    }

    // return count;
    if(deg(o2v(lo[j-1]))==T){
    return lo[j - 1];
    }else{
      printf("baka\n");
      exit(1);
    }
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
            printf("v.x[K2-1]=%d\n", v.x[K - 1]);
            // break;
            //
            exit(1);
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
        //B++;
    }
    if (count < T)
    {
        printf("error is too few\n");

        //AA++;
        // memcpy (zz, z1, sizeof (zz));
        /*
      printf("{");
      for (i = 0; i < D; i++)
        printf("%d,", z1[i]);
      printf("};\n");
      printpol(o2v(w));
      printf(" =========goppa\n");
      printsage(o2v(w));
      printf(" =========sage\n");
      printsage(o2v(f));
      printf(" =========syn\n");
      printpol(o2v(f));
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

    memset(z1, 0, sizeof(z1));

    while (j < num)
    {
        l = rand() % N;
        // printf ("l=%d\n", l);
        if (0 == z1[l])
        {
            z1[l] = 1;
            // printf("l=%d\n", l);
            j++;
        }
    }
}




int main(){
unsigned short i,count=0,c2=0;
unsigned short aaa[O]={0};

OP ff,k,uu1,uu2,vv1,vv2,s,l,u3,v3,u,ll,t;
unsigned short ss[K+1]={0};
vec x={0},v={0};
OP g={0};
unsigned short zz[1331]={0};

mkmf();
exit(1);

tas();



/*
van(K);
  printf("%d %d\n",fg[1111],fg[293]);
  printpol(o2v(hyoe[fg[1111]]));
  printf("\n");
  printf("plus=%d\n",plus(101,1111));
//exit(1);

for(i=10;i<15;i++)
{
  for(int j=10;j<K+10;j++){
    //printf("%d,",mat[i][j]);
  ss[j-10]=plus(ss[j-10],mat[i][j]); //(ss[j+1-10]+mat[i][j])%1331; 
//printf("%d,",mat[i][j]);
}
//printf("\n");
}
for(i=0;i<K+1;i++)
printf("|%d,",ss[i]);
printf("\n");
//exit(1);


//exit(1);
/*
mkerr(zz,T);
ff=synd(zz,K);
v=o2v(ff);
printpol(v);
memset(ss,0,sizeof(ss));
for(int i=0;i<K;i++)
ss[i+1]=v.x[i];
//exit(1);

g=bma(ss,K);
x=chen(g);
//for(i=0;i<1331;i++)
//if(x.x[i]>0)
//printf("%d,",i);
//printf("\n");
ero2(x);
for(i=0;i<1331;i++)
if(zz[i]>0)
printf("%d,",i);
printf("\n");
exit(1);
*/

//exit(1);
/*
for(i=0;i<27;i++){
  printpol(o2v(hyoe[i]));
  printf("\n");
}
//exit(1);
printf("unsigned short gf[27]={");
for(i=0;i<27;i++){
  
    printf("%d,",sue[i]);
}
printf("};\n");
printf("unsignd short fg[27]={\n");
for(i=0;i<27;i++){
   printf("%d,",sie[i]);
}
  printf("};\n");
  */
  /*
  printf("%d %d\n",gf[101]+gf[1111],fg[101]+fg[1111]);

  printf("%d\n",gf[mlt(fg[101],fg[1111])]);
  printf("%d\n",gf[mltn(100,fg[11])]);
  printf("%d %d %d\n",O,qinv(1111),gf[mlt(qinv(1111),fg[1111])]);
  exit(1);

  printf("hyou3=");
  printpol(o2v(hyoe[4]));
  printf("%d\n",xtrace(hyoe[4],P));
  printf("\n");
  printf("hyou4=");
  printpol(o2v(hyoe[5]));
  printf("\n");
  printf("%d\n",xtrace(hyoe[5],P));
  printf("\n");
  printpol(o2v(oadd(hyoe[4],hyoe[5])));
  printf("\n");

exit(1);


ff=setpol(f,K+1);
uu1=setpol(u1,K+1);
uu2=setpol(u2,K+1);
vv1=setpol(v1,K+1);
vv2=setpol(v2,K+1);

printpol(o2v(ff));
printf("\n");
printpol(o2v(uu1));
printf("\n");
printpol(o2v(uu2));
printf("\n");
printpol(o2v(vv1));
printf("\n");
printpol(o2v(vv2));
printf("\n");

ll=(omul(vv2,vv2));
printpol(o2v(ll));
printf("\n");
ll=(oadd(ff,minus(ll)));
printpol(o2v(ll));
printf("\n");
k=(odiv(ll,uu2));
printpol(o2v(k));
printf(" ('A`)\n");
ll=odiv(uu1,uu2);
printpol(o2v(ll));
printf("========div\n");
ll=omod(uu2,uu1);
printpol(o2v(ll));
printf("=======mod\n");
//exit(1);
OP tt={0};

tt=muri(uu2,uu1);
//tt=pinv(uu2,uu1);
printpol(o2v(tt));
//omod(omul(t,uu2),uu1);
printf(" ===inv\n");
//exit(1);
ll=oadd(vv1,minus(vv2));
printpol(o2v(ll));
printf("\n");
ll=omul(ll,tt);
s=omod(ll,uu1);
printpol(o2v(s));
printf("\n");
//exit(1);

l=omul(s,uu2);
printpol(o2v(l));
printf("\n");
u=odiv(oadd(k,minus(omul(s,oadd(l,omul(s,vv2))))),uu1);
printpol(o2v(u));
printf("\n");
u3=coeff(u);
printpol(o2v(u3));
printf(" =======u3\n");
v3=omod(minus(oadd(l,vv2)),u3);
printpol(o2v(v3));
printf(" =========v3\n");

//exit(1);
*/

/*
  mkmf();

  makefg(O);
tas();
printf("2*a^3=%d\n",hiku(12,12));

exit(1);

count=0;
c2=0;
for(i=0;i<O;i++){
  if(gf[i]>0)
  count++;
  if(fg[i]>0){
  c2++;
  }else{
   // printf("i=%d\n",i);
  }
}

printf("%d %d\n",count,c2);
//exit(1);
*/


return 0;
}


