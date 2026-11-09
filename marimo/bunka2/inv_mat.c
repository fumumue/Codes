//
// (over finite field) Gauss-Jordan法による逆行列
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//#include "global-p.h"
//#include "struct-p.h"

//#include "chash-p.c"

#define MAXN 4

void rp(unsigned short *a)
{
	int x, j;
	time_t t;

	srand(clock() + time(&t));

	for (int i = 0; i < N; i++)
	{
		a[i] = i;
	}
	for (int i = 0; i < N - 2; i++)
	{
		// rand from i+1 to F-1
		j = (rand() % (N - 1 - i)) + i + 1;

		// swap a[i] and a[j]
		x = a[j];
		a[j] = a[i];
		a[i] = x;
	}
	if (a[N - 1] == N - 1)
	{
		a[N - 1] = a[N - 2];
		a[N - 2] = N - 1;
	}
}

/*
int mlt(int x, int y){

	if(x==0||y==0)
		return 0;

  return ((x+y-2)%(N-1))+1;
}


int mltn(int n,int x){
  int i,j;

  if(n==0)
	return 1;
  i=x;
	for(j=0;j<n-1;j++)
	  i=mlt(i,x);

  return i;
}
*/

int Inv(unsigned short b)
{

	if (b == 0)
		return 0;

	for (int i = 0; i < N; i++)
	{
		if (gf[mlt(i, b)] == 1)
			return i;
	}

	return -1;
}

MTX gauss(MTX a)
{
	int buf = 0;
	unsigned short ff = 0, inv_a[F][F] = {0};
	MTX TT = {0}, b;
	b = a;
	// unsigned short a[F][F]={0};
	//\92P\88ʍs\97\F1\82\F0\8D\EC\82\E9
	for (int i = 0; i < F; i++)
	{
		for (int j = 0; j < F; j++)
		{
			inv_a[i][j] = (i == j) ? 1 : 0;
		}
	}
	//\91|\82\AB\8Fo\82\B5\96@
	for (int i = 0; i < F; i++)
	{
		buf = gf[Inv(a.x[i][i])];
		for (int j = 0; j < F; j++)
		{
			a.x[i][j] = gf[mlt(fg[a.x[i][j]], fg[buf])];
			inv_a[i][j] = gf[mlt(fg[inv_a[i][j]], fg[buf])];
		}
		for (int j = 0; j < F; j++)
		{
			if (i != j)
			{
				buf = a.x[j][i];
				for (int k = 0; k < F; k++)
				{
					a.x[j][k] ^= gf[mlt(fg[a.x[i][k]], fg[buf])];
					inv_a[j][k] ^= gf[mlt(fg[inv_a[i][k]], fg[buf])];
				}
			}
		}
	}
	//\8Bt\8Ds\97\F1\82\F0\8Fo\97\CD
	for (int i = 0; i < F; i++)
	{
		printf("{");
		for (int j = 0; j < F; j++)
		{
			printf(" %d,", inv_a[i][j]);
			TT.x[i][j] = inv_a[i][j];
		}
		printf("},\n");
	}

	for (int i = 0; i < F; i++)
	{
		for (int j = 0; j < F; j++)
		{
			for (int k = 0; k < F; k++)
				ff ^= TT.x[i][k] & b.x[k][j];
			TT.x[i][j] = ff;
		}
	}
	for (int i = 0; i < F; i++)
	{
		for (int j = 0; j < F; j++)
			printf("%d,", TT.x[i][j]);
		printf("\n");
	}

	return TT;
}

// inverse matrix
MTX matinv(MTX a, int n)
{

	// unsigned short a[F][F];     //={{1,2,0,1},{1,1,2,0},{2,0,1,1},{1,2,1,1}}; //入力用の配列
	unsigned short inv_a[N][N];	  // ここに逆行列が入る
	unsigned short buf;			  // 一時的なデータを蓄える
	unsigned short b[N][N] = {0}; //, dd[N][N] = {0};
	// int i, j, k, count;           // カウンタ
	//  MTX a={0};
	unsigned short c[N][N] = {0};
	MTX z = {0};
	//  unsigned short cc[N][N] = {0};

lab:
	memset(b, 0, sizeof(b));
	memset(a.x, 0, sizeof(a.x));
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a.x[i][j] = rand() % 256;
			printf("%d,", a.x[i][j]);
		}
		printf("\n");
	}
	// exit(1);
	//  printf("\n");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			c[i][j] = a.x[i][j];
	}
	// 単位行列を作る
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			inv_a[i][j] = (i == j) ? 1 : 0;
		}
	}
	// 掃き出し法
	// #pragma omp parallel for num_threads(omp_get_max_threads()) //private(i,j,k)
	for (int i = 0; i < n; i++)
	{
		buf = gf[Inv(fg[a.x[i][i]])];
		for (int j = 0; j < n; j++)
		{
			a.x[i][j] = gf[mlt(fg[buf], fg[a.x[i][j]])];
			inv_a[i][j] = gf[mlt(fg[buf], fg[inv_a[i][j]])];
		}
		for (int j = 0; j < n; j++)
		{
			if (i != j)
			{
				buf = a.x[j][i];
				for (int k = 0; k < n; k++)
				{
					a.x[j][k] ^= gf[mlt(fg[a.x[i][k]], fg[buf])];
					inv_a[j][k] ^= gf[mlt(fg[inv_a[i][k]], fg[buf])];
				}
			}
		}
	}

	int count = 0;
	// printf("\n\n逆行列を出力\n");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
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
		// printf("\n");
	}
	// exit(1);

	printf("行列を出力\n ={\n");
	// #pragma omp parallel for num_threads(omp_get_max_threads()) //private(i,j,k)
	for (int i = 0; i < n; i++)
	{
		printf("{");
		for (int j = 0; j < n; j++)
		{
			// a[i][j]=rand()%N;
			printf("%3d,", a.x[i][j]);
		}
		printf("},\n");
	}
	printf("};");
	count = 0;
	printf("\n逆行列を出力\n ={\n");
	for (int i = 0; i < n; i++)
	{
		count = 0;
		printf("{");
		for (int j = 0; j < n; j++)
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
	//  exit(1);

	memset(b, 0, sizeof(b));
	// 検算
	//   #pragma omp parallel for num_threads(omp_get_max_threads()) //private(i,j,k)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < n; k++)
				b[i][j] ^= gf[mlt(fg[c[i][k]], fg[inv_a[k][j]])];

			printf("%d,", b[i][j]);
		}
		printf("\n");
	}

	int flg = 0;
	for (int i = 0; i < n; i++)
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

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (b[i][j] == 0 && i != j)
				count++;
		}
	}
	if (flg == n && n * n - n == count)
		return z;

	goto lab;
}

CTX kenzan(CTX a, CTX inv_a, int p)
{
	short tmp = 0;
	static CTX z = {0};
	// 検算
	if (p == 1)
	{
		for (int i = 0; i < F; i++)
		{
			for (int j = 0; j < F; j++)
			{
				tmp = 0;
				for (int k = 0; k < F; k++)
					tmp += ((a.x[i][k] % Pr) * (inv_a.x[k][j] % Pr)) % Pr;
				z.x[i][j] = tmp % Pr;
			}
		}
		/*
		//return z;
		printf("もとの行列\n");
	  int count = 0;
	  for (int i = 0; i < F; i++)
	  {
		for (int j = 0; j < F; j++)
		{
		  printf("%d ", a.x[i][j] % Pr);
		}
		printf("\n");
	  }
	  printf("逆行列\n");
	  for (int i = 0; i < F; i++)
	  {
		for (int j = 0; j < F; j++)
		{
		  printf("*%d ", inv_a.x[i][j] % Pr);
		}
		printf("\n");
	  }

	  printf("検算\n");
	  for (int i = 0; i < F; i++)
	  {
		for (int j = 0; j < F; j++)
		{
		  printf("z%d ", z.x[i][j] % Pr);
		}
		printf("\n");
	  }
	  */
		return z;
	}

	if (p == 2)
	{
		for (int i = 0; i < F; i++)
		{
			for (int j = 0; j < M; j++)
			{
				tmp = 0;
				for (int k = 0; k < F; k++)
					tmp += ((a.x[i][k] % Pr) * (inv_a.x[k][j] % Pr)) % Pr;
				z.x[i][j] = tmp % Pr;
				// printf("%d ", z.x[i][j]);
			}
			// printf("\n");
		}
		/*
	  // int count = 0;
	  printf("かけた結果\n");
	  for (int i = 0; i < F; i++)
	  {
		for (int j = 0; j < M; j++)
		{
		  printf("%d, ", z.x[i][j] % Pr);
		}
		printf("\n");
	  }
	  */
		return z;
	}
}

MTX mulmat(MTX A, MTX B, int flg)
{
	// int i, j, k;
	MTX tmp = {0};

	if (flg == 1)
	{
		// #pragma omp parallel for num_threads(omp_get_max_threads()) //private(i,j,k)
		for (int i = 0; i < F; i++)
		{
			for (int j = 0; j < M; j++)
			{
				unsigned short u = 0;
				for (int k = 0; k < F; k++)
				{
					// tmp.z[j][i] ^= gf[mlt(fg[A.w[i][k]], fg[B.z[j][k]])];
					u += (A.x[i][k] * B.x[k][j]) % Pr;
					// tmp.d[j][i].fugo.b1 += A.d[i][k].fugo.b1 * B.d[j][k].fugo.b1%Pr;
					// tmp.d[j][i].fugo.b2 += A.d[i][k].fugo.b2 * B.d[j][k].fugo.b2%Pr;
				}
				tmp.x[i][j] = u % Pr;
				// printf("%d,",tmp.z[j][i]);
			}
			// printf("\n");
		}
		return tmp;
		// printf(" =====tmp.z\n");
		// exit(1);
	}
	if (flg == 2)
	{
		// #pragma omp parallel for num_threads(omp_get_max_threads()) // private(i,j,k)
		for (int i = 0; i < E * (K / 2 + 1); i++)
		{
			for (int j = 0; j < M; j++)
			{
				for (int k = 0; k < E * (K / 2 + 1); k++)
				{
					// tmp.w[j][i] ^= gf[mlt(fg[A.w[i][k]], fg[B.z[j][k]])];
					tmp.x[j][i] ^= A.x[i][k] & B.x[j][k];
				}
			}
		}
	}
	if (flg == 3)
	{
		short tmp2 = 0;
		short z[F][F] = {0};
		// 検算
		for (int i = 0; i < F; i++)
		{
			for (int j = 0; j < F; j++)
			{
				tmp2 = 0;
				for (int k = 0; k < F; k++)
					tmp2 += ((A.x[i][k]) * (B.x[k][j])) % Pr;
				z[i][j] = tmp2 % Pr;
			}
		}
		int count = 0;
		for (int i = 0; i < F; i++)
		{
			for (int j = 0; j < F; j++)
			{
				printf("=%d,", A.x[i][j] % Pr);
			}
			printf("\n");
		}
		// int count = 0;
		for (int i = 0; i < F; i++)
		{
			for (int j = 0; j < F; j++)
			{
				printf("^%d,", z[i][j] % Pr);
			}
			printf("\n");
			if (z[i][i] % Pr != 1)
			{
				// count++;
				printf("baka\n");
				exit(1);
			}
		}
	}
	/*
	//#pragma omp parallel for num_threads(omp_get_max_threads()) // private(i,j,k)
		for (int i= 0; i < F; i++)
		{
		  for (int j = 0; j < F; j++)
		  {
		  int u=0;
			for (int k = 0; k < F; k++)
			{
			  u += A.x[i][k]*B.x[k][j]%Pr;
			}
			tmp.x[i][j] =u%Pr;
		  }
		}
	  }
	*/
	return tmp;
}

/*
void mmul(MTX A, MTX B)
{
  int i, j, k;
  MTX tmp = {0};

  for (int i= 0; i < A.col; i++)
  {
	for (int j = 0; j < B.col; j++)
	{
	  for (int k = 0; k < A.row; k++)
	  {
		tmp.x[i][j] ^= gf[mlt(fg[A.x[i][k]], fg[B.x[k][j]])];
	  }
	  printf("%d,", tmp.x[i][j]);
	}
	printf("\n");
  }
  printf("\n");
}
*/

// Q-matrix
void matmul()
{
	int i, j, k, tmp[N][N] = {0};
	unsigned short x0[N]; //={1,2,3,4,5,6,7,0};
	// unsigned short x1[N]; //={2,3,1,6,5,7,0,4};
	// unsigned short x2[N] = {0};
	// unsigned short c[N][N] = {0};

	unsigned short a[N][N] = {0}; //{{0,3,0,0,},{0,0,4,0},{0,0,0,5},{6,0,0,0}};
	//{{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,0,0,0}};
	unsigned short inv_a[N][N] = {0};
	//{{0,0,0,1},{1,0,0,0},{0,1,0,0},{0,0,1,0}};
	//={{1,2,0,1},{1,1,2,0},{2,0,1,1},{1,2,1,1}}; //入力用の配列
	// unsigned short cc[N][N] = {0};
	for (int i = 0; i < N; i++)
		x0[i] = i;
	random_shuffle(x0, SIZE_OF_ARRAY(x0));

	printf("置換配列を表示\n");
	for (int i = 0; i < N; i++)
	{
		a[i][x0[i]] = 1; // rand()%N;
		printf("%d,", x0[i]);
	}
	printf("\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			inv_a[i][j] = a[j][i]; // gf[Inv(fg[a[j][i]])];//*inv_a[k][j];
		}
	}

	printf("Q1-置換行列を表示\n ={\n");
	for (int i = 0; i < N; i++)
	{
		printf("{");
		for (int j = 0; j < N; j++)
			printf("%3d,", a[j][i]);
		printf("},\n");
	}
	printf("};\n");

	printf("逆置換行列\n ={");
	for (int i = 0; i < N; i++)
	{
		printf("{");
		for (int j = 0; j < N; j++)
		{
			printf("%3d,", inv_a[j][i]);
		}
		printf("},\n");
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				tmp[i][j] ^= gf[mlt(fg[a[i][k]], fg[inv_a[k][j]])];
			}
		}
	}
	printf("};\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%d,", tmp[j][i]);
			// printf("%d,",inv_a[i][j]);
		}
		printf("\n");
	}
}

// 有限体の元の逆数
unsigned short
qinv(unsigned short a)
{

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

	for (int i = 0; i < K; i++)
	{
		p = a.x[i][i];

		for (int j = 0; j < (K + 1); j++)
		{
			a.x[i][j] = gf[mlt(fg[a.x[i][j]], qinv(p))];
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

/*
int main(){

	int i,j;
	double b[4],k=0;

	srand(clock());


	//g2();
 lab:
	printf("%d\n",gf[Inv(fg[3])]);
	printf("%d\n",gf[Inv(fg[4])]);
	printf("%d\n",gf[Inv(fg[5])]);
	printf("%d\n",gf[Inv(fg[6])]);
	matmul();
	matinv();
	printf("1=%d\n",gf[mlt(fg[3],fg[244])]);
	//if(det()!=1.0)
	//goto lab;

	return 0;
}
*/
