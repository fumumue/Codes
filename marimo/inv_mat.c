//
// (over finite field) Gauss-Jordan法による逆行列
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "global-p.h"
//#define F 8
//#define N 8
//#define Pr 1129
//#define E 11
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

#include <immintrin.h> // SSE4.1命令を使用するためのヘッダファイル
#include <pthread.h>  // マルチスレッドを使用するためのヘッダファイル
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#define MATRIX_SIZE 8 // 行列のサイズ
#define NUM_THREADS 4   // 使用するスレッドの数

// 行列を表示する関数
void print_matrix(int A[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%d ", A[i][j]);
        }
        printf("\n");
    }
}

/*
// ガウスの消去法をSIMDで実行する関数
void gaussian_elimination_simd(int A[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int k = 0; k < MATRIX_SIZE; k++) {
        for (int i = k + 1; i < MATRIX_SIZE; i++) {
            __m128i factor = _mm_set1_epi32(mlt(fg[A[i][k]] , Inv(A[k][k])));
            for (int j = k; j < MATRIX_SIZE; j += 4) {
                __m128i row1 = _mm_loadu_si128((__m128i *)&A[i][j]);
                __m128i row2 = _mm_loadu_si128((__m128i *)&A[k][j]);
                row1 = _mm_sub_epi32(row1, _mm_mullo_epi32(factor, row2));
                _mm_storeu_si128((__m128i *)&A[i][j], row1);
            }
        }
    }
}
*/
// スレッドごとに処理する行の範囲を指定する構造体
typedef struct {
    int (*matrix)[MATRIX_SIZE];
    int start_row;
    int end_row;
} ThreadData;

/*
// ガウスの消去法を実行する関数（マルチスレッド用）
void *gaussian_elimination_thread(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    int (*A)[MATRIX_SIZE] = data->matrix;
    int start_row = data->start_row;
    int end_row = data->end_row;

    // ガウスの消去法をSIMDで実行
    gaussian_elimination_simd(A);

    pthread_exit(NULL);
}


int maina_gausu(MTX A){
    //int A[MATRIX_SIZE][MATRIX_SIZE]; // 入力行列
    // 行列 A を初期化（例：ランダムな値を設定）
    
	for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A.x[i][j] = rand() % 100; // 0から99までのランダムな整数
        }
    }
	
    printf("Input Matrix:\n");
    print_matrix(A.x);

    // マルチスレッドのセットアップ
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];
    int rows_per_thread = MATRIX_SIZE / NUM_THREADS;

    // 各スレッドに行列 A の一部を処理させる
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].matrix = A.x;
        thread_data[i].start_row = i * rows_per_thread;
        thread_data[i].end_row = (i + 1) * rows_per_thread;
        pthread_create(&threads[i], NULL, gaussian_elimination_thread, &thread_data[i]);
    }

    // 各スレッドの終了を待つ
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("Result Matrix:\n");
    print_matrix(A.x);
}
*/

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


// 行列の次元を設定
#define Uz 8

// 行列の逆行列を計算する関数
void inverseMatrix(double A[Uz][Uz], double A_inv[Uz][Uz]) {
    int i, j, k;
    double temp;

    // 単位行列を初期化
    for (i = 0; i < Uz; i++) {
        for (j = 0; j < Uz; j++) {
            A_inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // ガウス・ジョルダン法による逆行列の計算
    for (k = 0; k < Uz; k++) {
        temp = A[k][k];
        for (j = 0; j < Uz; j++) {
            A[k][j] /= temp;
            A_inv[k][j] /= temp;
        }
        for (i = 0; i < Uz; i++) {
            if (i != k) {
                temp = A[i][k];
                for (j = 0; j < Uz; j++) {
                    A[i][j] -= A[k][j] * temp;
                    A_inv[i][j] -= A_inv[k][j] * temp;
                }
            }
        }
    }
}

// 行列の逆行列を計算する関数 (SSE4.1 SIMD版)
void matrix_inverse_simd(int A[MATRIX_SIZE][MATRIX_SIZE], double result[MATRIX_SIZE][MATRIX_SIZE]) {
    int pivot[MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; i++) {
        pivot[i] = -1;
    }

    __m128d one = _mm_set1_pd(1.0);

    for (int col = 0; col < MATRIX_SIZE; col++) {
        int pivot_row = -1;
        double max_value = 0.0;

        for (int row = 0; row < MATRIX_SIZE; row++) {
            if (pivot[row] != -1) continue;

            double val = fabs(A[row][col]);
            if (val > max_value) {
                max_value = val;
                pivot_row = row;
            }
        }

        if (pivot_row == -1) {
            fprintf(stderr, "Matrix is singular.\n");
            return;
        }

        pivot[pivot_row] = col;

        // Scale the pivot row
        double pivot_value = A[pivot_row][col];
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[pivot_row][j] /= pivot_value;
            result[pivot_row][j] = A[pivot_row][j];
        }

        // Eliminate non-zero entries below the pivot
        for (int row = 0; row < MATRIX_SIZE; row++) {
            if (row == pivot_row) continue;

            __m128d scale = _mm_set1_pd(A[row][col]);
            for (int j = 0; j < MATRIX_SIZE; j += 2) {
                __m128d row_pivot = _mm_loadu_pd(result[pivot_row] + j);
                __m128d scaled = _mm_mul_pd(scale, row_pivot);
                __m128d row_target = _mm_loadu_pd(result[row] + j);
                row_target = _mm_sub_pd(row_target, scaled);
                _mm_storeu_pd(result[row] + j, row_target);
            }
        }
    }
}

// inverse matrix
CTX matinv(CTX a, int n)
{

	// unsigned short a[F][F];     //={{1,2,0,1},{1,1,2,0},{2,0,1,1},{1,2,1,1}}; //入力用の配列
	unsigned short inv_a[N][N];	  // ここに逆行列が入る
	unsigned short buf;			  // 一時的なデータを蓄える
	unsigned short b[N][N] = {0}; //, dd[N][N] = {0};
	// int i, j, k, count;           // カウンタ
	//  MTX a={0};
	unsigned short c[N][N] = {0};
	CTX z = {0};
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

// 行列掛け算関数
void matrix_multiply(int A[MATRIX_SIZE][MATRIX_SIZE], int B[MATRIX_SIZE][MATRIX_SIZE], int *C, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            int sum = 0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                sum += A[i][k] * B[k][j];
            }
            C[i * MATRIX_SIZE + j] = sum%1129;
        }
    }
}

/*
// 行列掛け算関数（SIMD版、SSE4.1）
void matrix_multiply_simd(int A[MATRIX_SIZE][MATRIX_SIZE], int B[MATRIX_SIZE][MATRIX_SIZE], int *C, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j += 4) // 4つの要素を同時に処理するために4で割る
        {
            __m128i sum = _mm_setzero_si128(); // ゼロで初期化された128ビット整数ベクトルを作成

            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                __m128i a = _mm_set1_epi32(A[i][k]); // 行列Aの要素をベクトル化
                __m128i b = _mm_loadu_si128((__m128i *)&B[k][j]); // 行列Bの要素をベクトル化

                sum = _mm_add_epi32(sum, _mm_mullo_epi32(a, b))%N; // ベクトル同士の掛け算と足し算
            }

            // ベクトルの要素を合算して結果を格納
            int result[4];
            _mm_storeu_si128((__m128i *)result, sum);
            for (int k = 0; k < 4; k++)
            {
                C[i * MATRIX_SIZE + j + k] = result[k];
				C[i*MATRIX_SIZE+j+k]%=N;
            }
        }
    }
}
*/

#define SHM_KEY 1234

int task()
{
    int A[MATRIX_SIZE][MATRIX_SIZE];
    int B[MATRIX_SIZE][MATRIX_SIZE];
    //int C[MATRIX_SIZE][MATRIX_SIZE];

    // 行列 A と行列 B を初期化

    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            A[i][j] = (1129+i + j)%1129;
            B[i][j] = (1129+i - j)%1129;
			printf("%d,",A[i][j]);
        }
		printf("\n");
    }

    int num_processes = 8; // 4つのプロセスを使用して行列掛け算を高速化

    
    // 共有メモリを作成
    int shmid = shmget(SHM_KEY, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, IPC_CREAT | 0666);
    if (shmid == -1)
    {
        perror("shmget");
        exit(1);
    }

    int *C = (int *)shmat(shmid, NULL, 0);
    if (C == (int *)-1)
    {
        perror("shmat");
        exit(1);
    }
    
    // 各プロセスで一部の行を計算
    for (int i = 0; i < num_processes; i++)
    {
        pid_t pid = fork();

        if (pid == 0) // 子プロセスの処理
        {
            int start_row = i * (MATRIX_SIZE / num_processes);
            int end_row = (i + 1) * (MATRIX_SIZE / num_processes);

            matrix_multiply(A, B, C, start_row, end_row);
            //matrix_multiply_simd(A, B, C,start_row,end_row);

            exit(0); // 子プロセスを終了
        }
        else if (pid < 0)
        {
            // エラーハンドリング
            perror("fork");
            exit(1);
        }
    }

    // 親プロセスが子プロセスの終了を待つ
    for (int i = 0; i < num_processes; i++)
    {
        int status;
        wait(&status);
    }
    
    // 結果を表示（確認のため）
    printf("Result Matrix:\n");
    
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("%d ", C[i * MATRIX_SIZE + j]);
        }
        printf("\n");
    }
    
    // 共有メモリを解放
    if (shmdt(C) == -1)
    {
        perror("shmdt");
        exit(1);
    }
    if (shmctl(shmid, IPC_RMID, NULL) == -1)
    {
        perror("shmctl");
        exit(1);
    }

    return 0;
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
				tmp.x[i][j] = tmp2 % Pr;
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
				printf("^%d,", tmp.x[i][j] % Pr);
			}
			printf("\n");
			if (tmp.x[i][i] % Pr != 1)
			{
				// count++;
				//printf("baka\n");
				//exit(1);
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
