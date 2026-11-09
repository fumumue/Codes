#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <x86intrin.h> // SIMD命令を使用するためのヘッダファイル
#include <math.h>

#define MATRIX_SIZE 20 // 行列のサイズ
#define SHM_KEY 1234     // 共有メモリのキー

// 行列掛け算関数
void matrix_multiply2(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double *C, double start_row, double end_row)
{
double sum;
double dv=1129.0;
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            C[i*MATRIX_SIZE+j] = 0.0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                C[i*MATRIX_SIZE+j] += A[i][k] * B[k][j];
            }
            
            //while(sum>=dv){
            //sum-=dv;
            //}
            //C[i * MATRIX_SIZE + j] = sum;
            
            printf("%lf,",C[i*MATRIX_SIZE+j]);
        }
        printf("\n");
    }
}

// 行列掛け算関数
void matrix_multiply(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double *C, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            //w[i][j] = 0;
            double sum=0.0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                sum += A[i][k] * B[k][j];
            }
            C[i*MATRIX_SIZE+j]=sum;
        printf("%lf,",C[i*MATRIX_SIZE+j]);
        }
        //printf("\n");
    }
}


// 行列掛け算関数（SIMD版）
void matrix_multiply_simd(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE],double *C, double start_row, double end_row)
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

                sum = _mm_add_epi32(sum, _mm_mullo_epi32(a, b)); // ベクトル同士の掛け算と足し算

            }

            _mm_storeu_si128((__m128i *)&C[i*MATRIX_SIZE+j], sum); // 結果を行列Cに格納
        }
    }
}

#include <immintrin.h> // AVX2を使用するためのヘッダファイル


// 行列掛け算関数（SIMD版、AVX2）
void matrix_multiply_avx2(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double *C, double start_row,double end_row)
{
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            __m256i sum = _mm256_setzero_si256(); // ゼロで初期化された256ビット整数ベクトルを作成

            for (int k = 0; k < MATRIX_SIZE; k += 8) // 8つの要素を同時に処理するために8で割る
            {
                __m256i a = _mm256_loadu_si256((__m256i *)&A[i][k]); // 行列Aの要素をベクトル化
                __m256i b = _mm256_loadu_si256((__m256i *)&B[k][j]); // 行列Bの要素をベクトル化

                __m256i prod = _mm256_mullo_epi32(a, b); // ベクトル同士の掛け算
                sum = _mm256_add_epi32(sum, prod);       // ベクトル同士の足し算
            }

            // ベクトルの要素を合算して結果を格納
            double result[8];
            _mm256_storeu_si256((__m256i *)result, sum);
            for (int k = 0; k < 8; k++)
            {
                C[i*MATRIX_SIZE+j] += result[k];
            }
        }
    }
}

// 行列の逆行列を計算する関数
void inverseMatrix(double A[MATRIX_SIZE][MATRIX_SIZE], double A_inv[MATRIX_SIZE][MATRIX_SIZE]) {
    int i, j, k;
    double temp;

    // 単位行列を初期化
    for (i = 0; i < MATRIX_SIZE; i++) {
        for (j = 0; j < MATRIX_SIZE; j++) {
            A_inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // ガウス・ジョルダン法による逆行列の計算
    for (k = 0; k < MATRIX_SIZE; k++) {
        temp = A[k][k];
        for (j = 0; j < MATRIX_SIZE; j++) {
            A[k][j] /= temp;
            A_inv[k][j] /= temp;
        }
        for (i = 0; i < MATRIX_SIZE; i++) {
            if (i != k) {
                temp = A[i][k];
                for (j = 0; j < MATRIX_SIZE; j++) {
                    A[i][j] -= A[k][j] * temp;
                    A_inv[i][j] -= A_inv[k][j] * temp;
                }
            }
        }
    }
}


// 行列の逆行列を計算する関数 (SSE4.1 SIMD版)
void matrix_inverse_simd(double A[MATRIX_SIZE][MATRIX_SIZE], double result[MATRIX_SIZE][MATRIX_SIZE]) {
    double pivot[MATRIX_SIZE];
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


double main()
{
    double A[MATRIX_SIZE][MATRIX_SIZE];
    //double B[MATRIX_SIZE][MATRIX_SIZE];
    double AA[MATRIX_SIZE][MATRIX_SIZE];
    double A_inv[MATRIX_SIZE][MATRIX_SIZE];

    // 行列 A と行列 B を初期化
    double sum=0.0;
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            //sum = 1.0;
            //B[i][j] = (19.0+i - j);
            A[i][j]=A[i][j]=1.0+rand()%5; //1.0+i*N+j;//sum;
            AA[i][j]=A[i][j];
        }
    }

    int num_processes = 2; // 4つのプロセスを使用して行列掛け算を高速化

    
    // 共有メモリを作成
    int shmid = shmget(SHM_KEY, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, IPC_CREAT | 0666);
    if (shmid == -1)
    {
        perror("shmget");
        exit(1);
    }

    double *C = (double *)shmat(shmid, NULL, 0);
    if (C == (double *)-1)
    {
        perror("shmat");
        exit(1);
    }
    
    for(int i=0;i<MATRIX_SIZE;i++){
    for(int j=0;j<MATRIX_SIZE;j++)
    printf("%lf,",A[i][j]);
    printf("\n");
    }
    printf("\n");
    //matrix_inverse_simd(A,A_inv);
    inverseMatrix(A,A_inv);

    for(int i=0;i<MATRIX_SIZE;i++){
    for(int j=0;j<MATRIX_SIZE;j++)
    printf("%lf,",A_inv[i][j]);
    printf("\n");
    }
    //exit(1);
    // 各プロセスで一部の行を計算
    for (int i = 0; i < num_processes; i++)
    {
        pid_t pid = fork();

        if (pid == 0) // 子プロセスの処理
        {
            int start_row = i * (MATRIX_SIZE / num_processes);
            int end_row = (i + 1) * (MATRIX_SIZE / num_processes);
    
            //matrix_multiply(AA, A_inv, C, start_row, end_row);
            matrix_multiply_simd(AA, A_inv, C,start_row,end_row);


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
    printf("%d:\n",i);
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("%lf ", C[i * MATRIX_SIZE + j]);
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

