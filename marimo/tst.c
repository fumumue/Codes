#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <immintrin.h> // SSE4.1命令を使用するためのヘッダファイル
#include <pthread.h>  // マルチスレッドを使用するためのヘッダファイル
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>

// 行列の次元を設定

#define MATRIX_SIZE 4

double w[MATRIX_SIZE][MATRIX_SIZE];

// 行列掛け算関数
void matrix_multiply_direct(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double e[MATRIX_SIZE][MATRIX_SIZE])
{
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            e[i][j]=0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                e[i][j] += A[i][k] * B[k][j];
            }
        }
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
        //printf("%lf,",C[i][j]);
        }
        //printf("\n");
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

// 行列掛け算関数（SIMD版、SSE4.1）
void matrix_multiply_simd(int A[MATRIX_SIZE][MATRIX_SIZE], int B[MATRIX_SIZE][MATRIX_SIZE], double *C, int start_row, int end_row)
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

            // ベクトルの要素を合算して結果を格納
            int result[4];
            _mm_storeu_si128((__m128i *)result, sum);
            for (int k = 0; k < 4; k++)
            {
                //C[i][j]=result[k];
                C[i * MATRIX_SIZE + j + k] = result[k];
            }
        }
    }
}


#define SHM_KEY 1234

int main() {
    double A[MATRIX_SIZE][MATRIX_SIZE];
    double A_inv[MATRIX_SIZE][MATRIX_SIZE];
    double AA[MATRIX_SIZE][MATRIX_SIZE];
    double e[MATRIX_SIZE][MATRIX_SIZE]={0};
    
    for(int i=0;i<MATRIX_SIZE;i++){
    for(int j=0;j<MATRIX_SIZE;j++)
    A[i][j]=1.0+rand()%1129; //1.0+i*N+j;
    }
    A[2][2]=10.0;
    for(int i=0;i<MATRIX_SIZE;i++){
    for(int j=0;j<MATRIX_SIZE;j++)
    AA[i][j]=A[i][j];
    }

    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%lf,", A[i][j]);
        }
        printf("\n");
    }
    
    inverseMatrix(A, A_inv);
    //matrix_inverse_simd(A,A_inv);

    printf("Inverse Matrix:\n");
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%f\t", A_inv[i][j]);
        }
        printf("\n");
    }
    //exit(1);

    int num_processes = 4; // 4つのプロセスを使用して行列掛け算を高速化

    
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

    // 各プロセスで一部の行を計算
    for (int i = 0; i < num_processes; i++)
    {
        pid_t pid = fork();

        if (pid == 0) // 子プロセスの処理
        {
            int start_row = i * (MATRIX_SIZE / num_processes);
            int end_row = (i + 1) * (MATRIX_SIZE / num_processes);

            matrix_multiply(AA, A_inv, C, start_row, end_row);

            // 結果を表示
            printf("Process %d: Rows %d to %d completed\n", i, start_row, end_row);

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
            printf("%lf ", C[i*MATRIX_SIZE+j]);
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
