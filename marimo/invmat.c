#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h> // SSE4.1命令を使用するためのヘッダファイル
#include <pthread.h>  // マルチスレッドを使用するためのヘッダファイル
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>

#define MATRIX_SIZE 4   // 行列のサイズ (4x4行列を想定)
#define NUM_THREADS 4   // 使用するスレッドの数
#define SHM_KEY 1234     // 共有メモリのキー

// 行列を表示する関数
void print_matrix(double A[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
}

// 行列の逆行列を計算する関数 (SSE4.1 SIMD版)
void matrix_inverse_simd(double A[MATRIX_SIZE][MATRIX_SIZE], double result[MATRIX_SIZE][MATRIX_SIZE]) {
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

// スレッドごとに処理する行の範囲を指定する構造体
typedef struct {
    double (*matrix)[MATRIX_SIZE];
    double (*result)[MATRIX_SIZE];
    int start_row;
    int end_row;
} ThreadData;

// 逆行列を計算する関数（マルチスレッド用）
void *matrix_inverse_thread(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    double (*A)[MATRIX_SIZE] = data->matrix;
    double (*result)[MATRIX_SIZE] = data->result;
    int start_row = data->start_row;
    int end_row = data->end_row;

    // 逆行列を計算 (SIMD版)
    for (int row = start_row; row < end_row; row++) {
        matrix_inverse_simd(A, result);
    }

    pthread_exit(NULL);
}

// 行列掛け算関数（SIMD版）
void matrix_multiply_simd(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE],double *C, int start_row, int end_row)
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

// 行列掛け算関数
void matrix_multiply(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double *C, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            C[i*MATRIX_SIZE+j] = 0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                C[i*MATRIX_SIZE+j] += A[i][k] * B[k][j];
            }
        printf("%lf,",C[i*MATRIX_SIZE+j]);
        }
        printf("\n");
    }
}


int main() {
    double A[MATRIX_SIZE][MATRIX_SIZE]; // 入力行列
    double result[MATRIX_SIZE][MATRIX_SIZE]; // 逆行列
    double AA[MATRIX_SIZE][MATRIX_SIZE];
    // 行列 A を初期化（例：ランダムな値を設定）
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i][j] = rand() % 100; // 0から99までのランダムな整数
            AA[i][j]=A[i][j];
        }
    }

    printf("Input Matrix:\n");
    print_matrix(A);

    // マルチスレッドのセットアップ
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];
    int rows_per_thread = MATRIX_SIZE / NUM_THREADS;

    // 各スレッドに行列 A の一部を処理させる
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].matrix = A;
        thread_data[i].result =
        thread_data[i].result = result;
        thread_data[i].start_row = i * rows_per_thread;
        thread_data[i].end_row = (i == NUM_THREADS - 1) ? MATRIX_SIZE : (i + 1) * rows_per_thread;

        // スレッドを生成
        int ret = pthread_create(&threads[i], NULL, matrix_inverse_thread, &thread_data[i]);
        if (ret != 0) {
            fprintf(stderr, "Error creating thread %d\n", i);
            return 1;
        }
    }
    //double C[MATRIX_SIZE][MATRIX_SIZE];
    // スレッドの終了を待つ
    for (int i = 0; i < NUM_THREADS; i++) {
        int ret = pthread_join(threads[i], NULL);
        if (ret != 0) {
            fprintf(stderr, "Error joining thread %d\n", i);
            return 1;
        }
    }

    printf("Inverse Matrix:\n");
    print_matrix(result);


    int num_processes = 8; // 4つのプロセスを使用して行列掛け算を高速化
    
    // 共有メモリを作成
    int shmid = shmget(SHM_KEY, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, IPC_CREAT | 0666);
    if (shmid == -1)
    {
        perror("shmget");
        exit(1);
    }

    double *C = (int *)shmat(shmid, NULL, 0);
    double B[MATRIX_SIZE][MATRIX_SIZE];

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
            //matrix_multiply_avx2(A, B, C, start_row, end_row);
            //matrix_multiply_simd(A,result,C,start_row,end_row);

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
