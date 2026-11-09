#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>

#define MATRIX_SIZE 4   // 行列のサイズ (4x4行列を想定)
#define NUM_THREADS 4   // 使用するスレッドの数

// 行列を表示する関数
void print_matrix(double A[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
}

// 行列の検算を行う関数
void verify_matrix(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE]) {
    double C[MATRIX_SIZE][MATRIX_SIZE];

    // 行列 A と行列 B を掛け合わせた結果を行列 C に格納
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < MATRIX_SIZE; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    // 行列 C と単位行列を比較して検算を行う
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            if (i == j) {
                if (fabs(C[i][j] - 1.0) > 1e-6) {
                    printf("検算エラー: 単位行列ではありません\n");
                    return;
                }
            } else {
                if (fabs(C[i][j]) > 1e-6) {
                    printf("検算エラー: 単位行列ではありません\n");
                    return;
                }
            }
        }
    }

    printf("検算成功: 行列 A と B は正しく掛け合わされました。\n");
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

int main() {
    double A[MATRIX_SIZE][MATRIX_SIZE]; // 入力行列
    double result[MATRIX_SIZE][MATRIX_SIZE]; // 逆行列
    double X[4][4]={{1.0,0.0,0.0,0.0},{-1.0,0.0,0.0,1.0},{-1.0,0.0,1.0,0.0},{3.0,0.0,-2.0,-1.0}};
    int num_processes = NUM_THREADS; // スレッドの数

    // 行列 A を初期化（例：ランダムな値を設定）
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i][j] = j+1; //rand() % 100; // 0から99までのランダムな整数
        }
    }
    double D[4][4]={0};
    printf("Input Matrix:\n");
    print_matrix(A);
    for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
    for(int k=0;k<4;k++)
    D[i][j]+=A[i][k]*X[k][j];
    }
    }
    print_matrix(D);

    // マルチスレッドのセットアップ
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];
    int rows_per_thread = MATRIX_SIZE / NUM_THREADS;

    // 逆行列計算用のスレッドを作成
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].matrix = A;
        thread_data[i].result = result;
        thread_data[i].start_row = i * rows_per_thread;
        thread_data[i].end_row = (i + 1) * rows_per_thread;
        pthread_create(&threads[i], NULL, matrix_inverse_thread, &thread_data[i]);
    }

    // スレッドの終了を待つ
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("Inverse Matrix (Result):\n");
    print_matrix(result);

    // 逆行列と元の行列を掛けて単位行列を検算
    double identity[MATRIX_SIZE][MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            identity[i][j] = 0.0;
            for (int k = 0; k < MATRIX_SIZE; k++) {
                identity[i][j] += result[i][k] * A[k][j];
            }
        }
    }

    printf("Identity Matrix (Verification):\n");
    print_matrix(identity);

    // 検算を行う
    verify_matrix(identity, result);

    return 0;
}
