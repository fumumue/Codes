#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h> // SSE4.1命令を使用するためのヘッダファイル
#include <pthread.h>  // マルチスレッドを使用するためのヘッダファイル
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <immintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <mmintrin.h>

#define MATRIX_SIZE 4 // 行列のサイズ
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

// ガウスの消去法をSIMDで実行する関数
void gaussian_elimination_simd(int A[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int k = 0; k < MATRIX_SIZE; k++) {
        for (int i = k + 1; i < MATRIX_SIZE; i++) {
            __m128i factor = _mm_set1_epi32(A[i][k] / A[k][k]);
            for (int j = k; j < MATRIX_SIZE; j += 4) {
                __m128i row1 = _mm_loadu_si128((__m128i *)&A[i][j]);
                __m128i row2 = _mm_loadu_si128((__m128i *)&A[k][j]);
                row1 = _mm_sub_epi32(row1, _mm_mullo_epi32(factor, row2));
                _mm_storeu_si128((__m128i *)&A[i][j], row1);
            }
        }
    }
}

// スレッドごとに処理する行の範囲を指定する構造体
typedef struct {
    int (*matrix)[MATRIX_SIZE];
    int start_row;
    int end_row;
} ThreadData;

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

int main() {
    int A[MATRIX_SIZE][MATRIX_SIZE]; // 入力行列
    // 行列 A を初期化（例：ランダムな値を設定）
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i][j] = rand() % 100; // 0から99までのランダムな整数
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
        thread_data[i].start_row = i * rows_per_thread;
        thread_data[i].end_row = (i + 1) * rows_per_thread;
        pthread_create(&threads[i], NULL, gaussian_elimination_thread, &thread_data[i]);
    }

    // 各スレッドの終了を待つ
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("Result Matrix:\n");
    print_matrix(A);

    return 0;
}
