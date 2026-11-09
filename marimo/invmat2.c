#include <stdio.h>
#include <immintrin.h> // SIMD命令を使用するためのヘッダファイル

#define MATRIX_SIZE 4

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

            double val = (A[row][col]);
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
int main() {
    double A[MATRIX_SIZE][MATRIX_SIZE];        // 入力行列
    double result[MATRIX_SIZE][MATRIX_SIZE];    // 逆行列

    // 行列 A を初期化（例：単位行列を設定）
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            if (i == j) {
                A[i][j] = rand()%100;
            } else {
                A[i][j] = rand()%100;
            }
        }
    }

    printf("Input Matrix A:\n");
    print_matrix(A);

    // 逆行列を計算 (SIMD版)
    matrix_inverse_simd(A, result);

    printf("Inverse Matrix B:\n");
    print_matrix(result);

    // 行列 A と逆行列 B を掛けて単位行列を検証
    double identity[MATRIX_SIZE][MATRIX_SIZE] = {0};
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            for (int k = 0; k < MATRIX_SIZE; k++) {
                identity[i][j] += A[i][k] * result[k][j];
            }
        }
    }

    printf("Identity Matrix (A * B):\n");
    print_matrix(identity);

    return 0;
}
