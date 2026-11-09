#include <stdio.h>
#include <immintrin.h> // SIMD命令を使用するためのヘッダ

#define MATRIX_SIZE 4

// 行列を表示する関数
void printMatrix(float matrix[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            printf("%.2f ", matrix[i][j]);
        }
        printf("\n");
    }
}

// 行列の掛け算関数（SIMD版）
void matrix_multiply_simd(float A[MATRIX_SIZE][MATRIX_SIZE], float B[MATRIX_SIZE][MATRIX_SIZE], float C[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            __m128 sum = _mm_setzero_ps(); // 4つの要素を持つベクトルを0で初期化
            for (int k = 0; k < MATRIX_SIZE; ++k) {
                __m128 a = _mm_load_ps(&A[i][k]); // Aのi行k列のベクトルをロード
                __m128 b = _mm_set1_ps(B[k][j]); // Bのk行j列の要素をすべて同じ値に設定
                sum = _mm_add_ps(sum, _mm_mul_ps(a, b)); // ベクトル同士の掛け算と加算
            }
            _mm_store_ps(&C[i][j], sum); // 結果をCにストア
        }
    }
}

int main() {
    // 行列のデータ
    float A[MATRIX_SIZE][MATRIX_SIZE] __attribute__((aligned(16))) = {
        {1.0, 2.0, 3.0, 4.0},
        {5.0, 6.0, 7.0, 8.0},
        {9.0, 10.0, 11.0, 12.0},
        {13.0, 14.0, 15.0, 16.0}
    };

    float B[MATRIX_SIZE][MATRIX_SIZE] __attribute__((aligned(16))) = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0}
    };

    float C[MATRIX_SIZE][MATRIX_SIZE] __attribute__((aligned(16)));

    // 行列の掛け算（SIMD版）
    matrix_multiply_simd(A, B, C);

    // 結果を出力
    printf("Result Matrix:\n");
    printMatrix(C);

    return 0;
}

