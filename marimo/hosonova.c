#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>

#define MATRIX_SIZE 5

// 行列を表示する関数
void printMatrix(double matrix[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            printf("%.2f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int matmul_simd(double matrixA[MATRIX_SIZE][MATRIX_SIZE], double matrixB[MATRIX_SIZE][MATRIX_SIZE], double *resultMatrix, int start_row, int end_row)
{
    // double resultMatrix[MATRIX_SIZE][MATRIX_SIZE] = {0};

    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            __m128d sum = _mm_setzero_pd();

            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                __m128d a = _mm_set1_pd(matrixA[i][k]);
                __m128d b = _mm_loadu_pd(&matrixB[k][j]);
                sum = _mm_add_pd(sum, _mm_mul_pd(a, b));
            }

            _mm_storeu_pd(&resultMatrix[i * MATRIX_SIZE + j], sum);
        }
    }

    return 0;
}

// 行列の掛け算関数
void matrix_multiply(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double *C, int start_row, int end_row) {
    for (int i = start_row; i < end_row; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            double sum = 0.0;
            for (int k = 0; k < MATRIX_SIZE; k++) {
                sum += A[i][k] * B[k][j];
            }
            C[i * MATRIX_SIZE + j] = sum;
        }
    }
}

// 行列の単位行列を初期化する関数
void initializeIdentityMatrix(double identityMatrix[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            identityMatrix[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

// 行列の逆行列を計算する関数
void calculateInverseMatrix(double matrix[MATRIX_SIZE][MATRIX_SIZE], double inverseMatrix[MATRIX_SIZE][MATRIX_SIZE]) {
    // ガウス・ジョルダン法による逆行列の計算
    double tempMatrix[MATRIX_SIZE][2 * MATRIX_SIZE];
    
    // 入力行列と単位行列を結合
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            tempMatrix[i][j] = matrix[i][j];
            tempMatrix[i][j + MATRIX_SIZE] = (i == j) ? 1.0 : 0.0;
        }
    }

    // 前進消去
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        double pivot = tempMatrix[i][i];
        if (pivot == 0.0) {
            printf("Matrix is singular. Inverse does not exist.\n");
            return;
        }
        for (int j = 0; j < 2 * MATRIX_SIZE; ++j) {
            tempMatrix[i][j] /= pivot;
        }
        for (int k = 0; k < MATRIX_SIZE; ++k) {
            if (k != i) {
                double factor = tempMatrix[k][i];
                for (int j = 0; j < 2 * MATRIX_SIZE; ++j) {
                    tempMatrix[k][j] -= factor * tempMatrix[i][j];
                }
            }
        }
    }

    // 逆行列を抽出
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            inverseMatrix[i][j] = tempMatrix[i][j + MATRIX_SIZE];
        }
    }
}

void matrix_inverse_simd(double A[MATRIX_SIZE][MATRIX_SIZE], double result[MATRIX_SIZE][MATRIX_SIZE])
{
    double pivot[MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        pivot[i] = -1.0;
    }

    // 1. データ型を__m128dに変更し、SIMDレジスタを使用する
    __m128d one = _mm_set1_pd(1.0);

    for (int col = 0; col < MATRIX_SIZE; col++)
    {
        int pivot_row = -1;
        double max_value = 0.0;

        for (int row = 0; row < MATRIX_SIZE; row++)
        {
            if (pivot[row] != -1.0)
                continue;

            double val = fabs(A[row][col]);
            if (val > max_value)
            {
                max_value = val;
                pivot_row = row;
            }
        }

        if (pivot_row == -1)
        {
            fprintf(stderr, "Matrix is singular.\n");
            return;
        }

        pivot[pivot_row] = col;

        // Scale the pivot row
        double pivot_value = A[pivot_row][col];
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            A[pivot_row][j] /= pivot_value;
            result[pivot_row][j] = A[pivot_row][j];
        }

        // Eliminate non-zero entries below the pivot
        for (int row = 0; row < MATRIX_SIZE; row++)
        {
            if (row == pivot_row)
                continue;

            // 2. SIMDを使用して計算
            __m128d scale = _mm_set1_pd(A[row][col]);
            for (int j = 0; j < MATRIX_SIZE; j += 2)
            {
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
    double inputMatrix[MATRIX_SIZE][MATRIX_SIZE];
    double inverseMatrix[MATRIX_SIZE][MATRIX_SIZE];
    double C[MATRIX_SIZE][MATRIX_SIZE];

    srand(clock());
    // 5次元の行列を初期化
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            inputMatrix[i][j] = rand() % 10; // 0から9のランダムな値で初期化
        }
    }

    printf("Input Matrix:\n");
    printMatrix(inputMatrix);

    // 逆行列を計算
    calculateInverseMatrix(inputMatrix, inverseMatrix);
    //matrix_inverse_simd(inputMatrix,inverseMatrix);
    printf("Inverse Matrix:\n");
    printMatrix(inverseMatrix);

    // 行列の掛け算
    //matrix_multiply(inputMatrix, inverseMatrix, (double *)C, 0, MATRIX_SIZE);
    matmul_simd(inputMatrix, inverseMatrix, (double *)C, 0, MATRIX_SIZE);

    // 結果を出力
    printf("Original Matrix multiplied by Inverse Matrix:\n");
    printMatrix(C);

    return 0;
}

