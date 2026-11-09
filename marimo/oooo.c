#include <stdio.h>
#include <pthread.h>

#define MATRIX_SIZE 5

// 行列を格納する構造体
typedef struct {
    double data[MATRIX_SIZE][MATRIX_SIZE];
} Matrix;

// スレッド用のデータ構造体
typedef struct {
    Matrix *inputMatrix;
    Matrix *outputMatrix;
} ThreadData;

// 行列を表示する関数
void printMatrix(Matrix *matrix) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            printf("%.2f ", matrix->data[i][j]);
        }
        printf("\n");
    }
}

// 行列式を計算する関数
double determinant(Matrix *matrix, int n) {
    double det = 0;
    Matrix submatrix;
    if (n == 2) {
        return matrix->data[0][0] * matrix->data[1][1] - matrix->data[0][1] * matrix->data[1][0];
    } else {
        for (int k = 0; k < n; k++) {
            int sub_i = 0;
            for (int i = 1; i < n; i++) {
                int sub_j = 0;
                for (int j = 0; j < n; j++) {
                    if (j == k)
                        continue;
                    submatrix.data[sub_i][sub_j] = matrix->data[i][j];
                    sub_j++;
                }
                sub_i++;
            }
            det = det + (k % 2 == 0 ? 1 : -1) * matrix->data[0][k] * determinant(&submatrix, n - 1);
        }
    }
    return det;
}

// 逆行列を計算する関数
void *calculateInverse(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    Matrix *inputMatrix = data->inputMatrix;
    Matrix *outputMatrix = data->outputMatrix;

    double det = determinant(inputMatrix, MATRIX_SIZE);
    if (det == 0) {
        printf("Matrix is singular, inverse does not exist.\n");
        return NULL;
    }

    double invDet = 1.0 / det;
    Matrix adjoint;

    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            int sub_i = 0;
            for (int m = 0; m < MATRIX_SIZE; m++) {
                if (m == i)
                    continue;
                int sub_j = 0;
                for (int n = 0; n < MATRIX_SIZE; n++) {
                    if (n == j)
                        continue;
                    adjoint.data[sub_i][sub_j] = inputMatrix->data[m][n];
                    sub_j++;
                }
                sub_i++;
            }
            outputMatrix->data[i][j] = (i + j) % 2 == 0 ? 1 : -1 * determinant(&adjoint, MATRIX_SIZE - 1);
        }
    }

    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = i; j < MATRIX_SIZE; j++) {
            double temp = outputMatrix->data[i][j];
            outputMatrix->data[i][j] = outputMatrix->data[j][i];
            outputMatrix->data[j][i] = temp;
        }
    }

    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            outputMatrix->data[i][j] *= invDet;
        }
    }

    return NULL;
}

// 行列の掛け算を行う関数
void multiplyMatrices(Matrix *matrix1, Matrix *matrix2, Matrix *result) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            result->data[i][j] = 0;
            for (int k = 0; k < MATRIX_SIZE; ++k) {
                result->data[i][j] += matrix1->data[i][k] * matrix2->data[k][j];
            }
        }
    }
}

// 単位行列かどうかをチェックする関数
int isIdentityMatrix(Matrix *matrix) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            if (i == j && matrix->data[i][j] != 1.0) {
                return 0;
            } else if (i != j && matrix->data[i][j] != 0.0) {
                return 0;
            }
        }
    }
    return 1;
}

int main() {
    Matrix inputMatrix, inverseMatrix, resultMatrix;
    pthread_t thread;
    ThreadData data;

    // 入力行列の初期化
    inputMatrix.data[0][0] = 4.00;
    inputMatrix.data[0][1] = 6.00;
    inputMatrix.data[0][2] = 5.00;
    inputMatrix.data[0][3] = 5.00;
    inputMatrix.data[0][4] = 4.00;
    inputMatrix.data[1][0] = 4.00;
    inputMatrix.data[1][1] = 9.00;
    inputMatrix.data[1][2] = 15.00;
    inputMatrix.data[1][3] = 13.00;
    inputMatrix.data[1][4] = 10.00;
    inputMatrix.data[2][0] = 7.00;
    inputMatrix.data[2][1] = 11.00;
    inputMatrix.data[2][2] = 1.00;
    inputMatrix.data[2][3] = 5.00;
    inputMatrix.data[2][4] = 10.00;
    inputMatrix.data[3][0] = 8.00;
    inputMatrix.data[3][1] = 8.00;
    inputMatrix.data[3][2] = 4.00;
    inputMatrix.data[3][3] = 9.00;
    inputMatrix.data[3][4] = 3.00;
    inputMatrix.data[4][0] = 14.00;
    inputMatrix.data[4][1] = 0.00;
    inputMatrix.data[4][2] = 11.00;
    inputMatrix.data[4][3] = 4.00;
    inputMatrix.data[4][4] = 10.00;

    printf("Input Matrix:\n");
    printMatrix(&inputMatrix);

    data.inputMatrix = &inputMatrix;
    data.outputMatrix = &inverseMatrix;

    // 逆行列の計算
    pthread_create(&thread, NULL, calculateInverse, (void *)&data);
    pthread_join(thread, NULL);

    // 逆行列の出力
    printf("Inverse Matrix:\n");
    printMatrix(&inverseMatrix);

    // もとの行列と逆行列の掛け算
    multiplyMatrices(&inputMatrix, &inverseMatrix, &resultMatrix);

    // 結果の出力
    printf("Original Matrix multiplied by Inverse Matrix:\n");
    printMatrix(&resultMatrix);

    // 単位行列かどうかのチェック
    printf("Is it an Identity Matrix? %s\n", isIdentityMatrix(&resultMatrix) ? "Yes" : "No");

    return 0;
}

