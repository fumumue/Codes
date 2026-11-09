#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define MATRIX_SIZE 5
#define NUM_THREADS 5

// 行列を表示する関数
void printMatrix(double matrix[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            printf("%.2f ", matrix[i][j]);
        }
        printf("\n");
    }
}

// スレッドのデータ構造
struct ThreadData {
    double (*A)[MATRIX_SIZE]; // 行列A
    double (*B)[MATRIX_SIZE]; // 行列B
    double (*C)[MATRIX_SIZE]; // 結果行列C
    int start_row;            // 計算を始める行
    int end_row;              // 計算を終える行
};

// 行列の掛け算関数（スレッド版）
void *matrix_multiply_thread(void *threadArg) {
    struct ThreadData *myData;
    myData = (struct ThreadData *) threadArg;

    for (int i = myData->start_row; i < myData->end_row; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            double sum = 0.0;
            for (int k = 0; k < MATRIX_SIZE; k++) {
                sum += myData->A[i][k] * myData->B[k][j];
            }
            myData->C[i][j] = sum;
        }
    }

    pthread_exit(NULL);
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

    printf("Inverse Matrix:\n");
    printMatrix(inverseMatrix);

    // スレッドデータ
    pthread_t threads[NUM_THREADS];
    struct ThreadData threadData[NUM_THREADS];

    // スレッドの生成と行列の掛け算
    int rows_per_thread = MATRIX_SIZE / NUM_THREADS;
    for (int i = 0; i < NUM_THREADS; i++) {
        threadData[i].A = inputMatrix;
        threadData[i].B = inverseMatrix;
        threadData[i].C = C;
        threadData[i].start_row = i * rows_per_thread;
        threadData[i].end_row = (i + 1) * rows_per_thread;
        if (i == NUM_THREADS - 1) {
            threadData[i].end_row = MATRIX_SIZE;
        }
        pthread_create(&threads[i], NULL, matrix_multiply_thread, (void *)&threadData[i]);
    }

    // スレッドの終了を待つ
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    // 結果を出力
    printf("Original Matrix multiplied by Inverse Matrix:\n");
    printMatrix(C);

    return 0;
}

