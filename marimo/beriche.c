#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>

#define MATRIX_SIZE 3
#define EPSILON 1e-10 // 0に非常に近い値を表すためのイプシロン

// 行列を格納する構造体
typedef struct {
    double data[MATRIX_SIZE][MATRIX_SIZE];
    int row;
} Matrix;

// スレッド用のデータ構造体
typedef struct {
    Matrix *matrix;
    int start_row;
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

// 逆行列を計算する関数
void *calculateInverse(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    Matrix *matrix = data->matrix;
    int start_row = data->start_row;

    // ガウス・ジョルダン法を用いて逆行列を計算する
    for (int k = 0; k < MATRIX_SIZE; ++k) {
        // ピボット要素のスレッドごとの処理
        if (start_row <= k && k < start_row + matrix->row) {
            int row = k - start_row;
            double pivot = matrix->data[row][k];
            // ピボット要素が0に近い場合は特異行列として扱う
            if (fabs(pivot) < EPSILON) {
                fprintf(stderr, "Matrix is singular, inverse does not exist.\n");
                pthread_exit(NULL);
            }
            // ピボット要素を1にする
            for (int j = 0; j < MATRIX_SIZE; ++j) {
                matrix->data[row][j] /= pivot;
            }
            // 同じ行以外の要素を0にする
            for (int i = 0; i < MATRIX_SIZE; ++i) {
                if (i != row) {
                    double factor = matrix->data[i][k];
                    for (int j = 0; j < MATRIX_SIZE; ++j) {
                        matrix->data[i][j] -= factor * matrix->data[row][j];
                    }
                }
            }
        }
    }

    pthread_exit(NULL);
}

int main() {
    Matrix inputMatrix;
    // 入力行列の初期化
    for(int i=0;i<MATRIX_SIZE;i++){
    for(int j=0;j<MATRIX_SIZE;j++)
    inputMatrix.data[i][j] = rand()%17;
    }
    inputMatrix.row = MATRIX_SIZE;

    printf("Input Matrix:\n");
    printMatrix(&inputMatrix);

    // スレッド数
    const int num_threads = MATRIX_SIZE;

    // スレッドの宣言
    pthread_t threads[num_threads];
    ThreadData thread_data[num_threads];

    // スレッドの作成
    for (int i = 0; i < num_threads; ++i) {
        thread_data[i].matrix = &inputMatrix;
        thread_data[i].start_row = i;
        pthread_create(&threads[i], NULL, calculateInverse, (void *)&thread_data[i]);
    }

    // スレッドの終了を待つ
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }

    // 逆行列の出力
    printf("Inverse Matrix:\n");
    printMatrix(&inputMatrix);

    return 0;
}

