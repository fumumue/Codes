#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

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

// 行列の逆行列を計算する関数
void matrix_inverse(double A[MATRIX_SIZE][MATRIX_SIZE], double A_inv[MATRIX_SIZE][MATRIX_SIZE]) {
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

int main() {
    double A[MATRIX_SIZE][MATRIX_SIZE]; // 入力行列
    double result[MATRIX_SIZE][MATRIX_SIZE]; // 逆行列

    // 行列 A を初期化（例：ランダムな値を設定）
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i][j] = (double)(rand() % 100); // 0から99までのランダムな整数
        }
    }

    printf("Input Matrix:\n");
    print_matrix(A);

    // マルチプロセスのセットアップ
    int num_processes = 4; // 使用するプロセスの数
    pid_t child_pid;
    int status;

    for (int i = 0; i < num_processes; i++) {
        child_pid = fork();

        if (child_pid == 0) { // 子プロセスの処理
            int start_row = i * (MATRIX_SIZE / num_processes);
            int end_row = (i + 1) * (MATRIX_SIZE / num_processes);

            // 逆行列の計算を実行
            //matrix_inverse(A, result);
            inverseMatrix(A,result);

            // 逆行列を親プロセスに送信
            exit(0);
        } else if (child_pid < 0) {
            // エラーハンドリング
            perror("fork");
            exit(1);
        }
    }

    // 親プロセスが子プロセスの終了を待つ
    for (int i = 0; i < num_processes; i++) {
        wait(&status);
    }

    // 逆行列を表示
    printf("Inverse Matrix:\n");
    print_matrix(result);

    return 0;
}
