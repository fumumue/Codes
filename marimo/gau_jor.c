#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>

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

// 行列の逆行列を計算する関数（一部の行を処理）
void matrix_inverse_partial(double A[MATRIX_SIZE][MATRIX_SIZE], double result[MATRIX_SIZE][MATRIX_SIZE], int start_row, int end_row) {
    double temp;

    for (int k = 0; k < MATRIX_SIZE; k++) {
        if (fabs(A[k][k]) < 1e-15) {
            fprintf(stderr, "Matrix is singular.\n");
            exit(1); // 行列が特異行列の場合、逆行列は存在しない
        }

        if (k >= start_row && k < end_row) {
            // 部分行列の対角要素を1にする
            temp = A[k][k];
            for (int j = 0; j < MATRIX_SIZE; j++) {
                A[k][j] /= temp;
                result[k][j] /= temp;
            }
        }

        for (int i = start_row; i < end_row; i++) {
            if (i != k) {
                temp = A[i][k];
                for (int j = 0; j < MATRIX_SIZE; j++) {
                    A[i][j] -= A[k][j] * temp;
                    result[i][j] -= result[k][j] * temp;
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

    int num_processes = 4; // マルチプロセスの数
    int rows_per_process = MATRIX_SIZE / num_processes;
    int shmid;
    pid_t pid;

    // 共有メモリを作成
    shmid = shmget(IPC_PRIVATE, sizeof(double) * MATRIX_SIZE * MATRIX_SIZE, IPC_CREAT | 0666);
    if (shmid == -1) {
        perror("shmget");
        exit(1);
    }

    double *shared_A = (double *)shmat(shmid, NULL, 0);
    if (shared_A == (double *)-1) {
        perror("shmat");
        exit(1);
    }

    // 共有メモリに行列 A をコピー
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            shared_A[i * MATRIX_SIZE + j] = A[i][j];
        }
    }

    // マルチプロセスで並列計算
    for (int i = 0; i < num_processes; i++) {
        pid = fork();

        if (pid == 0) { // 子プロセス
            int start_row = i * rows_per_process;
            int end_row = (i + 1) * rows_per_process;

            // 逆行列の一部を計算
            matrix_inverse_partial(shared_A, result, start_row, end_row);

            exit(0); // 子プロセスを終了
        } else if (pid < 0) {
            perror("fork");
            exit(1);
        }
    }

    // 親プロセスが子プロセスの終了を待つ
    for (int i = 0; i < num_processes; i++) {
        int status;
        wait(&status);
    }

    // 共有メモリから逆行列を取得
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i][j] = shared_A[i * MATRIX_SIZE + j];
        }
    }

    // 共有メモリを解放
    if (shmdt(shared_A) == -1) {
        perror("shmdt");
        exit(1);
    }
    if (shmctl(shmid, IPC_RMID, NULL) == -1) {
        perror("shmctl");
        exit(1);
    }

    printf("Inverse Matrix:\n");
    print_matrix(result);

    return 0;
}

