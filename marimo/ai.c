#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>

#define MATRIX_SIZE 16
#define SHM_KEY 1234

// 行列の逆行列を計算する関数
void inverseMatrix(double A[MATRIX_SIZE][MATRIX_SIZE], double A_inv[MATRIX_SIZE][MATRIX_SIZE], int start_row, int end_row) {
    int i, j, k;
    double temp;

    for (i = start_row; i < end_row; i++) {
        temp = A[i][i];
        for (j = 0; j < MATRIX_SIZE; j++) {
            A[i][j] /= temp;
            A_inv[i][j] /= temp;
        }
        for (k = 0; k < MATRIX_SIZE; k++) {
            if (k != i) {
                temp = A[k][i];
                for (j = 0; j < MATRIX_SIZE; j++) {
                    A[k][j] -= A[i][j] * temp;
                    A_inv[k][j] -= A_inv[i][j] * temp;
                }
            }
        }
    }
}

int main() {
    double A[MATRIX_SIZE][MATRIX_SIZE];
    double A_inv[MATRIX_SIZE][MATRIX_SIZE];
    double AA[MATRIX_SIZE][MATRIX_SIZE];

    // 行列 A を初期化
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i][j] = 1.0 + rand() % 7;
            printf("%lf,", A[i][j]);
            AA[i][j] = A[i][j];
        }
        printf("\n");
    }

    // 行列 A のコピーを作成
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            AA[i][j] = A[i][j];
        }
    }

    // マルチプロセスで行列掛け算を並列化
    int num_processes = 4;
    int rows_per_process = MATRIX_SIZE / num_processes;

    int shmid = shmget(SHM_KEY, sizeof(double) * MATRIX_SIZE * MATRIX_SIZE, IPC_CREAT | 0666);
    if (shmid == -1) {
        perror("shmget");
        exit(1);
    }

    double *shared_C = (double *)shmat(shmid, NULL, 0);
    if (shared_C == (double *)-1) {
        perror("shmat");
        exit(1);
    }

    // 各プロセスで一部の行を計算
    for (int i = 0; i < num_processes; i++) {
        pid_t pid = fork();

        if (pid == 0) {
            int start_row = i * rows_per_process;
            int end_row = (i + 1) * rows_per_process;
            if (i == num_processes - 1) {
                // 最後のプロセスは残りの行を計算
                end_row = MATRIX_SIZE;
            }

            inverseMatrix(AA, A_inv, start_row, end_row);

            // 結果を表示
            printf("Process %d: Rows %d to %d completed\n", i, start_row, end_row);

            exit(0);
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

    printf("Result Matrix:\n");
    for (int i = 0; i < MATRIX_SIZE; i++) {
        printf(":%d\n", i);
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%lf,", A_inv[i][j]);
        }
        printf("\n");
    }

    // 共有メモリを解放
    if (shmdt(shared_C) == -1) {
        perror("shmdt");
        exit(1);
    }
    if (shmctl(shmid, IPC_RMID, NULL) == -1) {
        perror("shmctl");
        exit(1);
    }

    return 0;
}
