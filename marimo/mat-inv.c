#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>

#define MATRIX_SIZE 8
#define SHM_KEY 1234

// 行列の逆行列を計算する関数
void inverseMatrix(double A[MATRIX_SIZE][MATRIX_SIZE], double A_inv[MATRIX_SIZE][MATRIX_SIZE],double *C, int start_row, int end_row) {
    int i, j, k;
    double temp;

    // ガウス・ジョルダン法による逆行列の計算
    for (k = 0; k < MATRIX_SIZE; k++) {
        temp = A[k][k];
        for (j = 0; j < MATRIX_SIZE; j++) {
            A[k][j] /= temp;
            A_inv[k][j] /= temp;
        }
        for (i = start_row; i < end_row; i++) {
            if (i != k) {
                temp = A[i][k];
                for (j = 0; j < MATRIX_SIZE; j++) {
                    A[i][j] -= A[k][j] * temp;
                    A_inv[i][j] -= A_inv[k][j] * temp;
                }
            }
        }
    }
for(int i=0;i<MATRIX_SIZE;i++){
for(int j=0;j<MATRIX_SIZE;j++)
C[i*MATRIX_SIZE+j]=A_inv[i][j];
}
}

int main() {
    double A[MATRIX_SIZE][MATRIX_SIZE];
    double A_inv[MATRIX_SIZE][MATRIX_SIZE];
    
    // 行列 A を初期化
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i][j] = 1.0 + rand() % 7;
        }
    }

    // 共有メモリを作成
    int shmid = shmget(SHM_KEY, sizeof(double) * MATRIX_SIZE * MATRIX_SIZE, IPC_CREAT | 0666);
    if (shmid == -1) {
        perror("shmget");
        exit(1);
    }

    double *shared_A = (double *)shmat(shmid, NULL, 0);
    if (shared_A == (double *)-1) {
        perror("shmat");
        exit(1);
    }

    // 子プロセスの数を設定
    int num_processes = 4;
    int rows_per_process = MATRIX_SIZE / num_processes;

    // 各プロセスで一部の行を計算
    for (int i = 0; i < num_processes; i++) {
        pid_t pid = fork();

        if (pid == 0) {
            int start_row = i * rows_per_process;
            int end_row = (i + 1) * rows_per_process;

            inverseMatrix(A, A_inv,shared_A, start_row, end_row);

            // 子プロセスを終了
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

    // 結果を表示
    printf("Inverse Matrix:\n");
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%f\t", shared_A[i*MATRIX_SIZE+j]);
        }
        printf("\n");
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

    return 0;
}
