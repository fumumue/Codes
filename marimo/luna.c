#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>

#define MATRIX_SIZE 4 // 行列のサイズ (4x4行列を想定)

// 行列を表示する関数
void print_matrix(double *A) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%lf ", A[i*MATRIX_SIZE+j]);
        }
        printf("\n");
    }
}


#define MATRIX_SIZE 4

// 行列の逆行列を計算する関数
int matrix_inverse(double *A, double result[MATRIX_SIZE][MATRIX_SIZE]) {
    double temp;

    // 単位行列を生成
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            result[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // ガウス・ジョルダン法による逆行列計算
    for (int k = 0; k < MATRIX_SIZE; k++) {
        if (fabs(A[k*MATRIX_SIZE+k]) < 1e-15) {
            fprintf(stderr, "Matrix is singular.\n");
            return -1; // 行列が特異行列の場合、逆行列は存在しない
        }

        // 対角要素を1にする
        temp = A[k*MATRIX_SIZE+k];
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[k*MATRIX_SIZE+j] /= temp;
            result[k][j] /= temp;
        }

        for (int i = 0; i < MATRIX_SIZE; i++) {
            if (i != k) {
                temp = A[i*MATRIX_SIZE+k];
                for (int j = 0; j < MATRIX_SIZE; j++) {
                    A[i*MATRIX_SIZE+j] -= A[k*MATRIX_SIZE+j] * temp;
                    result[i][j] -= result[k][j] * temp;
                }
            }
        }
    }

    return 0;
}


// 行列掛け算関数
void matrix_multiply(double *A, double B[MATRIX_SIZE][MATRIX_SIZE], double *C, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            //w[i][j] = 0;
            double sum=0.0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                sum += A[i*MATRIX_SIZE+j] * B[k][j];
            }
            C[i*MATRIX_SIZE+j]=sum;
        //printf("%lf,",C[i][j]);
        }
        //printf("\n");
    }
}


int main() {
    double A[MATRIX_SIZE][MATRIX_SIZE]; // 入力行列
    double result[MATRIX_SIZE][MATRIX_SIZE]; // 逆行列
    double D[MATRIX_SIZE][MATRIX_SIZE];

    // 行列 A を初期化（例：ランダムな値を設定）
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i][j] = (double)(rand() % 100); // 0から99までのランダムな整数
        }
    }


    printf("Input Matrix:\n");
    print_matrix(A);

    // 逆行列を計算
    if (matrix_inverse(A, result) == 0) {
        printf("Inverse Matrix:\n");
        print_matrix(result);
    }
exit(1);

    printf("Input Matrix:\n");
    print_matrix(*A);

    // 共有メモリを作成
    int shmid = shmget(IPC_PRIVATE, sizeof(double) * MATRIX_SIZE * MATRIX_SIZE, IPC_CREAT | 0666);
    if (shmid == -1) {
        perror("shmget");
        exit(1);
    }

    double *shared_matrix = (double *)shmat(shmid, NULL, 0);
    if (shared_matrix == (double *)-1) {
        perror("shmat");
        exit(1);
    }

    // マルチプロセスのセットアップ
    int num_processes = 4; // 使用するプロセスの数
    pid_t child_pid;
    int status;

    for (int i = 0; i < num_processes; i++) {
        child_pid = fork();

        if (child_pid == 0) { // 子プロセスの処理
            int start_row = i * (MATRIX_SIZE / num_processes);
            int end_row = (i + 1) * (MATRIX_SIZE / num_processes);

            // 共有メモリに行列をコピー
            for (int row = start_row; row < end_row; row++) {
                for (int col = 0; col < MATRIX_SIZE; col++) {
                    shared_matrix[row * MATRIX_SIZE + col] = A[row][col];
                }
            }

            // 逆行列の計算を実行
            matrix_inverse(shared_matrix, result);
            matrix_multiply(shared_matrix,result,*D,start_row,end_row);

            // 逆行列を共有メモリに格納
            for (int row = start_row; row < end_row; row++) {
                for (int col = 0; col < MATRIX_SIZE; col++) {
                    A[row][col] = shared_matrix[row * MATRIX_SIZE + col];
                }
            }

            shmdt(shared_matrix);
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
    print_matrix(D);

    // 共有メモリを解放
    if (shmctl(shmid, IPC_RMID, NULL) == -1) {
        perror("shmctl");
        exit(1);
    }

    return 0;
}
