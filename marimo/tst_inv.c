#include <stdio.h>

#include <immintrin.h> // SSE4.1命令を使用するためのヘッダファイル
#include <pthread.h>  // マルチスレッドを使用するためのヘッダファイル
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>

// 行列の次元を設定
#define N 4
#define MATRIX_SIZE 4

// 行列掛け算関数
void matrix_multiply_direct(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double *C)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            int sum = 0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                sum += A[i][k] * B[k][j];
            }
            C[i * MATRIX_SIZE + j] = sum%1129;
        }
    }
}

// 行列掛け算関数
void matrix_multiply(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double *C, double start_row, double end_row)
{
double sum;
double dv=1129.0;
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            C[i*MATRIX_SIZE+j] = 0.0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                C[i*MATRIX_SIZE+j] += A[i][k] * B[k][j];
            }
            /*
            while(sum>=dv){
            sum-=dv;
            }
            C[i * MATRIX_SIZE + j] = sum;
            */
            printf("%lf,",C[i*MATRIX_SIZE+j]);
        }
        printf("\n");
    }
}


// 行列の逆行列を計算する関数
void inverseMatrix(double A[N][N], double A_inv[N][N]) {
    int i, j, k;
    double temp;

    // 単位行列を初期化
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A_inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // ガウス・ジョルダン法による逆行列の計算
    for (k = 0; k < N; k++) {
        temp = A[k][k];
        for (j = 0; j < N; j++) {
            A[k][j] /= temp;
            A_inv[k][j] /= temp;
        }
        for (i = 0; i < N; i++) {
            if (i != k) {
                temp = A[i][k];
                for (j = 0; j < N; j++) {
                    A[i][j] -= A[k][j] * temp;
                    A_inv[i][j] -= A_inv[k][j] * temp;
                }
            }
        }
    }
}

#define SHM_KEY 1234

int main() {
double A[N][N];
double B[N][N]={0};
double A_inv[N][N];
double C[N][N];
double AA[N][N];

for(int i=0;i<N;i++){
for(int j=0;j<N;j++){
A[i][j]=rand()%19*1.0;
AA[i][j]=A[i][j];
}
}

    printf("Matrix:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f\t", A[i][j]);
        }
        printf("\n");
    }

    inverseMatrix(A, A_inv);

    printf("Inverse Matrix:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f\t", A_inv[i][j]);
        }
        printf("\n");
    }
    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
    B[i][j]=0.0;
    for(int k=0;k<N;k++)
    B[i][j]+=AA[i][k]*A_inv[k][j];
    printf("%lf,",B[i][j]);
    //while(B[i][j]>1129.0)
    //B[i][j]-=1129.0;
    }
    printf("\n");
    }
    //exit(1);
    //matrix_multiply(A,A_inv,*B);
    
        int num_processes = 8; // 4つのプロセスを使用して行列掛け算を高速化

    
    // 共有メモリを作成
    int shmid = shmget(SHM_KEY, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, IPC_CREAT | 0666);
    if (shmid == -1)
    {
        perror("shmget");
        exit(1);
    }

    int *C = (int *)shmat(shmid, NULL, 0);
    if (C == (int *)-1)
    {
        perror("shmat");
        exit(1);
    }
    
    // 各プロセスで一部の行を計算
    for (int i = 0; i < num_processes; i++)
    {
        pid_t pid = fork();

        if (pid == 0) // 子プロセスの処理
        {
            int start_row = i * (MATRIX_SIZE / num_processes);
            int end_row = (i + 1) * (MATRIX_SIZE / num_processes);

            matrix_multiply(A, A_inv, C, start_row, end_row);
            //matrix_multiply_simd(A, B, C,start_row,end_row);

            exit(0); // 子プロセスを終了
        }
        else if (pid < 0)
        {
            // エラーハンドリング
            perror("fork");
            exit(1);
        }
    }

    // 親プロセスが子プロセスの終了を待つ
    for (int i = 0; i < num_processes; i++)
    {
        int status;
        wait(&status);
    }
  
    // 結果を表示（確認のため）
    printf("Result Matrix:\n");
    
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("%lf ", B[i * MATRIX_SIZE + j]);
        }
        printf("\n");
    }

   
    // 共有メモリを解放
    if (shmdt(C) == -1)
    {
        perror("shmdt");
        exit(1);
    }
    if (shmctl(shmid, IPC_RMID, NULL) == -1)
    {
        perror("shmctl");
        exit(1);
    }


    return 0;
}
