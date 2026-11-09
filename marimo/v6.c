#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>
#include <x86intrin.h> // SIMD命令を使用するためのヘッダファイル

#define MATRIX_SIZE 8
#define SHM_KEY 1234

// 行列を表示する関数
void print_matrix(double A[MATRIX_SIZE][MATRIX_SIZE])
{
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
}

// 行列の掛け算関数
void matrix_multiply(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double *C, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            double sum = 0.0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                sum += A[i][k] * B[k][j];
            }
            C[i * MATRIX_SIZE + j] = sum;
        }
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

// 行列の逆行列を計算する関数
void inverseMatrix(double A[MATRIX_SIZE][MATRIX_SIZE], double A_inv[MATRIX_SIZE][MATRIX_SIZE], int start_row, int end_row)
{
    int i, j, k;
    double temp;

    // 単位行列を初期化
    for (i = start_row; i < end_row; i++)
    {
        for (j = 0; j < MATRIX_SIZE; j++)
        {
            A_inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // ガウス・ジョルダン法による逆行列の計算
    for (k = start_row; k < end_row; k++)
    {
        temp = A[k][k];
        for (j = 0; j < MATRIX_SIZE; j++)
        {
            A[k][j] /= temp;
            A_inv[k][j] /= temp;
        }
        for (i = start_row; i < end_row; i++)
        {
            if (i != k)
            {
                temp = A[i][k];
                for (j = 0; j < MATRIX_SIZE; j++)
                {
                    A[i][j] -= A[k][j] * temp;
                    A_inv[i][j] -= A_inv[k][j] * temp;
                }
            }
        }
    }
}

// 行列の逆行列を計算する関数
void inverseMatrix2(double A[MATRIX_SIZE][MATRIX_SIZE], double A_inv[MATRIX_SIZE][MATRIX_SIZE])
{
    int i, j, k;
    double temp;

    // 単位行列を初期化
    for (i = 0; i < MATRIX_SIZE; i++)
    {
        for (j = 0; j < MATRIX_SIZE; j++)
        {
            A_inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // ガウス・ジョルダン法による逆行列の計算
    for (k = 0; k < MATRIX_SIZE; k++)
    {
        temp = A[k][k];
        for (j = 0; j < MATRIX_SIZE; j++)
        {
            A[k][j] /= temp;
            A_inv[k][j] /= temp;
        }
        for (i = 0; i < MATRIX_SIZE; i++)
        {
            if (i != k)
            {
                temp = A[i][k];
                for (j = 0; j < MATRIX_SIZE; j++)
                {
                    A[i][j] -= A[k][j] * temp;
                    A_inv[i][j] -= A_inv[k][j] * temp;
                }
            }
        }
    }
}

int main()
{
    double A[MATRIX_SIZE][MATRIX_SIZE];
    double A_inv[MATRIX_SIZE][MATRIX_SIZE];
    double C[MATRIX_SIZE][MATRIX_SIZE];
    double AA[MATRIX_SIZE][MATRIX_SIZE];

    // 行列 A を初期化
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            A[i][j] = 1.0 + rand() % 7;
            // printf("%lf,",A[i][j]);
            //  行列 A のコピーを作成
            AA[i][j] = A[i][j];
        }
        // printf("\n");
    }


     inverseMatrix2(A,A_inv);
    //matrix_inverse_simd(A,A_inv);
    //
    print_matrix(A);
        print_matrix(A_inv);
        matrix_multiply(A,A_inv,C,0,MATRIX_SIZE);
        //matmul_simd(A,A_inv,*C,0,MATRIX_SIZE);
        print_matrix(C);
    exit(1);
    // マルチプロセスで行列掛け算を並列化
    int num_processes = 8;
    int rows_per_process = MATRIX_SIZE / num_processes;

    int shmid = shmget(SHM_KEY, sizeof(double) * MATRIX_SIZE * MATRIX_SIZE, IPC_CREAT | 0666);
    if (shmid == -1)
    {
        perror("shmget");
        exit(1);
    }

    double *shared_C = (double *)shmat(shmid, NULL, 0);
    if (shared_C == (double *)-1)
    {
        perror("shmat");
        exit(1);
    }

    // 各プロセスで一部の行を計算
    for (int i = 0; i < num_processes; i++)
    {
        pid_t pid = fork();

        if (pid == 0)
        {
            int start_row = i * rows_per_process;
            int end_row = (i + 1) * rows_per_process;

            inverseMatrix(A, A_inv, start_row, end_row);
            //matmul_simd(AA, A_inv, shared_C, start_row, end_row);
            matrix_multiply(AA, A_inv, shared_C, start_row, end_row);

            // 結果を表示
            printf("Process %d: Rows %d to %d completed\n", i, start_row, end_row);

            exit(0);
        }
        else if (pid < 0)
        {
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

    // 結果を表示
    printf("Result Matrix:\n");
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("%lf ", shared_C[i * MATRIX_SIZE + j]);
        }
        printf("\n");
    }

    // 共有メモリを解放
    if (shmdt(shared_C) == -1)
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
