#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <immintrin.h> // SSE4.1を使用するためのヘッダファイル

#define MATRIX_SIZE 8 // 行列のサイズ
#define SHM_KEY 1234     // 共有メモリのキー

// 行列掛け算関数（SIMD版、SSE4.1）
void matrix_multiply_simd(int A[MATRIX_SIZE][MATRIX_SIZE], int B[MATRIX_SIZE][MATRIX_SIZE], int *C, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j += 4) // 4つの要素を同時に処理するために4で割る
        {
            __m128i sum = _mm_setzero_si128(); // ゼロで初期化された128ビット整数ベクトルを作成

            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                __m128i a = _mm_set1_epi32(A[i][k]); // 行列Aの要素をベクトル化
                __m128i b = _mm_loadu_si128((__m128i *)&B[k][j]); // 行列Bの要素をベクトル化

                sum = _mm_add_epi32(sum, _mm_mullo_epi32(a, b)); // ベクトル同士の掛け算と足し算
            }

            // ベクトルの要素を合算して結果を格納
            int result[4];
            _mm_storeu_si128((__m128i *)result, sum);
            for (int k = 0; k < 4; k++)
            {
                C[i * MATRIX_SIZE + j + k] = result[k];
            }
        }
    }
}

// 行列の逆行列を計算する関数
void inverseMatrix(int A[MATRIX_SIZE][MATRIX_SIZE], int A_inv[MATRIX_SIZE][MATRIX_SIZE]) {
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

int main()
{
    int A[MATRIX_SIZE][MATRIX_SIZE];
    int B[MATRIX_SIZE][MATRIX_SIZE];

    // 行列 A と行列 B を初期化
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            A[i][j] = (1129 + i + j) % 1129;
            B[i][j] = (1129 + i - j) % 1129;
        }
    }

    inverseMatrix(A,B);
    int num_processes = 8; // 16個のプロセスを使用して行列掛け算を高速化

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

            matrix_multiply_simd(A, B, C, start_row, end_row);

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

    // 共有メモリから結果を取得して表示
    printf("Result Matrix:\n");
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("%d ", A[i * MATRIX_SIZE + j]);
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
