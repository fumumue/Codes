/* S^-1 E T^-1 の平均重み(最小版)。gcc -O2 -o simple simple.c && ./simple */
#include <stdio.h>
#include <stdlib.h>
#define N 31
typedef int Mat[N][N];

/* A の逆行列を B に求める(掃き出し法, mod 2)。特異なら 0 */
int inverse(Mat A, Mat B) {
    Mat M; 
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) { M[i][j] = A[i][j]; B[i][j] = (i == j); }
    for (int c = 0; c < N; c++) {
        int p = c; while (p < N && !M[p][c]) p++;
        if (p == N) return 0;
        for (int j = 0; j < N; j++) { int t = M[c][j]; M[c][j] = M[p][j]; M[p][j] = t;
                                       t = B[c][j]; B[c][j] = B[p][j]; B[p][j] = t; }
        for (int r = 0; r < N; r++) if (r != c && M[r][c])
            for (int j = 0; j < N; j++) { M[r][j] ^= M[c][j]; B[r][j] ^= B[c][j]; }
    }
    return 1;
}
/* 密度 p のランダム正則行列(逆行列も返す) */
void random_invertible(Mat A, Mat Ainv, double p) {
    do for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) A[i][j] = rand() < p * RAND_MAX;
    while (!inverse(A, Ainv));
}
int main(void) {
    double p = 0.3; int w = 1, trials = 1000; long total = 0;
    for (int t = 0; t < trials; t++) {
        Mat S, Si, T, Ti, E = {{0}}, X, R;
        random_invertible(S, Si, p); random_invertible(T, Ti, p);
        for (int k = 0; k < w; k++) E[k][k] = 1;          /* 対角線上に w 個の 1 */
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {   /* X = Si*E */
            X[i][j] = 0; for (int k = 0; k < N; k++) X[i][j] ^= Si[i][k] & E[k][j]; }
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {   /* R = X*Ti */
            R[i][j] = 0; for (int k = 0; k < N; k++) R[i][j] ^= X[i][k] & Ti[k][j];
            total += R[i][j]; }
    }
    printf("p=%.2f w=%d  平均重み = %.1f (n^2/2 = %.0f)\n", p, w, (double)total / trials, N * N / 2.0);
}