/*
 * grow.c : 誤りを1個ずつ足していくと S^-1 E T^-1 の「1」がどう増えるかを表示する。
 *   gcc -O2 -o grow grow.c
 *   ./grow [-n 31] [-w 6] [-p 0.3] [-s 1] [-r] [-a]
 *     -n 行列サイズ   -w 誤りの重み   -p S,T の密度(中密度なら 0.2〜0.4)
 *     -s 乱数の種     -r 誤りを対角線でなくランダムな位置に置く   -a 全ステップの行列を表示
 * 誤り (i,j) を1個足すと、S^-1 の第i列 と T^-1 の第j行 の外積(ランク1)が XOR で加わる。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXN 64
typedef int Mat[MAXN][MAXN];
static int n = 31;

int inverse(Mat A, Mat B) {
    static Mat M;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { M[i][j] = A[i][j]; B[i][j] = (i == j); }
    for (int c = 0; c < n; c++) {
        int p = c; while (p < n && !M[p][c]) p++;
        if (p == n) return 0;
        for (int j = 0; j < n; j++) { int t = M[c][j]; M[c][j] = M[p][j]; M[p][j] = t;
                                       t = B[c][j]; B[c][j] = B[p][j]; B[p][j] = t; }
        for (int r = 0; r < n; r++) if (r != c && M[r][c])
            for (int j = 0; j < n; j++) { M[r][j] ^= M[c][j]; B[r][j] ^= B[c][j]; }
    }
    return 1;
}
void random_invertible(Mat A, Mat Ainv, double p) {
    do for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = rand() < p * RAND_MAX;
    while (!inverse(A, Ainv));
}
int weight(Mat R) { int s = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) s += R[i][j]; return s; }
void show(Mat R) {
    for (int i = 0; i < n; i++) { for (int j = 0; j < n; j++) putchar(R[i][j] ? '#' : '.'); putchar('\n'); }
}

int main(int argc, char **argv) {
    int w = 6, rnd = 0, all = 0; double p = 0.3; unsigned seed = 1;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-n") && i + 1 < argc) n = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-w") && i + 1 < argc) w = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-p") && i + 1 < argc) p = atof(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i + 1 < argc) seed = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-r")) rnd = 1;
        else if (!strcmp(argv[i], "-a")) all = 1;
    }
    if (n < 2 || n > MAXN || w < 1 || (!rnd && w > n)) { fprintf(stderr, "n は 2..%d、w は 1..n(対角のとき)\n", MAXN); return 1; }
    srand(seed);
    static Mat S, Si, T, Ti, R;
    random_invertible(S, Si, p); random_invertible(T, Ti, p);
    int ci = 0, cj = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { ci += Si[i][j]; cj += Ti[i][j]; }
    printf("n=%d  S,T の密度 p=%.2f  → 実際の S^-1 の密度 %.2f, T^-1 の密度 %.2f\n", n, p, ci / (double)(n * n), cj / (double)(n * n));
    printf("誤りの重み w=%d (%s)   n^2/2 = %.0f\n\n", w, rnd ? "ランダムな位置" : "対角線上", n * n / 2.0);

    int pos[MAXN * MAXN][2];                       /* 誤りの位置 */
    if (!rnd) {                                    /* 対角: 位置をシャッフルして先頭 w 個 */
        int d[MAXN]; for (int i = 0; i < n; i++) d[i] = i;
        for (int i = n - 1; i > 0; i--) { int j = rand() % (i + 1), t = d[i]; d[i] = d[j]; d[j] = t; }
        for (int k = 0; k < w; k++) pos[k][0] = pos[k][1] = d[k];
    } else {
        for (int k = 0; k < w; ) {
            int i = rand() % n, j = rand() % n, dup = 0;
            for (int m = 0; m < k; m++) if (pos[m][0] == i && pos[m][1] == j) dup = 1;
            if (!dup) { pos[k][0] = i; pos[k][1] = j; k++; }
        }
    }
    memset(R, 0, sizeof R);
    int trace[MAXN * MAXN];
    for (int k = 0; k < w; k++) {
        int i = pos[k][0], j = pos[k][1];
        for (int r = 0; r < n; r++) if (Si[r][i]) for (int c = 0; c < n; c++) R[r][c] ^= Ti[j][c];
        trace[k] = weight(R);
        if (all || k < 3 || k == w - 1) {
            printf("--- 誤り %d 個目 (%d,%d) を追加: 実効誤りの重み = %d  (%+d)\n", k + 1, i, j, trace[k], trace[k] - (k ? trace[k - 1] : 0));
            show(R); printf("\n");
        } else if (k == 3) printf("   (中略: -a で全ステップを表示)\n\n");
    }
    printf("重みの推移 (誤り個数 → 実効誤りの重み):\n");
    for (int k = 0; k < w; k++) {
        printf("%3d → %4d  ", k + 1, trace[k]);
        int bar = trace[k] * 50 / (n * n / 2 + 1); for (int b = 0; b < bar && b < 70; b++) putchar('|');
        putchar('\n');
    }
    return 0;
}