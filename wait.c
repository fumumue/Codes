/*
 * weight_amplification.c
 * C' = S C T + E を S^-1 C' T^-1 = C + S^-1 E T^-1 で戻すとき、
 * 対角型の誤り E(対角線上に w 個の 1)が実効誤り S^-1 E T^-1 でどれだけ重くなるかの統計。
 * S,T は密度 p のランダムなバイナリ行列を、正則になるまで引き直して作る(拒否サンプリング)。
 *
 * ビルド: gcc -O2 -o weight_amplification weight_amplification.c
 * 実行  : ./weight_amplification [-n 31] [-t 1000] [-s 1]      (n <= 64)
 *
 * 行列は 1 行 = uint64_t(bit j = 第 j 列)で持つ。
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef uint64_t u64;
static u64 rng_state = 88172645463325252ULL;

static u64 rng_next(void) {                 /* xorshift64* */
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    return rng_state * 2685821657736338717ULL;
}
static double rng_unif(void) { return (rng_next() >> 11) * (1.0 / 9007199254740992.0); }

/* A (n 行) の逆行列を Ai に書く。正則なら 1、特異なら 0 */
static int inv2(const u64 *A, u64 *Ai, int n) {
    u64 a[64], b[64];
    for (int r = 0; r < n; r++) { a[r] = A[r]; b[r] = 1ULL << r; }
    for (int c = 0; c < n; c++) {
        int p = -1;
        for (int r = c; r < n; r++) if ((a[r] >> c) & 1) { p = r; break; }
        if (p < 0) return 0;
        if (p != c) { u64 t = a[p]; a[p] = a[c]; a[c] = t; t = b[p]; b[p] = b[c]; b[c] = t; }
        for (int r = 0; r < n; r++)
            if (r != c && ((a[r] >> c) & 1)) { a[r] ^= a[c]; b[r] ^= b[c]; }
    }
    for (int r = 0; r < n; r++) Ai[r] = b[r];
    return 1;
}

/* 密度 p のランダム行列を正則になるまで引き直す。試行回数を返す */
static int rand_inv(u64 *A, u64 *Ai, int n, double p) {
    int tries = 0;
    for (;;) {
        tries++;
        for (int r = 0; r < n; r++) {
            u64 row = 0;
            for (int j = 0; j < n; j++) if (rng_unif() < p) row |= 1ULL << j;
            A[r] = row;
        }
        if (inv2(A, Ai, n)) return tries;
    }
}

static double density(const u64 *A, int n) {
    long s = 0;
    for (int r = 0; r < n; r++) s += __builtin_popcountll(A[r]);
    return (double)s / ((double)n * n);
}

/* S^-1 E T^-1 の重み。E = 対角線上の idx[0..w-1] に 1 */
static long err_weight(const u64 *Si, const u64 *Ti, const int *idx, int w, int n) {
    long wt = 0;
    for (int r = 0; r < n; r++) {
        u64 row = 0;
        for (int k = 0; k < w; k++)
            if ((Si[r] >> idx[k]) & 1) row ^= Ti[idx[k]];
        wt += __builtin_popcountll(row);
    }
    return wt;
}

static void pick_subset(int *idx, int w, int n) {     /* 部分 Fisher-Yates */
    int pool[64];
    for (int i = 0; i < n; i++) pool[i] = i;
    for (int k = 0; k < w; k++) {
        int j = k + (int)(rng_next() % (u64)(n - k));
        int t = pool[k]; pool[k] = pool[j]; pool[j] = t;
        idx[k] = pool[k];
    }
}

int main(int argc, char **argv) {
    int n = 31, trials = 1000; u64 seed = 1;
    for (int i = 1; i < argc - 1; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 'n') n = atoi(argv[++i]);
        else if (argv[i][0] == '-' && argv[i][1] == 't') trials = atoi(argv[++i]);
        else if (argv[i][0] == '-' && argv[i][1] == 's') seed = strtoull(argv[++i], NULL, 10);
    }
    if (n < 2 || n > 64 || trials < 1) { fprintf(stderr, "n は 2..64、trials >= 1\n"); return 1; }
    rng_state ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 20; i++) rng_next();

    int ws[6] = {1, 2, 4, 8, n / 2, n};
    double ps[6] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.50};
    printf("n=%d, 試行%d回。E = 対角線上に w 個の1(位置はランダム)。値は S^-1 E T^-1 の平均重み(ビット数)\n", n, trials);
    printf("%5s %8s %10s %9s |", "p", "dens(S)", "dens(S^-1)", "試行/採用");
    for (int k = 0; k < 6; k++) printf(" w=%-6d", ws[k]);
    printf("| 参考 n^2/2\n");

    for (int pi = 0; pi < 6; pi++) {
        double ds = 0, di = 0, tot[6] = {0};
        long tr = 0;
        for (int t = 0; t < trials; t++) {
            u64 S[64], Si[64], T[64], Ti[64]; int idx[64];
            tr += rand_inv(S, Si, n, ps[pi]);
            tr += rand_inv(T, Ti, n, ps[pi]);
            ds += density(S, n); di += density(Si, n);
            for (int k = 0; k < 6; k++) {
                pick_subset(idx, ws[k], n);
                tot[k] += (double)err_weight(Si, Ti, idx, ws[k], n);
            }
        }
        printf("%5.2f %8.3f %10.3f %9.1f |", ps[pi], ds / trials, di / trials, (double)tr / (2.0 * trials));
        for (int k = 0; k < 6; k++) printf(" %8.1f", tot[k] / trials);
        printf(" | %.0f\n", n * (double)n / 2.0);
    }
    return 0;
}