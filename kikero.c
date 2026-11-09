/*
 * 連接符号(正方形)
 *   内側: RM(1,4)[16,5,8] と 拡大BCH[16,7,6] の積符号 [256,35,48] (16x16)
 *         → 35 bit を 1 シンボルとして扱う
 *   外側: 7x7 の2次元巡回符号 [49,27,>=5] (F_8 上のゼロ点, Sakata 復号)
 *         49 個のシンボル(各35bit)を 7x7 に並べ、35 本のビット面それぞれに同じ符号を掛ける
 *   全体: (7*16) x (7*16) = 112 x 112 のビット配列, 情報 27*35 = 945 bit
 *
 * 復号
 *   1. 各 16x16 ブロックを積符号として復号 (行=RM最尤, 列=拡大BCH, 反復)
 *      失敗したブロックは「消失」, 誤訂正は「シンボル誤り」になる
 *   2. 各ビット面で errors-and-erasures 復号:
 *        消失位置の値を全通り(2^e)埋めて Sakata(2誤り訂正)で復号し、
 *        消失の外での不一致が最小の符号語を採用
 *      s 個のシンボル誤りと e 個の消失があるとき、2s+e <= 4 なら必ず訂正できる (d>=5)
 *
 * ビルド: gcc -O2 -Wall -o concat_pc_cyclic2d concat_pc_cyclic2d.c && ./concat_pc_cyclic2d
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t  u8;

#define R     16            /* 内側ブロックの一辺 */
#define K1    5             /* RM(1,4) の次元 */
#define K2    7             /* 拡大BCH の次元 */
#define KS    (K1 * K2)     /* シンボルのビット数 35 */
#define MAXIT 10

#define SB    7             /* 外側の一辺 */
#define NS    (SB * SB)     /* 外側の符号長(シンボル数) 49 */
#define N     (SB * R)      /* 全体の一辺 112 */
#define D     3             /* Sakata: 既知シンドロームの次数上限 */
#define TMAX  2
#define GM    8

/* =====================================================================
 *  内側: RM(1,4) x 拡大BCH[16,7,6]
 * ===================================================================== */
static int gexp[30], glog[16];
static void gf16_init(void) {
    gexp[0] = 1;
    for (int i = 1; i < 15; i++) {
        int v = gexp[i - 1] << 1;
        if (v & 16) v ^= 0x13;
        gexp[i] = v;
    }
    for (int i = 15; i < 30; i++) gexp[i] = gexp[i - 15];
    for (int i = 0; i < 15; i++) glog[gexp[i]] = i;
}

static u32 ecol[16];
static u16 erow[9];
static int epiv[9], erank, einfo[K2];
static u16 etab_pat[512];
static unsigned char etab_ok[512];

static void ebch_init(void) {
    for (int p = 0; p < 15; p++)
        ecol[p] = (u32)gexp[p] | ((u32)gexp[(3 * p) % 15] << 4) | (1u << 8);
    ecol[15] = 1u << 8;
    u16 row[9];
    for (int r = 0; r < 9; r++) {
        row[r] = 0;
        for (int p = 0; p < 16; p++) if ((ecol[p] >> r) & 1) row[r] |= (u16)1 << p;
    }
    erank = 0;
    for (int p = 0; p < 16 && erank < 9; p++) {
        int s = -1;
        for (int r = erank; r < 9; r++) if ((row[r] >> p) & 1) { s = r; break; }
        if (s < 0) continue;
        u16 t = row[s]; row[s] = row[erank]; row[erank] = t;
        for (int r = 0; r < 9; r++)
            if (r != erank && ((row[r] >> p) & 1)) row[r] ^= row[erank];
        epiv[erank] = p; erank++;
    }
    for (int r = 0; r < erank; r++) erow[r] = row[r];
    u16 pm = 0;
    for (int r = 0; r < erank; r++) pm |= (u16)1 << epiv[r];
    int k = 0;
    for (int p = 0; p < 16; p++) if (!((pm >> p) & 1)) einfo[k++] = p;
    if (erank != 9 || k != K2) { fprintf(stderr, "拡大BCHの次元が想定と違います\n"); exit(1); }
    memset(etab_ok, 0, sizeof etab_ok);
    for (int p = -1; p < 16; p++)
        for (int q = p + 1; q < 16; q++) {
            u32 s = (p >= 0 ? ecol[p] : 0) ^ ecol[q];
            u16 pat = (u16)(((p >= 0) ? (1u << p) : 0) | (1u << q));
            if (s == 0 || etab_ok[s]) { fprintf(stderr, "拡大BCH: 重み2以下が区別できません\n"); exit(1); }
            etab_ok[s] = 1; etab_pat[s] = pat;
        }
}
static u16 ebch_enc(unsigned m7) {
    u16 cw = 0;
    for (int t = 0; t < K2; t++) if ((m7 >> t) & 1) cw |= (u16)1 << einfo[t];
    for (int r = 0; r < erank; r++)
        if (__builtin_parity(erow[r] & cw)) cw |= (u16)1 << epiv[r];
    return cw;
}
static u32 esyn(u16 w) {
    u32 s = 0;
    for (int p = 0; p < 16; p++) if ((w >> p) & 1) s ^= ecol[p];
    return s;
}
static int ebch_dec(u16 w, u16 *fixed) {
    u32 s = esyn(w);
    if (s == 0) { *fixed = w; return 0; }
    if (!etab_ok[s]) return -1;
    *fixed = w ^ etab_pat[s];
    return 1;
}

static u16 rm_enc(unsigned m5) {
    int c0 = m5 & 1, a = m5 >> 1;
    u16 w = 0;
    for (int x = 0; x < 16; x++) if (c0 ^ __builtin_parity(a & x)) w |= (u16)1 << x;
    return w;
}
static unsigned rm_ext(u16 w) {
    int c0 = w & 1, a = 0;
    for (int i = 0; i < 4; i++) if (((w >> (1 << i)) & 1) ^ c0) a |= 1 << i;
    return c0 | (a << 1);
}
static int rm_dec(u16 w, u16 *cw) {
    int v[16];
    for (int x = 0; x < 16; x++) v[x] = 1 - 2 * ((w >> x) & 1);
    for (int len = 1; len < 16; len <<= 1)
        for (int i = 0; i < 16; i += 2 * len)
            for (int j = i; j < i + len; j++) {
                int a = v[j], b = v[j + len];
                v[j] = a + b; v[j + len] = a - b;
            }
    int best = 0, bi = 0, tie = 0;
    for (int a = 0; a < 16; a++) {
        int m = v[a] < 0 ? -v[a] : v[a];
        if (m > best) { best = m; bi = a; tie = 0; }
        else if (m == best) tie = 1;
    }
    if (tie) return -1;
    *cw = rm_enc((v[bi] > 0 ? 0 : 1) | (bi << 1));
    return (16 - best) / 2;
}

static void pc_encode(u64 msg, u16 arr[R]) {
    u16 rowinfo[K2];
    for (int i = 0; i < K2; i++) rowinfo[i] = rm_enc((msg >> (K1 * i)) & ((1u << K1) - 1));
    memset(arr, 0, sizeof(u16) * R);
    for (int c = 0; c < R; c++) {
        unsigned m7 = 0;
        for (int t = 0; t < K2; t++) m7 |= (unsigned)((rowinfo[t] >> c) & 1) << t;
        u16 col = ebch_enc(m7);
        for (int r = 0; r < R; r++) if ((col >> r) & 1) arr[r] |= (u16)1 << c;
    }
}
static u64 pc_extract(const u16 arr[R]) {
    u64 m = 0;
    for (int i = 0; i < K2; i++) m |= (u64)rm_ext(arr[einfo[i]]) << (K1 * i);
    return m;
}
static u16 getcol(const u16 a[R], int c) {
    u16 col = 0;
    for (int r = 0; r < R; r++) col |= (u16)(((a[r] >> c) & 1) << r);
    return col;
}
/* 0: 誤りなし, 1: 訂正, -1: 失敗(消失) */
static int pc_decode(const u16 in[R], u64 *msg) {
    u16 a[R];
    memcpy(a, in, sizeof a);
    for (int it = 0; it < MAXIT; it++) {
        for (int r = 0; r < R; r++) {
            u16 cw;
            int d = rm_dec(a[r], &cw);
            if (d > 0) a[r] = cw;
        }
        for (int c = 0; c < R; c++) {
            u16 fx;
            if (ebch_dec(getcol(a, c), &fx) == 1)
                for (int r = 0; r < R; r++)
                    a[r] = (u16)((a[r] & ~(1u << c)) | (((fx >> r) & 1u) << c));
        }
        int ok = 1;
        for (int r = 0; r < R && ok; r++) { u16 cw; if (rm_dec(a[r], &cw) != 0) ok = 0; }
        for (int c = 0; c < R && ok; c++) { u16 fx; if (ebch_dec(getcol(a, c), &fx) != 0) ok = 0; }
        if (ok) { *msg = pc_extract(a); return memcmp(a, in, sizeof a) ? 1 : 0; }
    }
    return -1;
}

/* =====================================================================
 *  外側: 7x7 の2次元巡回符号 [49,27,>=5] と Sakata 復号
 * ===================================================================== */
static int gf_exp[14], gf_log[8];
static void gf8_init(void) {
    gf_exp[0] = 1;
    for (int i = 1; i < 7; i++) {
        int v = gf_exp[i - 1] << 1;
        if (v & 8) v ^= 0xB;
        gf_exp[i] = v;
    }
    for (int i = 7; i < 14; i++) gf_exp[i] = gf_exp[i - 7];
    for (int i = 0; i < 7; i++) gf_log[gf_exp[i]] = i;
}
static inline u8 gmul(u8 a, u8 b) { return (a && b) ? (u8)gf_exp[gf_log[a] + gf_log[b]] : 0; }
static inline u8 ginv(u8 a) { return (u8)gf_exp[(7 - gf_log[a]) % 7]; }

static const int zero_reps[][2] = {{0,0},{1,0},{0,1},{1,1},{3,0},{0,3},{2,1},{1,2}};
#define NZ ((int)(sizeof(zero_reps) / sizeof(zero_reps[0])))

static u64 ocol[NS];
static int nsyn;
static u64 orrow[64];
static int opiv[64], orank, oinfo[64], KO;

typedef struct { u64 syn; u64 pat; } Ent;
static int cmp_ent(const void *a, const void *b) {
    u64 x = ((const Ent *)a)->syn, y = ((const Ent *)b)->syn;
    return (x > y) - (x < y);
}
static int in_Z(int i, int j) {
    for (int z = 0; z < NZ; z++) {
        int a = zero_reps[z][0], b = zero_reps[z][1];
        for (int t = 0; t < 3; t++) {
            if (a == i && b == j) return 1;
            a = (2 * a) % 7; b = (2 * b) % 7;
        }
    }
    return 0;
}
static void outer_init(void) {
    gf8_init();
    for (int i = 0; i <= D; i++)
        for (int j = 0; i + j <= D; j++)
            if (!in_Z(i, j)) { fprintf(stderr, "三角形がゼロ点集合に含まれていません\n"); exit(1); }
    nsyn = 0;
    for (int p = 0; p < NS; p++) ocol[p] = 0;
    for (int z = 0; z < NZ; z++) {
        int i = zero_reps[z][0], j = zero_reps[z][1];
        int width = (i == 0 && j == 0) ? 1 : 3;
        for (int p = 0; p < NS; p++) {
            int u = p / SB, v = p % SB;
            int val = (width == 1) ? 1 : gf_exp[(i * u + j * v) % 7];
            ocol[p] |= (u64)val << nsyn;
        }
        nsyn += width;
    }
    u64 row[64];
    for (int r = 0; r < nsyn; r++) {
        row[r] = 0;
        for (int p = 0; p < NS; p++) if ((ocol[p] >> r) & 1) row[r] |= (u64)1 << p;
    }
    orank = 0;
    for (int p = 0; p < NS && orank < nsyn; p++) {
        int s = -1;
        for (int r = orank; r < nsyn; r++) if ((row[r] >> p) & 1) { s = r; break; }
        if (s < 0) continue;
        u64 t = row[s]; row[s] = row[orank]; row[orank] = t;
        for (int r = 0; r < nsyn; r++)
            if (r != orank && ((row[r] >> p) & 1)) row[r] ^= row[orank];
        opiv[orank] = p; orank++;
    }
    for (int r = 0; r < orank; r++) orrow[r] = row[r];
    u64 pm = 0;
    for (int r = 0; r < orank; r++) pm |= (u64)1 << opiv[r];
    KO = 0;
    for (int p = 0; p < NS; p++) if (!((pm >> p) & 1)) oinfo[KO++] = p;

    /* d>=5 の確認: 重み2以下の誤りのシンドロームがすべて異なる */
    int nt = NS + NS * (NS - 1) / 2, e = 0;
    Ent *tab = malloc(sizeof(Ent) * nt);
    for (int p = 0; p < NS; p++) tab[e++] = (Ent){ocol[p], (u64)1 << p};
    for (int p = 0; p < NS; p++)
        for (int q = p + 1; q < NS; q++)
            tab[e++] = (Ent){ocol[p] ^ ocol[q], ((u64)1 << p) | ((u64)1 << q)};
    qsort(tab, nt, sizeof(Ent), cmp_ent);
    for (int t = 0; t < nt; t++)
        if (tab[t].syn == 0 || (t > 0 && tab[t].syn == tab[t - 1].syn)) {
            fprintf(stderr, "外側符号が d<5 です\n"); exit(1);
        }
    free(tab);
}
static u64 osyn(u64 r) {
    u64 s = 0;
    for (int p = 0; p < NS; p++) if ((r >> p) & 1) s ^= ocol[p];
    return s;
}

/* ---------- Berlekamp-Massey-Sakata ---------- */
typedef struct { u8 c[GM][GM]; int ex, ey; } Poly;
typedef struct { Poly g; int sx, sy; u8 d; } Rec;
static u8 SS[D + 1][D + 1];

static void compute_S(u64 r) {
    memset(SS, 0, sizeof SS);
    for (int p = 0; p < NS; p++) {
        if (!((r >> p) & 1)) continue;
        int u = p / SB, v = p % SB;
        for (int i = 0; i <= D; i++)
            for (int j = 0; i + j <= D; j++)
                SS[i][j] ^= (u8)gf_exp[(i * u + j * v) % 7];
    }
}
static int mcmp(int i, int j, int k, int l) {
    int d1 = i + j, d2 = k + l;
    if (d1 != d2) return d1 < d2 ? -1 : 1;
    if (j != l) return j < l ? -1 : 1;
    return 0;
}
static u8 disc(const Poly *f, int mx, int my) {
    u8 s = 0;
    int dx = mx - f->ex, dy = my - f->ey;
    for (int i = 0; i < GM; i++)
        for (int j = 0; j < GM; j++) {
            if (!f->c[i][j]) continue;
            if (i + dx > D || j + dy > D || i + j + dx + dy > D) {
                fprintf(stderr, "内部エラー: 未知のシンドロームを参照\n"); exit(2);
            }
            s ^= gmul(f->c[i][j], SS[i + dx][j + dy]);
        }
    return s;
}
static int addshift(Poly *dst, const Poly *src, int sx, int sy, u8 coef) {
    for (int i = 0; i < GM; i++)
        for (int j = 0; j < GM; j++) {
            if (!src->c[i][j]) continue;
            if (i + sx >= GM || j + sy >= GM) return -1;
            dst->c[i + sx][j + sy] ^= gmul(coef, src->c[i][j]);
        }
    return 0;
}
static int sakata(Poly *outF, int *outn, int *outd) {
    Poly cur[16], nf[16];
    int nF = 1;
    memset(cur, 0, sizeof cur);
    cur[0].c[0][0] = 1;
    u8 delta[GM][GM];
    memset(delta, 0, sizeof delta);
    int dsz = 0;
    Rec rec[64]; int nrec = 0;

    for (int deg = 0; deg <= D; deg++)
        for (int j = 0; j <= deg; j++) {
            int nx = deg - j, ny = j;
            u8 d[16]; int fail[16], any = 0;
            for (int f = 0; f < nF; f++) {
                if (nx >= cur[f].ex && ny >= cur[f].ey) {
                    d[f] = disc(&cur[f], nx, ny);
                    fail[f] = (d[f] != 0);
                } else { d[f] = 0; fail[f] = 0; }
                any |= fail[f];
            }
            if (!any) continue;
            u8 nd[GM][GM];
            memcpy(nd, delta, sizeof nd);
            for (int f = 0; f < nF; f++) {
                if (!fail[f]) continue;
                int rx = nx - cur[f].ex, ry = ny - cur[f].ey;
                for (int a = 0; a <= rx; a++)
                    for (int b = 0; b <= ry; b++) {
                        if (a >= GM || b >= GM) return -1;
                        nd[a][b] = 1;
                    }
            }
            int cnt = 0;
            for (int a = 0; a < GM; a++) for (int b = 0; b < GM; b++) cnt += nd[a][b];
            if (cnt > TMAX) return -1;

            int nn = 0;
            for (int a = 0; a < GM; a++)
                for (int b = 0; b < GM; b++) {
                    if (nd[a][b]) continue;
                    if (a > 0 && !nd[a - 1][b]) continue;
                    if (b > 0 && !nd[a][b - 1]) continue;
                    if (nn >= 16) return -1;
                    int oi = -1;
                    for (int f = 0; f < nF; f++) if (cur[f].ex == a && cur[f].ey == b) oi = f;
                    if (oi >= 0 && !fail[oi]) { nf[nn++] = cur[oi]; continue; }
                    int fi = oi;
                    if (fi < 0) {
                        for (int f = 0; f < nF; f++)
                            if (cur[f].ex <= a && cur[f].ey <= b)
                                if (fi < 0 || (fail[fi] && !fail[f])) fi = f;
                        if (fi < 0) return -1;
                    }
                    Poly np; memset(&np, 0, sizeof np);
                    np.ex = a; np.ey = b;
                    if (addshift(&np, &cur[fi], a - cur[fi].ex, b - cur[fi].ey, 1) < 0) return -1;
                    if (fail[fi] && a <= nx && b <= ny) {
                        int kx = nx - a, ky = ny - b, found = -1;
                        for (int r = nrec - 1; r >= 0; r--) {
                            int bx = rec[r].sx - kx - rec[r].g.ex;
                            int by = rec[r].sy - ky - rec[r].g.ey;
                            if (bx < 0 || by < 0) continue;
                            if (mcmp(rec[r].sx - kx, rec[r].sy - ky, a, b) >= 0) continue;
                            found = r; break;
                        }
                        if (found < 0) return -1;
                        int bx = rec[found].sx - kx - rec[found].g.ex;
                        int by = rec[found].sy - ky - rec[found].g.ey;
                        u8 lam = gmul(d[fi], ginv(rec[found].d));
                        if (addshift(&np, &rec[found].g, bx, by, lam) < 0) return -1;
                    }
                    nf[nn++] = np;
                }
            for (int f = 0; f < nF; f++)
                if (fail[f] && nrec < 64) { rec[nrec].g = cur[f]; rec[nrec].sx = nx; rec[nrec].sy = ny; rec[nrec].d = d[f]; nrec++; }
            memcpy(cur, nf, sizeof(Poly) * nn);
            nF = nn;
            memcpy(delta, nd, sizeof delta);
            dsz = cnt;
        }
    memcpy(outF, cur, sizeof(Poly) * nF);
    *outn = nF; *outd = dsz;
    return 0;
}

/* 2誤りまでの有界距離復号。成功なら補正後の符号語を *cw に返す (0: 誤りなし, 1: 訂正, -1: 失敗) */
static int od_sakata(u64 r, u64 *cw) {
    if (osyn(r) == 0) { *cw = r; return 0; }
    compute_S(r);
    Poly F[16]; int nF, dsz;
    if (sakata(F, &nF, &dsz) < 0) return -1;
    u64 err = 0; int cnt = 0;
    for (int u = 0; u < SB; u++)
        for (int v = 0; v < SB; v++) {
            int zero = 1;
            for (int f = 0; f < nF && zero; f++) {
                u8 s = 0;
                for (int i = 0; i < GM; i++)
                    for (int j = 0; j < GM; j++)
                        if (F[f].c[i][j]) s ^= gmul(F[f].c[i][j], (u8)gf_exp[(i * u + j * v) % 7]);
                if (s) zero = 0;
            }
            if (zero) { err |= (u64)1 << (u * SB + v); cnt++; }
        }
    if (cnt != dsz) return -1;
    u64 c = r ^ err;
    if (osyn(c) != 0) return -1;
    *cw = c;
    return cnt ? 1 : 0;
}

/* errors-and-erasures: emask の位置は無視。成功なら 0 */
static int outer_ee_decode(u64 r, u64 emask, u64 *cwout) {
    int e = __builtin_popcountll(emask);
    if (e > 4) return -1;
    int pos[4], n = 0;
    for (int p = 0; p < NS; p++) if ((emask >> p) & 1) pos[n++] = p;
    u64 base = r & ~emask, best = 0;
    int bestd = 99, found = 0;
    for (int f = 0; f < (1 << e); f++) {
        u64 w = base;
        for (int i = 0; i < e; i++) if ((f >> i) & 1) w |= (u64)1 << pos[i];
        u64 cw;
        if (od_sakata(w, &cw) < 0) continue;
        int dis = __builtin_popcountll((cw ^ r) & ~emask);
        if (dis < bestd) { bestd = dis; best = cw; found = 1; }
    }
    if (!found || 2 * bestd + e > 4) return -1;
    *cwout = best;
    return 0;
}

/* シンボル列 sym[49] の符号化 (各ビット面に同じ線形符号) */
static void outer_encode_sym(const u64 info[], u64 sym[NS]) {
    memset(sym, 0, sizeof(u64) * NS);
    for (int t = 0; t < KO; t++) sym[oinfo[t]] = info[t];
    for (int r = 0; r < orank; r++) {
        u64 x = 0;
        for (int j = 0; j < NS; j++)
            if (((orrow[r] >> j) & 1) && j != opiv[r]) x ^= sym[j];
        sym[opiv[r]] = x;
    }
}
/* シンボル列の復号。emask: 消失位置。成功なら 0 で sym を訂正後に書き換える */
static int outer_decode_sym(u64 sym[NS], u64 emask) {
    u64 out[NS];
    memset(out, 0, sizeof out);
    for (int b = 0; b < KS; b++) {
        u64 r = 0;
        for (int p = 0; p < NS; p++) r |= ((sym[p] >> b) & 1) << p;
        u64 cw;
        if (outer_ee_decode(r, emask, &cw) < 0) return -1;
        for (int p = 0; p < NS; p++) out[p] |= ((cw >> p) & 1) << b;
    }
    memcpy(sym, out, sizeof out);
    return 0;
}

/* =====================================================================
 *  全体
 * ===================================================================== */
typedef unsigned char Grid[N][N];

static void encode(const u64 info[], Grid g) {
    u64 sym[NS];
    outer_encode_sym(info, sym);
    for (int p = 0; p < NS; p++) {
        int I = p / SB, J = p % SB;
        u16 arr[R];
        pc_encode(sym[p], arr);
        for (int r = 0; r < R; r++)
            for (int c = 0; c < R; c++) g[I * R + r][J * R + c] = (arr[r] >> c) & 1;
    }
}
/* 戻り値 0: 成功, -1: 失敗。stat[0..2] = 内側が 誤りなし/訂正/消失 のブロック数 */
static int decode(Grid g, u64 info[], int stat[3]) {
    u64 sym[NS], emask = 0;
    stat[0] = stat[1] = stat[2] = 0;
    for (int p = 0; p < NS; p++) {
        int I = p / SB, J = p % SB;
        u16 arr[R];
        for (int r = 0; r < R; r++) {
            arr[r] = 0;
            for (int c = 0; c < R; c++) if (g[I * R + r][J * R + c]) arr[r] |= (u16)1 << c;
        }
        u64 m = 0;
        int st = pc_decode(arr, &m);
        sym[p] = (st < 0) ? 0 : m;
        if (st < 0) emask |= (u64)1 << p;
        stat[st < 0 ? 2 : st]++;
    }
    if (outer_decode_sym(sym, emask) < 0) return -1;
    for (int t = 0; t < KO; t++) info[t] = sym[oinfo[t]];
    return 0;
}

/* =====================================================================
 *  テスト
 * ===================================================================== */
static u64 rs = 88172645463325252ULL;
static u64 rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
#define MASK ((((u64)1) << KS) - 1)

static void flip_distinct(Grid g, int I, int J, int w) {
    /* ブロック (I,J) の異なる w 個のビットを反転 */
    u16 used[R]; memset(used, 0, sizeof used);
    for (int e = 0; e < w; ) {
        int p = rnd() % (R * R), r = p / R, c = p % R;
        if ((used[r] >> c) & 1) continue;
        used[r] |= (u16)1 << c;
        g[I * R + r][J * R + c] ^= 1;
        e++;
    }
}

int main(void) {
    gf16_init(); ebch_init(); outer_init();
    printf("内側: RM(1,4)[16,5,8] x 拡大BCH[16,7,6] = [256,%d,48] (16x16) → 1シンボル %d bit\n", KS, KS);
    printf("外側: 7x7 2次元巡回符号 [%d,%d,>=5] (F_8 上, Sakata, 検査 %d bit)\n", NS, KO, orank);
    printf("全体: %dx%d = %d bit, 情報 %d bit\n\n", N, N, N * N, KO * KS);

    /* --- 外側符号のみ: シンボル誤り s 個, 消失 e 個 (各 5000 回) --- */
    puts("外側のみ (s=シンボル誤り数, e=消失数, 条件 2s+e<=4 なら必ず訂正):");
    puts("  s  e  2s+e   成功   失敗   誤訂正");
    int cases[][2] = {{0,0},{1,0},{2,0},{0,1},{0,2},{0,3},{0,4},{1,1},{1,2},
                      {3,0},{2,1},{1,3},{0,5},{2,2}};
    for (unsigned ci = 0; ci < sizeof cases / sizeof cases[0]; ci++) {
        int s = cases[ci][0], e = cases[ci][1], ok = 0, fail = 0, mis = 0;
        for (int t = 0; t < 5000; t++) {
            u64 info[64], sym[NS], orig[NS];
            for (int i = 0; i < KO; i++) info[i] = rnd() & MASK;
            outer_encode_sym(info, sym);
            memcpy(orig, sym, sizeof sym);
            u64 used = 0, emask = 0;
            for (int i = 0; i < s; ) {
                int p = rnd() % NS;
                if ((used >> p) & 1) continue;
                used |= (u64)1 << p;
                u64 v; do v = rnd() & MASK; while (!v);
                sym[p] ^= v; i++;
            }
            for (int i = 0; i < e; ) {
                int p = rnd() % NS;
                if ((used >> p) & 1) continue;
                used |= (u64)1 << p; emask |= (u64)1 << p;
                sym[p] = rnd() & MASK;       /* 消失位置の中身は無意味 */
                i++;
            }
            if (outer_decode_sym(sym, emask) < 0) fail++;
            else if (memcmp(sym, orig, sizeof sym) == 0) ok++;
            else mis++;
        }
        printf(" %2d %2d  %3d   %5d  %5d  %6d%s\n", s, e, 2 * s + e, ok, fail, mis, (2 * s + e <= 4) ? "" : "   (保証外)");
    }

    /* --- 全体 --- */
    const int TR = 3000;
    int ok0 = 0;
    for (int t = 0; t < TR; t++) {
        u64 info[64], out[64];
        for (int i = 0; i < KO; i++) info[i] = rnd() & MASK;
        static Grid g; int st[3];
        encode(info, g);
        for (int I = 0; I < SB; I++) for (int J = 0; J < SB; J++) flip_distinct(g, I, J, rnd() % 25);
        if (decode(g, out, st) == 0 && memcmp(info, out, sizeof(u64) * KO) == 0) ok0++;
    }
    printf("\n全体 (各ブロック 0-24 誤り):                          成功 %d / %d\n", ok0, TR);

    for (int h = 1; h <= 6; h++) {
        int ok = 0, mis = 0;
        for (int t = 0; t < TR; t++) {
            u64 info[64], out[64];
            for (int i = 0; i < KO; i++) info[i] = rnd() & MASK;
            static Grid g; int st[3];
            encode(info, g);
            u64 heavy = 0;
            while (__builtin_popcountll(heavy) < h) heavy |= (u64)1 << (rnd() % NS);
            for (int p = 0; p < NS; p++) {
                int w = ((heavy >> p) & 1) ? 100 + rnd() % 60 : rnd() % 13;
                flip_distinct(g, p / SB, p % SB, w);
            }
            int r = decode(g, out, st);
            if (r == 0 && memcmp(info, out, sizeof(u64) * KO) == 0) ok++;
            else if (r == 0) mis++;
        }
        printf("全体 (%d ブロックが 100-159 誤り, 他 0-12 誤り):       成功 %d / %d (誤訂正 %d)\n", h, ok, TR, mis);
    }
    return 0;
}