// dh_kem_demo.c
// 記事「ライバルはHQC」の方式に、送信側の乱数と再暗号化チェックを加えた KEM の動作確認用プログラム。
// 単体でコンパイルできる(vmo.c などは不要)。
//
//   gcc -O2 -o dh_kem_demo dh_kem_demo.c
//   ./dh_kem_demo            # 重み W=3 で 1000 回(既定)
//   ./dh_kem_demo 3 2000     # 重み W=3 で 2000 回
//   ./dh_kem_demo 5 500      # 重みを大きくして、復号失敗が出る様子を見る
//
// 注意(重要):
//   * 研究・学習用のデモで、暗号として使ってはいけない。パラメータは非常に小さく、安全性はない。
//   * 乱数は /dev/urandom を使うが、疎ベクトルの生成は SplitMix64(暗号用ではない)で行っている。
//   * 符号は素体 F_P 上の一般化 Reed-Solomon 符号(Goppa 符号で m=1 の場合)。
//     本来の Goppa 符号を使うときは、syndrome() と decode_errors() を差し替える。
//
// 構成:
//   環 R = F_P[x]/(x^N - 1)(積は巡回畳み込み)、h は公開の一様乱数。
//   鍵生成(Bob):  s_B = x_B + h*y_B     公開鍵 (h, s_B)、秘密鍵 y_B
//   カプセル化(Alice): 乱数 m から (x', y') を決定的に作り、
//        s' = x' + h*y',  u = s_B*y',  a = H*u,  d = m xor Hash(u)
//        暗号文 (s', a, d),  K = Hash(m, 暗号文)
//   脱カプセル化(Bob): v = s'*y_B,  H*(u - v) = a - H*v  を復号して Delta = u - v を得る
//        (u - v = x_B*y' - x'*y_B: 共通項 h*y_B*y' が消える)
//        u = v + Delta,  m = d xor Hash(u),  再暗号化して一致を確かめ、不一致なら疑似乱数の鍵を返す。

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define P 257          // 素体の位数(素数)
#define N 251          // 環のサイズ・符号長(素数、N <= P)
#define TCAP 18        // 符号の訂正能力 t
#define NS (2 * TCAP)  // シンドロームの長さ(記号数)
#define WMAX 16

// ---------------------------------------------------------------- SHA-256
typedef struct { uint32_t h[8]; uint64_t len; uint8_t buf[64]; size_t n; } sha_t;
static const uint32_t K256[64] = {
 0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
 0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
 0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
 0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
 0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
 0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
 0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
 0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
#define ROR(x,n) (((x) >> (n)) | ((x) << (32 - (n))))
static void sha_block(sha_t *s, const uint8_t *p)
{
    uint32_t w[64], a, b, c, d, e, f, g, h;
    for (int i = 0; i < 16; i++)
        w[i] = (uint32_t)p[4*i] << 24 | (uint32_t)p[4*i+1] << 16 | (uint32_t)p[4*i+2] << 8 | p[4*i+3];
    for (int i = 16; i < 64; i++) {
        uint32_t s0 = ROR(w[i-15],7) ^ ROR(w[i-15],18) ^ (w[i-15] >> 3);
        uint32_t s1 = ROR(w[i-2],17) ^ ROR(w[i-2],19) ^ (w[i-2] >> 10);
        w[i] = w[i-16] + s0 + w[i-7] + s1;
    }
    a=s->h[0]; b=s->h[1]; c=s->h[2]; d=s->h[3]; e=s->h[4]; f=s->h[5]; g=s->h[6]; h=s->h[7];
    for (int i = 0; i < 64; i++) {
        uint32_t S1 = ROR(e,6) ^ ROR(e,11) ^ ROR(e,25);
        uint32_t ch = (e & f) ^ (~e & g);
        uint32_t t1 = h + S1 + ch + K256[i] + w[i];
        uint32_t S0 = ROR(a,2) ^ ROR(a,13) ^ ROR(a,22);
        uint32_t mj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t t2 = S0 + mj;
        h=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }
    s->h[0]+=a; s->h[1]+=b; s->h[2]+=c; s->h[3]+=d; s->h[4]+=e; s->h[5]+=f; s->h[6]+=g; s->h[7]+=h;
}
static void sha_init(sha_t *s)
{
    static const uint32_t iv[8] = {0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    memcpy(s->h, iv, sizeof iv); s->len = 0; s->n = 0;
}
static void sha_update(sha_t *s, const uint8_t *d, size_t n)
{
    s->len += n;
    while (n--) { s->buf[s->n++] = *d++; if (s->n == 64) { sha_block(s, s->buf); s->n = 0; } }
}
static void sha_final(sha_t *s, uint8_t out[32])
{
    uint64_t bits = s->len * 8;
    uint8_t pad = 0x80;
    sha_update(s, &pad, 1);
    pad = 0;
    while (s->n != 56) sha_update(s, &pad, 1);
    uint8_t lb[8];
    for (int i = 0; i < 8; i++) lb[i] = (uint8_t)(bits >> (56 - 8*i));
    sha_update(s, lb, 8);
    for (int i = 0; i < 8; i++) {
        out[4*i] = s->h[i] >> 24; out[4*i+1] = s->h[i] >> 16; out[4*i+2] = s->h[i] >> 8; out[4*i+3] = s->h[i];
    }
}

// ---------------------------------------------------------------- 乱数
static void os_random(uint8_t *buf, size_t n)
{
    FILE *f = fopen("/dev/urandom", "rb");
    if (!f || fread(buf, 1, n, f) != n) { fprintf(stderr, "cannot read /dev/urandom\n"); exit(1); }
    fclose(f);
}
typedef struct { uint64_t s; } rng_t;
static uint64_t rng_next(rng_t *r)
{
    uint64_t z = (r->s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static void rng_from_seed(rng_t *r, const char *tag, const uint8_t *seed, size_t n)
{
    sha_t s; uint8_t h[32];
    sha_init(&s); sha_update(&s, (const uint8_t *)tag, strlen(tag)); sha_update(&s, seed, n); sha_final(&s, h);
    r->s = 0; for (int i = 0; i < 8; i++) r->s = (r->s << 8) | h[i];
}

// ---------------------------------------------------------------- F_P と環 R = F_P[x]/(x^N - 1)
typedef struct { uint16_t c[N]; } vec;
static int md(int x) { x %= P; return x < 0 ? x + P : x; }
static int fpow(int a, int e) { long r = 1, b = a % P; while (e) { if (e & 1) r = r * b % P; b = b * b % P; e >>= 1; } return (int)r; }
static int finv(int a) { return fpow(a, P - 2); }

static void vec_add(vec *r, const vec *a, const vec *b) { for (int i = 0; i < N; i++) r->c[i] = (a->c[i] + b->c[i]) % P; }
static void vec_sub(vec *r, const vec *a, const vec *b) { for (int i = 0; i < N; i++) r->c[i] = (a->c[i] + P - b->c[i]) % P; }
static int vec_wt(const vec *a) { int w = 0; for (int i = 0; i < N; i++) w += a->c[i] != 0; return w; }
// r = sp * dn (巡回畳み込み)。sp は疎な側
static void ring_mul(vec *r, const vec *sp, const vec *dn)
{
    long acc[N]; memset(acc, 0, sizeof acc);
    for (int i = 0; i < N; i++) if (sp->c[i])
        for (int j = 0; j < N; j++) { int k = i + j; if (k >= N) k -= N; acc[k] += (long)sp->c[i] * dn->c[j]; }
    for (int i = 0; i < N; i++) r->c[i] = (uint16_t)(acc[i] % P);
}
static void sample_sparse(vec *v, int w, rng_t *r)
{
    memset(v, 0, sizeof *v);
    int cnt = 0;
    while (cnt < w) {
        int pos = (int)(rng_next(r) % N);
        if (v->c[pos]) continue;
        v->c[pos] = (uint16_t)(1 + rng_next(r) % (P - 1));
        cnt++;
    }
}
static void sample_dense(vec *v, rng_t *r) { for (int i = 0; i < N; i++) v->c[i] = (uint16_t)(rng_next(r) % P); }

// ---------------------------------------------------------------- 符号 H(一般化 Reed-Solomon / Goppa で m=1)
// 位置 i の評価点 alpha_i = i + 1。シンドローム S_j = sum_i e_i alpha_i^j (j = 0 .. NS-1)
static int pw_tab[NS][N];
static void init_code(void)
{
    for (int i = 0; i < N; i++) { int a = i + 1, x = 1; for (int j = 0; j < NS; j++) { pw_tab[j][i] = x; x = x * a % P; } }
}
static void syndrome(int S[NS], const vec *e)
{
    for (int j = 0; j < NS; j++) { long acc = 0; for (int i = 0; i < N; i++) if (e->c[i]) acc += (long)e->c[i] * pw_tab[j][i]; S[j] = (int)(acc % P); }
}
// Berlekamp-Massey + Chien 探索 + Forney。成功なら 1、失敗なら 0。
static int decode_errors(const int S[NS], vec *err)
{
    int C[NS + 2] = {1}, B[NS + 2] = {1}, T[NS + 2];
    int L = 0, m = 1, b = 1;
    for (int n = 0; n < NS; n++) {
        long d = S[n];
        for (int i = 1; i <= L; i++) d += (long)C[i] * S[n - i];
        d %= P;
        if (d == 0) { m++; continue; }
        int coef = (int)(d * finv(b) % P);
        if (2 * L <= n) {
            memcpy(T, C, sizeof T);
            for (int i = 0; i + m < NS + 2; i++) C[i + m] = md(C[i + m] - coef * B[i]);
            L = n + 1 - L; memcpy(B, T, sizeof B); b = (int)d; m = 1;
        } else {
            for (int i = 0; i + m < NS + 2; i++) C[i + m] = md(C[i + m] - coef * B[i]);
            m++;
        }
    }
    memset(err, 0, sizeof *err);
    if (L == 0) return 1;                 // 誤りなし
    if (L > TCAP) return 0;
    // Omega(x) = S(x) * Lambda(x) mod x^L
    int Om[NS] = {0};
    for (int j = 0; j < L; j++) { long a = 0; for (int k = 0; k <= j; k++) a += (long)C[k] * S[j - k]; Om[j] = (int)(a % P); }
    int roots = 0;
    for (int i = 0; i < N; i++) {
        int X = i + 1, Xi = finv(X);
        long lam = 0, xp = 1;
        for (int k = 0; k <= L; k++) { lam += xp * C[k]; xp = xp * Xi % P; }
        if (lam % P) continue;
        long om = 0; xp = 1;
        for (int k = 0; k < L; k++) { om += xp * Om[k]; xp = xp * Xi % P; }
        long dl = 0; xp = 1;
        for (int k = 1; k <= L; k++) { dl += (long)k * C[k] % P * xp; xp = xp * Xi % P; }
        om %= P; dl %= P;
        if (dl == 0) return 0;
        int val = md((int)(-(long)X * om % P * finv((int)dl) % P));
        err->c[i] = (uint16_t)val;
        roots++;
    }
    if (roots != L) return 0;
    // 検証: 復号した誤りのシンドロームが入力と一致するか
    int S2[NS]; syndrome(S2, err);
    for (int j = 0; j < NS; j++) if (S2[j] != S[j]) return 0;
    return 1;
}

// ---------------------------------------------------------------- KEM
static int W = 3;   // x, y の重み
typedef struct { vec h, sB, yB; uint8_t z[32]; } keypair_t;
typedef struct { vec sp; int a[NS]; uint8_t d[32]; } ct_t;

static void hash_u(uint8_t out[32], const vec *u)
{
    sha_t s; uint8_t b[2 * N]; sha_init(&s);
    sha_update(&s, (const uint8_t *)"m", 1);
    for (int i = 0; i < N; i++) { b[2*i] = u->c[i] & 255; b[2*i+1] = u->c[i] >> 8; }
    sha_update(&s, b, sizeof b); sha_final(&s, out);
}
static void hash_key(uint8_t out[32], const char *tag, const uint8_t *m, const ct_t *ct)
{
    sha_t s; uint8_t b[2 * N + 2 * NS]; sha_init(&s);
    sha_update(&s, (const uint8_t *)tag, strlen(tag)); sha_update(&s, m, 32);
    for (int i = 0; i < N; i++) { b[2*i] = ct->sp.c[i] & 255; b[2*i+1] = ct->sp.c[i] >> 8; }
    sha_update(&s, b, 2 * N);
    for (int i = 0; i < NS; i++) { b[2*i] = ct->a[i] & 255; b[2*i+1] = ct->a[i] >> 8; }
    sha_update(&s, b, 2 * NS); sha_update(&s, ct->d, 32); sha_final(&s, out);
}

static void keygen(keypair_t *kp)
{
    uint8_t seed[40]; rng_t r;
    os_random(seed, sizeof seed);
    rng_from_seed(&r, "keygen", seed, sizeof seed);
    os_random(kp->z, 32);
    vec xB, hy;
    sample_dense(&kp->h, &r);
    sample_sparse(&xB, W, &r);
    sample_sparse(&kp->yB, W, &r);
    ring_mul(&hy, &kp->yB, &kp->h);
    vec_add(&kp->sB, &xB, &hy);
}

// m から (x', y') を決定的に作り、暗号文と u を計算する
static void encaps_from_m(ct_t *ct, vec *u_out, const keypair_t *pk, const uint8_t m[32])
{
    rng_t r; vec xp, yp, hy, u; int S[NS];
    rng_from_seed(&r, "xy", m, 32);
    sample_sparse(&xp, W, &r);
    sample_sparse(&yp, W, &r);
    ring_mul(&hy, &yp, &pk->h);
    vec_add(&ct->sp, &xp, &hy);          // s' = x' + h y'
    ring_mul(&u, &yp, &pk->sB);          // u = s_B y'
    syndrome(S, &u);
    for (int j = 0; j < NS; j++) ct->a[j] = S[j];
    uint8_t hu[32]; hash_u(hu, &u);
    for (int i = 0; i < 32; i++) ct->d[i] = m[i] ^ hu[i];
    if (u_out) *u_out = u;
}
static void encaps(ct_t *ct, uint8_t K[32], vec *u_out, const keypair_t *pk)
{
    uint8_t m[32];
    os_random(m, 32);                    // 毎回新しい乱数
    encaps_from_m(ct, u_out, pk, m);
    hash_key(K, "K", m, ct);
}
// 戻り値: 1 = 正常に復号、0 = 復号失敗または再暗号化の不一致(疑似乱数の鍵を返す)
static int decaps(uint8_t K[32], const ct_t *ct, const keypair_t *kp, int *decode_ok)
{
    vec v, delta, u, d_vec; int sig[NS], Sv[NS];
    ring_mul(&v, &kp->yB, &ct->sp);      // v = s' y_B
    syndrome(Sv, &v);
    for (int j = 0; j < NS; j++) sig[j] = md(ct->a[j] - Sv[j]);   // H(u - v)
    int ok = decode_errors(sig, &delta);
    if (decode_ok) *decode_ok = ok;
    if (ok) {
        vec_add(&u, &v, &delta);         // u = v + (u - v)
        uint8_t hu[32], m[32]; hash_u(hu, &u);
        for (int i = 0; i < 32; i++) m[i] = ct->d[i] ^ hu[i];
        ct_t chk; encaps_from_m(&chk, NULL, kp, m);
        if (memcmp(&chk.sp, &ct->sp, sizeof chk.sp) == 0 && memcmp(chk.a, ct->a, sizeof chk.a) == 0 &&
            memcmp(chk.d, ct->d, 32) == 0) { hash_key(K, "K", m, ct); return 1; }
    }
    hash_key(K, "rej", kp->z, ct);       // 暗黙の拒否
    (void)d_vec;
    return 0;
}

int main(int argc, char **argv)
{
    // SHA-256 の自己テスト
    { sha_t s; uint8_t h[32]; sha_init(&s); sha_update(&s, (const uint8_t *)"abc", 3); sha_final(&s, h);
      static const uint8_t ex[4] = {0xba, 0x78, 0x16, 0xbf};
      if (memcmp(h, ex, 4)) { fprintf(stderr, "SHA-256 self-test failed\n"); return 1; } }

    int trials = 1000;
    if (argc > 1) W = atoi(argv[1]);
    if (argc > 2) trials = atoi(argv[2]);
    if (W < 1 || W > WMAX || trials < 1) { fprintf(stderr, "usage: %s [W 1..%d] [trials]\n", argv[0], WMAX); return 2; }
    init_code();
    if (!getenv("DEMO_QUIET"))
        printf("F_%d 上、環 R = F_%d[x]/(x^%d - 1)、符号の訂正能力 t = %d、重み W = %d\n"
               "差の重みの上界 2*W*W = %d  (t 以下なら必ず復号できる設計)\n\n", P, P, N, TCAP, W, 2 * W * W);

    int match = 0, dec_fail = 0, rejected = 0, max_dw = 0, distinct_ct = 0;
    uint8_t prev_s[2 * N]; memset(prev_s, 0, sizeof prev_s);
    keypair_t kp; keygen(&kp);
    for (int t = 0; t < trials; t++) {
        ct_t ct; uint8_t K1[32], K2[32]; vec u, v, delta;
        encaps(&ct, K1, &u, &kp);                       // Alice(毎回乱数を使う)
        int dec_ok = 0;
        int ok = decaps(K2, &ct, &kp, &dec_ok);         // Bob
        ring_mul(&v, &kp.yB, &ct.sp);
        vec_sub(&delta, &u, &v);
        int dw = vec_wt(&delta); if (dw > max_dw) max_dw = dw;
        if (!dec_ok) dec_fail++;
        if (!ok) rejected++;
        if (memcmp(K1, K2, 32) == 0) match++;
        if (memcmp(&ct.sp, prev_s, sizeof ct.sp) != 0) distinct_ct++;
        memcpy(prev_s, &ct.sp, sizeof ct.sp);
        // 暗号文が毎回変わること(乱数の効果)を、s' が前回と異なる回数で確認
    }
    printf("試行 %d 回:\n", trials);
    printf("  共通鍵が一致:        %d / %d\n", match, trials);
    printf("  差の復号に失敗:      %d\n", dec_fail);
    printf("  拒否(疑似乱数鍵):    %d\n", rejected);
    printf("  実測した差の最大重み: %d   (上界 %d, t = %d)\n", max_dw, 2 * W * W, TCAP);
    printf("  暗号文 s' が前回と異なる回数: %d / %d(乱数による)\n", distinct_ct, trials);

    // 改ざんした暗号文は、拒否されて別の鍵になる
    { ct_t ct; uint8_t K1[32], K2[32]; encaps(&ct, K1, NULL, &kp);
      ct.a[0] = (uint16_t)((ct.a[0] + 1) % P);
      int ok = decaps(K2, &ct, &kp, NULL);
      printf("  改ざんした暗号文: %s、鍵は %s\n", ok ? "受理(想定外)" : "拒否された",
             memcmp(K1, K2, 32) ? "一致しない(想定どおり)" : "一致(想定外)"); }
    return 0;
}
