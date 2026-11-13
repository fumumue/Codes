/*
 * F_17[x, y] の多項式演算（加算・減算・乗算・除算）
 *
 * コンパイル例:
 *   gcc -O2 -Wall -fopenmp poly_f17_omp.c -o poly
 *   （-fopenmp を付けなければ、並列化なしの通常の逐次プログラムになります）
 */
#include <stdio.h>
#include <string.h>

#define P 17
#define MAX_DEG 64
#define MAX_INPUT_DEG 32

/* 並列化を行う最小の作業量。これ未満ならスレッド起動コストの方が大きい */
#define PAR_MUL_MIN 1024
#define PAR_DIV_MIN 512

typedef struct {
    int c[MAX_DEG + 1][MAX_DEG + 1]; /* c[xの次数][yの次数] */
} Poly;

typedef struct {
    int i, j, c; /* 係数 c, xの次数 i, yの次数 j */
} Term;

static int mod(int a)
{
    a %= P;
    return a < 0 ? a + P : a;
}

static int mod_pow(int a, int n)
{
    int result = 1;
    a = mod(a);
    while (n > 0) {
        if (n & 1)
            result = mod(result * a);
        a = mod(a * a);
        n >>= 1;
    }
    return result;
}

static void add_term(Poly *p, int coeff, int xdeg, int ydeg)
{
    p->c[xdeg][ydeg] = mod(p->c[xdeg][ydeg] + coeff);
}

/* 項数のあとに「係数 xの次数 yの次数」を項数分入力する */
static int read_poly(const char *name, Poly *p)
{
    int n;
    memset(p, 0, sizeof(*p));

    printf("%sの項数を入力してください: ", name);
    if (scanf("%d", &n) != 1 || n < 0) {
        fprintf(stderr, "項数が不正です。\n");
        return 0;
    }

    printf("各項を「係数 xの次数 yの次数」の順に入力してください。\n");
    printf("例: 3*x^2*y は 3 2 1\n");

    for (int k = 0; k < n; ++k) {
        int coeff, xdeg, ydeg;
        printf("項%d: ", k + 1);
        if (scanf("%d %d %d", &coeff, &xdeg, &ydeg) != 3 ||
            xdeg < 0 || xdeg > MAX_INPUT_DEG ||
            ydeg < 0 || ydeg > MAX_INPUT_DEG) {
            fprintf(stderr, "入力が不正です。次数は0〜%dで指定してください。\n",
                    MAX_INPUT_DEG);
            return 0;
        }
        add_term(p, coeff, xdeg, ydeg);
    }
    return 1;
}

static void print_poly(const Poly *p)
{
    int first = 1;

    /* 辞書式順序（x優先）で大きい項から表示 */
    for (int i = MAX_DEG; i >= 0; --i) {
        for (int j = MAX_DEG; j >= 0; --j) {
            int coeff = p->c[i][j];
            if (coeff == 0)
                continue;

            if (!first)
                printf(" + ");
            first = 0;

            printf("%d", coeff);
            if (i > 0)
                printf("*x");
            if (i > 1)
                printf("^%d", i);
            if (j > 0)
                printf("*y");
            if (j > 1)
                printf("^%d", j);
        }
    }

    if (first)
        printf("0");
    printf(" (mod %d)\n", P);
}

/*
 * xの最大次数 *mx と yの最大次数 *my を求める（各変数について別々に最大）。
 * 零多項式のときは両方 -1。
 */
static void max_degs(const Poly *p, int *mx, int *my)
{
    *mx = -1;
    *my = -1;
    for (int i = 0; i <= MAX_DEG; ++i)
        for (int j = 0; j <= MAX_DEG; ++j)
            if (p->c[i][j] != 0) {
                if (i > *mx) *mx = i;
                if (j > *my) *my = j;
            }
}

/*
 * 加算・減算: 要素数が65*65しかなく、メモリ律速でもあるため
 * 並列化はせず、コンパイラの自動ベクトル化に任せる。
 */
static Poly add_poly(const Poly *a, const Poly *b)
{
    Poly r;
    for (int i = 0; i <= MAX_DEG; ++i)
        for (int j = 0; j <= MAX_DEG; ++j)
            r.c[i][j] = mod(a->c[i][j] + b->c[i][j]);
    return r;
}

static Poly sub_poly(const Poly *a, const Poly *b)
{
    Poly r;
    for (int i = 0; i <= MAX_DEG; ++i)
        for (int j = 0; j <= MAX_DEG; ++j)
            r.c[i][j] = mod(a->c[i][j] - b->c[i][j]);
    return r;
}

/*
 * 乗算: 出力の各項 (s, t) ごとに和を取る畳み込みの形にしてある。
 * 出力位置ごとに独立なので、書き込みの競合なしに並列化できる。
 * 係数は 0..16 なので、和は最大でも 65*65*16*16 程度で int に収まり、
 * mod は最後に1回だけ取れば足りる。
 * 戻り値: 1=成功, 0=結果の次数が上限を超える
 */
static int mul_poly(const Poly *a, const Poly *b, Poly *r)
{
    int ax, ay, bx, by;

    memset(r, 0, sizeof(*r));
    max_degs(a, &ax, &ay);
    max_degs(b, &bx, &by);

    if (ax < 0 || bx < 0)
        return 1; /* どちらかが0なら積も0 */

    if (ax + bx > MAX_DEG || ay + by > MAX_DEG) {
        fprintf(stderr, "結果の次数が上限(%d)を超えました。\n", MAX_DEG);
        return 0;
    }

    const int ns = ax + bx + 1;
    const int nt = ay + by + 1;

#pragma omp parallel for collapse(2) schedule(static) if (ns * nt >= PAR_MUL_MIN)
    for (int s = 0; s < ns; ++s) {
        for (int t = 0; t < nt; ++t) {
            const int i0 = s > bx ? s - bx : 0;
            const int i1 = s < ax ? s : ax;
            const int j0 = t > by ? t - by : 0;
            const int j1 = t < ay ? t : ay;
            int sum = 0;

            for (int i = i0; i <= i1; ++i)
                for (int j = j0; j <= j1; ++j)
                    sum += a->c[i][j] * b->c[s - i][t - j];

            r->c[s][t] = sum % P;
        }
    }
    return 1;
}

/*
 * 除算（辞書式順序・x優先）: a = q*b + r
 *
 * 改良点:
 *  - 先頭項の探索を毎回65*65走査するのをやめ、添字を大きい方から
 *    1回だけ下りながら処理する。割り算では、先頭項を消すと残りの項は
 *    すべて辞書式順序で小さくなるので、添字は後戻りしない。
 *  - 除数 b の非零項を最初に1回だけリスト化し、1ステップの減算は
 *    そのリストだけを回す。
 *  - 次数の上限チェックを、並列領域の外（ステップの先頭）で行う。
 *  - 項数が多いときだけ、ステップ内の減算を並列化する
 *    （書き込み先は項ごとに別位置なので競合しない）。
 *
 * 戻り値: 1=成功, 0=Bが0、または次数が上限を超えた
 */
static int div_poly(const Poly *a, const Poly *b, Poly *q, Poly *r)
{
    static Term terms[(MAX_DEG + 1) * (MAX_DEG + 1)];
    const int W = MAX_DEG + 1;
    Poly work = *a;
    int nt = 0;
    int bx = -1, by = -1, bc = 0; /* 先頭項 */
    int bxmax, bymax;

    memset(q, 0, sizeof(*q));
    memset(r, 0, sizeof(*r));

    /* b の先頭項（辞書式で最大）を探す */
    for (int i = MAX_DEG; i >= 0 && bx < 0; --i)
        for (int j = MAX_DEG; j >= 0; --j)
            if (b->c[i][j] != 0) {
                bx = i;
                by = j;
                bc = b->c[i][j];
                break;
            }
    if (bx < 0)
        return 0; /* 0で割ることはできない */

    max_degs(b, &bxmax, &bymax);

    for (int i = 0; i <= MAX_DEG; ++i)
        for (int j = 0; j <= MAX_DEG; ++j)
            if (b->c[i][j] != 0) {
                terms[nt].i = i;
                terms[nt].j = j;
                terms[nt].c = b->c[i][j];
                ++nt;
            }

    const int inv_bc = mod_pow(bc, P - 2);

    for (int k = W * W - 1; k >= 0; --k) {
        const int ax = k / W;
        const int ay = k % W;
        const int ac = work.c[ax][ay];

        if (ac == 0)
            continue;

        if (ax >= bx && ay >= by) {
            const int dx = ax - bx;
            const int dy = ay - by;
            const int factor = mod(ac * inv_bc);

            /* bx は x方向の最大次数でもあるので、x側は必ず収まる */
            if (bxmax + dx > MAX_DEG || bymax + dy > MAX_DEG)
                return 0;

            add_term(q, factor, dx, dy);

#pragma omp parallel for schedule(static) if (nt >= PAR_DIV_MIN)
            for (int t = 0; t < nt; ++t) {
                int *w = &work.c[terms[t].i + dx][terms[t].j + dy];
                *w = mod(*w - factor * terms[t].c);
            }
            /* 先頭項は (ax, ay) で打ち消されて 0 になる */
        } else {
            /* 先頭項が割れない場合は余りへ移す */
            add_term(r, ac, ax, ay);
            work.c[ax][ay] = 0;
        }
    }
    return 1;
}

int main(void)
{
    Poly a, b, result;

    printf("F_17[x, y] の多項式演算\n");
    printf("係数は mod 17 で計算します。各入力多項式の次数は最大%dです。\n",
           MAX_INPUT_DEG);

    if (!read_poly("A", &a) || !read_poly("B", &b))
        return 1;

    printf("\nA = ");
    print_poly(&a);
    printf("B = ");
    print_poly(&b);

    printf("\n演算を選んでください: 1=加算 2=減算 3=乗算 4=除算\n選択: ");
    int op;
    if (scanf("%d", &op) != 1 || op < 1 || op > 4) {
        fprintf(stderr, "演算の選択が不正です。\n");
        return 1;
    }

    switch (op) {
    case 1:
        result = add_poly(&a, &b);
        printf("A + B = ");
        print_poly(&result);
        break;
    case 2:
        result = sub_poly(&a, &b);
        printf("A - B = ");
        print_poly(&result);
        break;
    case 3:
        if (!mul_poly(&a, &b, &result))
            return 1;
        printf("A * B = ");
        print_poly(&result);
        break;
    case 4: {
        Poly q, r;
        if (!div_poly(&a, &b, &q, &r)) {
            fprintf(stderr, "除算に失敗しました（Bが0、または次数上限超過）。\n");
            return 1;
        }
        printf("商 Q = ");
        print_poly(&q);
        printf("余り R = ");
        print_poly(&r);
        break;
    }
    }
    return 0;
}