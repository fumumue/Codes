#include <stdio.h>
#include <string.h>

#define P 17
#define MAX_DEG 64
#define MAX_INPUT_DEG 32

typedef struct {
    int c[MAX_DEG + 1][MAX_DEG + 1]; /* c[xの次数][yの次数] */
} Poly;

static int mod(int a)
{
    a %= P;
    return a < 0 ? a + P : a;
}

static int mod_pow(int a, int n)
{
    int result = 1;
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

    /* 辞書式順序で大きい項から表示 */
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
                printf("*x%s", i == 1 ? "" : "");
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

static Poly add_poly(const Poly *a, const Poly *b)
{
    Poly r = {0};
    for (int i = 0; i <= MAX_DEG; ++i)
        for (int j = 0; j <= MAX_DEG; ++j)
            r.c[i][j] = mod(a->c[i][j] + b->c[i][j]);
    return r;
}

static Poly sub_poly(const Poly *a, const Poly *b)
{
    Poly r = {0};
    for (int i = 0; i <= MAX_DEG; ++i)
        for (int j = 0; j <= MAX_DEG; ++j)
            r.c[i][j] = mod(a->c[i][j] - b->c[i][j]);
    return r;
}

static Poly mul_poly(const Poly *a, const Poly *b)
{
    Poly r = {0};
    for (int i = 0; i <= MAX_DEG; ++i) {
        for (int j = 0; j <= MAX_DEG; ++j) {
            if (a->c[i][j] == 0)
                continue;
            for (int k = 0; k <= MAX_DEG; ++k) {
                for (int l = 0; l <= MAX_DEG; ++l) {
                    if (b->c[k][l] == 0)
                        continue;
                    if (i + k > MAX_DEG || j + l > MAX_DEG) {
                        fprintf(stderr, "結果の次数が上限を超えました。\n");
                        return r;
                    }
                    r.c[i + k][j + l] =
                        mod(r.c[i + k][j + l] + a->c[i][j] * b->c[k][l]);
                }
            }
        }
    }
    return r;
}

/* 先頭単項式を辞書式順序（x優先）で取得 */
static int leading_term(const Poly *p, int *xdeg, int *ydeg, int *coeff)
{
    for (int i = MAX_DEG; i >= 0; --i) {
        for (int j = MAX_DEG; j >= 0; --j) {
            if (p->c[i][j] != 0) {
                *xdeg = i;
                *ydeg = j;
                *coeff = p->c[i][j];
                return 1;
            }
        }
    }
    return 0;
}

/* 単項式除法による多項式の筆算。a = q*b + r */
static int div_poly(const Poly *a, const Poly *b, Poly *q, Poly *r)
{
    Poly work = *a;
    int bx, by, bc;

    memset(q, 0, sizeof(*q));
    memset(r, 0, sizeof(*r));

    if (!leading_term(b, &bx, &by, &bc))
        return 0; /* 0で割ることはできない */

    const int inv_bc = mod_pow(bc, P - 2);

    while (leading_term(&work, &(int){0}, &(int){0}, &(int){0})) {
        int ax, ay, ac;
        leading_term(&work, &ax, &ay, &ac);

        if (ax >= bx && ay >= by) {
            int dx = ax - bx;
            int dy = ay - by;
            int factor = mod(ac * inv_bc);

            add_term(q, factor, dx, dy);

            for (int i = 0; i <= MAX_DEG; ++i) {
                for (int j = 0; j <= MAX_DEG; ++j) {
                    if (b->c[i][j] == 0)
                        continue;
                    if (i + dx > MAX_DEG || j + dy > MAX_DEG)
                        return 0;
                    work.c[i + dx][j + dy] =
                        mod(work.c[i + dx][j + dy] -
                            factor * b->c[i][j]);
                }
            }
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
        result = mul_poly(&a, &b);
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