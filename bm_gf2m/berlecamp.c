#include <string.h>

// 多項式 f に体の元 c をかける(デバッグ出力なし)
static vec vscale(unsigned short c, vec f)
{
    int i, k = deg(f);
    vec h = {0};

    if (c == 0)
        return h;
    for (i = 0; i <= k; i++)
        if (f.x[i] > 0)
            h.x[i] = gf[mlt(fg[c], fg[f.x[i]])];
    return h;
}

// f * x^n
static vec vshift(vec f, int n)
{
    int i, k = deg(f);
    vec h = {0};

    for (i = 0; i <= k && i + n < DEG; i++)
        h.x[i + n] = f.x[i];
    return h;
}

// Berlekamp-Massey(0オリジン: s[0..2T-1])
vec bma(unsigned short s[])
{
    int L = 0, m = -1, k, i;
    unsigned short d, bd = 1, coef;
    vec f = {0}, g = {0}, h;

    f.x[0] = g.x[0] = 1;

    for (k = 0; k < 2 * T; k++)
    {
        // d = Σ_{i=0}^{L} f_i * s[k-i]
        d = 0;
        for (i = 0; i <= L; i++)
            d ^= gf[mlt(fg[f.x[i]], fg[s[k - i]])];

        if (d == 0)
            continue;

        h = f;
        coef = gf[mlt(fg[d], oinv(bd))]; // d / bd (oinv は対数の添字を返す)
        f = vadd(f, vshift(vscale(coef, g), k - m));

        if (2 * L <= k)
        {
            L = k + 1 - L;
            m = k;
            g = h;
            bd = d;
        }
    }

    if (deg(f) != T)
        fprintf(stderr, "bma: deg=%d (expected %d)\n", deg(f), T);
    return f;
}

// Berlekamp-Massey(1オリジン: s[1..2T])
vec bms(unsigned short s[], int kk)
{
    int i, j, L = 0;
    unsigned short d;
    vec C = {0}, B = {0}, old;

    (void)kk;
    C.x[0] = 1;
    B.x[0] = 1;

    for (j = 1; j <= 2 * T; j++)
    {
        d = s[j];
        for (i = 1; i <= L; i++)
            d ^= gf[mlt(fg[C.x[i]], fg[s[j - i]])];

        if (d == 0)
        {
            B = vshift(B, 1);
        }
        else if (2 * L > j - 1)
        {
            C = vadd(C, vshift(vscale(d, B), 1));
            B = vshift(B, 1);
        }
        else
        {
            old = C;
            C = vadd(C, vshift(vscale(d, B), 1));
            B = vscale(gf[oinv(d)], old); // B = old / d
            L = j - L;
        }
    }

    if (deg(C) != T)
        fprintf(stderr, "bms: deg=%d (expected %d)\n", deg(C), T);
    return C;
}