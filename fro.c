// ---------- p 乗根 mod g (g は無平方) ----------
// F_P 係数では (Σ y_i x^i)^P = Σ y_i x^(iP)。Frobenius 行列 F (列 i = x^(iP) mod g) の
// 逆行列を掛ければ P 乗根になる。g が無平方 <=> F が可逆。
static int Finv[MAXK + 2][MAXK + 2];
static int tdeg;

static int build_frobenius_inverse(poly g)
{
    int t = g.d;
    static int A[MAXK + 2][2 * (MAXK + 2)];
    poly xp = ppowmod(pshift(pconst(1), 1), P, g);   // x^P mod g
    poly col = pconst(1);                            // x^(iP) mod g, i = 0 から
    for (int i = 0; i < t; i++)
    {
        for (int j = 0; j < t; j++)
        {
            A[j][i] = (j <= col.d) ? col.c[j] : 0;   // 左半分: F の第 i 列
            A[j][t + i] = (i == j);                  // 右半分: 単位行列
        }
        col = pmulmod(col, xp, g);
    }
    for (int i = 0; i < t; i++)                      // ガウス・ジョルダン法で [F | I] -> [I | F^-1]
    {
        int p = -1;
        for (int j = i; j < t; j++)
            if (A[j][i]) { p = j; break; }
        if (p < 0)
            return 0; // 特異 = g が平方因子を持つ
        if (p != i)
            for (int k = 0; k < 2 * t; k++)
            {
                int tmp = A[i][k]; A[i][k] = A[p][k]; A[p][k] = tmp;
            }
        int iv = inv(A[i][i]);
        for (int k = 0; k < 2 * t; k++)
            A[i][k] = md((long)A[i][k] * iv);
        for (int j = 0; j < t; j++)
            if (j != i && A[j][i])
            {
                int f = A[j][i];
                for (int k = 0; k < 2 * t; k++)
                    A[j][k] = md(A[j][k] - (long)f * A[i][k]);
            }
    }
    for (int i = 0; i < t; i++)
        for (int j = 0; j < t; j++)
            Finv[i][j] = A[i][t + j];
    tdeg = t;
    return 1;
}

static poly proot(poly z)                            // z^(1/P) mod g = Finv * (z の係数ベクトル)
{
    poly r = pzero();
    for (int i = 0; i < tdeg; i++)
    {
        long s = 0;
        for (int j = 0; j <= z.d && j < tdeg; j++)
            s += (long)Finv[i][j] * z.c[j];
        r.c[i] = (int)(s % P);
    }
    r.d = tdeg - 1;
    ptrim(&r);
    return r;
}