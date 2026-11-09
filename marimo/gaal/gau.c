#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// #include "global-p.h"
// #include "struct-p.h"

// #define Pr 13
#define F K *E

short oinb(short o)
{
    if (ORD == 0)
        return 0;
    if (ORD == 1)
        return 1;
    for (int i = 0; i < Pr; i++)
    {
        for (int j = 0; j < Pr; j++)
        {
            if ((o * i) % Pr == 1)
                return i;
        }
    }
}

// short fu(short d)
//{
//     return (Pr - d) % Pr;
// }

void tet(CTX a, CTX inv_a)
{
    unsigned b[F][F] = {0};
    unsigned c[F][F] = {0};
    // CTX inv_a;   // ここに逆行列が入る
    unsigned buf; // 一時的なデータを蓄える

    printf("行列\n");
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            b[i][j] = a.x[i][j];
            printf("%d ", b[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    // 逆行列を出力
    printf("逆行列\n");
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            printf(" %d", inv_a.x[i][j]);
        }
        printf("\n");
    }
    // return inv_a;

    printf("検算\n");
    int tmp = 0;
    short z[F][F] = {0};
    // 検算
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            tmp = 0;
            for (int k = 0; k < F; k++)
                tmp += ((b[i][k] % Pr) * (inv_a.x[k][j] % Pr)) % Pr;
            z[i][j] = tmp % Pr;
            printf("%d ", z[i][j]);
        }
        printf("\n");
        // if(z[i][i]!=1)
        // exit(1);
    }
    exit(1);
    return inv_a;
    printf("もとの行列\n");
    int count = 0;
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            printf("%d ", b[i][j] % Pr);
        }
        printf("\n");
    }
    printf("\n");
    // int count = 0;
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            printf("%d", z[i][j] % Pr);
        }
        printf("\n");
        if (z[i][i] % Pr == 1)
            count++;
    }
}

CTX genS(CTX a)
{
    // short a.x[F][F] = {{0,2},{2,2}}; //{{2,1,1,1},{1,1,2,2},{2,0,1,1},{1,2,1,1}}; //入力用の配列
    CTX b = {0};
    unsigned c[F][F] = {0};
    CTX inv_a; // ここに逆行列が入る
    short buf; // 一時的なデータを蓄える
    // int i, j, k; // カウンタ
    // int n = F;   // 配列の次数

    srand(clock());
label:
    /*
    for (int i= 0; i < F; i++)
    {
        for (j = 0; j < F; j++)
            a.x[i][j] = rand() % Pr;
    }
    */
    b = a;
    // 単位行列を作る
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            inv_a.x[i][j] = (i == j) ? 1 : 0;
        }
    }
    // 掃き出し法
    for (int i = 0; i < F; i++)
    {
        printf("n=%d\n", i);

        if (a.x[i][i] == 0)
        {
            printf("aa.x[%d][%d]=%d\n", i, i, a.x[i][i]);
            int s = i;
            // short tmp[F] = {0}, tt[F] = {0};
            while (a.x[s][i] == 0)
            {
                s++;
                if (s == F)
                {
                    a.b = 0;
                    return a;
                }
            }
            for (int u = 0; u < F; u++)
            {
                a.x[i][u] = (a.x[i][u] + a.x[s][u]) % Pr;
                inv_a.x[i][u] = (inv_a.x[i][u] + inv_a.x[s][u]) % Pr;
            }
        }

        int s = 0;
        buf = oinb(a.x[i][i]);
        if (buf == 0)
        {
            printf("buf baka\n");
            exit(1);
        }
        printf("buf=%d\n", buf);
        // for(k=0;k<F;k++)
        if (buf > 0)
        {
            if (i == 0)
            {
                //buf = oinb(a.x[0][0]);
                for (int j = 0; j < F; j++)
                {
                    a.x[0][j] = (a.x[0][j] * buf) % Pr;
                    inv_a.x[0][j] = (inv_a.x[0][j] * buf) % Pr;
                }
                for (int jj = 1; jj < F; jj++)
                {
                    buf = a.x[jj][0];
                    for (int j = 0; j < F; j++)
                    {
                        a.x[jj][j] = (Pr + a.x[jj][j] - a.x[0][j] * buf) % Pr;
                        inv_a.x[jj][j] = (Pr + inv_a.x[jj][j] - inv_a.x[0][j] * buf) % Pr;
                    }
                    printf("\n");
                }
            }
            for (int j = 0; j < F; j++)
            {
                for (int i = 0; i < F; i++)
                    printf("%d,", a.x[j][i]);
                printf("\n");
            }
            // exit(1);
            if (i > 0)
            {   for(int k=0;k<F;k++){
                buf = oinb(a.x[k][k]);
                if(buf==0){
                int j=0;
                while(buf==0){
                buf=a.x[j++][k];
                }
                for(int i=0;i<F;i++)
                a.x[k][i]=a.x[k][i]+a.x[j][i];
                }
                buf=oinb(a.x[k][k]);
                for (int j = 0; j < F; j++)
                {
                    a.x[k][j] = a.x[k][j] * buf % Pr;
                    inv_a.x[k][j] = inv_a.x[k][j] * buf % Pr;
                     printf("%d ",a.x[k][j]);
                }
                printf("\n");
                }
                // exit(1);
                for (int j = 0; j < F; j++)
                {
                    if (i != j && i > 0)
                    {
                        printf("%d", a.x[i][i]);
                        buf = (a.x[j][i]);
                        for (int k = 0; k < F; k++)
                        {
                            a.x[j][k] = (Pr + a.x[j][k] - (a.x[i][k] * buf) % Pr) % Pr;
                            // if(inv_a.x[i][k]!=0)
                            inv_a.x[j][k] = (Pr + inv_a.x[j][k] - (inv_a.x[i][k] * buf) % Pr) % Pr;
                            // printf("%d,",inv_a.x[j][k]);
                        }
                        // printf("\n");
                    }
                }
                inv_a.b = 1;
            }
        }
    }
    unsigned cd = 0, zz[F][F] = {0};
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            cd = 0;
            for (int k = 0; k < F; k++)
                cd += (a.x[i][k] * inv_a.x[k][j]) % Pr;
            zz[i][j] = cd % Pr;
            printf("%d,", a.x[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    // exit(1);
    //  逆行列を出力
    printf("逆行列\n");
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            if (inv_a.x[i][j] < 0)
                inv_a.x[i][j] += Pr;
            // printf(" %d", inv_a.x[i][j]);
        }
        // printf("\n");
    }
    /*
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
            printf("%d ", a.x[i][j]);
        printf("\n");
    }
     exit(1);
    // return inv_a;
    */
    printf("検算\n");
    unsigned tmp = 0;
    unsigned z[F][F] = {0};
    int count = 0;
    int ca = 0;

    // 検算
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            tmp = 0;
            for (int k = 0; k < F; k++)
                tmp += ((b.x[i][k]) * (inv_a.x[k][j])) % Pr;
            z[i][j] = tmp % Pr;
            printf("%d ", z[i][j]);
            if (z[i][j] == 0 && i != j)
            {
                ca++;
            }
        }
        printf("\n");
        if (z[i][i] == 1)
        {
            count++;
        }
        else if (z[i][i] != 1)
        {
            inv_a.b = 0;
            printf("baka\n");
            // exit(1);
            return inv_a;
        }
    }
    exit(1);
    if (ca == F * F - F && count == F)
    {
        inv_a.b = 1;
        return inv_a;
    }
    else
    {
        exit(1);
        inv_a.b = 0;
        return inv_a;
    }

    printf("もとの行列\n");
    exit(1);
    count = 0;
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            printf("%d ", b.x[i][j] % Pr);
        }
        printf("\n");
    }
    printf("\n");
    // int count = 0;
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            printf("%d", z[i][j] % Pr);
        }
        printf("\n");
        if (z[i][i] % Pr == 1)
            count++;
    }

    if (count != F)
    {
        inv_a.b = 0;
        return inv_a;
    }

    exit(1);
    inv_a.b = 1;
    return inv_a;
}

CTX kenzan(CTX a, CTX inv_a, int p)
{
    unsigned tmp = 0;
    static CTX z = {0};
    // 検算
    if (p == 1)
    {
        for (int i = 0; i < F; i++)
        {
            for (int j = 0; j < F; j++)
            {
                tmp = 0;
                for (int k = 0; k < F; k++)
                    tmp += ((a.x[i][k] % Pr) * (inv_a.x[k][j] % Pr)) % Pr;
                z.x[i][j] = tmp % Pr;
            }
        }
        /*
        //return z;
        printf("もとの行列\n");
      int count = 0;
      for (int i = 0; i < F; i++)
      {
        for (int j = 0; j < F; j++)
        {
          printf("%d ", a.x[i][j] % Pr);
        }
        printf("\n");
      }
      printf("逆行列\n");
      for (int i = 0; i < F; i++)
      {
        for (int j = 0; j < F; j++)
        {
          printf("*%d ", inv_a.x[i][j] % Pr);
        }
        printf("\n");
      }

      printf("検算\n");
      for (int i = 0; i < F; i++)
      {
        for (int j = 0; j < F; j++)
        {
          printf("z%d ", z.x[i][j] % Pr);
        }
        printf("\n");
      }
      */
        return z;
    }

    if (p == 2)
    {
        for (int i = 0; i < F; i++)
        {
            for (int j = 0; j < M; j++)
            {
                tmp = 0;
                for (int k = 0; k < F; k++)
                    tmp += ((a.x[i][k] % Pr) * (inv_a.x[k][j] % Pr)) % Pr;
                z.x[i][j] = tmp % Pr;
                // printf("%d ", z.x[i][j]);
            }
            // printf("\n");
        }
        /*
      // int count = 0;
      printf("かけた結果\n");
      for (int i = 0; i < F; i++)
      {
        for (int j = 0; j < M; j++)
        {
          printf("%d, ", z.x[i][j] % Pr);
        }
        printf("\n");
      }
      */
        return z;
    }
}

/*
int main()
{
    MTX a = {0}; //{{0,2},{2,2}}; //{{2,1,1,1},{1,1,2,2},{2,0,1,1},{1,2,1,1}}; //入力用の配列
    short b[F][F] = {0};
    short c[F][F] = {0};
    MTX inv_a = {0}; // ここに逆行列が入る
    short buf;               // 一時的なデータを蓄える
    int i, j, k;             // カウンタ
    int n = F;               // 配列の次数
    MTX S, inv_S;

    srand(clock());
label:
do{
        for (int i= 0; i < F; i++)
        {
            for (j = 0; j < F; j++)
                a.x[i][j] = rand() % Pr;
        }
        inv_a=genS(a);
}while(inv_a.col!=1);
exit(1);

    // 単位行列を作る
    for (int i= 0; i < F; i++)
    {
        for (j = 0; j < F; j++)
        {
            b[i][j] = a[i][j];
            inv_a[i][j] = (i == j) ? 1 : 0;
        }
    }
    // 掃き出し法
    for (int i= 0; i < F; i++)
    {
        printf("n=%d\n", i);
    lulu:
        if (a[i][i] == 0)
        {
            printf("aa[%d][%d]=%d\n", i, i, a[i][i]);
            int s = i;
            short tmp[F] = {0}, tt[F] = {0};
            while (a[s][i] == 0)
            {
                s++;
                if (s == F)
                    goto label;
            }
            for (int u = 0; u < F; u++)
            {
                a[i][u] = (a[i][u] + a[s][u]) % Pr;
                inv_a[i][u] = (inv_a[i][u] + inv_a[s][u]) % Pr;
            }
            printf("ss=%d %d %d\n", i, s, a[i][i]);
            // i++;
            if (a[i][i] == 0)
                goto lulu;
            // exit(1);
        }
        // i++;
        buf = oinv(a[i][i]);
        if (buf == 0)
        {
            printf("a[%d][%d]=%d\n", i, i, a[i][i]);
            // buf=0;
            // a[i][i]=1;
            // printf("baka\n");
            // exit(1);
            goto label;
        }
        printf("buf=%d\n", buf);
        // for(k=0;k<F;k++)
        {
            for (j = 0; j < F; j++)
            {
                a[i][j] = a[i][j] * buf % Pr;
                inv_a[i][j] = inv_a[i][j] * buf % Pr;
                // printf("%d %d %d\n",i,j,a[i][j]);
            }
        }
        // exit(1);

        for (j = 0; j < F; j++)
        {
            if (i != j)
            {
                buf = a[j][i];
                for (k = 0; k < F; k++)
                {
                    a[j][k] = (a[j][k] - (a[i][k] * buf) % Pr) % Pr;
                    // if(inv_a[i][k]!=0)
                    inv_a[j][k] = (inv_a[j][k] - (inv_a[i][k] * buf) % Pr) % Pr;
                }
            }
        }
        int count = 0;
        for (k = 0; k < F; k++)
        {
            for (int t = 0; t < F; t++)
                printf("a%d ", a[k][t] % Pr);
            for (int t = 0; t < F; t++)
                printf("i%d ", (inv_a[k][t] + Pr) % Pr);
            printf("\n");
            if (a[k][k] % Pr == 1)
                count++;
        }
        printf("\n");
        if (count == F)
        {
            int cc = 0;
            for (int i = 0; i < F; i++)
            {
                for (int j = 0; j < F; j++)
                {
                    cc = 0;
                    for (int k = 0; k < F; k++)
                        cc += (b[i][k] * inv_a[k][j]) % Pr;
                    c[i][j] = (cc) % Pr;
                }
            }
            int c2 = 0;
            for (int i = 0; i < F; i++)
            {
                for (int j = 0; j < F; j++)
                    printf("%d ", (Pr + c[i][j]) % Pr);
                printf("\n");
                if (c[i][i] % Pr == 1)
                    c2++;
            }
            printf("\n");
        }
        for (j = 0; j < F; j++)
        {
            for (k = 0; k < F; k++)
            {
                if (a[j][k] < 0)
                    a[j][k] += Pr;
                // printf("%d ",a[j][k]);
            }
            for (k = 0; k < F; k++)
            {
                if (inv_a[j][k] < 0)
                    inv_a[j][k] += Pr;
                // printf(" %d,",inv_a[j][k]);
            }
            // printf("\n");
        }
        // exit(1);
    }

    // 逆行列を出力
    for (int i= 0; i < F; i++)
    {
        for (j = 0; j < F; j++)
        {
            printf(" %d", inv_a[i][j]);
        }
        printf("\n");
    }
    short tmp = 0;
    short z[F][F] = {0};
    // 検算
    for (int i= 0; i < F; i++)
    {
        for (j = 0; j < F; j++)
        {
            tmp = 0;
            for (k = 0; k < F; k++)
                tmp += ((b[i][k] % Pr) * (inv_a[k][j] % Pr)) % Pr;
            z[i][j] = tmp % Pr;
        }
    }
    int count = 0;
    for (int i= 0; i < F; i++)
    {
        for (j = 0; j < F; j++)
        {
            printf("%d,", b[i][j] % Pr);
        }
        printf("\n");
    }
    // int count = 0;
    for (int i= 0; i < F; i++)
    {
        for (j = 0; j < F; j++)
        {
            printf("%d,", z[i][j] % Pr);
        }
        printf("\n");
        if (z[i][i] % Pr == 1)
            count++;
    }

    if (count != F)
        goto label;

    return 0;
}
*/
