#include <stdio.h>
int main(void)
{
    int i, sum;    // 初期値無しで変数定義する
    printf("i = %d, sum = %d\n", i, sum);   // #1

    // それぞれに初期値を与える
    i = 999;
    sum = 0;
    printf("i = %d, sum = %d\n", i, sum);   // #2

    // for を使って 0～10の和（55になるはず）を計算させる
    for (i = 0; i <= 10; i++) {
        sum = sum + i;
    }

    // 結果を見る
    printf("i = %d, sum = %d\n", i, sum);   // #3
    return 0;
}