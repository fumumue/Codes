/*============================================================================
Project  : C Programing Tools
File     : mdump.c
Function : Memory Dump, debugging aid function
Revision : 1.00 2016/05/14 Created by A.Niida
Copyright(c) 2016 - 2020 A.Niida All rights reserved.
============================================================================*/
/*======================================================================
    プロジェクト名  ：C言語プログラミング
    ファイル名      ：mdump.c
    機能            ：メモリダンプ、デバッグ補助
    修正履歴        ：1.00    2016/05/14    新井田    作成
    Copyright(c) 2016 - 2020 A.Niida All Rights Reserved.
======================================================================*/

#include <stdio.h>
#include <ctype.h>      /* for isprint() */

/*============================================================================
Func Name   : void mdump(void *vp, int size, char *msg)
Function    : Dump memory
Param Input : void *vp   = memory address
              int size   = size to be dumped, 16 bytes align
              char *msg  = message to identify memory block
Param Output: None
Return      : None
Input Inf   : None
Output Inf  : None
Note        : How to compile with your 'prog.c',
                >  cl prog.c mdump.c
                $ gcc prog.c mdump.c
Revision    : 1.00 2016/05/14 Created by A.Niida
============================================================================*/
/*======================================================================
    関数名          ：void mdump(void *vp, int size, char *msg);
    機能            ：指定されたアドレスから、サイズ分、ダンプする
    入力引数説明    ：void *vp -- 変数（メモリ）のアドレス
                    ：int size -- ダンプする領域のバイト数。16バイト単位で表示
                    ：char *msg -- 先頭に表示する識別文字列。NULLなら空行のみ。
    出力引数説明    ：無し
    戻り値          ：無し
    入力情報        ：無し
    出力情報        ：無し
    特記事項        ：無し
    修正履歴        ：1.00 2016/05/14 新井田    作成
                    ：1.10 2020/11/01 dumpHexdec(), dumpAsciic() を独立
======================================================================*/
#define LINESIZE    16                    // 16バイトを一行に表示

static void dumpHexdec(unsigned char *mem, int ltop, int lnext);
static void dumpAsciic(unsigned char *mem, int ltop, int lnext);

void mdump(void *vp, int size, char *msg)
{
    unsigned char *mem = vp;  // 常に無符号char領域としてアクセス
    int ltop, lnext;    /* line-top and line-next */

    // 識別メッセージを表示
    if (msg != NULL) {        // つまり、NULLなら表示しない
        printf(" '%s'", msg);
    }
    printf("\n");             // でも最初に改行だけはする

    for (ltop = 0; ltop < size; ltop = lnext) {    // LINESIZE 毎に繰返す
        lnext = ltop + LINESIZE;            // 次の行の先頭位置

        printf("%p: ", &mem[ltop]);         // 今の行の先頭アドレス
        dumpHexdec(mem, ltop, lnext);       // メモリの値を16進数で表示
        dumpAsciic(mem, ltop, lnext);       // ASCII文字として
        printf("\n");
    }
}

static void dumpHexdec(unsigned char *mem, int ltop, int lnext)
{
    for (int i = ltop; i < lnext; i++)      // メモリの値を16進数で表示
        printf("%02x ", mem[i]);
}

static void dumpAsciic(unsigned char *mem, int ltop, int lnext)
{
    for (int i = ltop; i < lnext; i++) {      // 文字表示、ascii dump
        unsigned char c = mem[i];             // メモリの値
        printf("%c", isprint(c) ? c : '.');   // 表示不能なら '.' を
    }
}

#ifdef SAMPLE
/*
*    サンプルプログラムで使用法を示す
*        コンパイル方法 > cl mdump.c /D SAMPLE
*                      $ gcc mdump.c -D SAMPLE
*/
#define BUFFLEN    50        /* used be 40, 100 */
int main(void)
{
    int arr[6] = { 1, 10, 100, 1000, 10000 };
    char buff[BUFFLEN]  /* = { 0x31, 0xC3, 0x85, 0x61 }*/ ;
    char *cp;

    // mdump(stdin, sizeof(FILE) * 4, "_iob[0..3]"); // これは何???
    mdump("Mdump!\n", 3, "literal");
    do {
        mdump(arr, sizeof(arr), "int arr[6]");
        mdump(buff, BUFFLEN, NULL);
        mdump(&cp, sizeof(char*), "cp, char pointer");

        printf("文字列を入力してください: ");
        cp = fgets(buff, BUFFLEN, stdin);    // （正常終了が前提）
    } while (buff[0] != '\n');

    mdump(buff, 5, "");
    return 0;
}
#endif