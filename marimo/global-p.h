#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include "8192.h"

// 符号のパラーメータの指定。通常[N,K,T]として、
// Nは符号の長さ、Kが符号の次元、Tは訂正エラー数
// を表す。ここではDは符号長にしている。
#define N 2187 // 元の数
#define M  N // 符号長　M<=N 3072
#define K (16) // 符号の次元
#define DEG (K*2+1)
#define T (K/2) // エラーの数
#define E (7)     // 拡大体の拡大次数
// #define D (2187) //符号長（短縮符号）
#define F K *E // 2040
#define EXP E  // degree
#define Pr 3   // 基礎体
// #define O 2197 // 1331,2197,4913,6859,3125,2187,19683
#define ORD N // 1331,2197,4913,6859,3125,2187,19683,29791,16807,24389
#define O ORD


#define V 3 // 変数の数
//#define P 1024
//#define Q 5 // 基礎体
// #define N Q *Q *Q // 定義体
#define I Q + 1 // 曲線の次数
//#define J 3
// #define K I-2			//number h
//#define H (K + 1) * (K + 2) / 2 // シンドローム行列の横ベクトルの長さ
// #define F (J-K+1)*(J-K+2)/2	//シンドローム行列の縦ベクトル
#define U 26
// #define E 6

