#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"

//#include "8192.h"

unsigned char tmp[N][E * K] = {0};
//unsigned char pub[E * K][N] = {0};
//unsigned char BH[E * K][N] = {0};
unsigned short mat[N][K*E] = {0};
unsigned short ma[N][K*E] = {0};
unsigned short bm[N][K*E]={0};
unsigned short bm2[N][K*E]={0};

//unsigned short syn[K]={0};
unsigned short P[N] = {0};
unsigned short inv_P[N] = {0};
unsigned short uu = 0;
