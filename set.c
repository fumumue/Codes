#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main(){
    int i,j,k,count=0,a[64][64]={0},b[64][64]={0},c[64][64]={0};

    for(i=0;i<64;i++){
        for(j=0;j<64;j++){
            a[i][j]=rand()%2;
            b[i][j]=rand()%2;
        }
    }
    for(i=0;i<64;i++){
        for(j=0;j<64;j++){
            for(k=0;k<64;k++)
            c[i][j]^=a[i][k]&b[k][j];
            printf("%d ",c[i][j]);
        }
        printf("\n");
    }
return 0;
}