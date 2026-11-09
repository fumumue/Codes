#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include "global.h"
#include "struct.h"

void cyclic_convolution(
    const vec a,
    const vec b,
    vec *c,
    int n
) {
    for (int k = 0; k < n; k++) {
        c->x[k] = 0;

        for (int i = 0; i < n; i++) {
            c->x[k] ^= (a.x[i] & b.x[(k - i + n) % n]);
        }
    }
}

vec aa(int n){
    vec a = {0};
    int count = 0;

    while(count < n){
    {
        int i=rand()%8192;
        if(a.x[i] == 0){
                a.x[i] = 1;
                count++;
            }
            if(count == n)    break;

        }
    }

    return a;
}


int wt(vec a){
    int weight = 0;
    for (int i = 0; i < 8192; i++) {
        if (a.x[i] != 0) {
            weight++;
        }
    }
    return weight;
}


int main(void) {
    int n = 8192;

    vec a = {0};
    vec b = {0};
    vec c = {0};

    srand(time(NULL));
    a=aa(1024);
    b=aa(1024);
    cyclic_convolution(a, b, &c, n);
    int w=wt(c);
    printf("Weight of a, b, c: %d\n", wt(a), wt(b), w);
    vec v=aa(675);
    for (int i = 0; i < n; i++) {
        v.x[i]^=c.x[i];
    }
    printf("w= %d\n ",wt(v));
    exit(1);

    printf("Result: ");
    for (int i = 0; i < n; i++)
        printf("%d ", c.x[i]);

    printf("\n");
    return 0;
}