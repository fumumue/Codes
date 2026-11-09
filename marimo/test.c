#include "inv_mat.c"

int main(){
MTX A;
MTX B;
MTX C;
double D[8][8];
double W[8][8];

for(int i=0;i<8;i++){
for(int j=0;j<8;j++){
W[i][j]=(double)((1129+i+j)%1129);
B.x[i][j]=(1129+i-j)%1129;
}
}
double AA[8][8] = 
                 {{1.0, 2.0, 3.0,1.0,1.0, 2.0, 3.0,1.0},
                 {4.0, 5.0, 6.0,2.0,4.0, 5.0, 6.0,2.0},
                 {7.0, 8.0, 10.0,3.0,7.0, 8.0, 10.0,3.0},
                 {1.0,2.0,3.0,4.0,1.0,2.0,3.0,4.0},
                 {1.0, 2.0, 3.0,1.0,1.0, 2.0, 3.0,1.0},
                 {4.0, 5.0, 6.0,2.0,4.0, 5.0, 6.0,2.0},
                 {7.0, 8.0, 10.0,3.0,7.0, 8.0, 10.0,3.0},
                 {1.0,2.0,3.0,4.0,1.0,2.0,3.0,4.0}};
                
    task();
    inverseMatrix(AA,D);
    printf("Inverse Matrix:\n");
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            printf("%f\t", D[i][j]);
        }
        printf("\n");
    }
exit(1);



for(int i=0;i<8;i++){
for(int j=0;j<8;j++)
printf("%lf,",D[i][j]);
printf("\n");
}
//print_matrix((int)D);
//mulmat(A,B,3);

}