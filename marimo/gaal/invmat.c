//#include<iostream>
//#include<iomanip>
//#include<algorithm>
//#include<vector>
#include<math.h>
using namespace std;

#define D 3 //行列の次数を設定(今は３次元)

void Lu();
void inverse_L();
void inverse_U();

double a[D][D]; //逆行列を求めたい行列A
double l[D][D],u[D][D],q[D][D]; //LU分解で用いる行列
double l_[D][D],u_[D][D]; //下三角行列L,単位上三角行列Uの逆行列

int main(){
    /*
       a[D][D]の要素を決める
    */

    Lu();
    inverse_L();
    inverse_U();

    //Aの逆行列を計算
    vector<vector<double> > p(D,vector<double>(D,0));
    vector<vector<double> > a_(D,vector<double>(D,0)); //Aの逆行列を格納
    for(int i=0;i<D;i++){
        for(int j=0;j<D;j++){
            for(int k=0;k<D;k++){
                p[i][j]+=u_[i][k]*l_[k][j];
            }
        }
    }
    for(int i=0;i<D;i++){
        for(int j=0;j<D;j++){
            for(int k=0;k<D;k++){
                a_[i][j]+=p[i][k]*q[k][j];
            }
        }
    }
    //Aの逆行列を出力
    for(int i=0;i<D;i++){
        for(int j=0;j<D;j++){
            cout << a_[i][j] << ' ';
        }
        cout << endl;
    }

}

//LU分解する関数
void Lu(){
    // ピボット選択
    for(int i=0;i<D;i++){ for(int j=0;j<D;j++) u[i][j]=a[i][j];}
    for(int i=0;i<D;i++) q[i][i]=1.0;
    for(int k=0;k<D-1;k++){
        double hold_val=0;
        int hold_index=0;
        for(int i=k;i<D;i++){
            if(hold_val<abs(u[i][k])){
                hold_val=abs(u[i][k]);
                hold_index=i;
            }
        }

        if(hold_index!=k){
            for(int i=0;i<D;i++) swap(u[hold_index][i],u[k][i]);
            for(int i=0;i<D;i++) swap(l[hold_index][i],l[k][i]);
            for(int i=0;i<D;i++) swap(q[hold_index][i],q[k][i]);            
        }

        //上三角行列Lと単位下三角行列Uの構成
        for(int i=0;i<k;i++) l[i][k]=0;
        l[k][k]=1.0;

        for(int i=k+1;i<D;i++){
            double cat=u[i][k]/u[k][k];
            l[i][k]=cat;
            for(int j=0;j<D;j++) u[i][j]-=u[k][j]*cat;
        }
    }  
    l[D-1][D-1]=1.0;
}

//下三角行列Lの逆行列を求める
void inverse_L(){
    for(int k=0;k<D;k++){
        vector<double> x(D,0);
        x[k]=1.0;
        for(int i=k+1;i<D;i++){
            double tmp=0;
            for(int j=0;j<i;j++){
                tmp+=l[i][j]*x[j];
            }
            x[i]=-tmp;
        }  
        for(int i=0;i<D;i++) l_[i][k]=x[i];   
    }
}

//単位上三角行列Uの逆行列を求める
void inverse_U(){
    for(int k=D-1;k>=0;k--){
        vector<double> x(D,0);
        x[k]=1.0/u[k][k];
        for(int i=k-1;i>=0;i--){
            double tmp=0;
            for(int j=D-1;j>i;j--){
                tmp+=u[i][j]*x[j];
            }
            x[i]=-tmp/u[i][i];
        }  
        for(int i=0;i<D;i++) u_[i][k]=x[i];
    }
}
