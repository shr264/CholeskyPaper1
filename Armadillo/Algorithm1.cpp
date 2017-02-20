#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
    arma_rng::set_seed(12345);
    int p = 4;
    int n = 5;
    mat A = randu<mat>(n,p);
    mat S = A.t()*A;
    umat G = abs(S)>0;
    G(2,0) = 0; G(0,2) = 0; G(2,1) = 0; G(1,2) = 0; G(3,1) = 0; G(1,3) = 0;
    
    umat swaps; swaps.zeros((p-1),2);
    for(int j = 0; j<(p-1); j ++){
        uvec colsums; colsums.zeros(p-j);
        for(int i = j; i<p ; i++){
            colsums(i-j) = sum(G(j,j,size(p-j,p-j)).col(i-j));
            cout << colsums << endl;
        }
        swaps(j,0) = colsums.min();
        swaps(j,1) = j;
    
        if(swaps(j,1)<swaps(j,0)){
            cout << "Swapping " << swaps(j,0)+1 << " with " << swaps(j,1)+1 << endl;
            G.swap_rows(swaps(j,0),swaps(j,1));
            G.swap_cols(swaps(j,0),swaps(j,1));
            S.swap_rows(swaps(j,0),swaps(j,1));
            S.swap_cols(swaps(j,0),swaps(j,1));
        }
    }
    
    cout << swaps << endl;
    cout << S << endl;
    cout << G << endl;
    
    return 0;
}