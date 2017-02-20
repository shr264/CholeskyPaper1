#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

uvec setdiff(uvec x, uvec y){
    x = unique(x);
    y = unique(y);
    for (size_t j = 0; j < y.n_elem; j++) {
        uvec q1 = find(x == y[j]);
        if (!q1.empty()) {
            x.shed_row(q1(0));
        }
    }
    return x;
}

int main()
{
    arma_rng::set_seed(12345);
    int p = 6;
    int n = 7;
    mat A = randu<mat>(n,p);
    mat S = A.t()*A;
    umat G = abs(S)>0;
    G(2,0) = 0; G(0,2) = 0; G(2,1) = 0; G(1,2) = 0;
    G(3,1) = 0; G(1,3) = 0; G(3,2) = 0; G(2,3) = 0;
    G(4,2) = 0; G(2,4) = 0; G(4,3) = 0; G(3,4) = 0;
    G(5,1) = 0; G(1,5) = 0; G(5,3) = 0; G(3,5) = 0; G(5,0) = 0; G(0,5) = 0;
    
    uvec colsums; colsums.zeros(p);
    uvec order; order.zeros(p);
    
    int count = 0;
    
    for(uword i = 0; i<p ; i++){
        order(i) = i;
    }
    
    for(uword i = 0; i<p ; i++){
        colsums(i) = sum(G.col(i));
        cout << colsums << endl;
    }
    
    uvec neworder(p); neworder.fill(p);
    
    uword ii = colsums.index_min();
    uvec ngbrs = find(G.col(ii)>0);
    umat newG = G.submat(ngbrs, ngbrs);
    
    colsums.set_size(ngbrs.size());
    
    for(uword i = 0; i < ngbrs.size() ; i++){
        colsums(i) = sum(newG.col(i));
        cout << colsums << endl;
    }
    
    uword tempsize1 = 0;
    uword tempsize2 = ngbrs.size();
    uvec temporder = ngbrs(sort_index(colsums));
    
    for(uword i = tempsize1; i < tempsize1+tempsize2; i++){
        neworder(i+tempsize1) = temporder(i);
    }
    tempsize1 += tempsize2;
    
    uvec ndslft = setdiff(order,neworder);
    
    while(tempsize1<(p-1)){
        
        umat newGG = G.submat(ndslft,ndslft);
        
        cout << "Press Enter to Continue";
        cin.ignore();
        
        colsums.set_size(ndslft.size());
        for(uword i = 0; i<ndslft.size() ; i++){
            colsums(i) = sum(newGG.col(i));
            cout << colsums << endl;
        }
        
        ii = colsums.index_min();
        uword newsize = sum(newGG.col(ii)>0);
        ngbrs.set_size(newsize);
        ngbrs = find((G.submat(ndslft,ndslft))>0);
        newG.set_size(newsize,newsize);
        newG = G.submat(ngbrs, ngbrs);
        
        colsums.set_size(newsize);
        
        for(uword i = 0; i < newsize ; i++){
            colsums(i) = sum(newG.col(i));
            cout << colsums << endl;
        }
        
        tempsize2 += newsize;
        temporder.set_size(newsize);
        temporder = ngbrs(sort_index(colsums));
        
        for(uword i = tempsize1; i < tempsize1+tempsize2; i++){
            neworder(i+tempsize1) = temporder(i);
        }
        tempsize1 += tempsize2;
        
        uvec ndslft = setdiff(order,neworder);
    
    }
    
    cout << neworder << endl;
    cout << ngbrs << endl;
    cout << colsums << endl;
    
    cout << S << endl;
    cout << G << endl;
    
    return 0;
}


