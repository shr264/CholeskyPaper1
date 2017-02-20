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
    G(5,1) = 0; G(1,5) = 0; G(5,3) = 0; G(3,5) = 0; G(5,4) = 0; G(4,5) = 0;
    
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
    
    uvec neworder; neworder.zeros(p);
    
    uword ii = colsums.index_min();
    uvec ngbrs = find(G.col(ii)>0);
    uvec temporder = ngbrs(sort_index(colsums(ngbrs)));
    
    uword tempsize1 = temporder.size() ;
    
    for(uword i = 0; i<tempsize1; i++){
        neworder(i) = temporder(i);
    }
    count = tempsize1;
    uvec newngbrs;
    umat newG;
    uword tempsize2 = 0;
    
    cout << S << endl;
    cout << G << endl;
    
    while(count<p){
    
        newngbrs = setdiff(order,neworder);
        colsums.set_size(newngbrs.size());
        newG = G.submat(order,newngbrs);
        
        uvec tempngbrs(newngbrs.size());
        
        for(uword i = 0; i<newngbrs.size() ; i++){
            colsums(i) = sum(newG.col(i));
            cout << colsums << endl;
            tempngbrs(i) = i;
        }
    
        cout << ngbrs << endl;
        cout << colsums(tempngbrs) << endl;
        //cout << sort_index(colsums(tempngbrs)) << endl;
        //cout << ngbrs(sort_index(colsums(tempngbrs))) << endl;
        
        cout << "Press Enter to Continue";
        cin.ignore();
        
        
        ii = colsums.index_min();
        ngbrs = find(newG.col(ii)>0);
        
        uvec tempngbrs(ngbrs.size());
        
        for(uword i = 0; i<ngbrs.size() ; i++){
            colsums(i) = sum(newG.col(i));
            cout << colsums << endl;
            tempngbrs(i) = i;
        }
        
        temporder = ngbrs(sort_index(colsums(tempngbrs)));
        

    
        tempsize2 += tempsize1;
        tempsize1 = temporder.size();
    
        cout << "Press Enter to Continue";
        cin.ignore();
    
        for(uword i = tempsize2; i<tempsize2+tempsize1; i++){
            neworder(i) = temporder(i-tempsize2);
            cout << temporder(i-tempsize2) << endl;
        }
    
        cout << ngbrs << endl;
    
        cout << neworder <<endl;
        
        count = tempsize1 + tempsize2;
    }

    cout << S << endl;
    cout << G << endl;
    
    cout << S.submat(neworder,neworder) << endl;
    cout << G.submat(neworder,neworder) << endl;
    
    
    return 0;
}


