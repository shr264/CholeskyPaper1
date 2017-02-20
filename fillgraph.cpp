# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::export()]]
NumericMatrix gNM1(NumericMatrix Gfill, int p){
// create 4 by 6 matrix
// 0-valued entries
  for(int i=p;i>0;i--){
    for (int j=(i-1);j>0;j--){
      for(int k = 1;k<=j;k++){
                if((k<j)&&(j<i)){
                  if((Gfill(i-1,k-1)==1)&&(Gfill(j-1,k-1)==1)){
                    if(Gfill(i-1,j-1)!=1){
                      Gfill(i-1,j-1) = 1;
                      Gfill(j-1,i-1) = 1;
                        }
                    }
                }
            }
        }
    }
return(Gfill) ;
}
