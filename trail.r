    for(i in p:1){
        for (j in (i-1):1){
            for(k in 1:j){
                if((k<j)&&(j<i)){
                    if((Gfill[i,k]==1)&&(Gfill[j,k]==1)){
                        if(Gfill[i,j]!=1){
                            Gfill[i,j] = Gfill[j,i] = 1
                        }
                    }
                }
            }
        }
    }


sourceCpp("fillgraph.cpp")

gNM1(G,10) - Gfill
