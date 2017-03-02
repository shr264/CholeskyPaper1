rm(list=ls())

#par.vals - what you want to change it to. A matrix (use expand.grid)
#oat.nams - changed from. A matrix, same size
code.generate = function(par.vals, pat.nams, path="", path.opt=F,
                            templatename, out_name){ 
  count = 0
  for(i in 1:nrow(par.vals)){
    SourceCode=readLines(templatename) #read desired code

    for(j in 1:ncol(par.vals)){
      #replace with desired value
      SourceCode=gsub(pattern=pat.nams[i,j],replacement=par.vals[i,j],x=SourceCode) 
    }
    count = count + 1
    #output changed code
    path=ifelse(path==""|path.opt==F,"",paste(path,"/",sep=""))
    cat(SourceCode,file=paste(path, out_name,
                              paste(count,collapse=""),
                              '.R',sep=""),sep='\n')  
  }
}  

######################################################################
a = seq(0.6,60, length = 10)
AAseq = c(a)

numvals = expand.grid(list(AAseq, c(500,1500), 1:50))
nams = cbind(rep("AA",times=nrow(numvals)),
    rep("BB",times=nrow(numvals)),
    rep("CC",times=nrow(numvals)))

code.generate(numvals,nams,templatename='simulations_cscs3.r', out_name='simulations_cscs3_')


