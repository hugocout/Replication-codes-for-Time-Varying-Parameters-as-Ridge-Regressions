Zfun <- function(data){
    # PGC 18/05/14
    # Performs the regressors space expansion prescribded in GC (2018)
    #to obtain a closed-form solution for the TVP-VAR.
    #Inputs :
    #- original data
    #- the submodel: do the TVP, or just extract the PCs or KPCs (1,2 or 3)
    #- freq: at which frequency we allow parameters to vary
    #Ouput: the new regressors matrix.
    
    #Generate matrices
    X= cbind(matrix(1,(nrow(data)),1),data)
    Z=array(0,dim=c((nrow(data)),(nrow(data)-1),(ncol(data)+1)))
    X<-as.matrix(X)
    
    #The sequence of tau's

    #New regressors
    for(tt in 2:(dim(data)[1])){
                Z[tt,1:(tt-1),]= repmat(X[tt,],tt-1,1)
    }
    
    Zprime <- c()
    for(kk in 1:(ncol(data)+1)){
        Zprime = cbind(Zprime,Z[,,kk])
    }
    
Zprime <- cbind(Zprime,X)

    #New regressors for u's ***GARBAGE***
    #if(freq==1){
    #  for(tt in 1:(dim(data)[1])){
    #        for(j in 1:dim(X)[2])
    #        TX[tt,j]= sum(X[1:tt,j])
    #        }
    #Zprime <- cbind(Zprime,TX,X)
    #}
    #else{ #this also executes the other loop as a special but less efficiently 
    #for(tt in 1:(dim(data)[1])){
    #      tho <- max(which(append(1,param_seq) <= tt))
    #      for(j in 1:dim(X)[2])
    #      TX[tt,j]= sum(X[1:tho,j])
    #  }
    #  Zprime <- cbind(Zprime,TX,X)
    #}
    #}
    
    #Get rid of the (very) unlikely zero vector
    #print(dim(Zprime)[2])
   # getout <- c()
  #  ccc<- 0
  #  for(b in 1:dim(Zprime)[2]){
  #      if(is.element(sd(Zprime[,b]),0)){
  #          ccc <- ccc +1
  #          getout <- cbind(getout,b)
  #      }
  #  }
    #if(ccc>0){
    #Zprime = Zprime[,-getout]
    #  print('getout')
    #}
    
    return(Zprime)
}

