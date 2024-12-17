evaluator = function(truebeta,truesig,betahat,newstart=2,endofsample,crit=0.0001,visualizar=0){
  
  evalperiod = newstart:(dim(betahat)[2]-1-endofsample)
  
  mape <-mean(abs(betahat[,evalperiod]-truebeta[,evalperiod]))
  
  mspe <-mean((betahat[,evalperiod]-truebeta[,evalperiod])^2)
  
  #Discovery of Non-TVP
  umat = betahat[,evalperiod]-betahat[,(evalperiod-1)]
  sigmasq <- diag(umat%*%t(umat))
  
  sighat =  array(0,dim=dim(as.matrix(truesig)))
  for(ss in 1:length(sighat)){
    if((sigmasq[ss]/(mean(sigmasq)+0.0000001))>crit){sighat[ss]=1}
  }
  missrate=sum(abs(truesig-sighat))/length(truesig)
  
  if(visualizar==1){
    par(mfrow=c(3,3))
    for(j in 1:length(sighat)){ts.plot(cbind(betahat[j,evalperiod],truebeta[j,evalperiod]),col=c('blue','red'),
                                       lwd=3,ylim=c(min(cbind(betahat[j,evalperiod],truebeta[j,evalperiod]))-0.05,
                                                    max(cbind(betahat[j,evalperiod],truebeta[j,evalperiod]))+0.05))}
    par(mfrow=c(1,1))
  }
  return(c(mape,mspe,missrate))
}
