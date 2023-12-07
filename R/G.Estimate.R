G.Estimate<-function(Obs.N.all,All.Num.Dist,ML.curve=TRUE,ploidy=1,
										 col="black",col.est="darkcyan",col.ci="darkcyan"){
	
	if(!ploidy%in%c(1,2))stop("ploidy must 1 or 2")
	if(length(Obs.N.all)!=dim(All.Num.Dist)[2])stop("Obs.N.all vector needs to be of same length as the 2nd dimension in All.Num.dist array, i.e., the number of loci used")

  gmax<-dim(All.Num.Dist)[1]
  nloci<-length(Obs.N.all)
  P.array<-matrix(nrow=gmax,ncol=nloci)
  for(g in 1:gmax){
    for(l in 1:nloci){
			if(is.na(Obs.N.all[l])){P.array[g,l]<-NA}else{
      all.n.L<-as.numeric( names(table(All.Num.Dist[g,l,])))== Obs.N.all[l]
      if(any(all.n.L)){
        P.array[g,l]<-as.numeric(table(All.Num.Dist[g,l,])[all.n.L]/dim(All.Num.Dist)[3])
      }else{P.array[g,l]<-0}
			}
    }
  }
  ML.function<-apply(P.array,1,prod,na.rm=T)    # multiplying over different loci

 
  if(ploidy==1){
   g.estimate<-(1:gmax)[which.max(ML.function)]
  #CI estimation
  LogML<-rowSums(log(P.array),na.rm=T)
  cond<-LogML>=(LogML[g.estimate]-1.92)
  CI.range<-range((1:gmax)[cond])
  estimateM<-t(matrix( c( g.estimate,CI.range)))
  colnames(estimateM)<-c("G.hat","lower.95.CI","upper.95.CI")

  if(ML.curve){
    graphics::plot( LogML,type="l",xlab="Number of haploid genotypes",xlim=c(0,gmax),col=col)
    			#main=paste0("G.hat: ",estimateM$G.hat," 95.per.CI: ",estimateM$lower.95.CI,"-",estimateM$upper.95.CI ))
    lines(x=rep(g.estimate,2),y=c(-400,LogML[g.estimate]),col=col.est,lwd=2)
    lines(x=rep(CI.range[1],2),y=c(-400,LogML[CI.range[1]]),col=col.ci,lty=)
    lines(x=rep(CI.range[2],2),y=c(-400,LogML[CI.range[2]]),col=col.ci,lty=)
  }
  }else{ #diploids
  	g.estimate<-(1:gmax)[which.max(ML.function)]
  	g.estimate2N<-floor(((1:gmax)[which.max(ML.function)])/2)
  	#CI estimation
  	LogML<-rowSums(log(P.array),na.rm=T)
  	cond<-LogML>=(LogML[g.estimate]-1.92)
  	CI.range<-range((1:gmax)[cond])
  	CI.range2N<-floor(range((1:gmax)[cond])/2)
  	estimateM<-t(matrix( c( g.estimate2N,CI.range2N)))
  	colnames(estimateM)<-c("G.hat","lower.95.CI","upper.95.CI")
  	
  	if(ML.curve){
  		graphics::plot( LogML,type="l",xlab="Number of diploid genotypes",ylab="Log ML",xlim=c(1,gmax),col=col,axes=F)
  		graphics::axis(1,at=seq(2,gmax,2),labels=seq(1,gmax/2,1))
  		graphics::axis(2)
  		#main=paste0("G.hat: ",estimateM$G.hat," 95.per.CI: ",estimateM$lower.95.CI,"-",estimateM$upper.95.CI ))
  		lines(x=rep(g.estimate,2),y=c(-400,LogML[g.estimate]),col=col.est,lwd=2)
  		lines(x=rep(CI.range[1],2),y=c(-400,LogML[CI.range[1]]),col=col.ci,lty=1)
  		lines(x=rep(CI.range[2],2),y=c(-400,LogML[CI.range[2]]),col=col.ci,lty=1)
  	}
  	
  	
}
  
  return(list(estimateM=estimateM,LogML=LogML))
}
