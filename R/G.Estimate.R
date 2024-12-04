G.Estimate<-function(Obs.N.all,All.Num.Dist,ML.curve=TRUE,ploidy=1, col="black",col.est="darkcyan",col.ci="darkcyan",exclude=NULL){
	
	
	All.Freq.Table<-All.Num.Dist[[2]]
	Alleles.Per.Locus<-as.numeric(tapply(All.Freq.Table$Locus,factor(All.Freq.Table$Locus,levels=unique(All.Freq.Table$Locus)),length))
	All.Num.Dist<-All.Num.Dist[[1]]
	loc.names<-unique(All.Freq.Table$Locus)
	if(!is.null(exclude)){
	  if((length(loc.names)-length(exclude))<2)stop("You excluded too many loci leaving only one or fewer left")
	 }
	  
	if(!ploidy%in%c(1,2))stop("ploidy value must be 1 or 2")
	if(length(Obs.N.all)!=dim(All.Num.Dist)[2])stop("Obs.N.all vector needs to be of same length as the 2nd dimension in All.Num.Dist array, i.e., the number of loci used")

  gmax<-dim(All.Num.Dist)[1]
  nloci<-length(Obs.N.all)
  if(!is.numeric(exclude) & !is.null(exclude))stop("Parameter exclude needs to be numeric (the index of loci to remove if any")
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
  #P.array<-P.array+0.001  #Adding a small amount to avoid multiplying by zero
  if(!is.null(exclude)){
  ML.function<-apply(P.array[,-exclude],1,prod,na.rm=T)}else{  # multiplying over different loci, while excluding some loci, e.g., by observing high chisq values
  	ML.function<-apply(P.array,1,prod,na.rm=T)
  }    
  # excluding Locus that have observed number alleles equal to the number of alleles in that loci in the ref pop
  
  if(!is.null(exclude)){
    L1<-(Obs.N.all<Alleles.Per.Locus & !(1:nloci)%in%exclude)
    if(sum(L1,na.rm=T)==0){
      ML.function.Cor<-NA
      }else{
           if(sum(L1,na.rm=T)==1){
             ML.function.Cor<-P.array[,L1]
             }else{
                  ML.function.Cor<-apply( P.array[,L1],1,prod,na.rm=T)
             }
      }
    }else{
         L1<-Obs.N.all<Alleles.Per.Locus
         if(sum(L1,na.rm=T)==0){
            ML.function.Cor<-NA
            }else{
                if(sum(L1,na.rm=T)==1){
                     ML.function.Cor<-P.array[,L1]
                    }else{
                    ML.function.Cor<-apply( P.array[,L1],1,prod,na.rm=T)
                    }
            }
    }
  
  # Chisquare to detect problems with eDNA amplification
  Obs.Freq<-Obs.N.all
  LociExpFreq<-rowMeans(All.Num.Dist[dim(All.Num.Dist)[1],,])
  LociExpRelFreq<-LociExpFreq/sum(LociExpFreq)
  #LociExpRelFreq<-rowMeans(All.Num.Dist[dim(All.Num.Dist)[1],,])/sum(LociExpFreq)
  
  	NoNA.i<-which(!is.na(Obs.Freq)) #dealing with NA loci that did not amplify
  	Obs.Freq<-Obs.Freq[NoNA.i]
  	LociExpRelFreq<-LociExpRelFreq[NoNA.i]
  	LociExpRelFreq<-LociExpRelFreq/sum(LociExpRelFreq)
  
  Exp.Freq<-sum(Obs.Freq)*LociExpRelFreq
  Diff.O.minus.E<-Obs.Freq-Exp.Freq
  ChiSq.by.Loci<-((Diff.O.minus.E^2)/Exp.Freq)
 
  ChiSq.Stat<-sum(ChiSq.by.Loci)
  ChiSq.Pvalue<-pchisq(ChiSq.Stat,df=length(NoNA.i)-2,lower.tail=F)
  
  Amp.Qual<-matrix(as.numeric(c( round(ChiSq.by.Loci,5),round(ChiSq.Stat,5),round(ChiSq.Pvalue,6),
   round(Diff.O.minus.E,5),"","", Obs.Freq,"","", round(Exp.Freq,1),"","")),nrow=4,byrow=T)
  
  colnames(Amp.Qual)<-c(loc.names[NoNA.i],"ChiSq","P.value")
  row.names(Amp.Qual)<-c("Chi-square","O.minus.E","Observed.Freq","Expected.Freq")
  Amp.Qual<-data.frame(Amp.Qual)
  
 
  if(ploidy==1){
    if(all(is.na(ML.function))){
      g.estimate<-NA
      }else{
         g.estimate<-(1:gmax)[which.max(ML.function)]
      }
    if(all(is.na(ML.function.Cor))){
      g.estimate.Cor<-NA
      }else{
          g.estimate.Cor<-(1:gmax)[which.max(ML.function.Cor)]
      }
    
  #CI estimation
   if(!is.null(exclude)){
      LogML<-rowSums(log(P.array[,-exclude]),na.rm=T)
      }else{  LogML<-rowSums(log(P.array),na.rm=T) }
 
    #CI for corrected estimate  
   if(!is.null(exclude)){
      L2<-(Obs.N.all<Alleles.Per.Locus & !(1:nloci)%in%exclude)
      if(sum(L2,na.rm=T)==0){
         LogML.Cor<-NA
      }else{
          if(sum(L2,na.rm=T)==1){
            LogML.Cor<-log(P.array[,L2])
          }else{
            LogML.Cor<-rowSums(log(P.array[,L2]),na.rm=T)
          }
      }
    }else{
   	  L2<-Obs.N.all<Alleles.Per.Locus
   	  if(sum(L2,na.rm=T)==0){
   	    LogML.Cor<-NA
   	    }else{
   	       if(sum(L2,na.rm=T)==1){
   	         LogML.Cor<-log(P.array[,L2])
   	         }else{
                 LogML.Cor<-rowSums(log(P.array[,L2]),na.rm=T) }
   	           }
    }
    
  if(!all(is.na(LogML))){
                cond<-LogML>=(LogML[g.estimate]-1.92)
                CI.range<-range((1:gmax)[cond])
                }else{
                  CI.range<-rep(NA,2)
                }
    
    if(!all(is.na(LogML.Cor))){
               cond.Cor<-LogML.Cor>=(LogML.Cor[g.estimate.Cor]-1.92)
               CI.range.Cor<-range((1:gmax)[cond.Cor])
               }else{
                 CI.range.Cor<-rep(NA,2)
               }
    
  estimateM<-t(matrix( c( g.estimate,CI.range)))
  colnames(estimateM)<-c("G.hat","lower.95.CI","upper.95.CI")
  
  estimateM.Cor<-t(matrix( c( g.estimate.Cor,CI.range.Cor)))
  colnames(estimateM.Cor)<-c("G.hat.Cor","lower.95.CI","upper.95.CI")

  if(ML.curve & any(LogML!="-Inf") & !all(is.na(LogML))){
    graphics::plot( LogML,type="l",xlab="Number of haploid genotypes",xlim=c(0,gmax),col=col)
    			#main=paste0("G.hat: ",estimateM$G.hat," 95.per.CI: ",estimateM$lower.95.CI,"-",estimateM$upper.95.CI ))
    lines(x=rep(g.estimate,2),y=c(-400,LogML[g.estimate]),col=col.est,lwd=2)
    lines(x=rep(CI.range[1],2),y=c(-400,LogML[CI.range[1]]),col=col.ci,lty=)
    lines(x=rep(CI.range[2],2),y=c(-400,LogML[CI.range[2]]),col=col.ci,lty=)
  }
  
  }else{ #diploids
    
    
    if(all(is.na(ML.function))){
      g.estimate<-NA
      g.estimate2N<-NA
    }else{
      g.estimate<-(1:gmax)[which.max(ML.function)]
      g.estimate2N<-floor(((1:gmax)[which.max(ML.function)])/2)
    }
    
    if(all(is.na(ML.function.Cor))){
      g.estimate.Cor<-NA
      g.estimate2N.Cor<-NA
    }else{
      g.estimate.Cor<-(1:gmax)[which.max(ML.function.Cor)]
      g.estimate2N.Cor<-floor(((1:gmax)[which.max(ML.function.Cor)])/2) 
      }
    
  
  	#CI estimation
  	if(!is.null(exclude)){
  	LogML<-rowSums(log(P.array[,-exclude]),na.rm=T)
  	}else{
  		LogML<-rowSums(log(P.array),na.rm=T)
  	}
    
    
  	if(!is.null(exclude)){
  	  L2<-(Obs.N.all<Alleles.Per.Locus & !(1:nloci)%in%exclude)
  	  if(sum(L2,na.rm=T)==0){
  	    LogML.Cor<-NA}
  	    else{
  	      if(sum(L2,na.rm=T)==1){
  	        LogML.Cor<-log(P.array[,L2])
  	        }else{
  	          LogML.Cor<-rowSums(log(P.array[,L2]),na.rm=T)
  	        }
      }
  	  }else{
  	   L2<-Obs.N.all<Alleles.Per.Locus
  	   if(sum(L2,na.rm=T)==0){
  	     LogML.Cor<-NA
  	     }else{
  	       if(sum(L2,na.rm=T)==1){
  	         LogML.Cor<-log(P.array[,L2])
  	         }else{
                LogML.Cor<-rowSums(log(P.array[,L2]),na.rm=T)
  	         }
  	     } 
      }
    
    
    
    
    if(!all(is.na(LogML))){
      cond<-LogML>=(LogML[g.estimate]-1.92)
      CI.range<-range((1:gmax)[cond])
      CI.range2N<-floor(range((1:gmax)[cond])/2)
    }else{
      CI.range<-rep(NA,2)
      CI.range2N<-rep(NA,2)
    }
    
    if(!all(is.na(LogML.Cor))){
      cond.Cor<-LogML.Cor>=(LogML.Cor[g.estimate.Cor]-1.92)
      CI.range.Cor<-range((1:gmax)[cond.Cor])
      CI.range2N.Cor<-floor(range((1:gmax)[cond.Cor])/2)
    }else{
      CI.range.Cor<-rep(NA,2)
      CI.range2N.Cor<-rep(NA,2)
    }
    
    
  	estimateM<-t(matrix( c( g.estimate2N,CI.range2N)))
  	estimateM.Cor<-t(matrix( c( g.estimate2N.Cor,CI.range2N.Cor)))
  	
  	colnames(estimateM)<-c("G.hat","lower.95.CI","upper.95.CI")
  	colnames(estimateM.Cor)<-c("G.hat.Cor","lower.95.CI","upper.95.CI")
  	
  	
  	if(ML.curve & any(LogML!="-Inf") & !all(is.na(LogML))){
  		graphics::plot( LogML,type="l",xlab="Number of diploid genotypes",ylab="Log ML",xlim=c(1,gmax),col=col,axes=F)
  		graphics::axis(1,at=seq(2,gmax,2),labels=seq(1,gmax/2,1))
  		graphics::axis(2)
  		#main=paste0("G.hat: ",estimateM$G.hat," 95.per.CI: ",estimateM$lower.95.CI,"-",estimateM$upper.95.CI ))
  		lines(x=rep(g.estimate,2),y=c(-400,LogML[g.estimate]),col=col.est,lwd=2)
  		lines(x=rep(CI.range[1],2),y=c(-400,LogML[CI.range[1]]),col=col.ci,lty=1)
  		lines(x=rep(CI.range[2],2),y=c(-400,LogML[CI.range[2]]),col=col.ci,lty=1)
  	}
  	
  	
  }
  #correcting estimates when LogP are all Inf
  if(!all(is.na(LogML))){
    if(all(LogML=="-Inf")){estimateM[1:3]<-rep(NA,3)}
  }
  
  if(!all(is.na(LogML.Cor))){
    if(all(LogML.Cor=="-Inf")){estimateM.Cor[1:3]<-rep(NA,3)}
  }
  
  return(list(estimateM=estimateM,LogML=LogML,estimateM.Cor=estimateM.Cor,LogML.Cor=LogML.Cor, Amp.Qual=Amp.Qual))
}

