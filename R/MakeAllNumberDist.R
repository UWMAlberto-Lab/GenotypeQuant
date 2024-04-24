MakeAllNumberDist<-function(genotypes.df,ngametos=80,All.nchar=3,Fis=FALSE,Fis.treshold=0.1,npermutations=5000 ){
	requireNamespace("adegenet")
	if(!any(All.nchar%in%c(1,2,3)))stop("All.nchar, the number of digits to code each allele, needs to be 1, 2 or 3")
	#Converting to geneind class (adegenet)
	genotypes.GI<-adegenet::df2genind(
		X=genotypes.df,
		ncode = 3,
		ind.names = 1:nrow(genotypes.df),
		loc.names =names(genotypes.df) ,
		pop = NULL,
		NA.char = "0",
		ploidy = 2,
	)
	#Converting to genpop class (adegenet)
  genInd.pop<-adegenet::genind2genpop(genotypes.GI,quiet=T)
  #number of individuals per pop
  N<- dim(genotypes.GI@tab)[1]
  # getting loci names
  loc.names<-levels(genotypes.GI@loc.fac)
  nloci<-length(loc.names)
  #getting a table of allele frequencies
  genInd.AllFreq<-apply(genInd.pop@tab,MARGIN=2,function(x){ x/N } )
  nloci<-length(loc.names)

  genInd.AllFreq.DF<-data.frame(matrix(unlist(strsplit(names(genInd.AllFreq),split="\\.")),ncol=2,byrow=T))
  genInd.AllFreq.DF<-cbind(genInd.AllFreq.DF,genInd.AllFreq)
  names(genInd.AllFreq.DF)<-c("Locus","Allele","Frequency")

  Permutations.A3D.mod<-array(dim=c(ngametos,nloci, npermutations))

   if(Fis){      #For Fis correction
  	requireNamespace("genepop")
   	Ref.Gen<-genotypes.df
  	#export to Genepop from the data frame, my own function ----
  	
   	if(All.nchar==3){
  			Ref.Gen[Ref.Gen==0]<-"000000"
  						}else{if(All.nchar==2){Ref.Gen[Ref.Gen==0]<-"0000"
  									}else{if(All.nchar==1){Ref.Gen[Ref.Gen==0]<-"00"}
  												}
  												}
  			base::write(c("You can edit this line",names(Ref.Gen)),"GenepopIN",append=T)
  			base::write("Pop","GenepopIN",append=T)
  			N<-nrow(Ref.Gen)
  			utils::write.table(cbind(rep("Ref.Pop,",N),Ref.Gen),file="GenepopIN", append=T, quote=F, col.names=F,row.names=F )
  	
  			#using genepop to get FIS
  			genepop::genedivFis("GenepopIN",outputFile="FIS.out",verbose=F)
  			# 
  			FIS.lines<-readLines("FIS.out")
  			LociFisRows<-grep("All samples",FIS.lines)-2
  	
  			temp<-unlist(strsplit(FIS.lines[LociFisRows]," "))
  			FIS.DF<-data.frame(matrix(temp[!(nchar(temp)==0 | temp=="Ref.Pop")],ncol=3,byrow=T))
  			names(FIS.DF)<-c("1-Qintra","1-Qinter","Fis") #Use Qintra to sample correlated alleles when FIS is TRUE
  			system2("rm","GenepopIN")
  			system2("rm","FIS.out")
  			FIS.DF<-apply(FIS.DF,2,as.numeric)
  	  	Qintra<-(1-FIS.DF[,1])
  	
  	  	for(g in 1:ngametos){     #Creating the distributions using Qintra probability that two alleles are identical by descent (state really)
  	    	for(l in 1:nloci){
  		    	for(p in 1:npermutations){
  				
  			  	LocusDF<-genInd.AllFreq.DF[genInd.AllFreq.DF$Locus==loc.names[l],]
  					s.alleles<-1:g
  			  	if(	FIS.DF[l,3]<Fis.treshold){
  				     	s.alleles<-sample(LocusDF$Allele,size=g,replace=T,prob=LocusDF$Frequency)}else{
  							s.alleles[1]<-sample(LocusDF$Allele,size=1,replace=T,prob=LocusDF$Frequency)
  			        if(g>1){
  			        	 for(a in 2:g){
  		              	Same.L<-sample(c(1,0),size=1,prob=c(Qintra[l],(1-Qintra[l])))
  			              	if(Same.L){
  				    		        s.alleles[a]<-s.alleles[a-1]}else{
  					                   	s.alleles[a]<-sample(LocusDF$Allele,size=1,replace=T,prob=LocusDF$Frequency)
  					              	   }  				
  						            	}
  						         }
  					          }
 				        Permutations.A3D.mod[g,l,p]<-length(unique(s.alleles))
  				      }
  			      message(paste("N gametos:",g," Loci:",loc.names[l]))
  		       }
            }
  	  	}else{ #without FIS correction
  	for(g in 1:ngametos){
  		for(l in 1:nloci){
  			for(p in 1:npermutations){
  				
  				LocusDF<-genInd.AllFreq.DF[genInd.AllFreq.DF$Locus==loc.names[l],]
  				Permutations.A3D.mod[g,l,p]<-length(unique(sample(LocusDF$Allele,size=g,replace=T,prob=LocusDF$Frequency)))
  				
  			}
  			message(paste("N gametos:",g," Loci:",loc.names[l]))
  		}
  	}
  }
  
  return(Permutations.A3D.mod)
}
