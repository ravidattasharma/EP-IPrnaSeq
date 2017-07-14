#################################################################################
# Description
#############
# ExonPointer for finding cassette exon using RNA-Seq Data:
# require four matrices as input: 
# y=raw counts of reads corresponding to different terms
# (exons, junction and introns) X samples
# RMmeber= realtion of terms (exons, junction and introns) X exons
# where skipping junction marked via 0.5 instead of 1
# contrastM= contrast matrix that can be used by limma, contains contrast to be
# evaluated
# designM= designMatrix that can be used by limma, contains experiment information
# Groups= sample groups such as c(1,1,1,1,2,2,2,2) where first four samples belong 
# to first group and other four to second group
#
# Citation: Differential alternative splicing coupled to nonsense-mediated decay  
#           of mRNA ensures dietary restriction-induced longevity, Syed Shamsh 
#           Tabrez, Ravi Datta Sharma,Vaibhav Jain, Atif Siddiqui and Arnab 
#	          Mukhopadhyay, Nature communication, 2017        
#
# In case of any question or query please Contact: 
# Dr. Ravi D. Sharma
# ravidattasharma@gmail.com
# Amity Institute of Biotechnology,
# Amity University Haryana, Gurgaon, 122413, India
##################################################################################


ePrnaseqFunction <- function(counts=y,RMM= RMmeber,contrastM=contrastM,designM=designM, Groups=Groups){
#Main function#
#load required library
require(limma)
require(edgeR)
#checkpoint
      if (nrow(counts) != nrow(RMM) ) {
						      cat('Dimension missmatch?','\n',sep='');
						      return(NA)	     }
       if (nrow(RMM)< 3 ) { cat('Skipping short Gene','\n',sep=''); 
			    return(NA)  }
      if (is.null(rownames(counts)) || nrow(counts)< 3 ) {
						    cat('Skipping short Gene','\n',sep='');
					      return(NA) 	        }  
#converting expression
      ucounts <- counts
      counts <- DGEList(counts=counts,group=Groups)
      counts <- calcNormFactors(counts)
      try( fit<- voom(counts,designM),silent=TRUE)
      if(!exists("fit",inherits=FALSE))  return(NA) 
#removing zero counts false expression if any
      E <- ((ucounts !=0)*1)* fit$E
#Check low expressed reads #remove terms that have expression less than 6.32
      expression<- removeLECounts(E,designM,Groups)
#Running limma
      fit <- lmFit(expression, designM) 
      fit2 <- contrasts.fit(fit, contrastM)            
      fit2 <- eBayes(fit2) 
      index <- match(rownames(expression), rownames(RMM))
      RMM<- RMM[index,,drop=FALSE] 
      indexToremove<- match(colnames(RMM),rownames(RMM))
      RMM<- RMM[,!is.na(indexToremove),drop=FALSE]
      auxMatrix<- (RMM ==0.5 )* 0.5 
      RMMij<- (RMM > 1) *1
      RMM<- (RMM > 0.5) *1 
#Check the dimenstion of the matrices 
      if (!any(match(rownames(expression),rownames(RMM)),na.rm=TRUE)) {
						      cat('Dimension missmatch?','\n',sep='');
						      return(NA)	     }
#Running limma by contrast
      nbrOfcontrast<- ncol(contrastM)
      contrasts<- colnames(contrastM)
      FAP<- vector('list',length=nbrOfcontrast)
      FAP <- lapply(1:nbrOfcontrast, function(j) 
      {
	t <- data.frame(fit2$t[,j])
	pvalue2tail <- data.frame(fit2$p.value[,j])
	pvalue1tail <- pvalue2tail * .5 * (t>0) + (1 -pvalue2tail * .5) * (t<=0)
	SumOfPvalues <- t(as.matrix(pvalue1tail)) %*% ( (1 * (RMM > 0)) + ((-1)*(auxMatrix ==0.5))) + colSums(1*(auxMatrix ==0.5))
	NumberOfPvalues <- colSums((auxMatrix + RMM) !=0)
	EqPvalues<- as.matrix(unlist(lapply(1:length(SumOfPvalues), function(x) sumPvalsMethod(SumOfPvalues[x],NumberOfPvalues[x]))))
	constraints <- (colSums((RMMij > 0)* 1) >1)  &  (colSums((auxMatrix>0)*1) >0) 
	EqPvalues <- as.matrix(constraints * (( EqPvalues < 0.05 ) * EqPvalues ))
        rownames(EqPvalues)<- colnames(SumOfPvalues)
	EqPvalues <- as.matrix( EqPvalues[which(EqPvalues != 0),,drop=FALSE])
	EqPvalues <- cbind(EqPvalues,rep(contrasts[j],nrow(EqPvalues)))
	tGE<- as.matrix(t(as.matrix(t)) %*% ((RMM +auxMatrix)> 0)*1) # now it sums up the t-statistics of term that are used for summing up pvalues
	EqPvalues<- cbind(EqPvalues,tGE[,rownames(EqPvalues)])
	#Saving Pvalues
	if(nrow(EqPvalues)==0) return(NA) else {
						colnames(EqPvalues)<- c('EqPvalues','Contrast','t-statistics')
						return(EqPvalues) }	
	})
    return(FAP);
  } 

#Supportive function1
sumPvalsMethod <-  function(x,n) {
if (n < 10) {
psumunif(x,n)
} else {
pnorm(x,n/2,sqrt(n/12),lower=TRUE)
}}

#Supportive function 2
psumunif <-  function(x,n) 1/factorial(n) * sum(sapply(0:n, function(k) (-1)^k * choose(n,k) * ifelse(x > k,x-k,0)^(n)))


removeLECounts<-function(y,designM,Groups=NULL,threshold=6.32, ...)
{
#Supportive function3
#Function to remove low expressed terms
require(matrixStats)
y <- as.matrix(y)
if(is.null(Groups)) {
if(any(rowSums(designM) >1)) 
Groups <- rowSums(designM %*% t(designM))  else Groups<-which(designM!=0,arr.ind=TRUE)[,2]
} else Groups <- Groups
Gfactor<- as.factor(Groups)
n<- as.numeric(levels(Gfactor))
GIndex <- match(Gfactor, n)
rs<- lapply(1:length(n),function(x) rowMedians(y[,(Gfactor==n[x]),drop=FALSE]))
rs<-do.call('cbind',rs)
newy <- y[(rowMaxs(rs) >threshold),,drop=FALSE]
return(newy)
}

############################################################################
############################################################################