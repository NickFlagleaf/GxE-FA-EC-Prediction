

SNPprune<-function(marker.mat,#NxM marker matrix with row names and column names. Markers coded as 0,1 and 2 for hom, het and hom.
                   cutoff=.8, #Marker correlation cut of value to prune linked markers 
                   maf.range=.1, #Range of minor allele frequencies to first subset within before calculating correlations
                   verbose=T){ #Logical whether to print progress. Default is TRUE.
if(cutoff<1){
  if(verbose==T){print("Pruning:")}
  to.remove<-c() #create empty vector or marker names to remove
  maf<-apply(marker.mat,2,FUN = function(x) { #Make MAF function
    n0<-sum(na.omit(x=="0"))
    n2<-sum(na.omit(x=="2"))
    maf<-(min(n0,n2))/(n0+n2)        
    return(maf)
  })
  
  for(i in 1:ncol(marker.mat)){ #For each marker in order
    if(!colnames(marker.mat)[i]%in%to.remove){ #If the marker isn't already assigned for removal:
      pruned<-marker.mat[,!colnames(marker.mat)%in%to.remove] #Get the latest pruned subset of markers
      tocheck<-colnames(marker.mat)[(i+1):ncol(marker.mat)] #Get all the markers that haven't been checked yet
      tocheck<-tocheck[tocheck%in%colnames(pruned)] 
      in.maf.range<-abs(maf-maf[i])<maf.range #Make a logical vector of markers within the MAF range
      tocheck<-tocheck[tocheck%in%colnames(marker.mat)[in.maf.range]] #Subset the markers to check based on the MAF range
      if (length(tocheck)>1) {
        #Calculate correlations between the focal marker and all markers to check in the MAF range
        c.sub<-apply(pruned[,tocheck],2,function(x) cor(x,pruned[,colnames(marker.mat)[i]],
                                                        use = "pairwise.complete.obs"))
        linked<-names(c.sub)[abs(c.sub)>cutoff] #Get subset of linked markers above the cut off
        }else{linked=NULL}
      to.remove<-c(to.remove,linked[!linked==colnames(marker.mat)[i]]) #Assign linked markers for removal
    }
    #Print progress
    if(isTRUE(!round((1:ncol(marker.mat)/ncol(marker.mat))*100)[i]==round((1:ncol(marker.mat)/ncol(marker.mat))*100)[i-1])){
      if(verbose==T){cat("|",sep="")}
    }
  }
  if(verbose==T){cat("FINISHED! :)",sep="")}
  return(pruned)
}
}
