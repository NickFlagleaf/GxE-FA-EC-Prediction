


rfSNPimpute<-function(marker.mat,#NxM marker matrix with row names and column names. Markers coded as 0,1 and 2 for hom, het and hom.
                      n.cores=NULL, #number of cores to use 
                      maf.range=.1, #range of minor allele frequencies to subset within
                      nbest.markers=700, #The number of most closely linked markers for imputation models 
                      max.pred.error=.2, #The maximum out of bag prediction error to keep imputed markers
                      verbose=T) #Logical whether to print progress. Default is TRUE.

  pckgs<-c("parallel","doParallel","foreach","ranger")
  missingpckgs<-pckgs[!pckgs%in%installed.packages()[,"Package"]]
  if(length(missingpckgs)>0){print(paste("Install packages:",paste(missingpckgs,collapse = " ")))}

  {
  if(is.null(n.cores)){
    n.cores<-parallel::detectCores()
  }
  
  missing<-(apply(marker.mat,2,FUN = function(x) sum(is.na(x))))/nrow(marker.mat) #calculate NA freq
  is.mono<-apply(marker.mat,2,FUN = function(x) sum(na.omit(x=="0"))==0 | sum(na.omit(x=="2"))==0) #remove monomorphics
  marker.mat<-marker.mat[,!is.mono]
  maf<-apply(marker.mat,2,FUN = function(x) { #Calculate MAF
    n0<-sum(na.omit(x=="0"))
    n2<-sum(na.omit(x=="2"))
    maf<-(min(n0,n2))/(n0+n2)        
  return(maf)
    })

    require(foreach)
    rest.of.data<-marker.mat
    median.impute<-function(x){ #first impute with median
      x[is.na(x)]<-median(na.omit(x))
      return(x)
    }
    rest.of.data<-cbind(apply(rest.of.data,2,FUN = median.impute))
    file.remove("Imputing log.txt")
    if(verbose==T){print("Imputing log is output in directory")}
    #Set up parallel processing   
    cl <- parallel::makeCluster(n.cores,outfile="Imputing log.txt")
    doParallel::registerDoParallel(cl)
    imp.start<-Sys.time()
    imputed<-foreach::foreach(i=1:ncol(marker.mat),.combine = cbind,.multicombine = TRUE,.packages='ranger',.verbose = F) %dopar% { #for each marker:
      in.maf.range<-abs(maf-maf[i])<maf.range #subset based on maf range
      c.sub<-apply(marker.mat[,in.maf.range],2,function(x) cor(x,marker.mat[,i],use = "pairwise.complete.obs")) #subset based on correlation with focal marker
      rest.of.data.sub<-rest.of.data[,names(abs(c.sub))[order(abs(c.sub),decreasing=T)[1:min(nbest.markers,length(c.sub))]]]
      rest.of.data.sub<-rest.of.data.sub[,!colnames(rest.of.data.sub)==colnames(marker.mat)[i]]
      isna<-is.na(marker.mat[,i]) #take out the missing data
      imp.marker<-marker.mat[,i]
      if(sum(isna)>0){ #if there's any missing data, predict from all the linked markers using a RF model
        rfmod<-ranger::ranger(x = rest.of.data.sub[!isna,],
                              y = as.factor(marker.mat[!isna,i]),
                              num.trees = 100,importance = "none",
                              num.threads=1)
      if(rfmod$prediction.error<max.pred.error){ #only impute if accurate enough
      imp.marker[isna]<-as.numeric(as.character(predict(rfmod,data = rest.of.data.sub)$predictions[isna]))
      }
      }
      #print output
      if(isTRUE(!round((1:ncol(marker.mat)/ncol(marker.mat))*100)[i]==round((1:ncol(marker.mat)/ncol(marker.mat))*100)[i-1])){
        if(verbose==T){print(paste("|",round((1:ncol(marker.mat)/ncol(marker.mat))*100)[i],"%",sep=""))
          elapsed<-Sys.time()-imp.start
          percent.done<-((i/ncol(marker.mat))*100)
          print("Estimated time left:")
          print((100-percent.done)*(elapsed/percent.done))
        }
      }
      return(imp.marker)
    }
    parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
   dimnames(imputed)<-dimnames(marker.mat)
   return(imputed)
}

  