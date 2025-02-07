###############################################################################################################
###### Scripts for cross validation for prediction of environmental main effects and factor loadings###########
##This work © 2025 by Nick Fradgley, CSIRO is licensed under CC BY 4.0 ########################################

# 'all.fa.sums' is a list of FA model summaries from the '4. Multi year analysis' script
# 'all.mods is a list of fitted models from the '4. Multi year analysis' script
# 'WSmat' is a matrix of weather ans soil ECs for all environments from the '1. EC processing script' 

library(randomForest)
library(pls)
library(glmnet)

  #Get factor loadings from FA3 model
  Lam<-all.fa.sums$fa3fa1$Additive_G_effects$Lam
  #Get environmental main effects from FA3 model
  env.maineffs<-all.mods$fa3fa1$coefficients$fixed
  env.maineffs<-env.maineffs[1]+env.maineffs[2:length(env.maineffs)]
  names(env.maineffs)<-rownames(all.mods[[m+1]]$coefficients$fixed)[-1]
  
  #make empty matrices for predictions
  main.eff.OOF.preds<-matrix(NA,nrow = length(env.maineffs),ncol = 3,
                             dimnames=list(names(env.maineffs),c("RF","LASSO","PLS")))
  RF.OOF.predictd.lam<-matrix(NA,ncol=ncol(Lam),nrow=nrow(Lam),dimnames = dimnames(Lam))
  LASSO.OOF.predictd.lam<-matrix(NA,ncol=ncol(Lam),nrow=nrow(Lam),dimnames = dimnames(Lam))
  PLS.OOF.predictd.lam<-matrix(NA,ncol=ncol(Lam),nrow=nrow(Lam),dimnames = dimnames(Lam))
  
  #Leave one out cross validation 
  for(i in 1:nrow(Lam)){
    #Env main effects
    #RF model
    rfmod<-randomForest(x = WSmat[rownames(Lam),][-i,], y = env.maineffs[-i], ntree = 1000) #Train model
    main.eff.OOF.preds[i,"RF"]<-predict(rfmod,newdata = WSmat[rownames(Lam),][i,]) #Make predictions
    
    #LASSO
    cv.mod<-cv.glmnet(x = WSmat[rownames(Lam),][-i,],y = env.maineffs[-i], #CV LASSO model to find optimum lambda
                      nfolds = 10,lambda=seq(0.0001,0.6,length.out=100)^2)
    lassso.mod<-glmnet(x = WSmat[rownames(Lam),][-i,],y = env.maineffs[-i], #Fit LASSO model
                       lambda = cv.mod$lambda.min)
    main.eff.OOF.preds[i,"LASSO"]<-predict(lassso.mod,newx = WSmat[rownames(Lam),][i,]) #Make predictions
    
    #PLS
    plsdat<-cbind.data.frame(env.maineffs[-i],WSmat[rownames(Lam),][-i,])
    colnames(plsdat)[1]<-"y"
    cols.to.remove<-apply(plsdat,2,function(x) sum(duplicated(x)))>85
    plsmod<-plsr(y~., data=plsdat[,!cols.to.remove], scale=TRUE, validation="LOO",ncomp=50) #Find the optimum number of components to use
    opt.comps<-which.max(apply(plsmod$validation$pred[,1,],2,function(x) cor(x,plsmod$model$y))) 
    plsmod<-plsr(y~., data=plsdat[,!cols.to.remove], scale=TRUE, validation="LOO",ncomp=opt.comps) #Fit PLS model
    main.eff.OOF.preds[i,"PLS"]<-predict(object = plsmod,newdata = t(as.data.frame(WSmat[rownames(Lam),!cols.to.remove[-1]][i,])),ncomp = opt.comps) #Make predictions
    
    #Factor loading predictions
    for(f in 1:ncol(Lam)){ #For each factor separately
      cat("|")
      #RF
      rfmod<-randomForest(x = WSmat[rownames(Lam),][-i,],y = Lam[-i,f],ntree = 500) #Train model
      RF.OOF.predictd.lam[i,f]<-predict(rfmod,newdata = WSmat[rownames(Lam),][i,]) #Make predictions
      #LASSO
      cv.mod<-cv.glmnet(x = WSmat[rownames(Lam),][-i,],y = Lam[-i,f], #CV LASSO model to find optimum lambda
                        nfolds = 10,lambda=seq(0.0001,0.3,length.out=100)^2)
      lassso.mod<-glmnet(x = WSmat[rownames(Lam),][-i,],y = Lam[-i,f], #Fit LASSO model
                         lambda = cv.mod$lambda.min)
      LASSO.OOF.predictd.lam[i,f]<-predict(lassso.mod,newx = WSmat[rownames(Lam),][i,]) #Make predictions
      #PLS
      plsdat<-cbind.data.frame(Lam[-i,f],WSmat[rownames(Lam),][-i,])
      colnames(plsdat)[1]<-"y"
      cols.to.remove<-apply(plsdat,2,function(x) sum(duplicated(x)))>85
      plsmod<-plsr(y~., data=plsdat[,!cols.to.remove], scale=TRUE, validation="LOO",ncomp=50)#Find the optimum number of components to use
      opt.comps<-which.max(apply(plsmod$validation$pred[,1,],2,function(x) cor(x,plsmod$model$y)))
      plsmod<-plsr(y~., data=plsdat[,!cols.to.remove], scale=TRUE, validation="LOO",ncomp=opt.comps)#Fit PLS model
      PLS.OOF.predictd.lam[i,f]<-predict(object = plsmod,newdata = t(as.data.frame(WSmat[rownames(Lam),!cols.to.remove[-1]][i,])),ncomp = opt.comps) #Make predictions
    }
  }
  
  

