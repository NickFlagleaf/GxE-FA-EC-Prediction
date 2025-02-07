##########################################################################################
##################### Scripts for multi-year analysis#####################################
##This work © 2025 by Nick Fradgley, CSIRO is licensed under CC BY 4.0 ###################

#'plot.data' is a data frame with all trial plot data. Column names include:
#   'Env' Factor or environment names
#   'Grain.Yield..t.ha.' Grain yield of each plot numeric
#   'ROW and COL row and column coordinates factors
#   'Year' is a factor for trial years
#   'GID' is the ID factor for genotypes

#'bad.trials' is a vector of names of bad trials with low mean reliability
#'combined.A' is a genotype kinship matrix
#'All.best.mods' is a list of the best models for each environment identified by the '2. Single trial analysis' script
#'all.yr.DIAG.mods' is a list of fitted diagonal structured models from the single year analysis '3. Single year analysis' script
#'all.yr.FA1.mods'  is a list of fitted FA1 structured models from the single year analysis '3. Single year analysis' script

library(asreml)
library(stringr)
library(corpcor)
library(data.table)
#Source functions
source("FA.nxt.SVs function.R")
source("fa.sum function.R")

plot.data<-plot.data[!plot.data$Env%in%bad.trials,] #Remove bad trials
GIDs<-unique(plot.data$GID)

A<-combined.A[rownames(combined.A)%in%GIDs,rownames(combined.A)%in%GIDs]
envs<-as.character(unique(plot.data$Env))


#get a list of random and residual effect parameters from single trial best fitted models
envs<-unique(plot.data$Env)
env.ranefs<-lapply(All.best.mods,function(x) names(x$G.param)[!names(x$G.param)=="GID"])
env.ranefs<-env.ranefs[as.character(envs)]
env.ranefs<-env.ranefs[unlist(lapply(env.ranefs,function(x) length(x)>0))]

env.residefs<-unlist(lapply(All.best.mods,function(x) attributes(x$formulae$residual)$term.labels))
env.residefs<-env.residefs[as.character(envs)]
env.residefs<-env.residefs[envs]
env.residefs<-sapply(unique(env.residefs),function(x) names(which(env.residefs==x)))
names(env.residefs)[names(env.residefs)=="units"]<-"id(ROW):id(COL)"


#fill in missing row col coordinate plots with NAs
fill.ins<-as.data.frame(matrix(NA,ncol=ncol(plot.data),nrow=1,dimnames = list(1,colnames(plot.data))))
for (i in 1:length(envs)) {
  data.sub<-plot.data[plot.data$Env==envs[i],]
  rowmax<-max(as.numeric(data.sub$ROW))
  colmax<-max(as.numeric(data.sub$COL))
  all.combs<-apply(expand.grid(1:rowmax,1:colmax),1,function(x) paste(x,collapse = "_"))
  all.obs<-apply(data.frame(data.sub$ROW,data.sub$COL),1,function(x) paste(x,collapse = "_"))
  missing<-all.combs[!all.combs%in%all.obs]
  if(length(missing)>0){
  fill.ins<-as.data.frame(matrix(NA,ncol=ncol(plot.data),nrow=length(missing),dimnames = list(1:length(missing),colnames(plot.data))))
  fill.ins$Env<-envs[i]
  fill.ins$ROW<-unlist(lapply(strsplit(missing,split = "_"),function(x) x[1]))
  fill.ins$COL<-unlist(lapply(strsplit(missing,split = "_"),function(x) x[2]))
  fill.ins$Year<-unlist(lapply(strsplit(as.character(envs[i]),split = "_"),function(x) x[1]))
  fill.ins$Site<-unlist(lapply(strsplit(as.character(envs[i]),split = "_"),function(x) x[2]))
  plot.data<-rbind(plot.data,fill.ins)
  }
}
plot.data<-plot.data[order(plot.data$Env,plot.data$ROW,plot.data$COL),]

asreml.options(workspace="25028mb",trace=TRUE,maxit=300) #Set up ASREML workspace with plenty of memory

plot.data$Env<-factor(plot.data$Env)
plot.data$PLOT<-factor(plot.data$PLOT)

#Build the DIAG random effects formula
genetic.form<-"~diag(Env):vm(GID,A)+" 
if(length(env.ranefs)>0){
  DIAGranform<-paste(genetic.form,paste(unlist(lapply(1:length(env.ranefs),function(x) 
    paste("at(Env,'",names(env.ranefs)[x],"'):",env.ranefs[[x]],"+",sep=""))),collapse = ""))
  DIAGranform<-str_sub(string = DIAGranform,start = 1,end = nchar(DIAGranform)-1)
}else{
  DIAGranform<-str_sub(string = genetic.form,start = 1,end = nchar(genetic.form)-1)
}
#Build the residual covariance effects formula
resid.form<-paste("~dsum(",paste("~",paste(names(env.residefs),collapse = "+"),"|Env",sep=""),",levels=env.residefs)")  

#Get DIAG model default start values
diag.start<-asreml(fixed = Grain.Yield..t.ha.~Env,
                   random = as.formula(DIAGranform),
                   residual = as.formula(resid.form),
                  na.action = na.method(x='include'),
                  data=plot.data,start.values = T)

#Replace starting values with recycled fitted parameters from single year models
rownames(diag.start$vparameters.table)<-diag.start$vparameters.table$Component
for(e in 1:length(all.yr.DIAG.mods)){
  new.params<-all.yr.DIAG.mods[[e]]$vparameters
  names(new.params)<-gsub(pattern = "Asub","A",names(new.params))
  names(new.params)<-gsub(pattern = "ROW:COL",replacement = "",x = names(new.params))
  new.params<-new.params[names(new.params)%in%rownames(diag.start$vparameters.table)]
  diag.start$vparameters.table[names(new.params),"Value"]<-new.params
  }

#Fit DIAG model
lmmod.diag<-asreml(fixed = Grain.Yield..t.ha.~Env,
                   random = as.formula(DIAGranform),
                   residual = as.formula(resid.form),
                   G.param = diag.start$vparameters.table,
                   R.param = diag.start$vparameters.table,
                   na.action = na.method(x='include'),
                   data=plot.data)

if(lmmod.diag$converge==F){lmmod.diag<-update(lmmod.diag)} #make sure model converged ok
summary(lmmod.diag) #Check model summary
gc() #clean up memory

#Build the FA1 random effects formula
if(length(env.ranefs)>0){
  FA1ranform<-paste("~rr(Env,1):vm(GID,A)+diag(Env):vm(GID,A)+",
                    paste(unlist(lapply(1:length(env.ranefs),function(x) paste("at(Env,'",names(env.ranefs)[x],"'):",env.ranefs[[x]],"+",sep=""))),collapse = ""))
  FA1ranform<-str_sub(string = FA1ranform,start = 1,end = nchar(FA1ranform)-1)
}else{
  FA1ranform<-"~rr(Env,1):vm(GID,A)+diag(Env):vm(GID,A)"
}

#Get starting values for FA1 model----
fa1.start<-asreml(fixed = Grain.Yield..t.ha.~Env,
                  random = as.formula(FA1ranform),
                  residual = as.formula(resid.form),
                  na.action = na.method(x='include'),
                  data=plot.data,start.values = T)

sv.for.FA1<-fa1.start$vparameters.table
rownames(sv.for.FA1)<-sv.for.FA1$Component
#Replace parameters with the same name from the fitted DIAG model
sv.for.FA1[names(lmmod.diag$vparameters)[names(lmmod.diag$vparameters)%in%rownames(sv.for.FA1)],"Value"]<-lmmod.diag$vparameters[names(lmmod.diag$vparameters)%in%rownames(sv.for.FA1)]

#Get first factor starting values from single year models
for(e in 1:length(all.yr.FA1.mods)){      
  new.params<-all.yr.FA1.mods[[e]]$vparameters
  names(new.params)<-gsub(pattern = "Asub)","A)",names(new.params))
  new.FAparams<-new.params[grep(pattern = "\\!fa1",x = names(new.params))]
  sv.for.FA1[names(new.FAparams),"Value"]<-abs(new.FAparams)
}

#Fit FA1 model
lmmod.fa1.fa1<-asreml(fixed = Grain.Yield..t.ha.~Env,
                      random = as.formula(FA1ranform),
                      residual = as.formula(resid.form),
                      G.param = sv.for.FA1,
                      R.param = sv.for.FA1,
                      na.action = na.method(x='include'),
                      data=plot.data)
if(lmmod.fa1.fa1$converge==F){lmmod.fa1.fa1<-update(lmmod.fa1.fa1)} #make sure model converged ok
summary(lmmod.fa1.fa1) #Check model summary
gc() #clean up memory
rm(all.yr.DIAG.mods,all.yr.FA1.mods) #clean up unused large objects


#Build the FA2 random effects formula
if(length(env.ranefs)>0){
  FA2ranform<-paste("~rr(Env,2):vm(GID,A)+diag(Env):vm(GID,A)+",
                    paste(unlist(lapply(1:length(env.ranefs),function(x) paste("at(Env,'",names(env.ranefs)[x],"'):",env.ranefs[[x]],"+",sep=""))),collapse = ""))
  FA2ranform<-str_sub(string = FA2ranform,start = 1,end = nchar(FA2ranform)-1)
}else{
  FA2ranform<-"~rr(Env,2):vm(GID,A)+diag(Env):vm(GID,A)"
}
#Get starting values for FA2
fa2.start<-asreml(fixed = Grain.Yield..t.ha.~Env,
                  random = as.formula(FA2ranform),
                  residual = as.formula(resid.form),
                  na.action = na.method(x='include'),
                  data=plot.data,start.values = T)

#Update FA2 start vals from FA1 model 
sv.for.FA2<-FA.nxt.SVs(old.vparams = lmmod.fa1.fa1$vparameters,
                       new.default.SVs = fa2.start$vparameters.table,
                       new.add.FA.order = 2)

#Fit FA2 model
lmmod.fa2.fa1<-asreml(fixed = Grain.Yield..t.ha.~Env,
                      random = as.formula(FA2ranform),
                      residual = as.formula(resid.form),
                      G.param = sv.for.FA2,
                      R.param = sv.for.FA2,
                      na.action = na.method(x='include'),
                      data=plot.data)

if(lmmod.fa2.fa1$converge==F){lmmod.fa2.fa1<-update(lmmod.fa2.fa1)} #make sure model converged ok
summary(lmmod.fa2.fa1) #Check model summary
gc() #clean up memory

#Build the FA3 random effects formula
if(length(env.ranefs)>0){
  FA3ranform<-paste("~rr(Env,3):vm(GID,A)+diag(Env):vm(GID,A)+",
                    paste(unlist(lapply(1:length(env.ranefs),function(x) paste("at(Env,'",names(env.ranefs)[x],"'):",env.ranefs[[x]],"+",sep=""))),collapse = ""))
  FA3ranform<-str_sub(string = FA3ranform,start = 1,end = nchar(FA3ranform)-1)
}else{
  FA3ranform<-"~rr(Env,3):vm(GID,A)+diag(Env):vm(GID,A)"
}

#Get starting values for FA3
fa3.start<-asreml(fixed = Grain.Yield..t.ha.~Env,
                  random = as.formula(FA3ranform),
                  residual = as.formula(resid.form),
                  na.action = na.method(x='include'),
                  data=plot.data,start.values = T)

#Update FA3 start vals from FA2 model----
sv.for.FA3<-FA.nxt.SVs(old.vparams = lmmod.fa2.fa1$vparameters,
                       new.default.SVs = fa3.start$vparameters.table,
                       new.add.FA.order = 3)

#Fit FA3 model
lmmod.fa3.fa1<-asreml(fixed = Grain.Yield..t.ha.~Env,
                      random = as.formula(FA3ranform),
                      residual = as.formula(resid.form),
                      G.param = sv.for.FA3,
                      R.param = sv.for.FA3,
                      na.action = na.method(x='include'),
                      data=plot.data)
if(lmmod.fa3.fa1$converge==F){lmmod.fa3.fa1<-update(lmmod.fa3.fa1)} #make sure model converged ok
summary(lmmod.fa3.fa1) #Check model summary
gc() #clean up memory


#Process summaries for all FA models
all.mods<-list(DIAGDAIG=lmmod.diag,fa1fa1=lmmod.fa1.fa1,fa2fa1=lmmod.fa2.fa1,fa3fa1=lmmod.fa3.fa1)
all.fa.sums<-lapply(all.mods[-1],function(x) fa.sum(x,addGfac = "GID",non.add = F,add = T,Efac = "Env"))

#Make summary stats accross models
models.comp<-data.frame("nParams"=unlist(lapply(all.mods,function(x)length(x$vparameters))),
                        "Coefs"=unlist(lapply(all.mods,function(x)length(x$coefficients$random)+length(x$coefficients$fixed))),  
                        "AIC"=unlist(lapply(all.mods,function(x) summary(x)$aic[1])),
                        "BIC"=unlist(lapply(all.mods,function(x) summary(x)$bic[1])),
                        "Log likelihood"=unlist(lapply(all.mods,function(x) summary(x)$loglik[1])),
                        "Additive perc Var exp"=c(NA,unlist(lapply(all.fa.sums,function(x) x$Additive_G_effects$`total %vaf`))))








