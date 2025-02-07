##########################################################################################
##################### Scripts for single-year analysis####################################
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

library(asreml)
library(stringr)

plot.data<-plot.data[!plot.data$Env%in%bad.trials,] #Remove bad trials
GIDs<-unique(plot.data$GID)

A<-combined.A[rownames(combined.A)%in%GIDs,rownames(combined.A)%in%GIDs]
envs<-as.character(unique(plot.data$Env))


#fill in missing row col coordinate plots with NAs
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
    fill.ins$Year<-unlist(lapply(strsplit(envs[i],split = "_"),function(x) x[1]))
    fill.ins$Site<-unlist(lapply(strsplit(envs[i],split = "_"),function(x) x[2]))
    plot.data<-rbind(plot.data,fill.ins)
  }
}
plot.data<-plot.data[order(plot.data$Env,plot.data$ROW,plot.data$COL),]

asreml.options(workspace="10028mb",trace=TRUE,maxit=200,ai.sing = TRUE) #Set up ASREML workspace

plot.data$Env<-factor(plot.data$Env)
plot.data$PLOT<-factor(plot.data$PLOT)

#Set up empty lists
all.yr.DIAG.mods<-list()
all.yr.FA1.mods<-list()
all.yr.FA2.mods<-list()
all.yr.FA3.mods<-list()

for(i in 1:length(unique(plot.data$Year))){ #for all the trials in each year
  print(paste("Starting",unique(plot.data$Year)[i]))
  yr.sub.plot.data<-plot.data[plot.data$Year==unique(plot.data$Year)[i],] #subset data for trials only in one year
  
  envs<-unique(yr.sub.plot.data$Env) #Get subset of year environments
  #Extract the random effects to fit from the best single trial models
  env.ranefs<-lapply(All.best.mods,function(x) names(x$G.param)[!names(x$G.param)=="GID"]) 
  env.ranefs<-env.ranefs[as.character(envs)]
  env.ranefs<-env.ranefs[unlist(lapply(env.ranefs,function(x) length(x)>0))] 
  #Extract the residual covariance effects to fit from the best single trial models
  env.residefs<-unlist(lapply(All.best.mods,function(x) attributes(x$formulae$residual)$term.labels))
  env.residefs<-env.residefs[as.character(envs)]
  env.residefs<-sapply(unique(env.residefs),function(x) list(names(which(env.residefs==x)))) #restructure list for each combination of residual effects combinations
  names(env.residefs)[names(env.residefs)=="units"]<-"id(ROW):id(COL)"
  
  #Build the DIAG random effects formula
  genetic.form<-"~diag(Env):vm(GID,Asub)+" 
  if(length(env.ranefs)>0){
    DIAGranform<-paste(genetic.form,paste(unlist(lapply(1:length(env.ranefs),function(x) 
      paste("at(Env,'",names(env.ranefs)[x],"'):",env.ranefs[[x]],"+",sep=""))),collapse = ""))
    DIAGranform<-str_sub(string = DIAGranform,start = 1,end = nchar(DIAGranform)-1)
  }else{
    DIAGranform<-str_sub(string = genetic.form,start = 1,end = nchar(genetic.form)-1)
  }
  #Build the residual covariance effects formula
  resid.form<-paste("~dsum(",paste("~",paste(names(env.residefs),collapse = "+"),"|Env",sep=""),",levels=env.residefs)")
  
  Asub<-A[rownames(A)%in%yr.sub.plot.data$GID,colnames(A)%in%yr.sub.plot.data$GID] #Make subset of GRM
  
  #Fit DIAG model
  all.yr.DIAG.mods[[i]]<-asreml(fixed = Grain.Yield..t.ha.~Env,
                                random = as.formula(DIAGranform),
                                residual = as.formula(resid.form),
                                na.action = na.method(x='include'),
                                data=yr.sub.plot.data)
  
  #Build the FA1 random effects formula
  if(length(env.ranefs)>0){
    FA1ranform<-paste("~rr(Env,1):vm(GID,Asub)+diag(Env):vm(GID,Asub)+", #including FA and environment specific (diag) genotype effects
                       paste(unlist(lapply(1:length(env.ranefs),function(x) paste("at(Env,'",names(env.ranefs)[x],"'):",env.ranefs[[x]],"+",sep=""))),collapse = ""))
    FA1ranform<-str_sub(string = FA1ranform,start = 1,end = nchar(FA1ranform)-1)
  }else{
    FA1ranform<-"~rr(Env,1):vm(GID,Asub)+diag(Env):vm(GID,Asub)"
  }
  #Fit FA1 model
  all.yr.FA1.mods[[i]]<-asreml(fixed = Grain.Yield..t.ha.~Env,
                                random = as.formula(FA1ranform),
                                residual = as.formula(resid.form),
                                na.action = na.method(x='include'),
                                data=yr.sub.plot.data)
  #Only models up to FA1 are fitted because only the first factor is useful for starting values for multi-year full models
}
names(all.yr.DIAG.mods)<-unique(plot.data$Year)
names(all.yr.FA1.mods)<-unique(plot.data$Year)

#Check models converged ok
unlist(lapply(all.yr.DIAG.mods,function(x) x$converge))
unlist(lapply(all.yr.FA1.mods,function(x) x$converge))



