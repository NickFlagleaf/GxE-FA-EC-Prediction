

##This work © 2024 by Nick Fradgley, CSIRO is licensed under CC BY 4.0 
##'new.default.SVs' is the 'vparameters.table' form the starting values for the new model
##'old.vparams' is the old model 'vparameters'

FA.nxt.SVs<-function(old.vparams,
                     new.default.SVs,
                     add.fa.str="vm\\(",
                     non.add.fa.str="\\:GID!",
                     new.add.FA.order,
                     new.nonadd.FA.order){
  rownames(new.default.SVs)<-new.default.SVs$Component
  new.default.SVs[names(old.vparams)[names(old.vparams)%in%rownames(new.default.SVs)],"Value"]<-old.vparams[names(old.vparams)%in%rownames(new.default.SVs)]#Replace var params with the same names
  old.fas<-old.vparams[grepl("\\!fa",names(old.vparams))] #Get just the FA params for the old model
  old.fa.add.order<-max(as.numeric(sapply(names(old.fas)[grep(add.fa.str,names(old.fas))], function(x) stringr::str_sub(x,nchar(x),nchar(x))))) #work out the old model FA additive order
  old.fa.nonadd.order<-max(as.numeric(sapply(names(old.fas)[grep(non.add.fa.str,names(old.fas))], function(x) stringr::str_sub(x,nchar(x),nchar(x))))) #work out the old model FA non-additive order
  names(old.fas)<-gsub(pattern = paste(old.fa.add.order,"):",add.fa.str,sep=""),
                       replacement = paste(new.add.FA.order,"):",add.fa.str,sep=""), names(old.fas)) #update the additive fa order number to new one
  names(old.fas)<-gsub(pattern = paste(old.fa.nonadd.order,"):",non.add.fa.str,sep=""),
                       replacement = paste(new.nonadd.FA.order,"):",non.add.fa.str,sep=""), names(old.fas)) #update the non-additive fa order number to new one
  new.default.SVs[names(old.fas),"Value"]<-old.fas
  new.default.SVs$Value[new.default.SVs$Constraint=='P' & new.default.SVs$Value<1e-5]<-1e-4 #Move small components off the boundary
  return(new.default.SVs)
  }