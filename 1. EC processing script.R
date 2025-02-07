##########################################################################################
####### Scripts for deriving Environmental Covariates from weather and soil data #########
##This work © 2025 by Nick Fradgley, CSIRO is licensed under CC BY 4.0 ###################

# 'all.weather' list of weather data variables: ("Rain","Maxtemp","Mintemp","Solar radiation","Vapor pressure deficit")
# each with list of years
# each with 3D array values with dimensions of Lon, Lat and day of the year.


# 'site.info' is a data frame of info for each trial environment
#column names are: "Index","Environment","Year","Site","Lat","Long","Growing.region","Sowing.date","Min.heading.date",
#"Median.heading.date","Max.heading.date","Median.date.of.maturity","Harvest.date","Mean.Zadocks.score","Zadocks.score.date".

library(stringr)
library(chillR)
library(viridis)

#sort date formats
site.info$Sowing.date<-as.Date(site.info$Sowing.date,format = "%d/%m/%Y")
median.sowing.day.of.year<-as.numeric(median(na.omit(site.info$Sowing.date-as.Date(paste(site.info$Year,"-01-01",sep=""))))) #Get the median sowing day of the year
site.info$Sowing.date[is.na(site.info$Sowing.date)]<-as.Date(paste(site.info$Year[is.na(site.info$Sowing.date)],"-01-01",sep=""))+median.sowing.day.of.year
site.info$Min.heading.date<-as.Date(site.info$Min.heading.date,format = "%d/%m/%Y")
site.info$Median.heading.date<-as.Date(site.info$Median.heading.date,format = "%d/%m/%Y")
site.info$Harvest.date<-as.Date(site.info$Harvest.date,format = "%d/%m/%Y")
site.info$Zadocks.score.date<-as.Date(site.info$Zadocks.score.date,format = "%d/%m/%Y")

site.info$DTH<-as.numeric(site.info$Median.heading.date-site.info$Sowing.date) #Calculate days to heading

Envs<-site.info$Environment
rownames(site.info)<-Envs

#Work out x and y coordinate values to index against array
xcoords<-sapply(site.info$Long,function(x) which.min(abs(x-as.numeric(unlist(dimnames(all.weather$Maxtemp$'2015')[1])))))
ycoords<-sapply(site.info$Lat,function(x) which.min(abs(x-as.numeric(unlist(dimnames(all.weather$Maxtemp$'2015')[2])))))

#Work out start day of the year for each env
startday<-as.numeric(site.info$Sowing.date-as.Date(paste(site.info$Year,"-01-01",sep="")))
names(startday)<-Envs

#Get daily weather for each trial
all.env.weather<-list()
for(i in 1:length(Envs)){ #For each trial environment:
  cat("|")
  #Get min and max temps per day for the trial year
  max.temps<-all.weather$Maxtemp[[which(names(all.weather$Maxtemp)==site.info[Envs[i],"Year"])]][xcoords[i],ycoords[i],1:364]
  min.temps<-all.weather$Mintemp[[which(names(all.weather$Mintemp)==site.info[Envs[i],"Year"])]][xcoords[i],ycoords[i],1:364]
  #Get rainfall per day for the trial year
  rain<-all.weather$Rain[[which(names(all.weather$Rain)==site.info[Envs[i],"Year"])]][xcoords[i],ycoords[i],1:364]
  #Get solar radiation per day for the trial year
  solrad<-all.weather$`Solar radiation`[[which(names(all.weather$`Solar radiation`)==site.info[Envs[i],"Year"])]][xcoords[i],ycoords[i],1:364]
  #get day lengths
  dl<-daylength(latitude = site.info$Lat[site.info$Environment==Envs[i]]
                ,JDay = 1:364 , notimes.as.na = FALSE)
  #Get vapour pressure deficit
  VPD<-all.weather$`Vapor pressure deficit`[[which(names(all.weather$`Vapor pressure deficit`)==site.info[Envs[i],"Year"])]][xcoords[i],ycoords[i],1:364]
  
  all.env.weather[[i]]<-data.frame("Max_temp"=max.temps,
                 "Min_temp"=min.temps,
                 "Rain"=rain,
                 "Sol_rad"=solrad,
                 "Day_length"=dl$Daylength,
                 "VPD"=VPD)
  }
names(all.env.weather)<-paste("Yr_",Envs,sep="")


#Make TT degree days function
TTfun<-function(Tci){
  if(0<Tci & Tci<=26){out=Tci}
  if(26<Tci & Tci<=34){out=(26/8)*(34-Tci)}
  if(Tci<0 | Tci>34){out=0}
  return(out)
}

#Calculate daily thermal time for each day for each environment
for(i in 1:length(Envs)){
  Tc<-((all.env.weather[[i]]$Max_temp+all.env.weather[[i]]$Min_temp)/2)
  DailyTT<-sapply(Tc,TTfun)
  all.env.weather[[i]]$DailyTT<-DailyTT
}

#TT of heading growth stage by region....
aggregate(x = sapply(1:length(all.env.weather),function(x) cumsum(all.env.weather[[x]]$DailyTT[startday[x]:364])[site.info$DTH[x]]),
          by=list(site.info$Growing.region),FUN=function(x) mean(na.omit(x)))
reg.DTH.TT<-c("North"=1361,"South"=1301,"West"=1194)


#Estimate the dates of crop growth stages based on sowing dates and accumulated thermal time
stage.names<-c("Sowing","Emergence","End of juvenile","Heading","Flowering","Start of grain filling","End of grain filling","Maturity")
all.env.stages<-matrix(NA,nrow=length(Envs),ncol=length(stage.names),dimnames = list(Envs,stage.names))
for(i in 1:length(Envs)){ #For each environment:
 DailyTT<-all.env.weather[[i]]$DailyTT
 
 stages<-rep(NA,8)
 names(stages)<-stage.names
 stages[1]<-1
 wetsow<-sum(all.env.weather[[i]]$Rain[c(startday[i]-7):startday[i]])>3 #Was it a wet sow?
 if(wetsow){stages[2]<-14}else{stages[2]<-min(which(all.env.weather[[i]]$Rain[startday[i]:364]>0))+14} #Emergence two weeks after a wet sowing
 stages[3]<-which.min(abs(cumsum(DailyTT[startday[i]:364])- (cumsum(DailyTT[startday[i]:364])[stages[2]]+500)))
 stages[4]<-which.min(abs(cumsum(DailyTT[startday[i]:364])-reg.DTH.TT[site.info$Growing.region[i]]))  #different TT to DTH per region
 stages[5]<-which.min(abs(cumsum(DailyTT[startday[i]:364])- (cumsum(DailyTT[startday[i]:364])[stages[4]]+250)))
 stages[6]<-which.min(abs(cumsum(DailyTT[startday[i]:364])- (cumsum(DailyTT[startday[i]:364])[stages[5]]+250)))
 stages[7]<-which.min(abs(cumsum(DailyTT[startday[i]:364])- (cumsum(DailyTT[startday[i]:364])[stages[6]]+250)))
 stages[8]<-which.min(abs(cumsum(DailyTT[startday[i]:364])- (cumsum(DailyTT[startday[i]:364])[stages[7]]+400)))
 
 all.env.stages[i,]<-stages
 all.env.stages<-as.data.frame(all.env.stages)
 }



#Define stress covariates according to estimated growth intervals
interval.names<-c("Sow2Emer","Emer2Juv","Juv2He","He2Flw","Flw2Sgf","Sgf2Egf","Egf2mat")


#N days per stage
Ndays.per.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),function(s) all.env.stages[e,s]-all.env.stages[e,s-1])))
colnames(Ndays.per.stage)<-paste("Ndays_",interval.names,sep="")

#N days from sowing to flowering
Ndays_Sow2Flw<-all.env.stages[,5]
names(Ndays_Sow2Flw)<-Envs

#N days from flowering to end of grain fill
Ndays_Flw2Egf<-all.env.stages[,7]-all.env.stages[,5]
names(Ndays_Flw2Egf)<-Envs

#total rain per stage
Sum.rain.per.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                      function(s) sum(all.env.weather[[e]]$Rain[startday[e]:364][all.env.stages[e,s-1]:all.env.stages[e,s]]))))
colnames(Sum.rain.per.stage)<-paste("TotRain_",interval.names,sep="")

#total rain from sowing to flowering
TotRain_Sow2Flw<-sapply(1:length(Envs),function(e) sum(all.env.weather[[e]]$Rain[startday[e]:364][all.env.stages$Sowing[e]:all.env.stages$Flowering[e]]))
names(TotRain_Sow2Flw)<-Envs

#total rain from flowering to end of grain fill
TotRain_Flw2Egf<-sapply(1:length(Envs),function(e) sum(all.env.weather[[e]]$Rain[startday[e]:364][all.env.stages$Flowering[e]:all.env.stages$`End of grain filling`[e]]))
names(TotRain_Flw2Egf)<-Envs

#Total stored soil moisture before sowing
TotRain_priorSow<-sapply(1:length(Envs),function(e) sum(all.env.weather[[e]]$Rain[1:startday[e]]))
names(TotRain_priorSow)<-Envs

#mean temp per stage
mean.temp.per.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                                                               function(s) mean(c((all.env.weather[[e]]$Max_temp[startday[e]:364]+
                                                                                  all.env.weather[[e]]$Min_temp[startday[e]:364])/2)[all.env.stages[e,s-1]:all.env.stages[e,s]]))))
colnames(mean.temp.per.stage)<-paste("Avtemp_",interval.names,sep="")


#mean max temp per stage
mean.max.temp.per.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                                                                function(s) mean(c(all.env.weather[[e]]$Max_temp[startday[e]:364][all.env.stages[e,s-1]:all.env.stages[e,s]])))))
colnames(mean.max.temp.per.stage)<-paste("Avmaxtemp_",interval.names,sep="")

#mean min temp per stage
mean.min.temp.per.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                                                                    function(s) mean(c(all.env.weather[[e]]$Min_temp[startday[e]:364][all.env.stages[e,s-1]:all.env.stages[e,s]])))))
colnames(mean.min.temp.per.stage)<-paste("Avmintemp_",interval.names,sep="")

#Sum dry days per stage
Sum.drydays.per.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                                                             function(s) sum(c(all.env.weather[[e]]$Rain[startday[e]:364]==0)[all.env.stages[e,s-1]:all.env.stages[e,s]]))))
colnames(Sum.drydays.per.stage)<-paste("Ndd_",interval.names,sep="")

#Sum minT < 0C per stage
MinTbelow0per.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                                                                  function(s) sum(c(all.env.weather[[e]]$Min_temp[startday[e]:364]<0)[all.env.stages[e,s-1]:all.env.stages[e,s]]))))
colnames(MinTbelow0per.stage)<-paste("Ndays<0_",interval.names,sep="")

#Frosts at flowering
MinTbelow0FLW<-sapply(1:length(Envs),function(e) sum(c(all.env.weather[[e]]$Min_temp[startday[e]:364]<0)[(all.env.stages$Flowering[e]-7):(all.env.stages$Flowering[e]+7)]))
names(MinTbelow0FLW)<-Envs

#warm days
MaxToverr26per.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                                                                 function(s) sum(c(all.env.weather[[e]]$Max_temp[startday[e]:364]>26)[all.env.stages[e,s-1]:all.env.stages[e,s]]))))
colnames(MaxToverr26per.stage)<-paste("Ndays>26_",interval.names,sep="")

#Hot days
MaxToverr34per.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                                                                 function(s) sum(c(all.env.weather[[e]]$Max_temp[startday[e]:364]>34)[all.env.stages[e,s-1]:all.env.stages[e,s]]))))
colnames(MaxToverr34per.stage)<-paste("Ndays>34_",interval.names,sep="")

#mean day lengths per stage
MmeanDLper.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                                                                 function(s) mean(c(all.env.weather[[e]]$Day_length[startday[e]:364])[all.env.stages[e,s-1]:all.env.stages[e,s]]))))
colnames(MmeanDLper.stage)<-paste("AveDL_",interval.names,sep="")

#mean solar radiation per stage
Mmeansolradper.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                                                             function(s) mean(c(all.env.weather[[e]]$Sol_rad[startday[e]:364])[all.env.stages[e,s-1]:all.env.stages[e,s]]))))
colnames(Mmeansolradper.stage)<-paste("AveSR_",interval.names,sep="")

#mean solar radiation per stage
MeanVPDper.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                                                                 function(s) mean(c(all.env.weather[[e]]$VPD[startday[e]:364])[all.env.stages[e,s-1]:all.env.stages[e,s]]))))
colnames(MeanVPDper.stage)<-paste("AveVPD_",interval.names,sep="")

avePQ.per.stage<-t(sapply(1:length(Envs),function(e) sapply(2:ncol(all.env.stages),
                                          function(s) mean(c(all.env.weather[[e]]$Sol_rad[startday[e]:364]*.47)[all.env.stages[e,s-1]:all.env.stages[e,s]]/
                                                 c((all.env.weather[[e]]$Max_temp[startday[e]:364]+all.env.weather[[e]]$Min_temp[startday[e]:364])/2)[all.env.stages[e,s-1]:all.env.stages[e,s]]))))
colnames(avePQ.per.stage)<-paste("PQ_",interval.names,sep="")


#Combine all ECs into weather matrix
Wmat<-cbind(Ndays.per.stage,
            Ndays_Sow2Flw,
            Ndays_Flw2Egf,
            Sum.rain.per.stage,
            TotRain_Sow2Flw,
            TotRain_Flw2Egf,
            TotRain_priorSow,
            mean.temp.per.stage,
            mean.min.temp.per.stage,
            "Mintemp<0_Flw"=MinTbelow0FLW,
            mean.max.temp.per.stage,
            Sum.drydays.per.stage,
            MinTbelow0per.stage,
            MaxToverr26per.stage,
            MaxToverr34per.stage,
            MmeanDLper.stage,
            Mmeansolradper.stage,
            MeanVPDper.stage,
            avePQ.per.stage)
rownames(Wmat)<-Envs

Wmat<-Wmat[,!apply(Wmat,2,var)==0]




#Get soil covariates from SLGA database using the SLGACloud package
library(SLGACloud)
library(terra)

atts<-unlist(c(getParameterValues(Parameter = "Attribute"))) #Get possible soil attributes
atts<-atts[1:19]
atts<-atts[!atts%in%c("Depth of Regolith","Soil Organic Carbon (1\" resolution) ","Effective Cation Exchange Capacity")] #Don't use these
latlon<-site.info[,c("Long","Lat")]
all.env.soil<-data.frame("Env"=site.info$Environment)

for(a in 1:length(atts)){ #For each soil attribute
  print(paste("Starting",atts[a]))
  #Get list of raster options
  rasters <- getProductMetaData(Detail = 'High', Attribute = atts[a], 
                                Component = 'Modelled-Value',
                                isCurrentVersion = 1)
  depths<-paste(rasters$UpperDepth_m,"-",rasters$LowerDepth_m,"m",sep="") #Get available depths
  for(d in 1:length(depths)){
    print(paste("At",depths[d]))
    r<-NULL
    try(r<-rast(rasters$StagingPath[d])) #This part takes a while and sometimes doesn't work....
    if(!is.null(r)){
    vals<-extract(r,latlon) #Extract values for each environment from raster. This also takes a while... 
    colnames(vals)[2]<-paste(atts[a],depths[d],sep="_") #Name soil EC
    all.env.soil<-cbind.data.frame(all.env.soil,vals)
    }
    gc()
    }
}

all.env.soil<-all.env.soil[,!colnames(all.env.soil)=="ID"]
rownames(all.env.soil)<-all.env.soil$Env
all.env.soil<-all.env.soil[,!colnames(all.env.soil)=="Env"]
all.env.soil<-apply(all.env.soil,2,function(x){
  med<-median(na.omit(x))
  newx<-x
  newx[is.na(newx)]<-med
  return(newx)
  })

all.env.soil<-all.env.soil[,apply(all.env.soil,2,var)>0] #Remove ECs with zero variance
#Edit soil EC names
colnames(all.env.soil)[colnames(all.env.soil)%in%c("Soil Organic Carbon Fractions_0-0.05m",
                                                "Soil Organic Carbon Fractions_0.05-0.15m",
                                                "Soil Organic Carbon Fractions_0.15-0.3m")]<-c("Soil Organic Carbon Fractions_MAOC_0-0.05m",
                                                                                               "Soil Organic Carbon Fractions_MAOC0.05-0.15m",
                                                                                               "Soil Organic Carbon Fractions_MAOC0.15-0.3m")
colnames(all.env.soil)[colnames(all.env.soil)%in%c("Soil Organic Carbon Fractions_0-0.05m.1",
                                                   "Soil Organic Carbon Fractions_0.05-0.15m.1")]<-c("Soil Organic Carbon Fractions_POC_0-0.05m",
                                                                                                  "Soil Organic Carbon Fractions_POC0.05-0.15m")
colnames(all.env.soil)[colnames(all.env.soil)%in%c("Soil Organic Carbon Fractions_0-0.05m.2",
                                                   "Soil Organic Carbon Fractions_0.05-0.15m.2",
                                                   "Soil Organic Carbon Fractions_0.15-0.3m.2")]<-c("Soil Organic Carbon Fractions_PYOC_0-0.05m",
                                                                                                  "Soil Organic Carbon Fractions_PYOC_0.05-0.15m",
                                                                                                  "Soil Organic Carbon Fractions_PYOC_0.15-0.3m")


  #Combine weather and soil matrices
  WSmat<-cbind(Wmat,all.env.soil)

