#source functions used in code from github
library(devtools)
github = "https://raw.githubusercontent.com/kwiter/mckFUNCTIONS/master/mckFUNCTIONS.r"
source_url(github)

#set working directory
setwd('/home/mckwit/Dropbox/Guided')
#path to your data
path = 'Data/'

#uses custom function from mckFUNCTIONS to open zip into temp folder and pass paths to user
tmp = pathZip(paste0(path,'ld0814.zip'),fileName = character(),readFiles=T)

library(sas7bdat)
dat = read.sas7bdat(tmp$filePath) #read sas file 

#library(haven)
#dat = read_sas(tmp$filePath) # this is supposed to be a faster parser reuires haven

saveRDS(dat,file='Data/lyme.rds') 
dat = readRDS(file='Data/lyme.rds') #use saved r file

library(UScensus2010)
library(maptools)
library(mapproj)
library(rgeos)

library(rgdal)

#read in county level shapefile from zip
tmp = pathZip(paste0(path,'usCountyCensus.zip'),fileName = character(),readFiles=T)
cMap = readOGR(tmp$tmpPath,"gz_2010_us_050_00_20m")
plot(cMap) # check shape
cMap_fips = paste0(cMap@data[,'STATE'],cMap@data[,'COUNTY']) #creates county fips for shapefile

#read in state level shape file from zip
tmp = pathZip('/home/mckwit/Dropbox/Shapes/cb_2014_us_state_500k.zip',fileName = character(),readFiles=T)
sMap = readOGR(tmp$filePath[1],"cb_2014_us_state_500k")
plot(sMap)
sMap_fp =  sMap@data[,'STATEFP'] #create state fips id

#from data
lySt = table(dat[,'state1']) # creat count of incidence by state fip
whr = match(names(lySt),sMap_fp) # match data and map fip
sMapC= rep(0,length(sMap_fp )) # empty vector to hold counts for map
sMapC[whr] = lySt # add counts from data to map vector

cols = colorRange(sMapC/max(sMapC),colours = c('white','#027A40'),trans = 1) #custom function in mckFunctions
plot(sMap,col=cols,xlim = c(-125,-70),ylim=c(25,48))
title(paste ("Lyme Disease Incedence by State 2008-2014"))

cols = colorRange( (as.numeric(cut(sMapC/max(sMapC),c(1,quantile(sMapC/max(sMapC),c(.9,.75,.25)),-1)))-1)/3,colours = c('white','#027A40'),trans = 1)
plot(sMap,col=cols,xlim = c(-160,-40),ylim=c(25,71))

## county data 

#replace missing(999) county level indicators with counties based on incidence probability within state
#you could just eliminate these values 
stNa = unique(dat$stname)
for(i in 1:length(stNa)){
  whr = which(dat$stname == stNa[i])
  probs = table(dat[whr,'county1'])/sum(table(dat[whr,'county1']))
  if('999' %!in% names(probs)) next
  probs = probs[names(probs) != '999']
  for(j in 1:len(whr)) {
    if(a.n(dat[whr[j],'county1']) == 999) dat[whr[j],'fips'] = paste0(dat[whr[j],'state1'],sample(names(probs),1,prob=probs))
  }
  print(stNa[i])  
}

#checl fips data
lyCo = table(dat[,'fips'])
whr = match(names(lyCo),cMap_fips)

#still a few bad fips codes
#rplace fips with good within state codes based on probability
badFips = names(lyCo[which(is.na(whr))])
for(i in 1:len(badFips)){
  whrB = which(dat$fips == badFips[i])
  whr = which(dat$state1 == dat$state1[whrB])
  probs = table(dat[whr,'county1'])/sum(table(dat[whr,'county1']))
  probs = probs[names(probs) != dat$county1[whrB]]
  dat[whrB,'fips'] = paste0(dat[whrB,'state1'] , sample(names(probs),1,prob=probs))
}

saveRDS(dat,file='Data/lyme.rds') #save transfomed data

#match county map fips to data fips
lyCo = table(dat[,'fips'])
whr = match(names(lyCo),cMap_fips)

#map fips data to fips map
cMapC= rep(0,length(cMap_fips ))
cMapC[whr] = lyCo

#plot

cols = colorRange(cMapC/max(cMapC),colours = c('white','#027A40'),trans = 1) #make colors from counts must scale [0,1]
plot(cMap,col=cols,xlim = c(-160,-40),ylim=c(25,71),border=rgb(.5,.5,.5,.2))

par(mar=c(1,1,2,1))
cols = colorRange(cMapC/max(cMapC),colours = c('#ffffb2','#fecc5c','#fd8d3c','#f03b20','#bd0026'),trans = 1)
cols[cMapC == 0] = 'white' # make all places with no counts white
plot(cMap,col=cols,xlim = c(-125,-70),ylim=c(25,48),border=rgb(.5,.5,.5,.2)) #limits control bounding box

#you could color it make it catagorically
cols =colorRange( (as.numeric(cut(sMapC/max(sMapC),c(1,quantile(sMapC/max(sMapC),c(.9,.75,.25)),-1)))-1)/3,colours = c('white','#027A40'),trans = 1)
plot(sMap,col=cols,xlim = c(-160,-40),ylim=c(25,71))

#add a title and legend
title(paste ("Lyme Disease Incedence by County 2008-2014"))

myLegend(x = -75, y=27,h = 7,w = 3,colours = c('#ffffb2','#fecc5c','#fd8d3c','#f03b20','#bd0026'),
         highVal = max(cMapC),lowVal = min(cMapC))


