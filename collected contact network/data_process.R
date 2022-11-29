library(lubridate)
library(dplyr)
# library(tidyverse)
# Convert data format -----------------------------------------------------
convert_data_format <- function(data,dates_interval,saveTocsv=T,folder=NULL)
{
  if(dim(data)[1]==0|dim(data)[2]!=7) {stop("File data importation")}
  
  dat <- data[,-c(1,2,6,7)]
  colnames(dat) <- c('sender', 'receiver', 'contact_time')
  # format the variables class
  # dat$report_time <- as.POSIXct(stringr::str_replace_all(dat$report_time,c("'"),""),"%Y-%m-%d %H:%M:%S",tz='GMT')
  dat$sender <- stringr::str_replace_all(dat$sender,c("'"),"")
  dat$receiver <- stringr::str_replace_all(dat$receiver,c("'"),"")
  dat$contact_time <- as.POSIXct(stringr::str_replace_all(dat$contact_time,c("'"),""),"%Y-%m-%d %H:%M:%S",tz='GMT')
  # dat$battery_level <- stringr::str_replace_all(dat$battery_level,c("'"),"")
  dat <- dat%>%filter(contact_time%within%dates_interval)
  if(saveTocsv)
  {if(is.null(folder))
  {
    cat("The converted data will be saved into data folder!")
    write.csv(dat,"data/analysis_data_all.csv",row.names = F)
  }else{
    cat("The converted data will be saved into",folder,'folder!')
    write.csv(dat,paste0(folder,"/analysis_data_all.csv"),row.names = F)
  }
  }
  return(dat)
}
#dates_interval <-interval(as.Date('2021-02-13'),as.Date('2021-04-06'))#53 days
# sink(notefile)
cat('Study time window: 2021-02-13 to 2021-04-06. \n')

# Initial settings-------------------------------------------------------------------------------
# acacia = F
# cassia = T
if(acacia){
  notefile = 'notes_acacia.txt'
  foldername = 'working/acacia/'
  outputfolder = 'output/acacia/'
  survey_a <- read.csv('data/survey/dataentry.csv', as.is=T)
  cat('--------------survey forms details------------------',file=notefile,sep="\n",append=T)
  cat('Total number of survey forms in acacia:',nrow(survey_a),file=notefile,sep="\n",append=T)
  survey_a$tokenID <- paste0('0000',survey_a$QR)
  if(any(survey_a$tokenID=="0000")){
    survey_a <- survey_a[-which(survey_a$tokenID=="0000"),]
  }
  survey_a <- survey_a[,-which(colnames(survey_a)%in%c('Q4','Q6','DE.completed','QR'))]
  names(survey_a) <- c('seq_id','age','country','blk','room','Q5','tokenID')
  cat('Total number of scanable BP QR in acacia:',nrow(survey_a),file=notefile,sep="\n",append=T)
  possiblefiles = list.files('data/bluepass') #data/bluepass is the folder to save all extracted data
  filelist = na.omit(stringr::str_extract(possiblefiles,"BP\\s.+\\s(\\d+).txt.csv"))
  tt <- paste0('0000',substring(filelist,4,11),sep='')
  filelist <- filelist[which(tt%in%survey_a$tokenID)]
  tt1 <- paste0('0000',substring(filelist,4,11),sep='')
  survey_a$br <- paste0(survey_a$blk,'-',survey_a$room,sep='')
  rm_id_room <- which(is.na(str_extract(survey_a$br,paste0("[0-9]{3}",'\\-',"[0-9]{4}"))))#unique(survey_a$br)[c(6,20,30,41,45,54,75,93,96,118,119,128,129,144,146)]
  survey_id_allroom <- survey_a$tokenID[-rm_id_room]#[-which(survey_a$br%in%rm_id_room)]
  # survey_id_demo <- survey_a$tokenID[-which(survey_a$age==''|survey_a$country=='')]
  survey_id_keep <- survey_id_allroom[which(survey_id_allroom%in%tt1)]
  cat('Total number of valid BP with valid unit info in acacia:',length(survey_id_keep),file=notefile,sep="\n",append=T)
  rm(rm_id_room,tt,tt1,survey_id_allroom)
}
if(cassia){
  notefile = 'notes_cassia.txt'
  foldername = 'working/cassia/'
  outputfolder = 'output/cassia/'
  survey_c <- read.csv('data/survey/dataentry_cassia.csv', as.is=T)
  cat('Total number of survey forms in cassia:',nrow(survey_c),file= notefile,sep="\n",append=T)
  survey_c$tokenID <- paste0('0000',survey_c$QR)
  if(any(survey_c$tokenID=="0000")){
    survey_c <- survey_c[-which(survey_c$tokenID=="0000"),]
  }
  survey_c <- survey_c[,-which(colnames(survey_c)%in%c('Q4','Q6','DE.completed','QR'))]
  names(survey_c) <- c('seq_id','age','country','blk','room','Q5','tokenID')
  cat('Total number of scanable BP QR in cassia:',nrow(survey_c),file=notefile,sep="\n",append=T)
  possiblefiles = list.files('data/bluepass') #data/bluepass is the folder to save all extracted data
  filelist = na.omit(stringr::str_extract(possiblefiles,"BP\\s.+\\s(\\d+).txt.csv"))
  tt <- paste0('0000',substring(filelist,4,11),sep='')
  filelist <- filelist[which(tt%in%survey_c$tokenID)]
  tt1 <- paste0('0000',substring(filelist,4,11),sep='')
  survey_c$br <- paste0(survey_c$blk,'-',survey_c$room,sep='')
  rm_id_room <- which(is.na(str_extract(survey_c$br,paste0("[0-9]{2}",'\\-',"[0-9]{3}"))))#unique(survey_c$br)[c(6,19,22,29,31,36,66,77,81,83,84,100,108,121,128,135)]
  survey_id_allroom <- survey_c$tokenID[-rm_id_room]#[-which(survey_c$br%in%rm_id_room)]
  # survey_id_demo <- survey_c$tokenID[-which(survey_c$age==''|survey_c$country=='')]
  survey_id_keep <- survey_id_allroom[which(survey_id_allroom%in%tt1)]
  cat('Total number of valid BP with valid unit info in cassia:',
      length(survey_id_keep),file=notefile,sep="\n",append=T)
  rm(rm_id_room,tt,tt1,survey_id_allroom)
}

#-------------------------------------------------------------------------------------------
# sink()
if(!acacia & !cassia){
  fid_duplicated <- which(duplicated(substring(filelist,4,11)))
  if(length(fid_duplicated)>0)
  {
    filelist <- filelist[-fid_duplicated]
  };rm(fid_duplicated,possiblefiles)
  ## we only keep the BP data which is valid (i.e. consistent with survey_id_keep)
  tt <- paste0('0000',substring(filelist,4,11),sep='')
  keep_file_id <- which(tt%in%survey_id_keep)
  filelist <- filelist[keep_file_id]
  rm(tt,keep_file_id)
  URID=USID <- c()
  for(fid in 1:length(filelist))
  {
    single_BPdata <- read.csv(file.path('data','bluepass',filelist[fid]),as.is=T)
    single_BPdata <- convert_data_format(data = single_BPdata,dates_interval = dates_interval,saveTocsv = F)
    usid <- unique(single_BPdata$sender)
    USID <- unique(c(USID,usid))
  }
  USID <- stringr::str_replace_all(USID,c("'"),"") # all USID are located in the studing date interval.
  cat('Total number of sender in BP data pool with valid information:',length(USID),file = notefile,sep="\n",append = T)
  rm(usid,fid)
  tt <- paste0('0000',substring(filelist,4,11),sep='')
  fid_Nperiod <- which(!tt%in%USID) # remove sender which BP not in the studying date interval.
  if(length(fid_Nperiod)>0)
  {
    filelist <- filelist[-fid_Nperiod]
  };rm(fid_Nperiod)
  # get the URID based on the studying date interval.
  for(fid in 1:length(filelist))
  {
    single_BPdata <- read.csv(file.path('data','bluepass',filelist[fid]),as.is=T)
    single_BPdata <- convert_data_format(data = single_BPdata,dates_interval = dates_interval,saveTocsv = F)
    urid <- unique(single_BPdata$receiver)
    URID <- unique(c(URID,urid))
  }
  URID <- stringr::str_replace_all(URID,c("'"),"")
  save(USID,file=paste0(foldername,'USID.Rdata'))
  save(URID,file=paste0(foldername,'URID.Rdata'))
  cat('Total number of receiver in BP data pool:',length(URID),file = notefile,sep="\n",append = T)
}


# customized when acacia = T & cassia = T of published paper analysis
if(acacia&cassia){
  # notefile='notes_combine.txt'
  foldername = 'working/combine/'
  outputfolder = 'output/combine/'
  possiblefiles = list.files('data/bluepass') #data/bluepass is the folder to save all extracted data
  filelist = na.omit(stringr::str_extract(possiblefiles,"BP\\s.+\\s(\\d+).txt.csv"))
  
  ## about the survey data:
  tt1 <- paste0('0000',substring(filelist,4,11),sep='')
  survey_a <- read.csv('data/survey/dataentry.csv', as.is=T)
  survey_a$tokenID <- paste0('0000',survey_a$QR)
  if(any(survey_a$tokenID=="0000")){
    survey_a <- survey_a[-which(survey_a$tokenID=="0000"),]
  }
  survey_a <- survey_a[,-which(colnames(survey_a)%in%c('Q4','Q6','DE.completed','QR'))]
  names(survey_a) <- c('seq_id','age','country','blk','room','Q5','tokenID')
  survey_a$br <- paste0(survey_a$blk,'-',survey_a$room,sep='')
  rm_id_room <- which(is.na(str_extract(survey_a$br,paste0("[0-9]{3}",'\\-',"[0-9]{4}"))))#unique(survey_a$br)[c(6,20,30,41,45,54,75,93,96,118,119,128,129,144,146)]
  survey_id_allrooma <- survey_a$tokenID[-rm_id_room]#[-which(survey_a$br%in%rm_id_room)]
  survey_a$dorm <- 'A'
  
  survey_c <- read.csv('data/survey/dataentry_cassia.csv', as.is=T)
  survey_c$tokenID <- paste0('0000',survey_c$QR)
  if(any(survey_c$tokenID=="0000")){
    survey_c <- survey_c[-which(survey_c$tokenID=="0000"),]
  }
  survey_c <- survey_c[,-which(colnames(survey_c)%in%c('Q4','Q6','DE.completed','QR'))]
  names(survey_c) <- c('seq_id','age','country','blk','room','Q5','tokenID')
  survey_c$br <- paste0(survey_c$blk,'-',survey_c$room,sep='')
  rm_id_room <- which(is.na(str_extract(survey_c$br,paste0("[0-9]{2}",'\\-',"[0-9]{3}"))))#unique(survey_c$br)[c(6,19,22,29,31,36,66,77,81,83,84,100,108,121,128,135)]
  survey_id_allroomc <- survey_c$tokenID[-rm_id_room]#[-which(survey_c$br%in%rm_id_room)]
  survey_c$dorm <- "C"
  survey_id_allroom <- c(survey_id_allrooma,survey_id_allroomc)
  survey_id_keep <- survey_id_allroom[which(survey_id_allroom%in%tt1)]
  survey <- rbind(survey_a,survey_c)
  # print(nrow(survey))

  rm(survey_id_allrooma,survey_id_allroomc,rm_id_room,survey_id_allroom,tt1)
  
  ## about the URID and USID
  load('working/acacia/USID.Rdata')
  acacia_USID <- USID
  load('working/acacia/URID.Rdata')
  acacia_URID <- URID
  load('working/cassia/USID.Rdata')
  cassia_USID <- USID
  load('working/cassia/URID.Rdata')
  cassia_URID <- URID
  
  USID <- c(acacia_USID,cassia_USID)
  URID <- c(acacia_URID,cassia_URID)
  tt <- paste0('0000',substring(filelist,4,11),sep='')
  filelist <- filelist[-which(duplicated(tt))]
  tt <- paste0('0000',substring(filelist,4,11),sep='')
  
  keep_file <- which(tt %in% USID)
  filelist <- filelist[keep_file]
  ## keep all the valid survey data which is within the study period
  survey <- survey %>% filter(survey$tokenID %in% USID)
  survey <- survey[-which(duplicated(survey$tokenID)),]
  rm_sender1 <- survey$tokenID[which(survey$age=='')]
  rm_sender2 <- survey$tokenID[which(survey$country=='')]
  rm_sender <- c(rm_sender1,rm_sender2);rm(rm_sender1,rm_sender2)
  survey <- survey[-which(survey$tokenID%in%rm_sender),]
  cat('Total number of valid BP with valid unit info in two dorms:',
      nrow(survey))#,file=notefile,sep="\n",append=T)
  rm(survey_a,survey_c,acacia_URID,acacia_USID,cassia_URID,cassia_USID,tt)
}


