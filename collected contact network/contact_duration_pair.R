
library(lubridate)
contact_duration_pair <- function(receiver_id)
{
  tmp <- single_BPdata %>%
    filter(receiver==receiver_id) %>% 
    distinct()
  
  tmp_duration <- c(0,abs(diff(tmp$contact_time)))
  tmp_duration<- ifelse(tmp_duration<8,7,tmp_duration)
  tmp_index <- rep(c(1:length(rle(tmp_duration)$length)),rle(tmp_duration)$lengths)
  tmp$contact_index <- tmp_index;rm(tmp_index,tmp_duration)
  tmp$hour <- hour(tmp$contact_time)
  tmp <- tmp[,-which(colnames(tmp)=='contact_time')]
  t <- tmp %>% group_by(sender,receiver,contact_date,contact_wday,contact_index)%>%
    summarise(start_hour = min(hour))
  tt <- left_join(tmp,t,by=c("sender","receiver","contact_date","contact_wday","contact_index"));rm(t,tmp)
  subUniSenderBP <- tt %>% group_by(sender,receiver,contact_date,contact_wday,
                                    contact_index,start_hour) %>%
    summarise(contact_consecutive_duration = n()*7)
  subUniPairBP <- subUniSenderBP[,-which(colnames(subUniSenderBP)=='contact_index')]
  subUniPairBP
}
# calculate the consecutive contact duration 
# fid=1
library(foreach)
library(doParallel)
registerDoParallel(18)
allPairData <- NULL
a <- Sys.time()
# for(fid in 1:20)#length(filelist))
allPairData <- foreach(fid=1:length(filelist),
                       .combine = rbind,
                       .packages=c('dplyr','lubridate','foreach','doParallel')) %dopar%
{
  oneSenderData <- NULL
  single_BPdata <- read.csv(file.path('data','bluepass',filelist[fid]),as.is=T)
  single_BPdata <- convert_data_format(data = single_BPdata,dates_interval = dates_interval,saveTocsv = F)
  single_BPdata$contact_date <- as.Date(single_BPdata$contact_time)
  tmp_wday <- wday(single_BPdata$contact_date)-1
  tmp_wday[which(tmp_wday==0)] <- 7
  single_BPdata$contact_wday <- tmp_wday;rm(tmp_wday)
  uni_tmp_rid <- unique(single_BPdata$receiver)
  # for(rId in uni_tmp_rid){
  oneSenderData <- foreach(rId = 1:length(uni_tmp_rid),
                           .combine=rbind,
                           .packages=c('dplyr','lubridate')) %dopar%{
    pairContactDurationDat <- contact_duration_pair(receiver_id=uni_tmp_rid[rId])
    pairContactDurationDat
    # oneSenderData <- rbind(oneSenderData,pairContactDurationDat)
  }
  # allPairData <- rbind(allPairData,oneSenderData)
  save(oneSenderData,file = paste0(foldername,'contact_pair/sender_',USID[fid],'.Rdata'))
  oneSenderData
}

b=Sys.time()
cat("Time used to calculate the consecutive contact duration: ",round(as.numeric(difftime(b,a,units = 'mins')),1),'mins.\n',file = notefile,sep="\n",append = T)
allPairSumData <- allPairData %>% group_by(sender,receiver,contact_date,contact_wday) %>%
  summarise(mean_contact_consecutive_duration = mean(contact_consecutive_duration))
save(allPairData,file=paste0(foldername,'contact_pair/rawPairs.Rdata')) # to get the start_hour for the fig 1. add this saving on Oct10
save(allPairSumData,file = paste0(foldername,'contact_pair/toDefineContactType.Rdata'))
rm(allPairData,single_BPdata,fid,allPairSumData)


