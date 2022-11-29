# In the following analysis, we will combine two dormitory together.
# As for table 2, dataset was only considered by all completed paired data
# without demographic consideration.
library(dplyr)
foldername ='working/combine/'
load(paste0(foldername,'PairDataContactType.Rdata'))
# names(allPairSumData)
map_regular_index <- allPairSumData%>%select(sender,receiver,regular_contact,regular_index) 
map_regular_index <- unique(map_regular_index[,-1])
rm(allPairSumData)
# the all Pair Data should be used?
load(paste0('working/combine/','rawPairs.Rdata'))
load('working/map_regular_index.rdata')
temp <- left_join(allPairData,map_regular_index,by=c('sender','receiver'))
rm(allPairData,map_regular_index)
## Totally ------
# t <- table(allPairSumData$regular_index)#regular_contact)
# 0       1 
# 8869119  464295 # regular_contact
# 0        1 
# 8644964  688450 
table(temp$regular_index)
# 
# 0        1 
# 24675160  8658311 
temp%>%group_by(regular_index,contact_wday)%>%summarise(mean(n()))

table(temp$regular_index) # used in this table generation, to show

# 0        1 
# 31155110  8658311

# totNumR <- as.numeric(t[2])
# totNumT <- as.numeric(t[1])
# rm(t)
# 
# temp <- temp[,-which(colnames(temp)%in%'regular_contact')]
# temp <- unique(temp)


t_save_encounter <- temp %>% group_by(sender, receiver, contact_date, regular_index) %>% 
  summarise(contact_duration = mean(contact_consecutive_duration)) %>% ungroup %>%#sum(contact_consecutive_duration)) %>% ungroup %>%
  group_by(sender,receiver,regular_index) 
t_save_encounter %>%group_by(regular_index) %>% summarise(mean(contact_duration),sd(contact_duration))

t_save <- temp %>% group_by(sender, receiver, contact_date, regular_index) %>% 
  summarise(contact_duration = sum(contact_consecutive_duration)) %>% ungroup %>%#sum(contact_consecutive_duration)) %>% ungroup %>%
  group_by(sender,receiver,regular_index) 
t_save %>%group_by(regular_index) %>% summarise(mean(contact_duration),sd(contact_duration))

t_save_freq <- t_save%>%group_by(contact_date,regular_index)%>%summarise(uniq_pairs=n_distinct(paste0(sender,'-',receiver)))
t_save_freq %>% group_by(regular_index)%>%summarise(uniq_pairs_perD=mean(uniq_pairs))

t_save_dur <- t_save %>% summarise(mean_contact_duration_inday=sum(contact_duration)/52/60)
t_save_dur1 <- t_save_encounter%>% summarise(mean_contact_duration=mean(contact_duration))
t_save_dur1%>%group_by(regular_index)%>%summarise(mean(mean_contact_duration))
# in case, there are some outlier in the transient group, we need to clean it out
# summary(t_save_dur$mean_contact_duration[which(t_save$regular_index==0)])

t_save_encounter %>% summarise(mean_contact_duration_inday=sum(contact_duration)/52/60)%>% group_by(regular_index)%>%#group_by(regular_index)%>%#
  summarise(avg_contact_duration_inhour=mean(mean_contact_duration_inday),sd = sd(mean_contact_duration_inday))

t <- t_save_dur %>% group_by(regular_index)%>%#group_by(regular_index)%>%#
  summarise(avg_contact_duration_inhour=mean(mean_contact_duration_inday),sd = sd(mean_contact_duration_inday))
t

# t_save_wday <- temp %>% group_by(sender, receiver,contact_date,contact_wday,regular_index) %>% 
#   summarise(contact_duration = sum(contact_consecutive_duration)) %>% ungroup %>%
#   group_by(sender,receiver,contact_wday,regular_index) 
# 
# t_save_wday_freq <- t_save_wday%>%group_by(contact_date,contact_wday,regular_index)%>%summarise(uniq_pairs=n_distinct(paste0(sender,'-',receiver)))
# t_save_wday_freq %>% group_by(contact_wday,regular_index)%>%summarise(uniq_pairs_perD=mean(uniq_pairs))
# 
# new_t_save_wday_dur <- t_save_wday %>% group_by(sender,receiver,contact_wday,regular_index) %>% summarise(mean_contact_duration=sum(contact_duration)/8/60)
# 
# t_save_wday_dur <- t_save_wday%>% summarise(mean_contact_duration=sum(contact_duration)/8/60)
# # in case, there are some outlier in the transient group, we need to clean it out
# summary(t_save_wday_dur$mean_contact_duration[which(t_save_wday$regular_index==0)])
# 
# t <- t_save_wday_dur %>% group_by(contact_wday,regular_index)%>%#group_by(regular_index)%>%#
#   summarise(avg_contact_duration_inhour=mean(mean_contact_duration))
# t
# 
# rm(t_save_wday)

acacia = T
cassia = T
source('code/data_process.R') # to get survey
rm(URID,USID)
## Locations ------
### - combine with the survey data, then we will generate the elements by locations and dormitory
### remove empty age in the survey form
tt<- t_save[-which(t_save$sender%in%rm_sender),];# rm(t_save)
### combine the unit information of sender and receiver 
tt1<- left_join(tt,survey,by= c('sender'= 'tokenID'))
tt1<- left_join(tt1,survey,by= c('receiver'= 'tokenID'))
# allPairSumData1 <- as.data.frame(test)[complete.cases(test),]
# > length(unique(allPairSumData1$sender))
# [1] 1003
# > length(unique(allPairSumData1$receiver))
# [1] 1003
# rm(test)
tt1$unknown_location <- ifelse(is.na(tt1$br.y)|is.na(tt1$blk.y)|is.na(tt1$dorm.y)|is.na(tt1$room.y),'unknown','known')
# tk <- tt1 %>% group_by(regular_index,unknown_location) %>% #group_by(regular_contact,same_room) %>%
#   summarise(uniq_pairs=n_distinct(paste0(sender,'-',receiver)))#numPairs=n())
# tk_dur <- tt1 %>% group_by(regular_index,unknown_location) %>% #group_by(regular_contact,same_room) %>%
#   summarise(contact_duration=mean(mean_contact_duration))
tt1 <- as.data.frame(tt1)[complete.cases(tt1),]
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x)-n+1)
}
tt1$floor.x <- as.numeric(substrRight(tt1$room.x,3))
tt1$floor.y <- as.numeric(substrRight(tt1$room.y,3))

tt1$same_room <- ifelse(tt1$br.x==tt1$br.y & tt1$blk.x==tt1$blk.y,1,0)
tt1$same_floor <- ifelse(tt1$floor.x==tt1$floor.y & tt1$blk.x==tt1$blk.y,1,0)
tt1$same_blk <- ifelse(tt1$blk.x==tt1$blk.y,1,0)

tt1$same_floor_notroom <- ifelse(tt1$same_floor==1 & tt1$same_room==0,"sameFloorNotRoom",
                               ifelse(tt1$same_floor==1 & tt1$same_room==1,"sameFloorSameRoom",
                                      'nSameFloorNSameRoom'))
tt1$same_blk_notfloor <- ifelse(tt1$same_blk==1 & tt1$same_floor==0,"sameBlkNotFloor",
                                           ifelse(tt1$same_blk==1 & tt1$same_floor==1,"sameBlkSameFloor",
                                                  'nSameBlkNSameFloor'))
tt1$same_dorm <- ifelse(tt1$dorm.x==tt1$dorm.y,1,0)
tt1$same_dorm_notblk <- ifelse(tt1$same_dorm==1 & tt1$same_blk==0,"sameDormNotBlk",
                                           ifelse(tt1$same_dorm==1 & tt1$same_blk==1,"sameDormSameBlk",
                                                  'nSameDormNSameBlk'))
# rm_colnames <- c(7:22,33:34)
# tt1 <- tt1[,-rm_colnames]
tt11 <- unique(tt1[,c(1,2,4,25:31)])


tr <- tt1 %>% group_by(sender,contact_date,regular_index) %>%  #group_by(regular_contact,same_room) %>% 
summarise(uniq_contacts=n_distinct(receiver))
tr %>% group_by(regular_index)%>%summarise(sum(uniq_contacts)/1003/53)

# tr0 <- tr%>%filter(regular_index==0)
# tr1 <- tr%>%filter(regular_index==1)

# test_tr <- tr %>% group_by(sender,regular_index)%>%summarise(sum(uniq_contacts)/53)
# mean(test_tr$`sum(uniq_contacts)/53`[which(test_tr$regular_index==1)])
# mean(test_tr$`sum(uniq_contacts)/53`[which(test_tr$regular_index==0)])

# tr%>% group_by(regular_index)%>%summarise(mean(uniq_contacts)) #)numPairs=n()) %>% group_by(regular_index)%>%summarise(mean(uniq_contacts))
tall_dur <- tt1 %>% group_by(sender,receiver,regular_index)%>%summarise(contact_duration=sum(contact_duration)/53/60)#summarise(contact_duration=mean(mean_contact_duration))
tall_dur%>%group_by(regular_index)%>%summarise(mean(contact_duration))
# the one we used in the manuscript!!


# reply review,  to use boostrap to get the ratio and 95%CI ---------------
tData <- tall_dur%>%filter(regular_index==0)
rData <- tall_dur%>%filter(regular_index==1)
ratio = NULL
for(i in 1:5000){
  set.seed(100+i)
  # temp <- sample(tData$contact_duration,length(rData$contact_duration),replace = F)/rData$contact_duration
  temp <- tData$contact_duration/boot::boot(rData$contact_duration)
  ratio <- c(ratio,temp)
}


# reply reviewer, violin plot ---------------------------------------------

# please check in the violin_plot_reply_review.R





tsr <- tt1 %>% group_by(sender,contact_date,regular_index,same_room) %>% #group_by(regular_contact,same_room) %>% 
  summarise(uniq_pairs=n_distinct(paste0(sender,'-',receiver)))#)numPairs=n())
# tsr %>% group_by(regular_index,same_room)%>%summarise(mean(uniq_pairs))

# count unique contacts
tsr <- tt1 %>% group_by(sender,contact_date,regular_index,same_room) %>% #group_by(regular_contact,same_room) %>% 
  summarise(uniq_contacts=n_distinct(receiver))#)numPairs=n())
tsr %>% group_by(regular_index,same_room) %>% summarise(sum(uniq_contacts)/1003/53)

# test_tsr <- tsr %>% group_by(sender,regular_index,same_room)%>%summarise(sum(uniq_contacts)/53)
# roomNumR <- tsr %>% filter(regular_index==1) %>% #filter(regular_contact==1) %>% 
#   select(same_room,numPairs)
# roomNumT <- tsr %>% filter(regular_index==0) %>%#filter(regular_contact==0) %>% 
#   select(same_room,numPairs)
# rm(t)

tsr_dur <- tt1 %>% group_by(sender,receiver,same_room,regular_index)%>%summarise(contact_duration=sum(contact_duration)/53/60)#summarise(contact_duration=mean(mean_contact_duration))

# tsr_dur <- tt1 %>% group_by(contact_date,regular_index,same_room) %>% #group_by(regular_contact,same_room) %>% 
#   summarise(contact_duration=sum(contact_duration)/53/60)
tsr_dur%>%group_by(regular_index,same_room)%>%summarise(mean(contact_duration))


# count unique contacts
tsf <- tt1 %>% group_by(sender,contact_date,regular_index,same_floor_notroom) %>% #group_by(regular_contact,same_room) %>% 
  summarise(uniq_contacts=n_distinct(receiver))#)numPairs=n())
# tsf %>% group_by(regular_index,same_floor_notroom)%>%summarise(mean(uniq_contacts))
tsf %>% group_by(regular_index,same_floor_notroom)%>%summarise(sum(uniq_contacts)/1003/53)

# test_tsf <- tsf %>% group_by(sender,regular_index,same_floor_notroom)%>%summarise(sum(uniq_contacts)/53)
# mean(test_tsf$`sum(uniq_contacts)/53`[which(test_tsf$regular_index==1,test_tsf$same_floor_notroom==1)])
# mean(test_tsf$`sum(uniq_contacts)/53`[which(test_tsf$regular_index==0,test_tsf$same_floor_notroom==1)])
# 


tsf_dur <- tt1 %>% group_by(sender,receiver,same_floor_notroom,regular_index)%>%summarise(contact_duration=sum(contact_duration)/53/60)#summarise(contact_duration=mean(mean_contact_duration))

tsf_dur%>%group_by(regular_index,same_floor_notroom)%>%summarise(mean(contact_duration))
# tsf_dur%>%group_by(regular_index,same_floor_notroom)%>%summarise(sum(contact_duration))



# tsb <- tt1 %>% group_by(contact_date,regular_index,same_blk_notroom) %>% #group_by(regular_contact,same_blk_notroom) %>% 
#   summarise(uniq_pairs=n_distinct(paste0(sender,'-',receiver)))#numPairs=n())
# tsb %>% group_by(regular_index,same_blk_notroom)%>%summarise(mean(uniq_pairs))
# count unique contacts
tsb <- tt1 %>% group_by(sender,contact_date,regular_index,same_blk_notfloor) %>% #group_by(regular_contact,same_room) %>% 
  summarise(uniq_contacts=n_distinct(receiver))#)numPairs=n())
# tsb %>% group_by(regular_index,same_blk_notfloor)%>%summarise(mean(uniq_contacts))
tsb %>% group_by(regular_index,same_blk_notfloor)%>%summarise(sum(uniq_contacts)/1003/53)

# test_tsb <- tsb %>% group_by(sender,regular_index,same_blk_notfloor)%>%summarise(sum(uniq_contacts)/53)
# mean(test_tsb$`sum(uniq_contacts)/53`[which(test_tsb$regular_index==1,test_tsb$same_blk_notfloor==1)])
# mean(test_tsb$`sum(uniq_contacts)/53`[which(test_tsb$regular_index==0,test_tsb$same_blk_notfloor==1)])

tsb_dur <- tt1 %>% group_by(sender,receiver,same_blk_notfloor,regular_index)%>%summarise(contact_duration=sum(contact_duration)/53/60)#summarise(contact_duration=mean(mean_contact_duration))

# tsb_dur%>%group_by(regular_index,same_blk_notfloor)%>%summarise(mean(contact_duration))

# blkNumR <- t %>% filter(regular_index==1) %>% #filter(regular_contact==1) %>% 
#   select(same_blk_notroom,numPairs)
# blkNumT <- t %>%  filter(regular_index==0) %>% #filter(regular_contact==0) %>% 
#   select(same_blk_notroom,numPairs)
# rm(t)

# tsd <- tt1 %>% group_by(contact_date,regular_index,same_dorm_notblk) %>% #group_by(regular_contact,same_dorm_notblk) %>% 
#   summarise(uniq_pairs=n_distinct(paste0(sender,'-',receiver)))#numPairs=n())
# tsd %>% group_by(regular_index,same_dorm_notblk)%>%summarise(mean(uniq_pairs))
# count unique contacts
tsd <- tt1 %>% group_by(sender,contact_date,regular_index,same_dorm_notblk) %>% #group_by(regular_contact,same_room) %>% 
  summarise(uniq_contacts=n_distinct(receiver))#)numPairs=n())
# tsd %>% group_by(regular_index,same_dorm_notblk)%>%summarise(mean(uniq_contacts))
tsd %>% group_by(regular_index,same_dorm_notblk)%>%summarise(sum(uniq_contacts)/1003/53)

tsd_dur <- tt1 %>% group_by(sender,receiver,same_dorm_notblk,regular_index)%>%summarise(contact_duration=sum(contact_duration)/53/60)#summarise(contact_duration=mean(mean_contact_duration))
tsd_dur%>%group_by(regular_index,same_dorm_notblk)%>%summarise(mean(contact_duration))
# tsd_dur%>%group_by(regular_index,same_dorm_notblk)%>%summarise(sum(contact_duration)/1003/53)



# tsd_dur <- tt1 %>% group_by(regular_index,same_dorm_notblk) %>% #group_by(regular_contact,same_dorm_notblk) %>% 
#   summarise(contact_duration=mean(mean_contact_duration))
# dormNumR <- t %>% filter(regular_index==1) %>% #filter(regular_contact==1) %>% 
#   select(same_dorm_notblk,numPairs)
# dormNumT <- t %>% filter(regular_index==0) %>%#filter(regular_contact==0) %>%  
#   select(same_dorm_notblk,numPairs)
# rm(t)



# further check 0821
ttt<- temp[-which(temp$sender%in%rm_sender),];# rm(t_save)
### combine the unit information of sender and receiver 
ttt1<- left_join(ttt,survey,by= c('sender'= 'tokenID'))
ttt1<- left_join(ttt1,survey,by= c('receiver'= 'tokenID'))
# allPairSumData1 <- as.data.frame(test)[complete.cases(test),]
# > length(unique(allPairSumData1$sender))
# [1] 1003
# > length(unique(allPairSumData1$receiver))
# [1] 1003
# rm(test)
ttt1$unknown_location <- ifelse(is.na(ttt1$br.y)|is.na(ttt1$blk.y)|is.na(ttt1$dorm.y)|is.na(ttt1$room.y),'unknown','known')
# tk <- tt1 %>% group_by(regular_index,unknown_location) %>% #group_by(regular_contact,same_room) %>%
#   summarise(uniq_pairs=n_distinct(paste0(sender,'-',receiver)))#numPairs=n())
# tk_dur <- tt1 %>% group_by(regular_index,unknown_location) %>% #group_by(regular_contact,same_room) %>%
#   summarise(contact_duration=mean(mean_contact_duration))
ttt1 <- as.data.frame(ttt1)[complete.cases(ttt1),]
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x)-n+1)
}
ttt1$floor.x <- as.numeric(substrRight(ttt1$room.x,3))
ttt1$floor.y <- as.numeric(substrRight(ttt1$room.y,3))

ttt1$same_room <- ifelse(ttt1$br.x==ttt1$br.y & ttt1$blk.x==ttt1$blk.y,1,0)
# ttt1$mean_contact_consecutive_duration <- ttt1$contact_consecutive_duration/53
# tsr_dur1 <- ttt1 %>% group_by(sender,receiver,same_room,regular_index)%>%summarise(contact_duration=mean(mean_contact_consecutive_duration))
# 14mins

temp1 <- ttt1 %>% group_by(sender,receiver,same_room,regular_index,contact_date)%>%
  summarise(mean_contact_consecutive_duration=mean(contact_consecutive_duration))
tsr_dur1 <- temp1 %>% group_by(sender,receiver,same_room,regular_index)%>%summarise(contact_duration=sum(mean_contact_consecutive_duration)/52/60)
mean(tsr_dur1$contact_duration[which(tsr_dur1$same_room==1&tsr_dur1$regular_index==0)]) #9mins
mean(tsr_dur1$contact_duration[which(tsr_dur1$same_room==1&tsr_dur1$regular_index==1)]) #1hour7mins
