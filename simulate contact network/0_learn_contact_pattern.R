library(dplyr)
# Step 1: restrict to the completed survey data, get the contact matrix, and the proportion of Rc and Tc in different dormitory setting levels
load(paste0('working/combine/','rawPairs.Rdata')) # the cleaned data
load('working/map_regular_index.rdata')
temp <- left_join(allPairData,map_regular_index,by=c('sender','receiver'))
rm(allPairData,map_regular_index)
acacia = T
cassia = T
source('code/data_process.R') # to get survey
rm(URID,USID,filelist, foldername,keep_file,notefile,possiblefiles,survey_id_keep,rm_sender)
# temp0 <- temp
temp <- left_join(temp,survey,by=c('sender'='tokenID'))
head(temp)
temp <- as.data.frame(temp)[complete.cases(as.data.frame(temp)),]
temp <- left_join(temp,survey,by=c('receiver'='tokenID'))
temp <- as.data.frame(temp)[complete.cases(as.data.frame(temp)),]
head(temp)

# restrict to the unique contact pairs per day
temp <- temp[,-c(4,5,6,8)]
temp <- unique(temp)
temp$level.x <- substr(temp$room.x,1,nchar(temp$room.x)-2)
temp$level.x[which(temp$level.x==40)] <- 4
temp$level.y <- substr(temp$room.y,1,nchar(temp$room.y)-2)
temp$level.y[which(temp$level.y==40)] <- 4
temp <- temp[order(temp$dorm.x, temp$dorm.y,temp$br.x,temp$br.y,temp$level.x,temp$level.y,temp$blk.x,temp$blk.y),]
# Regular contact matrix ( adjcent matrix)
head(temp)
# rm(temp0)
unique_sender <- unique(temp$sender)
save(unique_sender,file='working/real contact matrix/contact_matrix_sender_index.Rdata')
unique_receiver <- unique_sender#unique(temp$receiver)
regular_mat <- matrix(0,length(unique_sender),length(unique_sender))
for(si in 1:length(unique_sender))
{
  sub_temp <-  temp %>% filter(sender==unique_sender[si])
  contact_index <- unique(sub_temp$receiver[which(sub_temp$regular_index==1)])
  col_id <- which(unique_receiver%in%contact_index)
  regular_mat[si,col_id] <- 1
}
save(regular_mat,file='working/real contact matrix/regular_contact_matrix.Rdata')
# Transient contact matrix

transient_mat <- matrix(0,length(unique_sender),length(unique_sender))
for(si in 1:length(unique_sender))
{
  sub_temp <-  temp %>% filter(sender==unique_sender[si])
  contact_index <- unique(sub_temp$receiver[which(sub_temp$regular_index==0)])
  col_id <- which(unique_receiver%in%contact_index)
  transient_mat[si,col_id] <- 1
}
save(transient_mat,file='working/real contact matrix/transient_contact_matrix.Rdata')
# running time: 1min in total

## Next, we will calculate the proportion of Rc and Tc in different domitory layer(including: room, level, blk, dorm)
## note that, should be consistent with the index of unique_sender. Notes: unique_sender, which is ordered by dorm, and unit number
# we found that, there is one room typo, should by 427 rather than 4027
survey$room[which(survey$room=='4027')] <- '427'
survey$br[which(survey$br=='15-4027')] <- '15-427'
t0 <- unique(survey[,c("room","blk","dorm")])
room <- t0$room
level <- substr(room,1,nchar(room)-2)
temp_table <-data.frame(dorm = t0$dorm,
                        blk = t0$blk,
                        level = level,
                        room = room)

temp_table<- temp_table[order(temp_table$dorm,temp_table$room,temp_table$level,temp_table$blk),]
temp_table$br <- paste0(temp_table$blk,'-',temp_table$room)
temp_table$number_regular_contact_room <- temp_table$number_transient_contact_room<- rep(0,nrow(temp_table))
temp$br.x[which(temp$room.x=='4027')] <- "15-427"
temp$br.y[which(temp$room.y=='4027')] <- "15-427"

unique_br <- unique(temp$br.x)
# in the following, we count number of unique regular pairs in the same room unit
for(i in 1:length(unique_br))
{
  sub_temp_room <- temp%>%filter(br.x==unique_br[i],br.y==unique_br[i])%>%distinct(sender,receiver,regular_index)

  if(nrow(sub_temp_room)>0){
    sub_temp_room$id <- c(1:nrow(sub_temp_room))
    t <- unique(sub_temp_room$receiver)
    for(j in 1:length(unique(t)))
    {
      if(any(sub_temp_room$sender%in%t[j]))
      {
        tt <- sub_temp_room[which(sub_temp_room$sender==t[j]),]
        rmid1 <- tt$id[which(tt$receiver %in% sub_temp_room$sender)]
        rmid <- rmid1[which(tt$receiver[which(tt$id%in%rmid1)]%in%sub_temp_room[which(sub_temp_room$receiver%in%tt$sender),]$sender)]
     
      if(length(rmid)>0){
        sub_temp_room<-sub_temp_room[-which(sub_temp_room$id%in%rmid),]
      } }
      # rm(rmid)
    }
 
    temp_table$number_regular_contact_room[which(temp_table$br==unique_br[i])] <- length(which(sub_temp_room$regular_index==1))
    temp_table$number_transient_contact_room[which(temp_table$br==unique_br[i])] <- length(which(sub_temp_room$regular_index==0))
    
  }else{
    temp_table$number_regular_contact_room[which(temp_table$br==unique_br[i])]<-0
    temp_table$number_transient_contact_room[which(temp_table$br==unique_br[i])]<-0
  }
}

str(temp_table)
# next step, we need to get the proportion of the Rc/Tc at different dormitory settings
# sub-step1, get the number of ppl in one room(ave) in our dataset, as well as floor-level, blk, dorm
head(survey)
numbr <- as.data.frame(table(survey$br))
mean(numbr$Freq)
colnames(numbr)[1] <- "br"
temp_table <- left_join(temp_table,numbr,by='br')
temp_table$total_pairs <- choose(temp_table$Freq,2)
temp_table$prop_Tc <- temp_table$number_transient_contact_room/choose(temp_table$Freq,2)
temp_table$prop_Rc <- temp_table$number_regular_contact_room/choose(temp_table$Freq,2)

write.csv(temp_table,'output/S_matrix_room_count_RandT_new.csv',row.names = F)
mean(temp_table$prop_Rc[complete.cases(temp_table$prop_Rc)])
# [1] 0.5472302  = 0.6
mean(temp_table$prop_Tc[complete.cases(temp_table$prop_Tc)])
# [1] 0.4260984 = 0.4

# the separation of the mean of A an dC were calculated from excel directly.

# in the following, we count number of unique regular pairs at the same floor-level (exclude same rooms)
survey$room[which(survey$room=='4027')] <- '427'
t0 <- unique(survey[,c("room","blk","dorm")])
room <- t0$room
level <- substr(room,1,nchar(room)-2)
temp_table1 <-unique(data.frame(dorm = t0$dorm,
                        blk = t0$blk,
                        level = level))
temp_table1$bl <- paste0(temp_table1$blk,'-',temp_table1$level)
temp$bl.x <- paste0(temp$blk.x,'-',temp$level.x)
temp$bl.y <- paste0(temp$blk.y,'-',temp$level.y)
temp_table1$number_regular_contact_level <- temp_table1$number_transient_contact_level<- rep(0,nrow(temp_table1))

unique_bl <- unique(temp$bl.x)
for(i in 1:length(unique_bl))
{
  sub_temp_level <- temp%>%filter(bl.x==unique_bl[i],bl.y==unique_bl[i])%>%distinct(sender,receiver,regular_index,br.x,br.y)
  sub_temp_level <- sub_temp_level[-which(sub_temp_level$br.x==sub_temp_level$br.y),]
  if(nrow(sub_temp_level)>0){
    sub_temp_level$id <- c(1:nrow(sub_temp_level))
    t <- unique(sub_temp_level$receiver)
    for(j in 1:length(unique(t)))
    {
      if(any(sub_temp_level$sender%in%t[j]))
      {
        tt <- sub_temp_level[which(sub_temp_level$sender==t[j]),]
        rmid1 <- tt$id[which(tt$receiver %in% sub_temp_level$sender)]
        rmid <- rmid1[which(tt$receiver[which(tt$id%in%rmid1)]%in%sub_temp_level[which(sub_temp_level$receiver%in%tt$sender),]$sender)]
     
      if(length(rmid)>0){
        sub_temp_level<-sub_temp_level[-which(sub_temp_level$id%in%rmid),]
      } }
    }
    temp_table1$number_regular_contact_level[which(temp_table1$bl==unique_bl[i])] <- length(which(sub_temp_level$regular_index==1))
    temp_table1$number_transient_contact_level[which(temp_table1$bl==unique_bl[i])] <- length(which(sub_temp_level$regular_index==0))
    
  }else{
    temp_table1$number_regular_contact_level[which(temp_table1$bl==unique_bl[i])]<-0
    temp_table1$number_transient_contact_level[which(temp_table1$bl==unique_bl[i])]<-0
  }
}
head(survey)
survey$bl <- paste0(survey$blk,'-',substr(survey$room,1,nchar(survey$room)-2))
numbr <- as.data.frame(table(survey$bl))
mean(numbr$Freq)
colnames(numbr)[1] <- "bl"
temp_table1 <- left_join(temp_table1,numbr,by='bl')
temp_table1$prop_Tc <- temp_table1$number_transient_contact_level/choose(temp_table1$Freq,2)
temp_table1$prop_Rc <- temp_table1$number_regular_contact_level/choose(temp_table1$Freq,2)
temp_table1$total_pairs <- choose(temp_table1$Freq,2)
write.csv(temp_table1,'output/S_matrix_level_count_RandT_new1.csv',row.names = F)
mean(temp_table1$prop_Rc[complete.cases(temp_table1$prop_Rc)])
# [1] 0.1712463  = 0.2
mean(temp_table1$prop_Tc[complete.cases(temp_table1$prop_Tc)])
# [1] 0.6589215 = 0.7



# in the following, we count number of unique regular pairs in the same blk
temp_table2 <-unique(data.frame(dorm = t0$dorm,
                               blk = t0$blk))
temp_table2$number_regular_contact_blk <- temp_table2$number_transient_contact_blk<- rep(0,nrow(temp_table2))

unique_blk <- unique(temp$blk.x)
for(i in 1:length(unique_blk))
{
  sub_temp_blk <- temp%>%filter(blk.x==unique_blk[i],blk.y==unique_blk[i])%>%distinct(sender,receiver,regular_index,bl.x,bl.y)
  sub_temp_blk <- sub_temp_blk[-which(sub_temp_blk$bl.x==sub_temp_blk$bl.y),]
  if(nrow(sub_temp_blk)>0){
    sub_temp_blk$id <- c(1:nrow(sub_temp_blk))
    t <- unique(sub_temp_blk$receiver)
    for(j in 1:length(unique(t)))
    {
      if(any(sub_temp_blk$sender%in%t[j]))
      {
        tt <- sub_temp_blk[which(sub_temp_blk$sender==t[j]),]
        rmid1 <- tt$id[which(tt$receiver %in% sub_temp_blk$sender)]
        rmid <- rmid1[which(tt$receiver[which(tt$id%in%rmid1)]%in%sub_temp_blk[which(sub_temp_blk$receiver%in%tt$sender),]$sender)]
      
      if(length(rmid)>0){
        sub_temp_blk<-sub_temp_blk[-which(sub_temp_blk$id%in%rmid),]
      }}
    }
    temp_table2$number_regular_contact_blk[which(temp_table2$blk==unique_blk[i])] <- length(which(sub_temp_blk$regular_index==1))
    temp_table2$number_transient_contact_blk[which(temp_table2$blk==unique_blk[i])] <- length(which(sub_temp_blk$regular_index==0))
    
  }else{
    temp_table2$number_regular_contact_blk[which(temp_table2$blk==unique_blk[i])]<-0
    temp_table2$number_transient_contact_blk[which(temp_table2$blk==unique_blk[i])]<-0
  }
}
numbr <- as.data.frame(table(survey$blk))
mean(numbr$Freq)
colnames(numbr)[1] <- "blk"
numbr$blk <- as.numeric(as.character(numbr$blk))
temp_table2 <- left_join(temp_table2,numbr,by='blk')
temp_table2$prop_Tc <- temp_table2$number_transient_contact_blk/choose(temp_table2$Freq,2)
temp_table2$prop_Rc <- temp_table2$number_regular_contact_blk/choose(temp_table2$Freq,2)
temp_table2$total_pairs <- choose(temp_table2$Freq,2)
write.csv(temp_table2,'output/S_matrix_blk_count_RandT_new1.csv',row.names = F)
mean(temp_table2$prop_Rc[complete.cases(temp_table2$prop_Rc)])
# [1] 0.07237977 = 0.1
mean(temp_table2$prop_Tc[complete.cases(temp_table2$prop_Tc)])
# [1] 0.6354207 = 0.6

# in the following, we count number of unique regular pairs in the same blk
temp_table3 <-unique(data.frame(dorm = t0$dorm))
temp_table3$number_regular_contact_dorm <- temp_table3$number_transient_contact_dorm<- rep(0,nrow(temp_table3))
unique_dorm <- unique(temp$dorm.x)
for(i in 1:2)
{
  sub_temp_dorm <- temp%>%filter(dorm.x==unique_dorm[i],dorm.y==unique_dorm[i])%>%distinct(sender,receiver,regular_index,blk.x,blk.y)
  sub_temp_dorm <- sub_temp_dorm[-which(sub_temp_dorm$blk.x==sub_temp_dorm$blk.y),]
  if(nrow(sub_temp_dorm)>0){
    sub_temp_dorm$id <- c(1:nrow(sub_temp_dorm))
    t <- unique(sub_temp_dorm$receiver)
    for(j in 1:length(unique(t)))
    {
      if(any(sub_temp_dorm$sender%in%t[j]))
      {
        tt <- sub_temp_dorm[which(sub_temp_dorm$sender==t[j]),]
        rmid1 <- tt$id[which(tt$receiver %in% sub_temp_dorm$sender)]
        rmid <- rmid1[which(tt$receiver[which(tt$id%in%rmid1)]%in%sub_temp_dorm[which(sub_temp_dorm$receiver%in%tt$sender),]$sender)]
      
      if(length(rmid)>0){
        sub_temp_dorm<-sub_temp_dorm[-which(sub_temp_dorm$id%in%rmid),]
      }}
    }
    temp_table3$number_regular_contact_dorm[which(temp_table3$dorm==unique_dorm[i])] <- length(which(sub_temp_dorm$regular_index==1))
    temp_table3$number_transient_contact_dorm[which(temp_table3$dorm==unique_dorm[i])] <- length(which(sub_temp_dorm$regular_index==0))
    
  }else{
    temp_table3$number_regular_contact_dorm[which(temp_table3$dorm==unique_dorm[i])]<-0
    temp_table3$number_transient_contact_dorm[which(temp_table3$dorm==unique_dorm[i])]<-0
  }
}
numbr <- as.data.frame(table(survey$dorm))
mean(numbr$Freq)
colnames(numbr)[1] <- "dorm"
temp_table3 <- left_join(temp_table3,numbr,by='dorm')
temp_table3$prop_Tc <- temp_table3$number_transient_contact_dorm/choose(temp_table3$Freq,2)
temp_table3$prop_Rc <- temp_table3$number_regular_contact_dorm/choose(temp_table3$Freq,2)
temp_table3$total_pairs <- choose(temp_table3$Freq,2)
mean(temp_table3$prop_Rc[complete.cases(temp_table3$prop_Rc)])
# .01127288=0.01 
mean(temp_table3$prop_Tc[complete.cases(temp_table3$prop_Tc)])
# [1] 0.4564791 = 0.46
write.csv(temp_table3,'output/S_matrix_dorm_count_RandT_new.csv',row.names = F)




# between dorm
sub_temp_between_dorm <- temp%>%distinct(sender,receiver,regular_index,dorm.x,dorm.y)
head(sub_temp_between_dorm)
sub_temp_between_dorm <- sub_temp_between_dorm[-which(sub_temp_between_dorm$dorm.x==sub_temp_between_dorm$dorm.y),]
table(sub_temp_between_dorm$regular_index)
270/choose(1003,2)
70275/choose(1003,2)
