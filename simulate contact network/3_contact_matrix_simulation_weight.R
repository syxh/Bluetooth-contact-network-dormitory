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


# we need to get the weight of the same room, same floor level, same blk and same dorm in A and C of Rc, Tc---------
#same room------------
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
temp_table$mean_Rcweight_room <- temp_table$mean_Tcweight_room<- rep(0,nrow(temp_table))
temp$br.x[which(temp$room.x=='4027')] <- "15-427"
temp$br.y[which(temp$room.y=='4027')] <- "15-427"

unique_br <- unique(temp$br.x)
for(i in 1:length(unique_br))
{
  sub_temp_room <- temp%>%filter(br.x==unique_br[i],br.y==unique_br[i])%>%distinct(sender,receiver,contact_date,regular_index)
  
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
      
      # print(length(rmid))
      # print(dim(sub_temp_room))
      if(length(rmid)>0){
        sub_temp_room<-sub_temp_room[-which(sub_temp_room$id%in%rmid),]
      }}
      # print(unique(sub_temp_room[,c('sender','receiver')]))
    }
    sub_temp_room$sr <- paste0(sub_temp_room$sender,'-',sub_temp_room$receiver)
    tt <- as.data.frame(table(sub_temp_room$sr,sub_temp_room$regular_index))
    colnames(tt)<-c('sender-receiver','regular_index','days')
    temp_table$mean_Rcweight_room[which(temp_table$br==unique_br[i])] <- mean(tt$days[which(tt$regular_index==1)]/53)
    temp_table$mean_Tcweight_room[which(temp_table$br==unique_br[i])] <- mean(tt$days[which(tt$regular_index==0)]/53)
    
  }else{
    temp_table$mean_Rcweight_room[which(temp_table$br==unique_br[i])]<-0
    temp_table$mean_Tcweight_room[which(temp_table$br==unique_br[i])]<-0
  }
}

write.csv(temp_table,'output/S_matrix_room_weight_RandT.csv',row.names = F)

# same floor different room-----------
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
temp_table1$mean_Rcweight_level = temp_table1$mean_Tcweight_level= rep(0,nrow(temp_table1))

unique_bl <- unique(temp$bl.x)
for(i in 1:length(unique_bl))
{
  sub_temp_level <- temp%>%filter(bl.x==unique_bl[i],bl.y==unique_bl[i])%>%distinct(sender,receiver,contact_date,regular_index,br.x,br.y)
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
    sub_temp_level$sr <- paste0(sub_temp_level$sender,'-',sub_temp_level$receiver)
    tt <- as.data.frame(table(sub_temp_level$sr,sub_temp_level$regular_index))
    colnames(tt)<-c('sender-receiver','regular_index','days')
    temp_table1$mean_Rcweight_level[which(temp_table1$bl==unique_bl[i])] <- mean(tt$days[which(tt$regular_index==1)]/53)
    temp_table1$mean_Tcweight_level[which(temp_table1$bl==unique_bl[i])] <- mean(tt$days[which(tt$regular_index==0)]/53)
    
  }else{
    temp_table1$mean_Rcweight_level[which(temp_table1$bl==unique_bl[i])]<-0
    temp_table1$mean_Tcweight_level[which(temp_table1$bl==unique_bl[i])]<-0
  }
}
write.csv(temp_table1,'output/S_matrix_level_weight_RandT.csv',row.names = F)

# same blk different floors---------
temp_table2 <-unique(data.frame(dorm = t0$dorm,
                                blk = t0$blk))
temp_table2$mean_Rcweight_blk <- temp_table2$mean_Tcweight_blk<- rep(0,nrow(temp_table2))

unique_blk <- unique(temp$blk.x)
for(i in 1:length(unique_blk))
{
  sub_temp_blk <- temp%>%filter(blk.x==unique_blk[i],blk.y==unique_blk[i])%>%distinct(sender,receiver,contact_date,regular_index,bl.x,bl.y)
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
    sub_temp_blk$sr <- paste0(sub_temp_blk$sender,'-',sub_temp_blk$receiver)
    tt <- as.data.frame(table(sub_temp_blk$sr,sub_temp_blk$regular_index))
    colnames(tt)<-c('sender-receiver','regular_index','days')
    temp_table2$mean_Rcweight_blk[which(temp_table2$blk==unique_blk[i])] <- mean(tt$days[which(tt$regular_index==1)]/53)
    temp_table2$mean_Tcweight_blk[which(temp_table2$blk==unique_blk[i])] <- mean(tt$days[which(tt$regular_index==0)]/53)
    
  }else{
    temp_table2$mean_Rcweight_blk[which(temp_table2$blk==unique_blk[i])]<-0
    temp_table2$mean_Tcweight_blk[which(temp_table2$blk==unique_blk[i])]<-0
  }
}
write.csv(temp_table2,'output/S_matrix_blk_weight_RandT.csv',row.names = F)


# same dorm different blks------------
temp_table3 <-unique(data.frame(dorm = t0$dorm))
temp_table3$mean_Rcweight_dorm <- temp_table3$mean_Tcweight_dorm<- rep(0,nrow(temp_table3))
unique_dorm <- unique(temp$dorm.x)
for(i in 1:2)
{
  sub_temp_dorm <- temp%>%filter(dorm.x==unique_dorm[i],dorm.y==unique_dorm[i])%>%distinct(sender,receiver,contact_date,regular_index,blk.x,blk.y)
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
    sub_temp_dorm$sr <- paste0(sub_temp_dorm$sender,'-',sub_temp_dorm$receiver)
    tt <- as.data.frame(table(sub_temp_dorm$sr,sub_temp_dorm$regular_index))
    colnames(tt)<-c('sender-receiver','regular_index','days')
    temp_table3$mean_Rcweight_dorm[which(temp_table3$dorm==unique_dorm[i])] <- mean(tt$days[which(tt$regular_index==1)]/53)
    temp_table3$mean_Tcweight_dorm[which(temp_table3$dorm==unique_dorm[i])] <- mean(tt$days[which(tt$regular_index==0)]/53)
    
  }else{
    temp_table3$mean_Rcweight_dorm[which(temp_table3$dorm==unique_dorm[i])]<-0
    temp_table3$mean_Tcweight_dorm[which(temp_table3$dorm==unique_dorm[i])]<-0
  }
}

write.csv(temp_table3,'output/S_matrix_dorm_weight_RandT.csv',row.names = F)

# between dormitories--------
sub_temp_between_dorm <- temp%>%distinct(sender,receiver,regular_index,contact_date,dorm.x,dorm.y)
dim(sub_temp_between_dorm)
head(sub_temp_between_dorm)
sub_temp_between_dorm <- sub_temp_between_dorm[-which(sub_temp_between_dorm$dorm.x==sub_temp_between_dorm$dorm.y),]
dim(sub_temp_between_dorm)

sub_temp_between_dorm$sr <- paste0(sub_temp_between_dorm$sender,'-',sub_temp_between_dorm$receiver)
tt <- as.data.frame(table(sub_temp_between_dorm$sr,sub_temp_between_dorm$regular_index))
colnames(tt)<-c('sender-receiver','regular_index','days')
mean(tt$days[which(tt$regular_index==1)]/53)
# [1] 0.002497269
median(tt$days[which(tt$regular_index==1)]/53)
# 0
mean(tt$days[which(tt$regular_index==0)]/53)
# 0.106248
median(tt$days[which(tt$regular_index==0)]/53)
# 0.03773585
rm(list = ls())


# create the 8256*8256 weight matrix of contact
ARRc <- 0.45/0.45
ARTc <- 0.4/0.45
AFRc <- 0.04/0.45
AFTc <- 0.13/0.45
ABRc <- 0.04/0.45
ABTc <- 0.11/0.45
ADRc <- 0.0003/0.45
ADTc <- 0.05/0.45

CRRc <- 0.68/0.68
CRTc <- 0.36/0.68
CFRc <- 0.05/0.68
CFTc <- 0.25/0.68
CBRc <- 0.03/0.68
CBTc <- 0.1/0.68
CDRc <- 0.003/0.68
CDTc <- 0.06/0.68

betweenRc <- 0.002/0.68
betweenTc <- 0.04/0.68

load('C:/Users/ephsy/Documents/In Hand Project (PhD)/Dorm/spread_between_dorm/working/simulate contact matrix/simulate_contact.Rdata')
sim_regular_contact <- 1*as.matrix(sim_contact==1)
sim_transient_contact <- 1*as.matrix(sim_contact==2)
rm(sim_contact)
n <- 4800+3456
weight_sim_regular_contact <- matrix(0,n,n)
weight_sim_transient_contact <- matrix(0,n,n)
#dorm A
di=1
for(bi in 1:6)
{
  for(fi in 1:5)
  {
    for(ri in 1:8)
    {
      for(i in (1+(bi*fi-1)*20*8 +(ri-1)*20):(20+(bi*fi-1)*20*8 +(ri-1)*20))
        for(j in i:(20+(bi*fi-1)*20*8 +(ri-1)*20))
        {
          weight_sim_regular_contact[i,j] <- ARRc*sim_regular_contact[i,j]
          weight_sim_transient_contact[i,j] <- ARTc*sim_transient_contact[i,j]
        }
      if(ri<8){
        for(restr_row in (1+(bi*fi-1)*160):(20+(bi*fi-1)*20*8 +(ri-1)*20))#160:n_bedder*n_room on each floor
        {
          for(restr_col in min((141+(bi*fi-1)*20*8),(20+1+(ri-1)*20+(bi*fi-1)*20*8)):min((160+(bi*fi-1)*20*8),(ri*20+20+(bi*fi-1)*20*8)))
          {
            weight_sim_regular_contact[restr_row,restr_col] <- AFRc*sim_regular_contact[restr_row,restr_col]
            weight_sim_transient_contact[restr_row,restr_col] <- AFTc*sim_transient_contact[restr_row,restr_col]
          }
        }
      }
    }
    # To fill in the same blk but different room cells
    if(fi<5)
    {
      for(restf_row in (1+800*(bi-1)):(160+160*(fi-1)+800*(bi-1)))
      {
        for(restf_col in (161+160*(fi-1)+800*(bi-1)):(320+160*(fi-1)+800*(bi-1)))
        {
          weight_sim_regular_contact[restf_row,restf_col] <- ABRc*sim_regular_contact[restf_row,restf_col]
          weight_sim_transient_contact[restf_row,restf_col] <- ABTc*sim_transient_contact[restf_row,restf_col]
        }
      }
    }
  }
  # To fill in the same dorm but different blk cells
  if(bi<6)
  {
    for(restb_row in (1+4800*(di-1)):(800+800*(bi-1)+4800*(di-1)))
    {
      for(restb_col in (801+800*(bi-1)+4800*(di-1)):(1600+800*(bi-1)+4800*(di-1)))
      {
        weight_sim_regular_contact[restb_row,restb_col] <- ADRc*sim_regular_contact[restb_row,restb_col]
        weight_sim_transient_contact[restb_row,restb_col] <- ADTc*sim_transient_contact[restb_row,restb_col]
      }
    }
  }
}
table(weight_sim_regular_contact[1:20,1:20])
table(weight_sim_transient_contact[1:20,1:20])
table(weight_sim_regular_contact[1:160,1:160])
table(weight_sim_transient_contact[1:160,1:160])
table(weight_sim_regular_contact[1:800,1:800])
table(weight_sim_transient_contact[1:800,1:800])
table(weight_sim_regular_contact[1:4800,1:4800])
table(weight_sim_transient_contact[1:4800,1:4800])


di=2
n_bedder <- 16
n_room <- 12
n_level = 6
n_blk = 3
for(bi in 1:3)# 6 blk in acacia
{ 
  for(fi in 1:6) #each blk 5 floor-levels
  {
    for(ri in 1:n_room) # each floor 8 rooms
    {
      
      for(i in (4800+1+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder):(4800+n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))
      {
        for(j in i:(4800+n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))#(i+(fi-1)*8*20):(20+(fi-1)*8*20+(ri-1)*20))
        {
          if(i!=j)
          {
            weight_sim_regular_contact[i,j] <- CRRc*sim_regular_contact[i,j]
            weight_sim_transient_contact[i,j] <- CRTc*sim_transient_contact[i,j]
          }
        }
      }
      
      # To fill in the same floor but different room contacts
      if(ri<12){
        for(restr_row in (4800+1+(bi*fi-1)*n_bedder*n_room):(4800+n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))#160:n_bedder*n_room on each floor
        {
          for(restr_col in min((4800+n_bedder*n_room-(n_bedder)+1+(bi*fi-1)*n_bedder*n_room),(4800+n_bedder+1+(ri-1)*n_bedder+(bi*fi-1)*n_bedder*n_room )):min((4800+n_bedder*n_room+(bi*fi-1)*n_bedder*n_room ),(4800+ri*n_bedder+n_bedder+(bi*fi-1)*n_bedder*n_room)))
          {                    
            weight_sim_regular_contact[restr_row,restr_col] <- CFRc*sim_regular_contact[restr_row,restr_col]
            weight_sim_transient_contact[restr_row,restr_col] <- CFTc*sim_transient_contact[restr_row,restr_col]
          }
        }
      }
    }
    # To fill in the same blk but different room contacts
    if(fi<6)
    {
      for(restf_row in (4800+1+n_bedder*n_room*n_level*(bi-1)):(4800+n_bedder*n_room+n_bedder*n_room*(fi-1)+1152*(bi-1)))
      {
        for(restf_col in (4800+n_bedder*n_room+1+n_bedder*n_room*(fi-1)+n_bedder*n_room*n_level*(bi-1)):(4800+2*n_bedder*n_room+n_bedder*n_room*(fi-1)+n_bedder*n_room*n_level*(bi-1)))
        {
          weight_sim_regular_contact[restf_row,restf_col] <- CBRc*sim_regular_contact[restf_row,restf_col]
          weight_sim_transient_contact[restf_row,restf_col] <- CBTc*sim_transient_contact[restf_row,restf_col]
        }
      }
    }
  }
  # To fill in the same dorm but different blk contacts
  if(bi<3)
  {
    for(restb_row in (1+4800*(di-1)):(n_bedder*n_room*n_level+n_bedder*n_room*n_level*(bi-1)+4800*(di-1)))
    {
      for(restb_col in (n_bedder*n_room*n_level+1+n_bedder*n_room*n_level*(bi-1)+4800*(di-1)):(2*n_bedder*n_room*n_level+n_bedder*n_room*n_level*(bi-1)+4800*(di-1)))
      {
        weight_sim_regular_contact[restb_row,restb_col] <- CDRc*sim_regular_contact[restb_row,restb_col]
        weight_sim_transient_contact[restb_row,restb_col] <- CDTc*sim_transient_contact[restb_row,restb_col]
      }
    }
  }
}
table(weight_sim_regular_contact[4801:(4800+16),4801:(4800+16)])
table(weight_sim_transient_contact[4801:(4800+16),4801:(4800+16)])

table(weight_sim_regular_contact[4801:(4800+192),4801:(4800+192)])
table(weight_sim_transient_contact[4801:(4800+192),4801:(4800+192)])

table(weight_sim_regular_contact[4801:(4800+1152),4801:(4800+1152)])
table(weight_sim_transient_contact[4801:(4800+1152),4801:(4800+1152)])

table(weight_sim_regular_contact[4801:n,4801:n])
table(weight_sim_transient_contact[4801:n,4801:n])

#betwenn dormitories
for(rest_row in 1:4800)
{
  for(rest_col in 4801:n)
  {
    weight_sim_regular_contact[rest_row,rest_col] <- betweenRc*sim_regular_contact[rest_row,rest_col]
    weight_sim_transient_contact[rest_row,rest_col] <- betweenTc*sim_transient_contact[rest_row,rest_col]
    
  }
}
write.csv(weight_sim_regular_contact,file='working/simulate contact matrix/simulate_regular_contact_weight.csv',row.names = F)
write.csv(weight_sim_transient_contact,file='working/simulate contact matrix/simulate_transient_contact_weight.csv',row.names = F)
save(weight_sim_regular_contact,file='working/simulate contact matrix/simulate_regular_contact_weight.Rdata')
save(weight_sim_transient_contact,file='working/simulate contact matrix/simulate_transient_contact_weight.Rdata')



