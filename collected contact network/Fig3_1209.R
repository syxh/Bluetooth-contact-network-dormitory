library(grid)
library(dplyr)
library(lubridate)
colors <- c("#ff7373","#6897bb") # regular (red) and transient (blue)
foldername <- 'working/combine/'
outputfoler <- 'output/combine/'
# source('code/functionsOfPlot.R')
# aim to plot the proportion of contacts and contact duration 

load(paste0(foldername,'PairDataContactType.Rdata'))
# proportion of the duration paneller -------
proportionDurationPaneller <- function(numGrp=5,data,yType,fontsize)
{
  if(yType=='pairs')
  {
    ylab <- "Proportion of contact events"
    # data <- data %>% group_by(sender,receiver,contact_date,contact_consecutive_duration,regular_contact) %>% summarise(mean_contact_consecutive_duration = mean(contact_consecutive_duration))
    data$durationCat <- ifelse(data$mean_contact_consecutive_duration<=7,1,
                               ifelse(data$mean_contact_consecutive_duration<=20 & data$mean_contact_consecutive_duration>7,2,
                                      ifelse(data$mean_contact_consecutive_duration<=60 & data$mean_contact_consecutive_duration>20,3,
                                             ifelse(data$mean_contact_consecutive_duration<=120 & data$mean_contact_consecutive_duration>60,4,5))))
    # data$durationCat <- ifelse(data$contact_consecutive_duration<=7,1,
    #                            ifelse(data$contact_consecutive_duration<=20 & data$contact_consecutive_duration>7,2,
    #                                   ifelse(data$contact_consecutive_duration<=60 & data$contact_consecutive_duration>20,3,
    #                                          ifelse(data$contact_consecutive_duration<=120 & data$contact_consecutive_duration>60,4,5))))
    tmp <- as.data.frame(prop.table(table(data$regular_index,data$durationCat),margin=2))
    names(tmp) <- c('regular_index','duration_cat','prop_freq')
    tmp$prop_freq <- round(tmp$prop_freq*100,1)
  }
  # if(yType='duration')
  # {
  #   ylab <- "Proportion of contact duration"
  #   tmp <- allPairSumData %>% group_by(sender, receiver, regular_index) %>% 
  #     summarise(tn = n())
  #   daily <- tmp[which(tmp$tn==52),]
  #   t_daily <- as.data.frame(round(prop.table(table(daily$regular_index)),2))
  #   weekly <- tmp[which(tmp$tn>=7),]
  #   t_weekly <- as.data.frame(round(prop.table(table(weekly$regular_index)),2))
  #   monthly <- tmp[which(tmp$tn>=2),]
  #   t_monthly <- as.data.frame(round(prop.table(table(monthly$regular_index)),2))
  # }
  
  xlm <- c(0,numGrp)
  ylm <- c(0,100)
  grid.newpage()
  pushViewport(plotViewport(c(3.5,3.5,1,1),xscale=xlm,yscale=ylm))
  for(bar_id in 1:numGrp)
  {
    grid.polygon(x = unit(c((bar_id-1),bar_id,bar_id,(bar_id-1),(bar_id-1)),'native'),
                 y = unit(c(0,0,
                            tmp$prop_freq[which(tmp$duration_cat==bar_id&tmp$regular_index==1)],
                            tmp$prop_freq[which(tmp$duration_cat==bar_id&tmp$regular_index==1)],
                            0),'native'),
                 gp = gpar(fill=colors[1], col='white'))
    
    grid.polygon(x = unit(c((bar_id-1),bar_id,bar_id,(bar_id-1),(bar_id-1)),'native'),
                 y = unit(c(tmp$prop_freq[which(tmp$duration_cat==bar_id&tmp$regular_index==1)],
                            tmp$prop_freq[which(tmp$duration_cat==bar_id&tmp$regular_index==1)],
                            100,100,
                            tmp$prop_freq[which(tmp$duration_cat==bar_id&tmp$regular_index==1)]),'native'),
                 gp = gpar(fill=colors[2],col='white'))
  }
  # grid.polygon(x = unit(c(0,1,1,0,0),'native'),
  #              y = unit(c(0,0,tmp$prop_freq[which(tmp$duration_cat==i,tmp$regular_index==1)],
  #                         t_daily$Freq[which(t_daily$Var1==1)],0),'native'),
  #              gp = gpar(fill=colors[1], col='white'))
  # grid.polygon(x = unit(c(0,1,1,0,0),'native'),
  #              y = unit(c(t_daily$Freq[which(t_daily$Var1==1)],t_daily$Freq[which(t_daily$Var1==1)],1,1,
  #                         t_daily$Freq[which(t_daily$Var1==1)]),'native'),
  #              gp = gpar(fill=colors[2],col='white'))
  
  grid.yaxis()
  grid.text(label = c("7 mins","21 mins","1 hr","2 hrs","4+ hrs"),
            x = unit(c(1,2,3,4,4.9),'native'),
            y = unit(-1,'lines'),gp=gpar(fontsize=fontsize))
  grid.text(label = ylab,x = unit(-3,'lines'),
            y = unit(0.5,'npc'), rot=90, gp=gpar(fontface='bold',fontsize=fontsize))
  grid.lines(unit(c(0,1,1,0,0),'npc'),unit(c(0,0,1,1,0),'npc'),gp=gpar(lwd=1))
  
  popViewport()
}

# draw the png ------------------------------
png(filename = paste0(outputfolder,'Fig3_proportion_durationCat_1231.png'),
    height = 12,width = 12,units = 'cm',pointsize = 12,res = 900)
proportionDurationPaneller(data=allPairSumData,yType='pairs',fontsize=14)
dev.off()

rm(allPairSumData)
# -------------------------------------------


load(paste0('working/combine/','rawPairs.Rdata'))
load('working/map_regular_index.rdata')
temp <- left_join(allPairData,map_regular_index,by=c('sender','receiver'))
temp$week <- week(temp$contact_date)
temp$month <- month(temp$contact_date)

rm(allPairData,map_regular_index)
# proportion of the duration paneller -------
proportionTimePaneller <- function(numGrp=3,data,yType,fontsize=14)
{
  if(yType=='pairs')
  {
    ylab <- "Proportion of contact events"
    tmp_daily <- data %>% group_by(sender, receiver, regular_index) %>%
      summarise(tn = length(unique(contact_date)))
    daily <- tmp_daily[which(tmp_daily$tn>=52),]
    t_daily <- as.data.frame(round(prop.table(table(daily$regular_index)),2))
    tmp_weekly <- data %>% group_by(sender, receiver, regular_index) %>%
      summarise(tn = length(unique(week)))
    weekly <- tmp_weekly[which(tmp_weekly$tn>7),]
    t_weekly <- as.data.frame(round(prop.table(table(weekly$regular_index)),2))
    tmp_monthly <- data %>% group_by(sender, receiver, regular_index) %>%
      summarise(tn = length(unique(month)))
    monthly <- tmp_monthly[which(tmp_monthly$tn>=2),]
    t_monthly <- as.data.frame(round(prop.table(table(monthly$regular_index)),2))
  }

  xlm <- c(0,numGrp)
  ylm <- c(0,1)
  grid.newpage()
  pushViewport(plotViewport(c(3.5,3.5,1,1),xscale=xlm,yscale=ylm))

  grid.polygon(x = unit(c(0,1,1,0,0),'native'),
               y = unit(c(0,0,t_daily$Freq[which(t_daily$Var1==1)],
                          t_daily$Freq[which(t_daily$Var1==1)],0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  grid.polygon(x = unit(c(0,1,1,0,0),'native'),
               y = unit(c(t_daily$Freq[which(t_daily$Var1==1)],t_daily$Freq[which(t_daily$Var1==1)],1,1,
                          t_daily$Freq[which(t_daily$Var1==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.polygon(x = unit(c(1,2,2,1,1),'native'),
               y = unit(c(0,0,t_weekly$Freq[which(t_weekly$Var1==1)],
                          t_weekly$Freq[which(t_weekly$Var1==1)],0),'native'),
               gp = gpar(fill=colors[1],col='white'))
  grid.polygon(x = unit(c(1,2,2,1,1),'native'),
               y = unit(c(t_weekly$Freq[which(t_weekly$Var1==1)],t_weekly$Freq[which(t_weekly$Var1==1)],
                          1,1,t_weekly$Freq[which(t_weekly$Var1==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))

  grid.polygon(x = unit(c(2,3,3,2,2),'native'),
               y = unit(c(0,0,t_monthly$Freq[which(t_monthly$Var1==1)],
                          t_monthly$Freq[which(t_monthly$Var1==1)],0),'native'),
               gp = gpar(fill=colors[1],col='white'))
  grid.polygon(x = unit(c(2,3,3,2,2),'native'),
               y = unit(c(t_monthly$Freq[which(t_monthly$Var1==1)],t_monthly$Freq[which(t_monthly$Var1==1)],
                          1,1, t_monthly$Freq[which(t_monthly$Var1==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  grid.yaxis()
  grid.text(label = c("Daily","Weekly","Monthly"),x = unit(c(0.5,1.5,2.5),'native'),
            y = unit(-1,'lines'),gp=gpar(fontsize=fontsize))
  grid.text(label = ylab,x = unit(-3,'lines'),
            y = unit(0.5,'npc'), rot=90, gp = gpar(fontface='bold',fontsize=fontsize))
  grid.lines(unit(c(0,1,1,0,0),'npc'),unit(c(0,0,1,1,0),'npc'),gp=gpar(lwd=1))
  
  popViewport()
}
# -------------------------------------------

# draw the png ------------------------------
png(filename = paste0(outputfolder,'Fig3_proportion_DWMly_upd_1231.png'),
    height = 12,width = 12,units = 'cm',pointsize = 12,res = 900)
proportionTimePaneller(data=temp,fontsize = 14,yType='pairs')
dev.off()
# -------------------------------------------


load(paste0(foldername,'PairDataContactType.Rdata'))
names(allPairSumData)
allPairSumData <- allPairSumData[,-9]
# i also dont need to know the contact date, freq and duration. only need to know the unique pairs?
unique_contact_pairs <- unique(allPairSumData[,c(1,2,9)])
acacia = T
cassia = T
source('code/data_process.R') # to get survey
rm(URID,USID)


unique_contact_pairs <- unique_contact_pairs[-which(unique_contact_pairs$sender%in%rm_sender),]
# keep all receivers with demographic data
# ### combine the unit information of sender and receiver 
t<- left_join(unique_contact_pairs,survey,by= c('sender'= 'tokenID'))
t<- left_join(t,survey,by= c('receiver'= 'tokenID'))
head(t)
# remove the unknown records
t <- as.data.frame(t)
t <- t[complete.cases(t),]
# proportion of the duration paneller -------
proportionLocPaneller <- function(numGrp=3,data,yType,fontsize)
{
  if(yType=='pairs')
  {
    ylab <- "Proportion of unique contact pairs"
    data$same_room <- ifelse(data$br.x==data$br.y,1,0)
    data$same_blk <- ifelse(data$blk.x==data$blk.y,1,0)
    data$same_blk_notroom <- ifelse(data$same_blk==1 & data$same_room==0,"sameBlkNotRoom",
                                               ifelse(data$same_blk==1 & data$same_room==1,"sameBlkSameRoom",
                                                      'nSameBlkNSameRoom'))
    data$same_dorm <- ifelse(data$dorm.x==data$dorm.y,1,0)
    data$same_dorm_notblk <- ifelse(data$same_dorm==1 & data$same_blk==0,"sameDormNotBlk",
                                               ifelse(data$same_dorm==1 & data$same_blk==1,"sameDormSameBlk",
                                                      'nSameDormNSameBlk'))
    tmp0 <- as.data.frame(prop.table(table(data$regular_index,data$same_room),margin=2))
    names(tmp0) <- c('regular_index','same_room','prop_freq_sameroom')
    tmp0$prop_freq_sameroom <- round(tmp0$prop_freq*100,1)
    
    tmp1 <- as.data.frame(prop.table(table(data$regular_index,data$same_blk_notroom),margin=2))
    names(tmp1) <- c('regular_index','same_blk','prop_freq_sameblk')
    tmp1$prop_freq_sameblk <- round(tmp1$prop_freq*100,1)

    tmp2 <- as.data.frame(prop.table(table(data$regular_index,data$same_dorm_notblk),margin=2))
    names(tmp2) <- c('regular_index','same_dorm','prop_freq_samedorm')
    tmp2$prop_freq_samedorm <- round(tmp2$prop_freq*100,1)
  }
  # if(yType='duration')
  # {
  #   ylab <- "Proportion of contact duration"
  #   tmp <- allPairSumData %>% group_by(sender, receiver, regular_index) %>% 
  #     summarise(tn = n())
  #   daily <- tmp[which(tmp$tn==52),]
  #   t_daily <- as.data.frame(round(prop.table(table(daily$regular_index)),2))
  #   weekly <- tmp[which(tmp$tn>=7),]
  #   t_weekly <- as.data.frame(round(prop.table(table(weekly$regular_index)),2))
  #   monthly <- tmp[which(tmp$tn>=2),]
  #   t_monthly <- as.data.frame(round(prop.table(table(monthly$regular_index)),2))
  # }
  
  xlm <- c(0,numGrp)
  ylm <- c(0,100)
  grid.newpage()
  pushViewport(plotViewport(c(3.5,3.5,1,1),xscale=xlm,yscale=ylm))
  
  grid.polygon(x = unit(c(0,1,1,0,0),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq_sameroom[which(tmp0$same_room==1&tmp0$regular_index==1)],
                          tmp0$prop_freq_sameroom[which(tmp0$same_room==1&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(0,1,1,0,0),'native'),
               y = unit(c(tmp0$prop_freq_sameroom[which(tmp0$same_room==1&tmp0$regular_index==1)],
                          tmp0$prop_freq_sameroom[which(tmp0$same_room==1&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq_sameroom[which(tmp0$same_room==1&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.polygon(x = unit(c(1,2,2,1,1),'native'),
               y = unit(c(0,0,
                          tmp1$prop_freq_sameblk[which(tmp1$same_blk=="sameBlkSameRoom"&tmp1$regular_index==1)],
                          tmp1$prop_freq_sameblk[which(tmp1$same_blk=="sameBlkSameRoom"&tmp1$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(1,2,2,1,1),'native'),
               y = unit(c(tmp1$prop_freq_sameblk[which(tmp1$same_blk=="sameBlkNotRoom"&tmp1$regular_index==1)],
                          tmp1$prop_freq_sameblk[which(tmp1$same_blk=="sameBlkNotRoom"&tmp1$regular_index==1)],
                          100,100,
                          tmp1$prop_freq_sameblk[which(tmp1$same_blk=="sameBlkSameRoom"&tmp1$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.polygon(x = unit(c(2,3,3,2,2),'native'),
               y = unit(c(0,0,
                          tmp2$prop_freq_samedorm[which(tmp2$same_dorm=="sameDormNotBlk"&tmp1$regular_index==1)],
                          tmp2$prop_freq_samedorm[which(tmp2$same_dorm=="sameDormNotBlk"&tmp1$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(2,3,3,2,2),'native'),
               y = unit(c(tmp2$prop_freq_samedorm[which(tmp2$same_dorm=="sameDormNotBlk"&tmp2$regular_index==1)],
                          tmp2$prop_freq_samedorm[which(tmp2$same_dorm=="sameDormNotBlk"&tmp2$regular_index==1)],
                          100,100,
                          tmp2$prop_freq_samedorm[which(tmp2$same_dorm=="sameDormNotBlk"&tmp2$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.yaxis()
  grid.text(label = c("Same Room","Same Block","Same Dorm"),
            x = unit(c(0.5,1.5,2.5),'native'),
            y = unit(-1,'lines'),gp=gpar(fontsize))
  grid.text(label = ylab,x = unit(-3,'lines'),
            y = unit(0.5,'npc'), rot=90, gp=gpar(fontface='bold',fontsize))
  grid.lines(unit(c(0,1,1,0,0),'npc'),unit(c(0,0,1,1,0),'npc'),gp=gpar(lwd=1))
  
  popViewport()
}

# draw the png ------------------------------
png(filename = paste0(outputfolder,'Fig3_proportion_loc_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionLocPaneller(yType='pairs',data=t,fontsize = 12)
dev.off()
# -------------------------------------------

# proportion of the Age paneller -------
proportionAgePaneller <- function(numGrp=6,sender_ageGrp="<24",data,yType,fontsize)
{
  if(yType=='pairs')
  {
    ylab <- "Proportion of unique contact pairs"
    subdata <- data %>% filter(age.x==sender_ageGrp)
    subdata$age.y[which(subdata$age.y=='>45')] <- '45+'
    tmp0 <- as.data.frame(prop.table(table(subdata$regular_index,subdata$age.y),margin=1)) 
    names(tmp0) <- c('regular_index','receiver_ageGrp','prop_freq')
    tmp0$prop_freq <- round(tmp0$prop_freq*100,1)

  }
  # if(yType='duration')
  # {
  #   ylab <- "Proportion of contact duration"
  #   tmp <- allPairSumData %>% group_by(sender, receiver, regular_index) %>% 
  #     summarise(tn = n())
  #   daily <- tmp[which(tmp$tn==52),]
  #   t_daily <- as.data.frame(round(prop.table(table(daily$regular_index)),2))
  #   weekly <- tmp[which(tmp$tn>=7),]
  #   t_weekly <- as.data.frame(round(prop.table(table(weekly$regular_index)),2))
  #   monthly <- tmp[which(tmp$tn>=2),]
  #   t_monthly <- as.data.frame(round(prop.table(table(monthly$regular_index)),2))
  # }
  
  xlm <- c(0,numGrp)
  ylm <- c(0,100)
  grid.newpage()
  pushViewport(plotViewport(c(3.5,3.5,1,1),xscale=xlm,yscale=ylm))
  
  grid.polygon(x = unit(c(0,1,1,0,0),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='<24'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='<24'&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(0,1,1,0,0),'native'),
               y = unit(c(tmp0$prop_freq[which(tmp0$receiver_ageGrp=='<24'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='<24'&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='<24'&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.polygon(x = unit(c(1,2,2,1,1),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='25-29'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='25-29'&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(1,2,2,1,1),'native'),
               y = unit(c(tmp0$prop_freq[which(tmp0$receiver_ageGrp=='25-29'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='25-29'&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='25-29'&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.polygon(x = unit(c(2,3,3,2,2),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='30-34'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='30-34'&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(2,3,3,2,2),'native'),
               y = unit(c(tmp0$prop_freq[which(tmp0$receiver_ageGrp=='30-34'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='30-34'&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='30-34'&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.polygon(x = unit(c(3,4,4,3,3),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='35-39'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='35-39'&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(3,4,4,3,3),'native'),
               y = unit(c(tmp0$prop_freq[which(tmp0$receiver_ageGrp=='35-39'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='35-39'&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='35-39'&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.polygon(x = unit(c(4,5,5,4,4),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='40-44'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='40-44'&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(4,5,5,4,4),'native'),
               y = unit(c(tmp0$prop_freq[which(tmp0$receiver_ageGrp=='40-44'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='40-44'&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='40-44'&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  grid.polygon(x = unit(c(5,6,6,5,5),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='45+'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='45+'&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(5,6,6,5,5),'native'),
               y = unit(c(tmp0$prop_freq[which(tmp0$receiver_ageGrp=='45+'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='45+'&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq[which(tmp0$receiver_ageGrp=='45+'&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  grid.yaxis()
  grid.text(label = c("<24","25-30","30-34",'35-39','40-44','45+'),
            x = unit(c(0.5,1.5,2.5,3.5,4.5,5.5),'native'),
            y = unit(-1,'lines'),gp=gpar(fontsize))
  grid.text(label = ylab,x = unit(-3,'lines'),
            y = unit(0.5,'npc'), rot=90, gp=gpar(fontface='bold',fontsize))
  grid.lines(unit(c(0,1,1,0,0),'npc'),unit(c(0,0,1,1,0),'npc'),gp=gpar(lwd=1))
  
  popViewport()
}

# draw the png ------------------------------
png(filename = paste0(outputfolder,'Fig3_proportion_ageless25_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionAgePaneller(yType='pairs',data=t,fontsize = 12)
dev.off()

png(filename = paste0(outputfolder,'Fig3_proportion_age2529_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionAgePaneller(yType='pairs',sender_ageGrp="25-29",data=t,fontsize = 12)
dev.off()

png(filename = paste0(outputfolder,'Fig3_proportion_age3034_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionAgePaneller(yType='pairs',sender_ageGrp="30-34",data=t,fontsize = 12)
dev.off()

png(filename = paste0(outputfolder,'Fig3_proportion_age3539_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionAgePaneller(yType='pairs',sender_ageGrp="35-39",data=t,fontsize = 12)
dev.off()

png(filename = paste0(outputfolder,'Fig3_proportion_age4044_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionAgePaneller(yType='pairs',sender_ageGrp="40-44",data=t,fontsize = 12)
dev.off()


png(filename = paste0(outputfolder,'Fig3_proportion_age45plus_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionAgePaneller(yType='pairs',sender_ageGrp=">45",data=t,fontsize = 12)
dev.off()
# -------------------------------------------

# proportion of the original country paneller -------
proportionCountryPaneller <- function(numGrp=5,sender_country="Bangladesh",data,yType='pairs',fontsize)
{
  if(yType=='pairs')
  {
    ylab <- "Proportion of unique contact pairs"
    subdata <- data %>% filter(country.x==sender_country)
    # subdata$age.y[which(subdata$country.y=='>45')] <- '45+'
    tmp0 <- as.data.frame(prop.table(table(subdata$regular_index,subdata$country.y),margin=2))
    names(tmp0) <- c('regular_index','receiver_country','prop_freq')
    tmp0$prop_freq <- round(tmp0$prop_freq*100,1)
    
  }
  # if(yType='duration')
  # {
  #   ylab <- "Proportion of contact duration"
  #   tmp <- allPairSumData %>% group_by(sender, receiver, regular_index) %>% 
  #     summarise(tn = n())
  #   daily <- tmp[which(tmp$tn==52),]
  #   t_daily <- as.data.frame(round(prop.table(table(daily$regular_index)),2))
  #   weekly <- tmp[which(tmp$tn>=7),]
  #   t_weekly <- as.data.frame(round(prop.table(table(weekly$regular_index)),2))
  #   monthly <- tmp[which(tmp$tn>=2),]
  #   t_monthly <- as.data.frame(round(prop.table(table(monthly$regular_index)),2))
  # }
  
  xlm <- c(0,numGrp)
  ylm <- c(0,100)
  grid.newpage()
  pushViewport(plotViewport(c(3.5,3.5,1,1),xscale=xlm,yscale=ylm))
  
  grid.polygon(x = unit(c(0,1,1,0,0),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq[which(tmp0$receiver_country=='Bangladesh'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_country=='Bangladesh'&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(0,1,1,0,0),'native'),
               y = unit(c(tmp0$prop_freq[which(tmp0$receiver_country=='Bangladesh'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_country=='Bangladesh'&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq[which(tmp0$receiver_country=='Bangladesh'&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.polygon(x = unit(c(1,2,2,1,1),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq[which(tmp0$receiver_country=='India'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_country=='India'&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(1,2,2,1,1),'native'),
               y = unit(c(tmp0$prop_freq[which(tmp0$receiver_country=='India'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_country=='India'&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq[which(tmp0$receiver_country=='India'&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.polygon(x = unit(c(2,3,3,2,2),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq[which(tmp0$receiver_country=='China'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_country=='China'&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(2,3,3,2,2),'native'),
               y = unit(c(tmp0$prop_freq[which(tmp0$receiver_country=='China'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_country=='China'&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq[which(tmp0$receiver_country=='China'&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.polygon(x = unit(c(3,4,4,3,3),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq[which(tmp0$receiver_country=='Myanmar'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_country=='Myanmar'&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(3,4,4,3,3),'native'),
               y = unit(c(tmp0$prop_freq[which(tmp0$receiver_country=='Myanmar'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_country=='Myanmar'&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq[which(tmp0$receiver_country=='Myanmar'&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.polygon(x = unit(c(4,5,5,4,4),'native'),
               y = unit(c(0,0,
                          tmp0$prop_freq[which(tmp0$receiver_country=='Others'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_country=='Others'&tmp0$regular_index==1)],
                          0),'native'),
               gp = gpar(fill=colors[1], col='white'))
  
  grid.polygon(x = unit(c(4,5,5,4,4),'native'),
               y = unit(c(tmp0$prop_freq[which(tmp0$receiver_country=='Others'&tmp0$regular_index==1)],
                          tmp0$prop_freq[which(tmp0$receiver_country=='Others'&tmp0$regular_index==1)],
                          100,100,
                          tmp0$prop_freq[which(tmp0$receiver_country=='Others'&tmp0$regular_index==1)]),'native'),
               gp = gpar(fill=colors[2],col='white'))
  
  grid.yaxis()
  grid.text(label = c("Bangladesh","India","China",'Myanmar','Others'),
            x = unit(c(0.5,1.5,2.5,3.5,4.5),'native'),
            y = unit(-1,'lines'),gp=gpar(fontsize))
  grid.text(label = ylab,x = unit(-3,'lines'),
            y = unit(0.5,'npc'), rot=90, gp=gpar(fontface='bold',fontsize))
  grid.lines(unit(c(0,1,1,0,0),'npc'),unit(c(0,0,1,1,0),'npc'),gp=gpar(lwd=1))
  
  popViewport()
}

# draw the png ------------------------------
png(filename = paste0(outputfolder,'Fig3_proportion_countryBangladesh_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionCountryPaneller(yType='pairs',data=t,fontsize = 12)
dev.off()

png(filename = paste0(outputfolder,'Fig3_proportion_countryIndia_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionCountryPaneller(yType='pairs',sender_country="India",data=t,fontsize = 12)
dev.off()

png(filename = paste0(outputfolder,'Fig3_proportion_countryChina_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionCountryPaneller(yType='pairs',sender_country="China",data=t,fontsize = 12)
dev.off()

png(filename = paste0(outputfolder,'Fig3_proportion_countryMyanmar_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionCountryPaneller(yType='pairs',sender_country="Myanmar",data=t,fontsize = 12)
dev.off()

png(filename = paste0(outputfolder,'Fig3_proportion_countryOthers_upd_1209.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
proportionCountryPaneller(yType='pairs',sender_country="Others",data=t,fontsize = 12)
dev.off()


# -------------------------------------------


# reply reviewers (two by two table) --------------------------------------
proportionAgeTabeller <- function(sender_ageGrp="<24",data)
{
  # if(yType=='pairs')
  # {
    # ylab <- "Proportion of unique contact pairs"
    subdata <- data %>% filter(age.x==sender_ageGrp)
    subdata$age.y[which(subdata$age.y=='>45')] <- '45+'
    tmp0 <- as.data.frame(prop.table(table(subdata$regular_index,subdata$age.y),margin=2)) # separate the regular and transient
    names(tmp0) <- c('regular_index','receiver_ageGrp','prop_freq')
    tmp0$prop_freq <- round(tmp0$prop_freq*100,1)
    
  # }
temp_age <- unique(tmp0$receiver_ageGrp)
# create the showing table
age_table = NULL
# for(j in 1:length(temp_age)){
  tab <- NULL
  for(i in 1:length(temp_age)){
    temp_tab <- data.frame(
                           # Age = c(as.character(temp_age[j]),rep(NA,length(temp_age))),
                           Receiver_Age = temp_age[i],
                           Regular = tmp0$prop_freq[which(tmp0$receiver_ageGrp==temp_age[i]&tmp0$regular_index==1)],
                           Transit= tmp0$prop_freq[which(tmp0$receiver_ageGrp==temp_age[i]&tmp0$regular_index==0)]
    )
    
    tab <- rbind(tab,temp_tab)
  }
  tab <- cbind(Age = c(as.character(temp_age[j]),rep(NA,length(temp_age)-1)),tab)
  # age_table <- rbind(age_table,tab)
# }
  tab
}
age_j <- unique(t$age.x)
age_table <- NULL
for(j in 1:length(age_j)){

  temp_age_table <- proportionAgeTabeller(sender_ageGrp = as.character(age_j[j]),data = t)
  # update after discuss with Alex on Nov 22
  temp_age_table <- temp_age_table[,-4]
  temp_age_table$p_value <- c(round(t.test(temp_age_table$Regular-mean(temp_age_table$Regular))$p.value,3),
                              rep(NA,(length(temp_age_table$Receiver_Age)-1)))
  
  # temp_age_table$p_value <- do.call("c",lapply(1:nrow(temp_age_table),function(i){
  #   round(t.test(c(temp_age_table$Regular[i],temp_age_table$Transit[i]))$p.value,3)}))
  age_table <- rbind(age_table,temp_age_table)
}

write.csv(age_table,"review_16_table_1_updated.csv",row.names = F,na = "")



proportionCountryTabeller <- function(sender_Grp,data)
{
  # if(yType=='pairs')
  # {
  # ylab <- "Proportion of unique contact pairs"
  subdata <- data %>% filter(country.x==sender_Grp)
  # subdata$age.y[which(subdata$age.y=='>45')] <- '45+'
  tmp0 <- as.data.frame(prop.table(table(subdata$regular_index,subdata$country.y),margin=2)) # separate the regular and transient
  names(tmp0) <- c('regular_index','receiver_Grp','prop_freq')
  tmp0$prop_freq <- round(tmp0$prop_freq*100,1)
  
  # }
  temp_country <- unique(data$country.x)
  # create the showing table
  country_table = NULL
  # for(j in 1:length(temp_age)){
  tab <- NULL
  for(i in 1:length(temp_country)){
    temp_tab <- data.frame(
      # Age = c(as.character(temp_age[j]),rep(NA,length(temp_age))),
      Receiver_country = temp_country[i],
      Regular = tmp0$prop_freq[which(tmp0$receiver_Grp==temp_country[i]&tmp0$regular_index==1)],
      Transit= tmp0$prop_freq[which(tmp0$receiver_Grp==temp_country[i]&tmp0$regular_index==0)]
    )
    
    tab <- rbind(tab,temp_tab)
  }
  tab <- cbind(country = c(as.character(temp_country)[j],rep(NA,length(temp_country)-1)),tab)
  # age_table <- rbind(age_table,tab)
  # }
  tab
}
country_j <- unique(t$country.x)
country_table <- NULL
t$country.x <- as.character(t$country.x)
for(j in 1:length(country_j)){
  
  temp_country_table <- proportionCountryTabeller(sender_Grp = as.character(country_j[j]),data = t)
  temp_country_table <- temp_country_table[,-4]
  temp_country_table$p_value <- c(round(t.test(temp_country_table$Regular-mean(temp_country_table$Regular))$p.value,3),
                              rep(NA,(length(temp_country_table$Receiver_country)-1)))
  
  # temp_country_table$p_value <- do.call("c",lapply(1:nrow(temp_country_table),function(i){
  #   round(t.test(c(temp_country_table$Regular[i],temp_country_table$Transit[i]))$p.value,3)}))
  country_table <- rbind(country_table,temp_country_table)
}

write.csv(country_table,"review_16_table_2_update.csv",row.names = F,na = "")

