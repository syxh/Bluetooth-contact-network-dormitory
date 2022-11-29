library(grid)

source("code/functionsOfPlot.R")
load(paste0(foldername,'PairDataContactType.Rdata'))
load(paste0(foldername,'rawPairs.Rdata'))
## Merge the data above, to get the data with the contact start_hour, and information of contact type.
tt <- unique(allPairSumData[,c("sender","receiver","regular_index")])
tmp <- left_join(allPairData,tt,by=c("sender","receiver"))
#wavk test Doy of Week
numPairWdayHour <- tmp %>% group_by(contact_date,contact_wday,start_hour,regular_index) %>% 
  summarise(numPairs=n())
t_wd <- numPairWdayHour%>%group_by(contact_date)%>%summarise(sum_pairs =sum(numPairs))
t_wd_r <- t_wd%>%filter(regular_index==1)
t_wd_t <- t_wd%>%filter(regular_index==0)
wavk_test(t_wd_r$sum_pairs~t_wd_r$contact_date)
wavk_test(t_wd_r$sum_pairs~t_wd_r$contact_date)

numPairWdayHour <- tmp %>% group_by(contact_date,contact_wday,start_hour,regular_index) %>% 
  summarise(numPairs=n()) %>% ungroup() %>%
  group_by(contact_wday,start_hour,regular_index) %>% 
  summarise(meanNumPairs = mean(numPairs))

numPairDailyHour <- tmp %>% group_by(contact_date,start_hour,regular_index) %>% 
  summarise(numPairs=n()) %>% ungroup() %>%
  group_by(start_hour,regular_index) %>% 
  summarise(meanNumPairs = mean(numPairs))

durationPairWdayHour <- tmp %>% group_by(contact_wday,start_hour,regular_index) %>% 
  summarise(durationPairs=mean(contact_consecutive_duration)/60)

rm(tt,tmp,allPairData,allPairSumData)

timeLinePaneller <- function(dt,yType)
{
  xlm <- c(0,length(paste0(rep(1:7,each=24),'-',c(1:24))))
  # colors <- c("#ff7373","#6897bb") # regular (red) and transient (blue)
  if(yType=='pairs')
  {
    ylm <- c(0,max(dt$meanNumPairs[which(dt$regular_index==1)]+dt$meanNumPairs[which(dt$regular_index==0)])*1.1)
    ylab <- c("Number of interactions (on the hour scale)")
  }
  if(yType=='duration')
  {
    ylm <- c(0,max(dt$durationPairs[which(dt$regular_index==1)]+dt$durationPairs[which(dt$regular_index==0)])*1.1)
    ylab <- c("Contact duration (hrs)")
  }
  
  grid.newpage()
  pushViewport(plotViewport(c(3.5,3.5,1,1),xscale=xlm,yscale=ylm))
  xaxes2(xlm)
  tmp_r <- as.data.frame(dt[which(dt$regular_index==1),])
  tmp_r <- rbind(tmp_r,tmp_r[168,])
    tmp_r <- tmp_r[order(tmp_r$contact_wday,tmp_r$start_hour),]
    
  tmp_r$seq <- c(1:nrow(tmp_r))
  grid.polygon(x = unit(c(tmp_r$seq-1,rev(tmp_r$seq-1)),'native'), 
               y = unit(c(rep(0,length(tmp_r[,4])),rev(tmp_r[,4])),'native'),
               gp = gpar(fill=colors[1]))
  tmp_t <- as.data.frame(dt[which(dt$regular_index==0),])
  tmp_t <- rbind(tmp_t,tmp_t[168,])
  
    tmp_t <- tmp_t[order(tmp_t$contact_wday,tmp_t$start_hour),]
    
  tmp_t$seq <- c(1:nrow(tmp_t))
  grid.polygon(x = unit(c(tmp_t$seq-1,rev(tmp_t$seq-1)),'native'), 
               y = unit(c(tmp_r[,4],rev(tmp_t[,4]+tmp_r[,4])),'native'),
               gp = gpar(fill=colors[2]))
  for(i in 1:7)
  {
    grid.lines(x=unit(unit(0 +0 +(i-1)*24,'native'),'native'),y = unit(c(0,1),'npc'),
               gp=gpar(lty='dotted',col='gray'))
    grid.lines(x=unit(unit(0 +6 +(i-1)*24,'native'),'native'),y = unit(c(0,1),'npc'),
               gp=gpar(lty='dotted',col='gray'))
    grid.lines(x=unit(unit(0 +12 +(i-1)*24,'native'),'native'),y = unit(c(0,1),'npc'),
               gp=gpar(lty='dotted',col='gray'))
    grid.lines(x=unit(unit(0 +18 +(i-1)*24,'native'),'native'),y = unit(c(0,1),'npc'),
               gp=gpar(lty='dotted',col='gray'))
  }
  grid.lines(unit(c(0,168,168,0,0),'native'),unit(c(0,0,1,1,0),'npc'),gp=gpar(lwd=1.5))
  yaxis(ylm)
  grid.text(ylab, x=unit(-3,'lines'),
            rot=90,gp=gpar(fontface='bold'))
  grid.rect(x=unit(0.82,'npc'),y=unit(0.98,'npc'),
            width = unit(0.02,'npc'),height = unit(0.01,'npc'),
            gp = gpar(col=colors[2],fill=colors[2]))
  grid.text('Transient',x=unit(0.85,'npc'),y=unit(0.98,'npc'),just='left')
  grid.rect(x=unit(0.82,'npc'),y=unit(0.92,'npc'),
            width = unit(0.02,'npc'),height = unit(0.01,'npc'),
            gp = gpar(col=colors[1],fill=colors[1]))
  grid.text('Regular',x=unit(0.85,'npc'),y=unit(0.92,'npc'),just = 'left')
  popViewport()
}

# png(filename = paste0(outputfolder,'Fig1_numberPairs_1014.png'),
#     height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
png(filename = paste0(outputfolder,'Fig1_numberPairs_1231.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
timeLinePaneller(numPairWdayHour,yType='pairs')
dev.off()

png(filename = paste0(outputfolder,'Fig1_duration_0124.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
timeLinePaneller(durationPairWdayHour,yType='duration')
dev.off()
rm(durationPairWdayHour,numPairWdayHour)
# rm(ylab,ylm,yType,xlm,xaxes,xaxes2,yaxis,tmp_r,tmp_t,dt)




png(filename = paste0(outputfolder,'/FigS5a.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
timeLinePaneller(numPairWdayHour,yType='pairs')
dev.off()

timeLinePanellerDaily <- function(dt,yType)
{
  # xlm <- c(0,length(paste0(rep(1:7,each=24),'-',c(1:24))))
  xlm <- c(0,24)
  # colors <- c("#ff7373","#6897bb") # regular (red) and transient (blue)
  if(yType=='pairs')
  {
    ylm <- c(0,max(dt$meanNumPairs[which(dt$regular_index==1)]+dt$meanNumPairs[which(dt$regular_index==0)])*1.1)
    ylab <- c("Number of interactions (on the hour scale)")
  }
  if(yType=='duration')
  {
    ylm <- c(0,max(dt$durationPairs[which(dt$regular_index==1)]+dt$durationPairs[which(dt$regular_index==0)])*1.1)
    ylab <- c("Contact duration (hrs)")
  }
  
  grid.newpage()
  pushViewport(plotViewport(c(3.5,3.5,1,1),xscale=xlm,yscale=ylm))
  xaxes3(xlm)
  tmp_r <- as.data.frame(dt[which(dt$regular_index==1),])
  tmp_r <- rbind(tmp_r,tmp_r[24,])
  # tmp_r <- tmp_r[order(tmp_r$contact_wday,tmp_r$start_hour),]
  
  tmp_r$seq <- c(1:nrow(tmp_r))
  grid.polygon(x = unit(c(tmp_r$seq-1,rev(tmp_r$seq-1)),'native'), 
               y = unit(c(rep(0,length(tmp_r[,3])),rev(tmp_r[,3])),'native'),
               gp = gpar(fill=colors[1]))
  tmp_t <- as.data.frame(dt[which(dt$regular_index==0),])
  tmp_t <- rbind(tmp_t,tmp_t[24,])
  
  # tmp_t <- tmp_t[order(tmp_t$contact_wday,tmp_t$start_hour),]
  
  tmp_t$seq <- c(1:nrow(tmp_t))
  grid.polygon(x = unit(c(tmp_t$seq-1,rev(tmp_t$seq-1)),'native'), 
               y = unit(c(tmp_r[,3],rev(tmp_t[,3]+tmp_r[,3])),'native'),
               gp = gpar(fill=colors[2]))
  # for(i in 1:7)
  # {
    grid.lines(x=unit(unit(0 +0 +(i-1)*24,'native'),'native'),y = unit(c(0,1),'npc'),
               gp=gpar(lty='dotted',col='gray'))
    grid.lines(x=unit(unit(0 +6 +(i-1)*24,'native'),'native'),y = unit(c(0,1),'npc'),
               gp=gpar(lty='dotted',col='gray'))
    grid.lines(x=unit(unit(0 +12 +(i-1)*24,'native'),'native'),y = unit(c(0,1),'npc'),
               gp=gpar(lty='dotted',col='gray'))
    grid.lines(x=unit(unit(0 +18 +(i-1)*24,'native'),'native'),y = unit(c(0,1),'npc'),
               gp=gpar(lty='dotted',col='gray'))
  # }
  grid.lines(unit(c(0,24,24,0,0),'native'),unit(c(0,0,1,1,0),'npc'),gp=gpar(lwd=1.5))
  yaxis(ylm)
  grid.text(ylab, x=unit(-3,'lines'),
            rot=90,gp=gpar(fontface='bold'))
  grid.rect(x=unit(0.82,'npc'),y=unit(0.98,'npc'),
            width = unit(0.02,'npc'),height = unit(0.01,'npc'),
            gp = gpar(col=colors[2],fill=colors[2]))
  grid.text('Transient',x=unit(0.85,'npc'),y=unit(0.98,'npc'),just='left')
  grid.rect(x=unit(0.82,'npc'),y=unit(0.92,'npc'),
            width = unit(0.02,'npc'),height = unit(0.01,'npc'),
            gp = gpar(col=colors[1],fill=colors[1]))
  grid.text('Regular',x=unit(0.85,'npc'),y=unit(0.92,'npc'),just = 'left')
  popViewport()
}

png(filename = paste0(outputfolder,'/FigS5b.png'),
    height = 12,width = 12,units = 'cm',pointsize = 10,res = 900)
timeLinePanellerDaily(numPairDailyHour,yType='pairs')
dev.off()
