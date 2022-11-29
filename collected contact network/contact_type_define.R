load(paste0(foldername,'contact_pair/toDefineContactType.Rdata'))
# Define the regular contact and transit contact ------------------------------
## Close Contact and meet 4 times a week will be defined as Regular Contact
allPairSumData$nweek <- week(allPairSumData$contact_date)
### Close Contact (MOH) definition: contact more than 15 mins, here we put threshold as 20 mins
allPairSumData$close_contact_index <- ifelse(allPairSumData$mean_contact_consecutive_duration>20,1,0)
### number of contact-days per week 
tmp <- allPairSumData %>% group_by(sender, receiver, nweek) %>%
  summarise(cloContNumPerWeek = length(which(close_contact_index==1)))
tmp$regular_contact <- ifelse(tmp$cloContNumPerWeek>=4,1,0)
tmp <- unique(tmp)
tmp <- tmp %>% group_by(sender,receiver) %>% summarise(regular_index=sum(unique(regular_contact)))
allPairSumData <- left_join(allPairSumData,tmp,
                            by=c("sender","receiver","nweek"))


rm(tmp)
save(allPairSumData,file=paste0(foldername,'PairDataContactType.Rdata'))
rm(allPairSumData)