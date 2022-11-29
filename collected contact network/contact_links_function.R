# - Explore the cliques according to dataset. 
# load(paste0(foldername,'PairDataContactType.Rdata'))
order1_cliq_network <- function(data)
{
  network_D1 <- list()
  # s <- unique(data$sender)
  s <- unique(data$index_ID_sender)
  
  for(i in 1:length(s))
  {
    # i=1
    subs <- s[i]
    subdata <- data %>% filter(index_ID_sender%in%subs)
    dirct_r <- unique(subdata$index_ID_receiver) # the receiver (d=1)
    
    # data <- anti_join(data,subdata)
    d1 <- list()
    cliq_d1 <- list()
    if(length(which(dirct_r%in%s))>0){
      d1S <- dirct_r[which(dirct_r%in%s)]
      
      for(j in 1:length(d1S))
      {
        d1Sub <- data %>% filter(index_ID_sender%in%d1S[j])
        d1_r <- unique(d1Sub$index_ID_receiver)
        tmp_cliq_d1 <- d1_r[which(d1_r%in%dirct_r)]
        cliq_d1[[j]] <- list(s=d1S[j],r=tmp_cliq_d1)
      }
      # interestnodes<- c(subs,unique(unlist(cliq_d1)))
    }
    d1 <- list(sender_index=subs,
               initial_receiver = dirct_r,
               kids = cliq_d1)
    network_D1[[i]] <- d1
    names(network_D1)[[i]] <-s[i]
  }
  network_D1
}

generate_plot_data <- function(networkList)
{
  plotData <- list()
  for(n in 1:length(networkList))
  {
    sub_data <- networkList[[n]]
    if(length(sub_data$kids)>0){
      sub_g_nodes_critical <- unique(c(sub_data$sender_index,unlist(lapply(1:length(sub_data$kids),function(i) sub_data$kids[[i]]$s)),unlist(lapply(1:length(sub_data$kids),function(i) sub_data$kids[[i]]$r))))
      sub_g_nodes_others <- setdiff(sub_data$initial_receiver,unique(c(unlist(lapply(1:length(sub_data$kids),function(i) sub_data$kids[[i]]$s)),unlist(lapply(1:length(sub_data$kids),function(i) sub_data$kids[[i]]$r)))))
      sub_g_edges <- data.frame(from = c(rep(sub_data$sender_index,length(sub_data$initial_receiver)),unlist(lapply(1:length(sub_data$kids),function(i){rep(sub_data$kids[[i]]$s,length(sub_data$kids[[i]]$r))}))),
                                to = c(sub_data$initial_receiver,unlist(lapply(1:length(sub_data$kids),function(i) sub_data$kids[[i]]$r))))
    }else{
      sub_g_nodes_critical <- NULL
      sub_g_nodes_others <- sub_data$initial_receiver
      sub_g_edges <- data.frame(from = rep(sub_data$sender_index,length(sub_data$initial_receiver)),
                                to = sub_data$initial_receiver)
    }
    
    plotData[[n]] <- list(nodes_critical = sub_g_nodes_critical,
                          nodes_others = sub_g_nodes_others,
                          edges = sub_g_edges)
  }
  plotData
  
}
contact_links <- function(selected_date)
{
  # cat('calculation----',selected_date,'------\n')
  # Aim: find the triangle cliques on the selected date
  dailyRawRegData <- allPairSumData %>% filter(contact_date==selected_date,
                                               regular_index==1)
  # head(dailyRawRegData)
  # temp_reg_send <- unique(dailyRawRegData$sender)
  # temp_reg_rec <- unique(dailyRawRegData$receiver)
  # rec_also_send_reg <- temp_reg_rec[which(temp_reg_rec%in%temp_reg_send)]
  # on March 10, there are 847 receivers in the senders list (the total number of sender on that day is 930)
  
  dailyRawTransData <- allPairSumData %>% filter(contact_date==selected_date,
                                                 regular_index==0)
  # temp_trans_send <- unique(dailyRawTransData$sender)
  # temp_trans_rec <- unique(dailyRawTransData$receiver)
  # rec_also_send_trans <- temp_trans_rec[which(temp_trans_rec%in%temp_trans_send)]
  # on March 10, there are 943 receivers in the senders list (the total number of sender on that day is 957)
  
  # allPairSumData <- allPairSumData[order(allPairSumData$sender,allPairSumData$receiver),]
  # head(allPairSumData)
  
  
  
  # trace all the edges from each individual (order = 1)
  # going to find the sub-triangle among the total network, aim to find the strongest connection within the network.
  # actually, it is only the part one, need another round screening to get the only cliques members. that would be subset of the current network.(cliques alrdy)
  
  ## replace the token ID to the index ID (i.e. 1 to xxx)
  
  map_all <- unique(c(dailyRawRegData$sender,dailyRawRegData$receiver,dailyRawTransData$sender,dailyRawTransData$receiver))
  map_all <- data.frame(token_ID = map_all, index_ID = c(1:length(map_all)))
  dailyRawRegData <- left_join(dailyRawRegData,map_all,by=c("sender"="token_ID"))
  dailyRawRegData <- left_join(dailyRawRegData,map_all,by=c("receiver"="token_ID"))
  names(dailyRawRegData)[11:12] <- c("index_ID_sender","index_ID_receiver")
  dailyRawTransData <- left_join(dailyRawTransData,map_all,by=c("sender"="token_ID"))
  dailyRawTransData <- left_join(dailyRawTransData,map_all,by=c("receiver"="token_ID"))
  names(dailyRawTransData)[11:12] <- c("index_ID_sender","index_ID_receiver")
  # write.csv(map_all,paste0('working/links/map_tokenID_indexID_',selected_date,".csv"),row.names = F)

  # a = Sys.time()
  # order1_cliq_network_reg <- order1_cliq_network(dailyRawRegData)
  order1_cliq_network_trans <- order1_cliq_network(dailyRawTransData)
  
  # save(order1_cliq_network_reg,file='working/order1_cliq_network_reg.Rdata')
  # b=Sys.time()
  # 4mins of the regular contact on March 10
  
  # a = Sys.time()
  # order1_cliq_network_trans <- order1_cliq_network(dailyRawTransData)
  # b=Sys.time()
  
  ## try to draw the cliques plot

  # plot_data <- generate_plot_data(networkList = order1_cliq_network_reg)
  plot_data <- generate_plot_data(networkList = order1_cliq_network_trans)
  
  # try to combine the plot data (i.e. edges into one data set, remove the duplicated-direction link-nodes)
  combind_edge <- do.call('rbind',lapply(1:length(plot_data),function(i) plot_data[[i]]$edges))
  check_from_id <- which(combind_edge$from %in% combind_edge$to)
  temp_check_from <- unique(combind_edge$from[check_from_id])
  rm_ids <- c()
  for(n in 1:length(temp_check_from))
  {
    tmpid <- which(combind_edge$to==temp_check_from[n])
    t <- which(combind_edge$from[tmpid]%in%combind_edge$to[which(combind_edge$from==temp_check_from[n])])
    rm_id <- tmpid[t]
    rm_ids <- c(rm_ids,rm_id)
  }
  combind_edge_unique <- combind_edge[-rm_ids,]
  write.csv(combind_edge_unique,paste0('working/links/contact_on_',selected_date,'.csv'),row.names=F)
  return(combind_edge_unique)
  rm(combind_edge,check_from_id,temp_check_from,rm_ids,dailyRawRegData,dailyRawTransData,plot_data)
  
}

