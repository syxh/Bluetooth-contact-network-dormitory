# process the results for better plot present
# step 1: to combine the scenarios results into 

fn <- list.files(paste0(getwd(),'/policies/working'))
split_fn_acacia <- list()
split_fn_cassia <- list()
for(i in 1:7)
{
  split_fn_acacia[[i]] <- grep(paste0('Scenario_',i,'_theta_0.4'),fn,value=TRUE)#fn[substr(fn,10,10)==i & substr(fn,20,21)==40]
  split_fn_cassia[[i]] <-  grep(paste0('Scenario_',i,'_theta_0.87'),fn,value=TRUE)#fn[substr(fn,10,10)==i & substr(fn,20,21)==87]
}

gen_CI <- function(data)
{
  CI=matrix(0,dim(data)[1],3)
  for(j in 1:dim(data)[1])CI[j,] = quantile(as.numeric(data[j,-1]),c(25,500,975)/1000)
  CI = as.data.frame(CI)
  names(CI) = c("lower",'central','upper')
  CI
}

for(scenario_id in 1:7){
  combine_acacia <- data.frame("Day"=c(1:120))
  combine_cassia <- data.frame("Day"=c(1:120))
  
  for(i in 1:length(split_fn_acacia[[1]]))
  {
    tmp1 <- read.csv(paste0('policies/working/',split_fn_acacia[[scenario_id]][[i]]))
    tmp1[1,1] <- 20
    # tmp1 <- cumsum(tmp1)
    combine_acacia <- cbind(combine_acacia,tmp1)
    
    tmp2 <- read.csv(paste0('policies/working/',split_fn_cassia[[scenario_id]][[i]]))
    tmp2[1,1] <- 20
    # tmp2 <- cumsum(tmp2)
    combine_cassia <- cbind(combine_cassia,tmp2)
  }
  write.csv(combine_acacia,file=paste0('policies/working/combine_Scenario_',scenario_id,'_theta_acacia.csv'),row.names = F)
  CI <- gen_CI(combine_acacia)
  # CI$doy <- c(1:nrow(CI))
  write.csv(CI,paste0('policies/output/data/Scenario_',scenario_id,'_acacia_daily.csv'),row.names = FALSE)
  
  write.csv(combine_cassia,file=paste0('policies/working/combine_Scenario_',scenario_id,'_theta_cassia_daily.csv'),row.names = F)
  CI <- gen_CI(combine_cassia)
  # CI$doy <- c(1:nrow(CI))
  write.csv(CI,paste0('policies/output/data/Scenario_',scenario_id,'_cassia_daily.csv'),row.names = FALSE)
}

rm(tmp1,tmp2,CI,i,j,fn,scenario_id)

