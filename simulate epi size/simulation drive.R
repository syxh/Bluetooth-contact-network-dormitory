# set.seed(1234)
alrdy_id <- c(1083,1277,1516,1884,3873,3904,4248,5508,5775,6163,6556,6946,
              7336,7757,8110,8224,8558,8714,8726,8903,9124,9746)

set.seed(9876)
start_time = Sys.time()
# source('code/branching_function.r')
# source('code/branching_setup.r')
source('code/branching_function_sim_contact.r')
source('code/branching_setup_sim_contact.r')

# simulation--------------------------------------------------------



# load the simulated contact matrix
load('working/simulate contact matrix/simulate_contact.Rdata')
sim_regular_contact <- 1*as.matrix(sim_contact==1)
sim_regular_contact<- as.data.frame(sim_regular_contact + t(sim_regular_contact))
sim_transient_contact <- 1*as.matrix(sim_contact==2)
sim_transient_contact <- as.data.frame(sim_transient_contact + t(sim_transient_contact))


contact_ind_index <- 1:nrow(sim_regular_contact)
rm(sim_contact)

#Load the simulated transmission tree (i.e weighted contact matrix)
load('working/simulate contact matrix/simulate_regular_contact_weight.Rdata')
load('working/simulate contact matrix/simulate_transient_contact_weight.Rdata')
weight_sim_regular_contact <- weight_sim_regular_contact+t(weight_sim_regular_contact)
weight_sim_transient_contact <- weight_sim_transient_contact+t(weight_sim_transient_contact)
# to reduce the 50% of layout above room, edit on March 30.
weight_sim_regular_contact[weight_sim_regular_contact<1] <- weight_sim_regular_contact[weight_sim_regular_contact<1]*0.5
weight_sim_transient_contact[weight_sim_transient_contact<0.8] <- weight_sim_transient_contact[weight_sim_transient_contact<0.8]*0.5

set_ids <- sample(c(1000:10000),size=10,replace=F)
for(seed_id in set_ids)
{
  if(!seed_id%in%alrdy_id){
    source('code/simulate_theta_0106.R') # here is the simulation drive to trigger the branching process, input the contact matrices in and get the epidemic size
    storage_epi <- as.data.frame(storage_epi)
    # write.csv(storage_epi,file=paste0('working/simulation results 0112/sim_epi_',seed_id,'.csv'),row.names = F) # original
    
    # to investigate the reduce 50% of contact outside the room unit. Edit on March 30
    write.csv(storage_epi,file=paste0('working/simulation results reduce50/sim_epi_',seed_id,'.csv'),row.names = F) # reduce 50% of the contact rate between blks and dormitories
  }

}
end_time=Sys.time()



