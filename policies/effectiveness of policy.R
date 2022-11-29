# Aim of this step:
# 1) we need to do the validation
# 2) simulation under different scenarios



# functions sourcing ------------------------------------------------------
# source('code/branching_function_sim_contact_0112.r')
source('code/draw_distribution.R')
source('code/generate_new_infections.R')
source('code/create_sim_data.R')
source('code/step_sim.R')
source('code/simulation_set_up.R')
library(Matrix)
simulation_time = 1#commandArgs(trailingOnly=TRUE)[1]
# settings ----------------------------------------------------------------
# load the simulated contact matrix
regular_contact <- readMM("data/sim_regular_contact.mtx")
transit_contact <- readMM("data/sim_transient_contact.mtx")


trans_tree_regular <- readMM("data/weight_sim_regular_contact.mtx")
trans_tree_transient<- readMM("data/weight_sim_transient_contact.mtx")

# simulation --------------------------------------------------------------

# one simulation without any further intervention policy, the running time is 2.6mins
# try 1000 simulation to get the new cases
scenario <- read.csv('data/different intervention upd.csv')
simulation_length = 120#nrow(dat)
total_sim_times=100
seed_id = 9876
a = Sys.time()
for(scenario_id in 5)
{
  theta <- c(0.4,0.87)#c(0.59,0.91)#scenario$theta[scenario_id]
  mean_iso_delay_projection = scenario$mean_iso_delay_projection[scenario_id]
  iso_rc_only =  scenario$only_rc_iso[scenario_id]
  iso_length = scenario$iso_length[scenario_id]
  iso = scenario$isolation[scenario_id]
  set.seed(seed_id+as.numeric(simulation_time))
  ####  no intensive intervention policy within the dormitory setting, they still go to work
  new_cases <- matrix(NA,ncol = simulation_length+2, nrow=1)
  
   for(theta_id in 2)
   {
    for(simulation_time in 10)#:total_sim_times)
    {
      cat("simulation order:",simulation_time,"\n")
      cat("theta value:",theta[theta_id],"\n")
      nsims_extinct<-0 # all cases went extinct
      nsims_epidemic<-0  # more than [max_cases] total cases
      sim_params <- initialize_sim_params(p_infect = theta[theta_id],#infect_rate, 
                                          incub_params = list(shape=(5.5*5.5)/2.8,scale=2.8/5.5),#list(shape=(6.4*6.4)/2.3,scale=2.3/6.4),
                                          iso_delay_params=list(shape=20,scale=mean_iso_delay_projection), 
                                          dt=1,
                                          iso = iso,
                                          seed_id = as.numeric(simulation_time),
                                          iso_rc_only = iso_rc_only,
                                          iso_length = iso_length)
      initial_cases <- 20
      start_case_ids <- initialize_case(initial_cases)
      sim_status <- initialize_sim_status(start_time = 0,start_case_ids = start_case_ids)
      occupied_id <- sim_status$alrdy_infected
      state_df   <- create_state_df(case_id=start_case_ids ,sim_params,sim_status,initialize=T)
      record_df <- create_record_df(state_df, sim_status, initialize = T)
      # add the isolation contact data here
      isolation_df <- isolation_id_record(state_df,sim_params)
      
      timemax=simulation_length#nrow(dat)#ceiling(12*7/dt) # 12 weeks after initial cases (Fig 3 caption)
      n_new_cases=rep(0,timemax) # number of new cases per step
      n_cum=rep(0,timemax) # number of new cases per step
      
      early_exit <- FALSE
      for (ii in 1:timemax){
        set.seed(seed_id)
        # occupied ids: 2 parts. 
        ## - already infected, who will always be zero in contact matrix and transmission matrix
        ## - isolation ones, who will be returned back to contact matrix and transmission matrix
        if(length(sim_status$isolation_contact)>0|length(sim_status$alrdy_infected)>0){
          occupied_id <- unique(c(sim_status$alrdy_infected,sim_status$isolation_contact))
        }else{
          occupied_id <- NULL
        }
        occupied_id <- occupied_id
        
        out <- step_simulation(sim_status, state_df, record_df, iso_contact_df = isolation_df, sim_params)
        sim_status <- out$status
        state_df <- out$state
        record_df <- out$record
        isolation_df <- out$iso_contact
        n_new_cases[ii] <- out$new_sec_cases
        if (nrow(state_df)==0){
          nsims_extinct <- nsims_extinct + 1
          early_exit <- TRUE
          break
        } 
        if(ii%%30==0) print(paste0('simulation days:',ii))
      }
      
      new_cases[1,2:(simulation_length+1)] <- n_new_cases
      new_cases[1,1] <- theta[theta_id]
      new_cases[1,(simulation_length+2)] <- nrow(record_df)
      new_cases <- as.data.frame(new_cases)
      write.csv(n_new_cases,file = paste0('newworking/Scenario_',scenario_id,'_theta_',theta[theta_id], '+', as.numeric(simulation_time),'.csv'),row.names = F)
      write.csv(record_df$case_id,file = paste0('newworking/caseid_',scenario_id,'_theta_',theta[theta_id], '+', as.numeric(simulation_time),'.csv'),row.names = F)
      write.csv(record_df$contact_number_regular,file = paste0('newworking/regularCNum_',scenario_id,'_theta_',theta[theta_id], '+', as.numeric(simulation_time),'.csv'),row.names = F)
      write.csv(record_df$contact_number_transient,file = paste0('newworking/transientCNum_',scenario_id,'_theta_',theta[theta_id], '+', as.numeric(simulation_time),'.csv'),row.names = F)
      
      }

    
  }
}

b = Sys.time()


#paper number checking
# dat31 <- read.csv(paste0('working/effectIntervention/Scenario_',3,'_theta_',theta[1],'.csv'))
# summary(apply(dat31,MARGIN = 1,sum))/8256
# dat32 <- read.csv(paste0('working/effectIntervention/Scenario_',3,'_theta_',theta[2],'.csv'))
# summary(apply(dat32,MARGIN = 1,sum))/8256
# 
# 
# dat21 <- read.csv(paste0('working/effectIntervention/Scenario_',2,'_theta_',theta[1],'.csv'))
# summary(apply(dat21,MARGIN = 1,sum))/8256
# dat22 <- read.csv(paste0('working/effectIntervention/Scenario_',2,'_theta_',theta[2],'.csv'))
# summary(apply(dat22,MARGIN = 1,sum))/8256
# 
# 
# dat11 <- read.csv(paste0('working/effectIntervention/Scenario_',1,'_theta_',theta[1],'.csv'))
# summary(apply(dat11,MARGIN = 1,sum))/8256
# dat12 <- read.csv(paste0('working/effectIntervention/Scenario_',1,'_theta_',theta[2],'.csv'))
# summary(apply(dat12,MARGIN = 1,sum))/8256





