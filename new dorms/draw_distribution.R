# Draw useful functions to create secondary infections when new cases are generated

# The function to draw incubation period of secondary infections, based on the gamma distribution
draw_incubation_interval <- function(n, sim_params){
  incubation_ints <- rgamma(n, shape=sim_params$incub_params$shape, scale  = sim_params$incub_param$scale)
  return(incubation_ints)
}

# The function to draw symptomatic status of secondary infections, based on the uniform distribution
draw_symptomatic_status <- function(n, sim_params){
  return(runif(n) < sim_params$p_sym)
}

# to get the 0-1 vector to show whether the secondary infections of new case are traced by BP device
draw_traced_status <- function(n, sim_params){
  if (sim_params$vary_trace){
    if(n<10){
      p_trace<-sim_params$p_trace_vary[1] # there is unaware to trace, p_trace_vary = 0.1
    }
    else if(n<20){
      p_trace<-sim_params$p_trace_vary[2] # some could not be traced, p_trace_vary =0.5
    } else{
      p_trace<-sim_params$p_trace_vary[3] # some could not be traced, p_trace_vary = 0.7 # 30% of the BP broken
    }
    return(runif(n) < p_trace)
  } else{ # otherwise, use fixed tracing value
    return(runif(n) < sim_params$p_trace)
  }
}

# generation interval: from infected time point to next infected time point
draw_generation_interval <- function(n, sim_params)
{
  generation_ints <- rgamma(n, shape=sim_params$generation_int_params$shape, 
                            rate=sim_params$generation_int_params$rate)
  return(generation_ints)
}

# The function to draw the delay days of all cases to be isolated
## The delay days to be isolated follows a weibull distribution
## If the infection is traced, then it delay to be isolated = onset-sampled delay to be isolated
## Need to note that days until isolated (d_until_iso) are calculated
draw_isolation_delay_period <- function(state_df, sim_params,
                                        primary_state_df=NULL,
                                        primary_case_ids=NULL){
  iso_delay_params <- sim_params$iso_delay_params
  # iso_delay_params$dist=='weibull'
  n <- nrow(state_df)
  # Draw from Weibull distribution initially
  iso_delay<-rweibull(
    n,
    shape=iso_delay_params$shape, # the isolation rate??
    scale=iso_delay_params$scale # mean/around mean of isolation lengths
  )
  # Modify actual delay time based on tracing and status of primary case
  # Only when primary_state_df provided (i.e. a secondary case)
  if (!is.null(primary_state_df)){
    # Get rows of primary_state_df matching primary cases
    infector_df <- subset(primary_state_df,case_id %in% primary_case_ids)
    # Calculate how many days between now and infector being isolated
    infector_df$d_until_iso <- infector_df$incubation_length +
      infector_df$isolation_delay -
      infector_df$days_infected
    # For each new case, determine the actual isolation delay
    iso_delay2 <- vapply(X=seq_along(iso_delay), FUN=function(ii){
      if (state_df[ii,]$is_traced){ # If traced, then check infector isolation time
        d_until_onset <- state_df$incubation_length[ii]
        d_until_pri_iso <- subset(infector_df,case_id==primary_case_ids[ii])$d_until_iso
        # If secondary case onset is *after* primary's isolation time, delay set to 0
        if (d_until_onset > d_until_pri_iso){
          return(0) # NB: return is to vapply
        } else{ # otherwise, delay is until primary case is isolated
          return(d_until_pri_iso-d_until_onset) # NB: return is to vapply
        }
      } else{ # If not traced, then use drawn iso_delay value
        return(iso_delay[ii]) # NB: return is to vapply
      }
    }, FUN.VALUE=999)
    # Update iso_delay vector to account for tracing
    iso_delay <- iso_delay2
  }
  
  return(iso_delay)
}

# Based on known contact matrix of new cases, we will find the potential/susceptible infections 
# (preparing for secondary infections) on it.
# Return a list, including the susceptible IDs of each new cases, 
# and number of susceptible cases of each new cases
find_contact<- function(case_id, contact_matrix, sim_status,sim_params) 
{
  
  # if(sim_params$iso){
  #   # if isolate yes, then remove all alrdy infected and isolated contacts from contact matrix (i.e. 
  #   # suspectible population), to make sure no infection will cause by them, and they will not in the
  #   # new_infections pool in the following procedure.
  #   # Need to note that, isolated contacts will be returned back to the suspectible population. 
  #   remove_id <- na.omit(c(sim_status$alrdy_infected,sim_status$isolation_contact))
  # }else{
  #   # if isolate is no, then only remove all infected individuals from suspectibel population.
  #   remove_id <- sim_status$alrdy_infected
  # }
  
  id <- lapply(seq_along(case_id),function(i) {
    # look for the index of the case id
    # new_dorm_ids are the real individual ids
    # look for the contact relationship from the matrix, should get the correct row/column number of the matrix
    index <- which(new_dorm_ids%in%case_id[i])
    allcontact_index <- which(contact_matrix[index,]==1)
    allcontact <- new_dorm_ids[allcontact_index]
    if(length(which(allcontact%in%occupied_id))>0){ # occupied_id is updating outside the function. pls find in the time-for-loop layer
      allcontact = allcontact[-which(allcontact%in%occupied_id)]
    }
    # print(allcontact)
    allcontact
  }
)
  number = do.call('c',lapply(id,length))
  return(list(id=id,number=number))
}

sanity_check_new_case <- function(new_case_id) # to make sure case will be infected from unique sources
{
  temp_new_id <- new_case_id[[1]]
  for(ii in seq_along(new_case_id))
  {
    if(ii>=2){
      rm <- which(unlist(new_case_id[[ii]]) %in% temp_new_id)
      if(length(rm)>0){
        new_case_id[[ii]] <- new_case_id[[ii]][-rm]
      }else{
        new_case_id[[ii]] <- new_case_id[[ii]]
      }
    }else{
      new_case_id[[ii]] <- new_case_id[[ii]]
    }
    temp_new_id <- unique(c(temp_new_id,unlist(new_case_id[[ii]])))
  }
  return(new_case_id)
}

draw_sec_infects_df <- function(state_df, sim_params, sim_status, import=FALSE){
  n_cases = nrow(state_df)
  # Following binomial distribution to determine
  # number of sec. infections of each infector 
  # and "generation interval"? of each sec. infection
  col_names <- c('n_infect', 'incubation_int', 'serial_int')
  n_cols <- length(col_names)
  infect_rate <- exp(-0.5*8256/(8256-length(sim_status$alrdy_infected)))*sim_params$p_infect

  if(n_cases>=1)
  {
    n_infect_reg <- state_df$contact_number_regular*infect_rate
  }
  set.seed(1111)
  new_case_id_reg <- list()
  for(i in seq_along(state_df$case_id)){
    row_index <- which(new_dorm_ids%in%state_df$case_id[i])
    col_index <- which(new_dorm_ids%in%state_df$contact_id_regular[[i]])
    if(n_infect_reg[i]<1|all(trans_tree_regular[row_index,col_index]<0.00000001)) {new_case_id_reg[[i]]=integer(0)
    }else{
      new_case_id_reg[[i]] <- sample(state_df$contact_id_regular[[i]],#contact_ind_index[-case_id],
                                     n_infect_reg[i],replace=FALSE, prob = trans_tree_regular[row_index,col_index])

    }
  }

  if(n_cases>=1)
    # if(!is.null(n_transient))
  {
    # n_infect_trans <- unlist(lapply(1:nrow(n_transient),
    #                                 function(i)
    #                                 {sum(n_transient[i,])*infect_rate*0.009}))
    n_infect_trans <- state_df$contact_number_transient*infect_rate*0.009
  }
  
  # if(length(sim_status$alrdy_infected)>1000)
  # n_infect_reg <- rnbinom(n_cases, mu = state_df$contact_number_regular*infect_rate, size = 0.7)
  # n_infect_trans <- rnbinom(n_cases, mu = state_df$contact_number_transient*infect_rate*0.01, size = 0.7)
  # number of secondary cases from infector i, here number of i =n_cases
  # to assign the infectees
  set.seed(1111)
   new_case_id_trans <- list()
   for(i in seq_along(state_df$case_id)){
     row_index <- which(new_dorm_ids%in%state_df$case_id[i])
     col_index <- which(new_dorm_ids%in%state_df$contact_id_transient[[i]])
    if(n_infect_trans[i]<1|all(trans_tree_transient[row_index,col_index]<0.00000001)) {new_case_id_trans[[i]]=integer(0)
    }else{
      new_case_id_trans[[i]] <- sample(state_df$contact_id_transient[[i]],#contact_ind_index[-case_id],
                                     n_infect_trans[i],replace=FALSE, prob = trans_tree_transient[row_index,col_index])
      # if(length(new_case_id_trans[[i]])>0){
      #   update_contact_ind_index <- update_contact_ind_index[-new_case_id_trans[[i]]]
      # }
    }
    
  }

  new_case_id_reg <- sanity_check_new_case(new_case_id_reg)
  new_case_id_trans <- sanity_check_new_case(new_case_id_trans)
  n_infect_reg <- unlist(lapply(new_case_id_reg,length))# number of infections from each infectee
  n_infect_trans <- unlist(lapply(new_case_id_trans,length)) # number of infections from each infectee
  n_infect <- n_infect_reg + n_infect_trans

  new_case_id <- lapply(1:length(state_df$case_id), function(i) {
    unique(c(new_case_id_reg[[i]],new_case_id_trans[[i]]))})
  
  #[April update] contact to be isolated
  # regular only
  
  # all contact to be isolated
  
  set.seed(sim_params$seed_id)
  # Determine incubation period and serial interval of each secondary infections from each infector source case
  incubation_period <- mapply(draw_incubation_interval, # FUN to be called on each X
                              n_infect, # X to loop over, n_infect is the secondary cases
                              MoreArgs=list(sim_params=sim_params), # additional required input for FUN
                              SIMPLIFY = FALSE) # force return as list
  
  # lapply(seq_along(incubation_length),function(i) unlist(incubation_length[i])+unlist(serial_int[i])) # not run on Nov 24, since both of them are list
  
  first_day_contagious <- rep(0,n_cases)
  # Infect-er
  # Switches for each contagious scenario
  is_T_and_S <- state_df$is_traced * state_df$is_symptomatic # traced and symptomatic
  is_T_and_nS <- state_df$is_traced * (1-state_df$is_symptomatic) # traced and not sympt
  is_nT_and_S <- (1-state_df$is_traced) * state_df$is_symptomatic # not traced and sympt
  is_nT_and_nS <- (1-state_df$is_traced) * (1-state_df$is_symptomatic) # not traced and not sympt
  # Length of contagious period for each scenario
  T_and_S_time <- state_df$incubation_length + state_df$isolation_delay # Contagious period. Delay for primary not yet isolated folded into isolation delay: which is contagiouslength
  T_and_nS_time <- state_df$infection_length # no isolation if no symptoms, infection_length is the constant = 14.
  nT_and_S_time <- state_df$incubation_length + state_df$isolation_delay # isolated some time after symptoms
  nT_and_nS_time <- state_df$infection_length  # no isolation if no symptoms
  last_day_contagious <- is_T_and_S * T_and_S_time   +
    is_T_and_nS * T_and_nS_time  +
    is_nT_and_S * nT_and_S_time  +
    is_nT_and_nS * nT_and_nS_time
  # Ensure no "last day contagious" is larger than infection_length days (would be removed as inactive)
  longer_than_contag_len <- last_day_contagious > state_df$infection_length
  last_day_contagious[longer_than_contag_len] <- state_df$infection_length[longer_than_contag_len]
  # Also, to avoid weird negative values in case somehow the isolated infector still caused
  # this case to happen, make all negative values be zero (would not allow any sec infects)
  neg_last_day <- last_day_contagious < 0
  last_day_contagious[neg_last_day] <- 0
  
  # if the first day of contagious is later than generation interval (from infector), then it would be kept in the infectee list. 
  # genertion interval is generated from gamma distribution.
  generation_int <- mapply(draw_generation_interval, # FUN to be called on each X
                           n_infect, # X to loop over
                           MoreArgs=list(sim_params=sim_params), # additional required input for FUN
                           SIMPLIFY = FALSE) # force return as list
  generation_int <- lapply(generation_int,sort)
  # first_day_contagious <- lapply(seq_along(generation_int),function(i) rep(first_day_contagious[i],length(generation_int[[i]])))
  # last_day_contagious <- lapply(seq_along(generation_int),function(i) rep(last_day_contagious[i],length(generation_int[[i]])))
  # Split the generation interval list.
  # Keep the valid infections for model, store the rest for record keeping
  generation_keep <- lapply(
    seq_along(generation_int),
    function(ii, generation, first_day, last_day,iso_rc_only){
      index_to_keep <- (generation[[ii]]  > first_day[ii]) &
        (generation[[ii]]< last_day[ii])
      return(generation[[ii]][index_to_keep])
    },
    generation = generation_int,
    first_day = first_day_contagious,
    last_day = last_day_contagious
    # iso_rc_only = iso_rc_only
  )
  
  generation_reject <- lapply(
    seq_along(generation_int),
    function(ii, generation, first_day, last_day,iso_rc_only){
      index_to_keep <- (generation[[ii]]  > first_day[ii]) &
        (generation[[ii]]< last_day[ii]) #& (!new_case_id[[ii]]%in%new_case_id_trans[[ii]])
      
      return(generation[[ii]][!index_to_keep])
    },
    generation=generation_int,
    first_day=first_day_contagious,
    last_day=last_day_contagious
    # iso_rc_only
  )
  
  infectee_keep <- lapply(
    seq_along(new_case_id),
    function(ii, infectee, generation, first_day, last_day,iso_rc_only){
      index_to_keep <- (generation[[ii]] > first_day[ii]) &
        (generation[[ii]]  < last_day[ii]) 
      return(infectee[[ii]][index_to_keep])
    },
    infectee = new_case_id,
    generation = generation_int,
    first_day = first_day_contagious,
    last_day = last_day_contagious
    #Jan 12
    # iso_rc_only=sim_params$iso_rc_only
  )
  infectee_reject <- lapply(
    seq_along(generation_int),
    function(ii, infectee,generation, first_day, last_day,iso_rc_only){
      
      index_to_keep <- (generation[[ii]] > first_day[ii]) &
        (generation[[ii]]  < last_day[ii])
      return(infectee[[ii]][!index_to_keep])
    },
    infectee = new_case_id,
    generation=generation_int,
    first_day=first_day_contagious,
    last_day=last_day_contagious
    #Jan 12
    # iso_rc_only=sim_params$iso_rc_only
  )
  
  
  # state_df$infectee_id <- infectee_keep
  
  # Get number of actual infections
  n_infect <- sapply(generation_keep,length)
  # Will probably move this outside of the if-block when other methods added
  return(list(n_reg=n_infect_reg,n_trans=n_infect_trans,n=n_infect,generation=generation_keep,
              infectee=infectee_keep,non_infects=generation_reject,non_infectee=infectee_reject))
  
}

