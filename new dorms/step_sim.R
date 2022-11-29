
step_simulation <- function(sim_status, state_df, 
                            rec_df, iso_contact_df,
                            sim_params){
  # Increment time
  sim_status$t <- sim_status$t + sim_params$dt
  # Update timestamp of current infections
  state_df$days_infected <- state_df$days_infected + sim_params$dt
  iso_contact_df$days_isolated <- iso_contact_df$days_isolated + sim_params$dt
  # update the isolation contact IDs when needed.
  if(sim_params$iso)
    {
    # Determine whether the isolated contact will be released to the community (to susceptible population) or not
    iso_contact_df$isolated[which(iso_contact_df$days_isolated>=0 & iso_contact_df$days_isolated<=sim_params$iso_length)] <- 1 
    contact_be_isolated_row <- which(iso_contact_df$isolated==1)
    contact_not_isolated_row <- which(iso_contact_df$isolated==0)
      
    # Need to know that, isolation_contact in sim_status is dynamic!!
    if(length(contact_not_isolated_row)>0){
      release_contact <- NULL
      for(m in contact_not_isolated_row)
      {
        temp <- iso_contact_df$contact_id[[m]]
        release_contact <- unique(c(release_contact,temp))
      }
      if(length(which(sim_status$isolation_contact%in%release_contact))>0){
        sim_status$isolation_contact <- sim_status$isolation_contact[-which(sim_status$isolation_contact%in%release_contact)] # remove the ones who can be released
      }
    }
      
    if(length(contact_be_isolated_row)>0){ # add isolation contact IDs into sim_status
      sim_status$isolation_contact <- unique(c(sim_status$isolation_contact,unlist(do.call('c',lapply(contact_be_isolated_row,function(ii){
        iso_contact_df$contact_id[ii]
      })))))
    }
  }

  # only update isolation contact IDs in sim_status
  
  ## prep remove the individual who have alrdy been infected
  sim_status$alrdy_infected <- unique(c(sim_status$alrdy_infected,state_df$case_id))
  
  # Determine which cases will cause infections this step
  gen_sec_infect_out <- generate_secondary_infections(state_df, sim_params) #update the new infections, based on updated days_infect
  state_df <- gen_sec_infect_out$updated_state_df
  sec_infs_source <- gen_sec_infect_out$sec_infection_sources
  infectee_case_id <- gen_sec_infect_out$infectee_case_id
  n_sec_infections <- length(sec_infs_source)

  # If there are some new infections, create a state data and record data for the new cases
  # update sim_status for already infected IDs
  # update isolation information
  if (n_sec_infections > 0){
    # # Create new state_df for new cases
    new_case_id <- infectee_case_id
    if(any(new_case_id%in%sim_status$alrdy_infected)){
      rm <- which(new_case_id%in%sim_status$alrdy_infected)
      new_case_id <- new_case_id[-rm]
      sec_infs_source <- sec_infs_source[-rm]
    }
    sim_status$alrdy_infected <- unique(c(sim_status$alrdy_infected,new_case_id))
    sim_status$new_case_id <- new_case_id
    sim_status$sourced_id <- c(sim_status$sourced_id,sec_infs_source)
    
    
    sec_cases_state_df <- create_state_df(case_id = new_case_id,#n_sec_infections,
                                          sim_params=sim_params,
                                          sim_status=sim_status
    )
    
    if(length(which(duplicated(sec_cases_state_df$case_id)))>0)
    {
      rm <- which(duplicated(sec_cases_state_df$case_id))
      sec_cases_state_df <-  sec_cases_state_df[-rm,]
      sec_infs_source <- sec_infs_source[-rm]
    }
    if(nrow(sec_cases_state_df)>0)
    {
      sec_cases_iso_contact_df <- isolation_id_record(state_df = sec_cases_state_df,sim_params = sim_params)
    }else{
      sec_cases_iso_contact_df <- NULL
    }
    # Update last_case_id
    # sim_status$new_case_id <- sec_cases_state_df$case_id
    # Create new rec_df for new cases
    sec_cases_rec_df <- create_record_df(state_df = sec_cases_state_df, sim_status, infection_source=sec_infs_source)
    # TODO: Update rec_df of the origin cases (so no tracking of this yet)
    n_sec_infections <- nrow(sec_cases_state_df)
    
  }
  
  set.seed(seed_id)
  # First, identify and mark all cases where infection has ended, regardless of status
  cases_to_inact <- state_df$days_infected > state_df$infection_length
  n_cases_to_inact <- sum(cases_to_inact)
  if (n_cases_to_inact > 0){
    # Update state_df
    state_df$status[cases_to_inact] <- "inactive"
    # Update rec_df
    ids_to_inact <- state_df$case_id[cases_to_inact]
    rec_df$s_status[rec_df$case_id %in% ids_to_inact] <- "inactive"
  }
  
  # Identify all remaining cases that are eligible for advancement to next stage
  # Doing it now prevents checking advanced cases for another advancement
  cases_adv_inc <- (state_df$status=="incubation") &
    (
      (state_df$days_infected > state_df$incubation_length) |
        (state_df$days_infected > state_df$incubation_length + state_df$isolation_delay)
    )
  cases_adv_symp <- (state_df$status=="symptomatic") &
    (state_df$days_infected > (state_df$incubation_length + state_df$isolation_delay))
  cases_adv_asymp <- (state_df$status=="asymptomatic") &
    (state_df$days_infected > (state_df$incubation_length + state_df$isolation_delay))
  
  # Advance current infections from incubation to next stage
  # Goes to symptomatic or asymptomatic based on $is_symptomatic
  # Symptomatic cases that are traced goes right to isolation
  n_cases_adv_inc <- sum(cases_adv_inc)
  if (n_cases_adv_inc > 0){
    # Check whether advancing right to isolation (i.e. delay short enough)
    past_iso_delay <- state_df$days_infected > (state_df$incubation_length + state_df$isolation_delay)
    # Update state_df
    to_symp <- cases_adv_inc & state_df$is_symptomatic & !past_iso_delay
    to_iso  <- cases_adv_inc & past_iso_delay
    to_asymp <- cases_adv_inc & !state_df$is_symptomatic & !past_iso_delay
    state_df$status[to_symp] <- "symptomatic"
    state_df$status[to_iso] <- "isolated"
    state_df$status[to_asymp] <- "asymptomatic"
    # Update rec_df
    # New symptomatic cases
    ids_to_symp <- state_df$case_id[to_symp]
    # rec_df$t_symp[rec_df$case_id %in% ids_to_symp] <- sim_status$t
    rec_df$s_status[rec_df$case_id %in% ids_to_symp] <- "symptomatic"
    # New isolated cases
    ids_to_iso <- state_df$case_id[to_iso]
    # rec_df$t_symp[rec_df$case_id %in% ids_to_iso] <- sim_status$t
    # rec_df$t_iso[rec_df$case_id %in% ids_to_iso] <- sim_status$t
    rec_df$s_status[rec_df$case_id %in% ids_to_iso] <- "isolated"
    # New asymptomatic cases
    ids_to_asymp <- state_df$case_id[to_asymp]
    rec_df$s_status[rec_df$case_id %in% ids_to_asymp] <- "asymptomatic"
  }
  
  # Advance current infections from symptomatic to isolated after delay
  n_cases_adv_symp <- sum(cases_adv_symp)
  if (n_cases_adv_symp > 0){
    # Update state_df
    state_df$status[cases_adv_symp] <- "isolated"
    # Update rec_df
    ids_to_iso <- state_df$case_id[cases_adv_symp]
    # rec_df$t_iso[rec_df$case_id %in% ids_to_iso] <- sim_status$t
    rec_df$s_status[rec_df$case_id %in% ids_to_iso] <- "isolated"
  }
  
  # Advance current infections from asymptomatic to isolated after delay
  n_cases_adv_asymp <- sum(cases_adv_asymp)
  if (n_cases_adv_asymp > 0){
    # Update state_df
    state_df$status[cases_adv_asymp] <- "isolated"
    # Update rec_df
    ids_to_iso <- state_df$case_id[cases_adv_asymp]
    # rec_df$t_iso[rec_df$case_id %in% ids_to_iso] <- sim_status$t
    rec_df$s_status[rec_df$case_id %in% ids_to_iso] <- "isolated"
  }
  
  # record the contact individuals
  # iso_contact_df_new <- isolation_id_record(state_df,sim_params)
  
  # Remove isolated and inactive cases from state_df
  trunc_state_df <- state_df[(state_df$status != 'isolated') & (state_df$status != 'inactive'),]
  # trunc_state_df<-state_df[(state_df$status == 'isolated') | (state_df$status == 'inactive'),]
  
  # Add secondary infections to state_df and rec_df
  if (n_sec_infections > 0){
    out_state_df <- rbind(trunc_state_df, sec_cases_state_df)
    out_rec_df <- rbind(rec_df, sec_cases_rec_df)
    out_iso_df <- rbind(iso_contact_df,sec_cases_iso_contact_df)
  } else{
    out_state_df <- trunc_state_df
    out_rec_df <- rec_df
    out_iso_df <- iso_contact_df
  }
  
  # Add imported infections to output state_df and rec_df
  # if (n_imp_infections > 0){
  #   out_state_df <- rbind(out_state_df, imp_cases_state_df)
  #   # out_rec_df <- rbind(out_rec_df, imp_cases_rec_df)
  # }
  # 
  # Return updated inputs
  return(list(status=sim_status, state=out_state_df, 
              record=out_rec_df, iso_contact = out_iso_df,
              new_sec_cases=length(sec_infs_source)#,new_imp_cases=n_imp_infections
  ))
}