# generate the potential secondary cases of infectees
generate_secondary_infections <- function(state_df, sim_params){
  sec_infection_sources <- c()
  infectee_case_id <- c()
  if (nrow(state_df)>0){
    # TODO: Can we replace this for-loop?
    for (row in 1:nrow(state_df)){
      generation_int <- state_df$generation_intervals[[row]]
      infectee_id <- state_df$infectee_id[[row]]
      # Find which entries in generation_int are causing new secondary infections from this case
      sec_inf_ind <- state_df$days_infected[row] > generation_int
      n_sec_inf <- sum(sec_inf_ind)
      
      if (n_sec_inf < 1){
        next # no new infections caused by this case this step
      }
      
      # Add case_id to output list once for each new secondary infection case
      source_case_id <- state_df$case_id[row]
      sec_infection_sources <- c(sec_infection_sources,rep(source_case_id,n_sec_inf))
      infectee_case_id <- c(infectee_case_id,infectee_id[sec_inf_ind])
      # Decrease infection counter and remove from generation interval list
      state_df$n_sec_infects[row] <- state_df$n_sec_infects[row] - n_sec_inf
      state_df$generation_intervals[[row]] <- generation_int[!sec_inf_ind]
      state_df$infectee_id[[row]] <- infectee_id[!sec_inf_ind]
      
    }
    if(length(which(duplicated(state_df$case_id)))>0){
      rm <- which(duplicated(state_df$case_id))
      state_df <- state_df[-rm,]
      sec_infection_sources <- sec_infection_sources[-rm]
      infectee_case_id <- infectee_case_id[-rm]
    }
    
    return(list(updated_state_df=state_df, sec_infection_sources=sec_infection_sources, infectee_case_id = infectee_case_id))
  }
  else{
    return(list(updated_state_df=state_df, sec_infection_sources=NULL, infectee_case_id=NULL))
  }
}


