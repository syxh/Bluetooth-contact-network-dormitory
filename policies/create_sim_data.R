## to create isolation contact data

isolation_id_record <- function(state_df,sim_params)
{
  col_names <- c("days_isolated", "contact_id",'isolated') 
  n_cols <- length(col_names)
  # case_id <- unique(case_id)
  n_days <- nrow(state_df)#length(unique(state_df$days_infected))
  # Create data frame
  iso_contact_df <- data.frame(matrix(nrow=n_days,ncol=n_cols, dimnames=list(NULL,col_names)))
  iso_contact_df$days_isolated = 0-state_df$isolation_delay
  # iso_contact_df$isolated <- 1 # Yes, isolated
  # iso_contact_df$isolated[which(iso_contact_df$days_isolated<0|iso_contact_df$days_isolated>sim_params$iso_length)] <- 0 # No, not be isolated yet
  iso_contact_df$isolated<-0
  # [April] if infectee were isolated, then infectee_case_id will be 
  # for (row in 1:nrow(iso_contact_df)){
  if(!sim_params$iso)
  {
    iso_contact_df$contact_id <- NA
  }else{
    if(sim_params$iso_rc_only)
    { # only isolate regular contacts
      # look for the contacts in the isolation status
      # remove the ones alrdy been infected

      contact_id <- lapply(seq_along(state_df$case_id),function(i){
        tmp <- unlist(state_df$contact_id_regular[i])
        if(length(which(tmp%in%unlist(state_df$infectee_id[i])))>0){
          tmp <-  tmp[-which(tmp%in%unlist(state_df$infectee_id[i]))]
        }
        return(tmp)
      })
      iso_contact_df$contact_id <- contact_id
      #record time stamp
      # iso_contact_df$days_isolated <- iso_contact_df$days_isolated + state_df$days_infected#sim_params$dt
    }else{
      contact_id <- lapply(seq_along(state_df$case_id),function(i){
        tmp <- c(unlist(state_df$contact_id_regular[i]),unlist(state_df$contact_id_transient[i]))
        if(length(which(tmp%in%unlist(state_df$infectee_id[i])))>0){
          tmp <- tmp[-which(tmp%in%unlist(state_df$infectee_id[i]))]
        }
        return(tmp)
      })
      iso_contact_df$contact_id <- contact_id
      #record time stamp
      # iso_contact_df$days_isolated <- iso_contact_df$days_isolated + state_df$days_infected#sim_params$dt
    }
  }
  
  
  return(iso_contact_df)
}



### to track active cases in the simulation.
create_state_df <- function(case_id,#n_cases, 
                            sim_params, sim_status,
                            initialize=FALSE, 
                            import=FALSE,
                            primary_state_df=NULL, 
                            primary_case_ids=NULL){
  # List of columns in state df
  col_names <- c("case_id", "status", 
                 "is_traced","is_symptomatic", 
                 "days_infected", "incubation_length",
                 "isolation_delay", "infection_length", #'generation_interval',
                 "contact_number_regular","contact_number_transient",
                 "n_sec_infects") # infection_length here means valid period to generate secondary cases.
  n_cols <- length(col_names)
  # case_id <- unique(case_id)
  n_cases <- length(case_id)
  # Create data frame
  state_df <- data.frame(matrix(nrow=n_cases,ncol=n_cols, dimnames=list(NULL,col_names)))
  
  # If no cases being added, then return empty data frame
  # This should only happen when initializing
  if (n_cases==0){
    return(state_df)
  }
  
  # Fill in start values
  if (initialize){
    state_df$case_id <- case_id#1:n_cases # special case for starting out
  } else{
    state_df$case_id <- case_id#sim_status$new_case_id#c(case_id,sim_status$new_case_id)#sim_status$last_case_id + 1:n_cases # n_cases is the secondary cases
  }
  
  state_df$status <- rep("incubation", n_cases)
  if (initialize || import){
    state_df$is_traced <- rep(FALSE, n_cases)
  } else{
    state_df$is_traced <- draw_traced_status(n_cases,sim_params)
  }
  
  # if (!is.null(primary_state_df)){
  #   state_df$is_traced <- draw_traced_status(n = ,sim_params = sim_params)
  # } else {
  #   state_df$is_traced <- rep(FALSE, n_cases)
  # }
  
  state_df$is_symptomatic <- draw_symptomatic_status(n_cases,sim_params)
  state_df$days_infected <- rep(0, n_cases)
  state_df$incubation_length <- draw_incubation_interval(n_cases,sim_params)
  # state_df$generation_interval <- draw_generation_interval(n_cases,sim_params)
  state_df$infection_length <- 14
  
  state_df$isolation_delay <- draw_isolation_delay_period(state_df,sim_params
                                                          # primary_state_df=primary_state_df,
                                                          # primary_case_ids=primary_case_ids
  )
  
  # state_df$contact_number_regular <- find_contact_number(case_id, contact_ind_index=contact_ind_index, 
  #                                                        contact_matrix = sim_params$regular_contact_matrix,sim_status)
  # state_df$contact_number_transient <- find_contact_number(case_id, contact_ind_index=contact_ind_index, 
  #                                                          contact_matrix = sim_params$transit_contact_matrix,sim_status)
  # April edited
  state_df$contact_id_regular <- find_contact(case_id,contact_matrix =regular_contact,sim_status,sim_params)$id
  state_df$contact_number_regular <- find_contact(case_id, contact_matrix = regular_contact,sim_status,sim_params)$number
  
  state_df$contact_id_transient <- find_contact(case_id,contact_matrix = transit_contact, sim_status,sim_params)$id
  state_df$contact_number_transient <- find_contact(case_id, contact_matrix = transit_contact, sim_status,sim_params)$number
  
  
  
  sec_infect_out <- draw_sec_infects_df(state_df, sim_params, sim_status, import=F)
  # Set up the transmission tree, combine with the transmission tree, to define who infect whom
  
  state_df$n_sec_infects<- sec_infect_out$n
  # state_df$n_sec_infects_reg<- sec_infect_out$n_reg
  # state_df$n_sec_infects_trans<- sec_infect_out$n_trans
  state_df$generation_intervals<- sec_infect_out$generation
  state_df$non_generation <- sec_infect_out$non_infects
  state_df$infectee_id <- sec_infect_out$infectee
  state_df$non_infectee_id <- sec_infect_out$non_infectee
  
  return(state_df)
}



### to record active cases in the simulation.
create_record_df <- function(state_df, sim_status, initialize=FALSE, infection_source=NULL){
  # List of columns to record
  col_names <- c("case_id", "source", "is_traced", "is_symptomatic", 
                 "d_incub", "d_iso_delay","d_infection", "num_of_contacts", 
                 "n_sec_infects", "d_generation_ints",
                 "t_inf", #"t_symp", "t_iso", "t_inact",
                 "s_status", "non_infect_generations")
  n_cols <- length(col_names)
  # Get number of rows
  n_rows <- nrow(state_df)
  # Create data frame
  rec_df <- data.frame(matrix(nrow=n_rows,ncol=n_cols, dimnames=list(NULL,col_names)))
  
  # If no cases being added, then return empty data frame
  # This should only happen when initializing
  if (n_rows==0){
    return(rec_df)
  }
  
  # Populate data frame
  rec_df$case_id <- state_df$case_id
  if (initialize){
    rec_df$source <- rep("initial", n_rows) # initial cases
  } else{
    rec_df$source <- infection_source # vector of source case_ids or import source provided
  }
  rec_df$is_traced <- state_df$is_traced
  rec_df$is_symptomatic <- state_df$is_symptomatic
  rec_df$d_incub <- state_df$incubation_length
  rec_df$d_iso_delay <- state_df$isolation_delay
  rec_df$d_infection <- state_df$infection_length
  # rec_df$num_of_contacts <- state_df$contact_number
  rec_df$contact_number_regular <- state_df$contact_number_regular
  rec_df$contact_number_transient <- state_df$contact_number_transient
  rec_df$n_sec_infects <- state_df$n_sec_infects
  rec_df$d_generation_ints <- state_df$generation_intervals
  rec_df$infectee_case_ids <- state_df$infectee_id
  rec_df$t_inf <- rep(sim_status$t, n_rows)
  rec_df$s_status <- state_df$status
  rec_df$non_infect_generations <- state_df$non_infect_generations
  
  return(rec_df)
}

