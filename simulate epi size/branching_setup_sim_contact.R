
####################### Simulation setup #################
## In this session, we include the parameters setting, and the initial individual cases data
initialize_sim_params <- function(p_infect,
                                  # infect_rate,
                                  generation_int_params = list(shape=3.5,rate=0.8),  #refer: doi: 10.2807/1560-7917.ES.2020.25.17.2000257
                                  incub_params = list(shape=9.65,scale=0.54),
                                  p_sym = 0.6,
                                  iso_delay_params = list(shape=5,scale=3),
                                  trans_tree_reg,
                                  trans_tree_trans,
                                  regular_contact_matrix, #related to mcmc_contact_rate_recults
                                  transit_contact_matrix,
                                  vary_trace=T, p_trace_vary = list(0.1,0.5,0.7),
                                  sec_infect_params = list(disp=0.7), #nbinom parameter, mean was defined by contact number, disp only
                                  # refer of dispersion value: doi: 0.1186/s12889-020-09624-2
                                  dt, # Advance one day at a time
                                  seed_id
                                  ){
  sim_params <- list(
    # infect_rate = infect_rate,
    p_infect = p_infect, 
    generation_int_params = generation_int_params, # parameters for distribution of serial interval
    incub_params = incub_params,           # parameters for distribution of incubation length
    p_sym = p_sym,
    iso_delay_params = iso_delay_params,
    trans_tree_reg =  trans_tree_reg,
    trans_tree_trans =  trans_tree_trans,
    regular_contact_matrix = regular_contact_matrix, # parameters for calculation of mean number of contact in Poisson distribution
    transit_contact_matrix = transit_contact_matrix, # parameters for calculation of mean number of contact in Poisson distribution
    vary_trace = vary_trace,
    p_trace_vary = p_trace_vary,
    sec_infect_params = sec_infect_params,
    dt = dt,
    seed_id=seed_id

  )
  return(sim_params)
}

initialize_case <- function(n_cases,same_room=F)
{
  if(same_room)
  {
    case_id <- sample(1:16,n_cases,replace = F)
  }else{
    case_id <- sample(1:8256,n_cases,replace=F)
  }
  case_id
}

initialize_sim_status <- function(start_time, start_case_ids
                                  # regular_contact = sim_regular_contact,
                                  # transient_contact = sim_transient_contact
                                  ){
  sim_status <- list(
    t = start_time,
    sourced_id = start_case_ids,
    alrdy_infected = start_case_ids,
    new_case_id = NULL
    # regular_contact = regular_contact,
    # transient_contact = transient_contact
  )
  return(sim_status)
}






