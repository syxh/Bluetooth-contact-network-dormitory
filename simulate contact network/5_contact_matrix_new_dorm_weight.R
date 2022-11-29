# create the 8256*8256 weight matrix of contact
ARRc <- 0.45/0.45
ARTc <- 0.4/0.45
AFRc <- 0.04/0.45
AFTc <- 0.13/0.45
ABRc <- 0.04/0.45
ABTc <- 0.11/0.45
ADRc <- 0.0003/0.45
ADTc <- 0.05/0.45

CRRc <- 0.68/0.68
CRTc <- 0.36/0.68
CFRc <- 0.05/0.68
CFTc <- 0.25/0.68
CBRc <- 0.03/0.68
CBTc <- 0.1/0.68
CDRc <- 0.003/0.68
CDTc <- 0.06/0.68

betweenRc <- 0.002/0.68
betweenTc <- 0.04/0.68

load('C:/Users/ephsy/Documents/In Hand Project (PhD)/Dorm/spread_between_dorm/working/simulate contact matrix/simulate_regular_contact_weight.Rdata')

load('C:/Users/ephsy/Documents/In Hand Project (PhD)/Dorm/spread_between_dorm/working/simulate contact matrix/simulate_transient_contact_weight.Rdata')

contact_ind_index <- 1:nrow(weight_sim_regular_contact)

# extract the individuals who will be kept in the existed dorms (need to consistent with the IDs in simulate_contact_3dorm_30%.Rdata)
set.seed(1)
n <- 4800+3456
n_new_dorm = 4*5*8*16
# select ids who will be moved to the new dormitory
move_new_id <- sample(contact_ind_index,n_new_dorm) #random selected
n_bedder <- 16
n_room <- 8
n_level = 5
n_blk = 4

stay_contact_index <- contact_ind_index[-move_new_id]

stay_contact_weight_network_reg <- weight_sim_regular_contact[-move_new_id,-move_new_id]; rm(weight_sim_regular_contact)
stay_contact_weight_network_trans <- weight_sim_transient_contact[-move_new_id,-move_new_id]; rm(weight_sim_transient_contact)

new_contact_weight_newtork_reg <- matrix(0,n,n)
new_contact_weight_newtork_trans <- matrix(0,n,n)
new_contact_weight_newtork_reg[1:length(stay_contact_index),1:length(stay_contact_index)] <- stay_contact_weight_network_reg
new_contact_weight_newtork_trans[1:length(stay_contact_index),1:length(stay_contact_index)] <- stay_contact_weight_network_trans

rm(stay_contact_weight_network_reg, stay_contact_weight_network_trans)


load('working/simulate contact matrix/simulate_contact_3dorm.Rdata')
sim_regular <- 1*as.matrix(new_sim_contact==1)
sim_transient <- 1*as.matrix(new_sim_contact==2)
sim_regular_contact <- sim_regular[(n-n_new_dorm+1):n,(n-n_new_dorm+1):n]
sim_transient_contact <- sim_transient[(n-n_new_dorm+1):n,(n-n_new_dorm+1):n]

temp_regular_contact <- matrix(0,n_new_dorm,n_new_dorm)
temp_transient_contact <- matrix(0,n_new_dorm,n_new_dorm)


# di=3
for(bi in 1:n_blk)
{
  for(fi in 1:n_level)
  {
    for(ri in 1:n_room)
      {
        for(i in (1+(bi*fi-1)*n_room*n_bedder +(ri-1)*n_bedder):(n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))
          for(j in i:(n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))
          {
            temp_regular_contact[i,j] <- ARRc*sim_regular_contact[i,j]
            temp_transient_contact[i,j] <- ARTc*sim_transient_contact[i,j]
          }
        if(ri<n_room){
          for(restr_row in (1+(bi*fi-1)*n_room*n_bedder):(n_bedder+(bi*fi-1)*n_room*n_bedder +(ri-1)*n_bedder))#160:n_bedder*n_room on each floor
          {
            for(restr_col in min((n_bedder*(n_room-1)+1+(bi*fi-1)*n_room*n_bedder),
                                 (n_bedder+1+(ri-1)*n_bedder+(bi*fi-1)*n_room*n_bedder)):min((n_room*n_bedder+(bi*fi-1)*n_room*n_bedder),
                                                                                 (ri*n_bedder+n_bedder+(bi*fi-1)*n_room*n_bedder)))
            {
              temp_regular_contact[restr_row,restr_col] <- AFRc*sim_regular_contact[restr_row,restr_col]
              temp_transient_contact[restr_row,restr_col] <- AFTc*sim_transient_contact[restr_row,restr_col]
            }
          }
        }
      }
      # To fill in the same blk but different room cells
      if(fi<n_level)
      {
        for(restf_row in (1+n_bedder*n_room*n_level*(bi-1)):(n_bedder*n_room+n_bedder*n_room*(fi-1)+n_bedder*n_room*n_level*(bi-1)))
        {
          for(restf_col in (n_bedder*n_room+1+n_bedder*n_room*(fi-1)+n_bedder*n_room*n_level*(bi-1)):
            (2*n_bedder*n_room+n_bedder*n_room*(fi-1)+n_bedder*n_room*n_level*(bi-1)))
          {
            temp_regular_contact[restf_row,restf_col] <- ABRc*sim_regular_contact[restf_row,restf_col]
            temp_transient_contact[restf_row,restf_col] <- ABTc*sim_transient_contact[restf_row,restf_col]
          }
        }
      }
    }
    # To fill in the same dorm but different blk cells
    if(bi<n_blk)
    {
      for(restb_row in 1:(n_bedder*n_room*n_level+n_bedder*n_room*n_level*(bi-1)))
      {
        for(restb_col in (n_bedder*n_room*n_level+1+n_bedder*n_room*n_level*(bi-1)):
            (2*n_bedder*n_room*n_level+n_bedder*n_room*n_level*(bi-1)))
        {
          temp_regular_contact[restb_row,restb_col] <- ADRc*sim_regular_contact[restb_row,restb_col]
          temp_transient_contact[restb_row,restb_col] <- ADTc*sim_transient_contact[restb_row,restb_col]
        }
      }
    }
}

new_contact_weight_newtork_reg[(n-n_new_dorm+1):n,(n-n_new_dorm+1):n] <- 
  temp_regular_contact; rm(temp_regular_contact)
new_contact_weight_newtork_trans[(n-n_new_dorm+1):n,(n-n_new_dorm+1):n] <- 
  temp_transient_contact; rm(temp_transient_contact)


#fill up between dorms
# upper part
for(rest_row in 1:(n-n_new_dorm))
{
  for(rest_col in 5697:n)
    {
    new_contact_weight_newtork_reg[rest_row,rest_col] <- betweenRc*sim_regular[rest_row,rest_col]
    new_contact_weight_newtork_trans[rest_row,rest_col] <- betweenTc*sim_transient[rest_row,rest_col]
    }
}
for(rest_row in 5697:n)
{
  for(rest_col in 1:(n-n_new_dorm))
  {
    new_contact_weight_newtork_reg[rest_row,rest_col] <- betweenRc*sim_regular[rest_row,rest_col]
    new_contact_weight_newtork_trans[rest_row,rest_col] <- betweenTc*sim_transient[rest_row,rest_col]
  }
}
colnames(new_contact_weight_newtork_reg) <- c(stay_contact_index,move_new_id)
rownames(new_contact_weight_newtork_reg) <- c(stay_contact_index,move_new_id)
colnames(new_contact_weight_newtork_trans) <- c(stay_contact_index,move_new_id)
rownames(new_contact_weight_newtork_trans) <- c(stay_contact_index,move_new_id)

save(new_contact_weight_newtork_reg,file='working/simulate contact matrix/simulate_regular_contact_weight_new_dorm.Rdata')
save(new_contact_weight_newtork_trans,file='working/simulate contact matrix/simulate_transient_contact_weight_new_dorm.Rdata')





  
