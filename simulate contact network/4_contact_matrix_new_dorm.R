# Add one more dorm
# import current contact matrix
load('C:/Users/ephsy/Documents/In Hand Project (PhD)/Dorm/spread_between_dorm/working/simulate contact matrix/simulate_contact.Rdata')
# sim_contact <- sim_contact+t(sim_contact)
# sim_regular_contact <- 1*as.matrix(sim_contact==1)
# sim_regular_contact<- as.data.frame(sim_regular_contact + t(sim_regular_contact))
# sim_transient_contact <- 1*as.matrix(sim_contact==2)
# sim_transient_contact <- as.data.frame(sim_transient_contact + t(sim_transient_contact))
set.seed(1)
contact_ind_index <- 1:nrow(sim_contact)
n <- 4800+3456
# add one column and one row for the contact matrix
# the blks in the new dorm is same as acacia, but only 16 person in one room
n_new_dorm = 4*5*8*16

# select ids who will be moved to the new dormitory
move_new_id <- sample(contact_ind_index,n_new_dorm) #random selected
# partial_regular_contact <- sim_regular_contact[-move_new_id,-move_new_id]
# partial_transient_contact <- sim_transient_contact[-move_new_id,-move_new_id]
partial_contact <- sim_contact[-move_new_id,-move_new_id]
contact_stay_index <- contact_ind_index[-move_new_id]
rm(sim_contact)
# the new id will update the simulate contact
new_sim_contact <- matrix(0,n,n)

new_sim_contact[1:length(contact_stay_index),1:length(contact_stay_index)] <- partial_contact
rm(partial_contact)

# in the following part, the connection in the new dorm will be simulated.
# the setting in the new dormitory will be similar to acacia with 6 blks, 5 floors/blk, 8 rooms/floor, 20 beds/room
seed_newdorm <- sample(100:1000,240)
di=3 # New dorm
para_arr <- read.csv('working/simulate contact matrix/Acacia_regular_Room_BetaMCMC.csv')
p_arr <- rbeta(8*5*4,mean(para_arr$alpha),mean(para_arr$beta))
para_alr <- read.csv('working/simulate contact matrix/Acacia_regular_Level_BetaMCMC.csv')
p_alr <- rbeta(5*4,mean(para_alr$alpha),mean(para_alr$beta))
para_alt <- read.csv('working/simulate contact matrix/Acacia_transient_Level_BetaMCMC.csv')
p_alt <- rbeta(5*4,mean(para_alt$alpha),mean(para_alt$beta))
# para_abr <- read.csv('working/simulate contact matrix/Acacia_regular_Blk_BetaMCMC.csv')
# p_abr <- rbeta(6,mean(para_abr$alpha),mean(para_abr$beta))
# para_abt <- read.csv('working/simulate contact matrix/Acacia_transient_Blk_BetaMCMC.csv')
# p_abt <- rbeta(6,mean(para_abt$alpha),mean(para_abt$beta))

sim_contact <- matrix(0,n_new_dorm,n_new_dorm)
h=0
hl=0
#acacia, there are 6 blks
n_bedder <- 16
n_room <- 8
n_level = 5
n_blk = 4
set.seed(9876)
# sample_dormcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level*n_blk,2)-choose(n_bedder*n_room*n_level,2)*n_blk,replace=T,
#                                   prob=c(1-0.0001-0.279,0.0001,0.279))
sample_dormcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level*n_blk,2)-choose(n_bedder*n_room*n_level,2)*n_blk,replace=T,
                                  prob=c(1-0.0002*4/53-0.279/20,0.0002*4/53,0.279/20))
d = 0
cat('dorms',di,"-rest blks",table(sample_dormcontact_info),'\n')

# To fill in the same room Rc and Tc
count=0 # to get the seed of the sample below
for(bi in 1:n_blk)# 6 blk in acacia
{ 
  count = count+1
  set.seed(seed_newdorm[count])
  # hb=hb+1
  # sample_blkcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level,2)-choose(n_bedder*n_room,2)*n_level,replace=T,prob=c(1-0.005-0.56,0.005,0.56))
  # sample_blkcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level,2)-choose(n_bedder*n_room,2)*n_level,replace=T,prob=c(1-p_abr[hb]-p_abt[hb],p_abr[hb],p_abt[hb]))
  sample_blkcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level,2)-choose(n_bedder*n_room,2)*n_level,
                                   replace=T,prob=c(1-0.04*14/53-0.56/10,0.04*14/53,0.56/10))
  
  c = 0
  cat('blk',bi,"-rest floors",table(sample_blkcontact_info),'\n')
  
  for(fi in 1:n_level) #each blk 5 floor-levels
  {
    sample_floorcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room,2)-choose(n_bedder,2)*n_room,replace=T,prob=c(1-0.01-0.58,0.01,0.58))
    hl=hl+1
    # sample_floorcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room,2)-choose(n_bedder,2)*n_room,replace=T,prob=c(1-p_alt[hl]-p_alr[hl],p_alr[hl],p_alt[hl]))
    # sample_floorcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room,2)-choose(n_bedder,2)*n_room,
    #                                    replace=T,prob=c(1-p_alt[hl]/53-p_alr[hl]*28/53,p_alr[hl]*28/53,p_alt[hl]/53))
    
    b=0
    cat('blk',bi,"-floor",fi,'rest rooms:',table(sample_floorcontact_info),'\n')
    
    for(ri in 1:n_room) # each floor 8 rooms
    {
      h=h+1
      # sample_floorcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room,2)-choose(n_bedder,2)*n_room,replace=T,prob=c(1-0.01-0.58,0.01,0.58))
      sample_roomcontact_info <- sample(c(1,2),choose(n_bedder,2),replace=T,prob=c(p_arr[h],1-p_arr[h]))
      
      # cat('blk',bi,"-floor",fi,"-room",ri,':',table(sample_roomcontact_info),'\n')
      a = 0
      for(i in (1+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder):(n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))
      {
        for(j in i:(n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))#(i+(fi-1)*8*20):(20+(fi-1)*8*20+(ri-1)*20))
        {
          if(i!=j)
          {
            a=a+1
            sim_contact[i,j] <- sample_roomcontact_info[a]
          }
        }
      }
      
      # To fill in the same floor but different room cells
      if(ri<8){
        for(restr_row in (1+(bi*fi-1)*n_bedder*n_room):(n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))#160:n_bedder*n_room on each floor
        {
          for(restr_col in min((n_bedder*n_room-n_bedder+1+(bi*fi-1)*n_bedder*n_room),(n_bedder+1+(ri-1)*n_bedder+(bi*fi-1)*n_bedder*n_room )):min((n_bedder*n_room+(bi*fi-1)*n_bedder*n_room ),(ri*n_bedder+n_bedder+(bi*fi-1)*n_bedder*n_room)))
          {
            b=b+1
            sim_contact[restr_row,restr_col] <- sample_floorcontact_info[b]
          }
        }
      }
    }
    # To fill in the same blk but different room cells
    if(fi<5)
    {
      for(restf_row in (1+n_bedder*n_room*n_level*(bi-1)):(n_bedder*n_room+n_bedder*n_room*(fi-1)+n_bedder*n_room*n_level*(bi-1)))
      {
        for(restf_col in (n_bedder*n_room+1+n_bedder*n_room*(fi-1)+n_bedder*n_room*n_level*(bi-1)):(n_bedder*n_room*2+n_bedder*n_room*(fi-1)+n_bedder*n_room*n_level*(bi-1)))
        {
          c=c+1
          sim_contact[restf_row,restf_col] <- sample_blkcontact_info[c]
        }
      }
    }
  }
  # To fill in the same dorm but different blk cells
  if(bi<4)
  {
    for(restb_row in (1+n_bedder*n_room*n_level*n_blk*(di-3)):(n_bedder*n_room*n_level+n_bedder*n_room*n_level*(bi-1)+n_bedder*n_room*n_level*n_blk*(di-3)))
    {
      for(restb_col in (1+n_bedder*n_room*n_level+n_bedder*n_room*n_level*(bi-1)+n_bedder*n_room*n_level*n_blk*(di-3)):(n_bedder*n_room*n_level*2+n_bedder*n_room*n_level*(bi-1)+n_bedder*n_room*n_level*n_blk*(di-3)))
      {
        d=d+1
        sim_contact[restb_row,restb_col] <- sample_dormcontact_info[d]
      }
    }
  }
}

new_sim_contact[(length(contact_stay_index)+1):n,(length(contact_stay_index)+1):n] <- sim_contact

sample_between_dormcontact_info <- sample(c(0,1,2),
                                          choose(n,2)-choose(length(contact_stay_index),2)-choose(n_new_dorm,2),replace=T,
                                          prob=c(1-0.0005*4/53-0.1398/20,0.0005*4/53,0.1398/20))
e = 0
for(rest_row in 1:length(contact_stay_index))
{
  for(rest_col in (length(contact_stay_index)+1):n)
  {
    e=e+1
    new_sim_contact[rest_row,rest_col] <- sample_between_dormcontact_info[e]
  }
}
e = 0
for(rest_row in (length(contact_stay_index)+1):n)
{
  for(rest_col in 1:length(contact_stay_index))
  {
    e=e+1
    new_sim_contact[rest_row,rest_col] <- sample_between_dormcontact_info[e]
  }
}




t <- 1*as.matrix(new_sim_contact==1)

a<-c();for(k in 1:n){temp<-sum(t[k,]);a<-c(a,temp)}
mean(a)

t <- 1*as.matrix(new_sim_contact==2)

a<-c();for(k in 1:n){temp<-sum(t[k,]);a<-c(a,temp)}
mean(a)

#save(new_sim_contact,file='working/simulate contact matrix/simulate_contact_3dorm.Rdata')
write.csv(new_sim_contact,file='working/simulate contact matrix/simulate_contact_3dorm.csv',row.names = F) # for Joel

new_sim_contact <- as.data.frame(new_sim_contact)
colnames(new_sim_contact) <- c(contact_stay_index,move_new_id)
rownames(new_sim_contact) <- c(contact_stay_index,move_new_id)
save(new_sim_contact,file='working/simulate contact matrix/simulate_contact_3dorm.Rdata')

# check the original dorms of individual who move to new dormitory
length(which(move_new_id%in%c(1:4800)))/length(move_new_id) 
#57.5 from acacia

