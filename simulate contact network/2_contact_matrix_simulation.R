# According to the pre-calculation of the regular and transient contact in the simualted dormitory
# named as "data/simulation_2whole_dormitory.csv"
# We will construct the contact matrix, as dimension: (4800+3456)^2
n <- 4800+3456
sim_contact <- matrix(0,n,n)
set.seed(1)
seed_acacia <- sample(100:1000,240)
di=1 # Acacia
para_arr <- read.csv('working/simulate contact matrix/Acacia_regular_Room_BetaMCMC.csv')
p_arr <- rbeta(8*5*6,mean(para_arr$alpha),mean(para_arr$beta))
para_alr <- read.csv('working/simulate contact matrix/Acacia_regular_Level_BetaMCMC.csv')
p_alr <- rbeta(5*6,mean(para_alr$alpha),mean(para_alr$beta))
para_alt <- read.csv('working/simulate contact matrix/Acacia_transient_Level_BetaMCMC.csv')
p_alt <- rbeta(5*6,mean(para_alt$alpha),mean(para_alt$beta))
# para_abr <- read.csv('working/simulate contact matrix/Acacia_regular_Blk_BetaMCMC.csv')
# p_abr <- rbeta(6,mean(para_abr$alpha),mean(para_abr$beta))
# para_abt <- read.csv('working/simulate contact matrix/Acacia_transient_Blk_BetaMCMC.csv')
# p_abt <- rbeta(6,mean(para_abt$alpha),mean(para_abt$beta))


h=0
hl=0
# hb=0
if(di==1)
{
  #acacia, there are 6 blks
  n_bedder <- 20
  n_room <- 8
  n_level = 5
  n_blk = 6
  set.seed(9876)
  # sample_dormcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level*n_blk,2)-choose(n_bedder*n_room*n_level,2)*n_blk,replace=T,
  #                                   prob=c(1-0.0001-0.279,0.0001,0.279))
  sample_dormcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level*n_blk,2)-choose(n_bedder*n_room*n_level,2)*n_blk,replace=T,
                                    prob=c(1-0.0002*4/53-0.279/20,0.0002*4/53,0.279/20))
  d = 0
  cat('dorms',di,"-rest blks",table(sample_dormcontact_info),'\n')
  
  # To fill in the same room Rc and Tc
  count=0 # to get the seed of the sample below
  for(bi in 1:6)# 6 blk in acacia
  { 
    count = count+1
    set.seed(seed_acacia[count])
    # hb=hb+1
    # sample_blkcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level,2)-choose(n_bedder*n_room,2)*n_level,replace=T,prob=c(1-0.005-0.56,0.005,0.56))
    # sample_blkcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level,2)-choose(n_bedder*n_room,2)*n_level,replace=T,prob=c(1-p_abr[hb]-p_abt[hb],p_abr[hb],p_abt[hb]))
    sample_blkcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level,2)-choose(n_bedder*n_room,2)*n_level,
                                     replace=T,prob=c(1-0.04*14/53-0.56/10,0.04*14/53,0.56/10))
    
    c = 0
    cat('blk',bi,"-rest floors",table(sample_blkcontact_info),'\n')
    
    for(fi in 1:5) #each blk 5 floor-levels
    {
      sample_floorcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room,2)-choose(n_bedder,2)*n_room,replace=T,prob=c(1-0.01-0.58,0.01,0.58))
      hl=hl+1
      # sample_floorcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room,2)-choose(n_bedder,2)*n_room,replace=T,prob=c(1-p_alt[hl]-p_alr[hl],p_alr[hl],p_alt[hl]))
      # sample_floorcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room,2)-choose(n_bedder,2)*n_room,
      #                                    replace=T,prob=c(1-p_alt[hl]/53-p_alr[hl]*28/53,p_alr[hl]*28/53,p_alt[hl]/53))
      
      b=0
      cat('blk',bi,"-floor",fi,'rest rooms:',table(sample_floorcontact_info),'\n')

      for(ri in 1:8) # each floor 8 rooms
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
          for(restr_row in (1+(bi*fi-1)*160):(n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))#160:n_bedder*n_room on each floor
          {
            for(restr_col in min((141+(bi*fi-1)*n_bedder*n_room),(n_bedder+1+(ri-1)*n_bedder+(bi*fi-1)*n_bedder*n_room )):min((160+(bi*fi-1)*n_bedder*n_room ),(ri*n_bedder+n_bedder+(bi*fi-1)*n_bedder*n_room)))
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
        for(restf_row in (1+800*(bi-1)):(160+160*(fi-1)+800*(bi-1)))
        {
          for(restf_col in (161+160*(fi-1)+800*(bi-1)):(320+160*(fi-1)+800*(bi-1)))
          {
            c=c+1
            sim_contact[restf_row,restf_col] <- sample_blkcontact_info[c]
          }
        }
      }
    }
    # To fill in the same dorm but different blk cells
    if(bi<6)
    {
      for(restb_row in (1+4800*(di-1)):(800+800*(bi-1)+4800*(di-1)))
      {
        for(restb_col in (801+800*(bi-1)+4800*(di-1)):(1600+800*(bi-1)+4800*(di-1)))
        {
          d=d+1
          sim_contact[restb_row,restb_col] <- sample_dormcontact_info[d]
        }
      }
    }
  }
}
# save(sim_contact,file='working/temp_sim_contact_Acacia.Rdata')
# table(sim_contact[1:20,1:20])
# table(sim_contact[1:160,1:160]) # 805-87*8 
# table(sim_contact[1:800,1:800]) #14321-805*5
# table(sim_contact[1:4800,1:4800]) #76671-14321-13559-11854-12584-11860-11007
# > table(sim_contact[801:1600,801:1600])
# 
# 0      1      2 
# 453752  13559 172689 
# > table(sim_contact[1601:2400,1601:2400])
# 
# 0      1      2 
# 470298  11854 157848 
# > table(sim_contact[2401:3200,2401:3200])
# 
# 0      1      2 
# 462170  12584 165246 
# > table(sim_contact[3201:4000,3201:4000])
# 
# 0      1      2 
# 470205  11860 157935 
# > table(sim_contact[4001:4800,4001:4800])
# 
# 0      1      2 
# 478438  11007 150555 

 
di=2 #caccia
set.seed(1)
para_crr <- read.csv('working/simulate contact matrix/Cassia_regular_Room_BetaMCMC.csv')
p_crr <- rbeta(12*6*3,mean(para_crr$alpha),mean(para_crr$beta))
para_clr <- read.csv('working/simulate contact matrix/Cassia_regular_Level_BetaMCMC.csv')
p_clr <- rbeta(6*3,mean(para_clr$alpha),mean(para_clr$beta))
para_clt <- read.csv('working/simulate contact matrix/Cassia_transient_Level_BetaMCMC.csv')
p_clt <- rbeta(6*3,mean(para_clt$alpha),mean(para_clt$beta))
# para_cbr <- read.csv('working/simulate contact matrix/Cassia_regular_Blk_BetaMCMC.csv')
# p_cbr <- rbeta(6,mean(para_cbr$alpha),mean(para_cbr$beta))
# para_cbt <- read.csv('working/simulate contact matrix/Cassia_transient_Blk_BetaMCMC.csv')
# p_cbt <- rbeta(6,mean(para_cbt$alpha),mean(para_cbt$beta))
h=0
hl=0
# hb=0
if(di==2)
{
  n_bedder <- 16
  n_room <- 12
  n_level = 6
  n_blk = 3
  di=2
  set.seed(9876)
  # sample_dormcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level*n_blk,2)-choose(n_bedder*n_room*n_level,2)*n_blk,replace=T,prob=c(1-0.0001-0.28,
  sample_dormcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level*n_blk,2)-choose(n_bedder*n_room*n_level,2)*n_blk,
                                    replace=T,prob=c(1-0.001*4/53-0.28/20, 0.001*4/53,0.28/20))
  d = 0
  cat('dorms',di,"-rest blks",table(sample_dormcontact_info),'\n')
  
  # To fill in the same room Rc and Tc
  count=0 # to get the seed of the sample below
  for(bi in 1:3)# 6 blk in acacia
  { 
    count = count+1
    set.seed(seed_acacia[count])
    # sample_blkcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level,2)-choose(n_bedder*n_room,2)*n_level,
    #                                  replace=T,prob=c(1-0.005-0.37,0.005,0.37))
    # sample_blkcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level,2)-choose(n_bedder*n_room,2)*n_level,replace=T,prob=c(1-p_cbr[hb]-p_cbt[hb],p_cbr[hb],p_cbt[hb]))
    sample_blkcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room*n_level,2)-choose(n_bedder*n_room,2)*n_level,
                                     replace=T,prob=c(1-0.01*14/53-0.37/10,0.01*14/53,0.37/10))
    c = 0
    cat('blk',bi,"-rest floors",table(sample_blkcontact_info),'\n')
    
    for(fi in 1:6) #each blk 5 floor-levels
    {
      # sample_floorcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room,2)-choose(n_bedder,2)*n_room,replace=T,prob=c(1-0.008-0.53,0.01,0.53))
      hl=hl+1
      sample_floorcontact_info <- sample(c(0,1,2),choose(n_bedder*n_room,2)-choose(n_bedder,2)*n_room,
                                         replace=T,prob=c(1-p_clt[hl]-p_clr[hl],p_clr[hl],p_clt[hl]))
      
      b=0
      cat('blk',bi,"-floor",fi,'rest rooms:',table(sample_floorcontact_info),'\n')
      
      for(ri in 1:n_room) # each floor 8 rooms
      {
        
        # sample_roomcontact_info <- sample(c(0,1,2),choose(n_bedder,2),replace=T,prob=c(1-0.62-0.34,0.62,0.34))
        h=h+1
        sample_roomcontact_info <- sample(c(1,2),choose(n_bedder,2),replace=T,prob=c(p_crr[h],1-p_crr[h]))
        # cat('blk',bi,"-floor",fi,"-room",ri,':',table(sample_roomcontact_info),'\n')
        a = 0
        for(i in (4800+1+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder):(4800+n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))
        {
          for(j in i:(4800+n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))#(i+(fi-1)*8*20):(20+(fi-1)*8*20+(ri-1)*20))
          {
            if(i!=j)
            {
              a=a+1
              sim_contact[i,j] <- sample_roomcontact_info[a]
            }
          }
        }
        
        # To fill in the same floor but different room contacts
        if(ri<12){
          for(restr_row in (4800+1+(bi*fi-1)*n_bedder*n_room):(4800+n_bedder+(bi*fi-1)*n_bedder*n_room +(ri-1)*n_bedder))#160:n_bedder*n_room on each floor
          {
            for(restr_col in min((4800+n_bedder*n_room-(n_bedder)+1+(bi*fi-1)*n_bedder*n_room),(4800+n_bedder+1+(ri-1)*n_bedder+(bi*fi-1)*n_bedder*n_room )):min((4800+n_bedder*n_room+(bi*fi-1)*n_bedder*n_room ),(4800+ri*n_bedder+n_bedder+(bi*fi-1)*n_bedder*n_room)))
            {
              b=b+1
              sim_contact[restr_row,restr_col] <- sample_floorcontact_info[b]
            }
          }
        }
      }
      # To fill in the same blk but different room contacts
      if(fi<6)
      {
        for(restf_row in (4800+1+n_bedder*n_room*n_level*(bi-1)):(4800+n_bedder*n_room+n_bedder*n_room*(fi-1)+1152*(bi-1)))
        {
          for(restf_col in (4800+n_bedder*n_room+1+n_bedder*n_room*(fi-1)+n_bedder*n_room*n_level*(bi-1)):(4800+2*n_bedder*n_room+n_bedder*n_room*(fi-1)+n_bedder*n_room*n_level*(bi-1)))
          {
            c=c+1
            sim_contact[restf_row,restf_col] <- sample_blkcontact_info[c]
          }
        }
      }
    }
    # To fill in the same dorm but different blk contacts
    if(bi<3)
    {
      for(restb_row in (1+4800*(di-1)):(n_bedder*n_room*n_level+n_bedder*n_room*n_level*(bi-1)+4800*(di-1)))
      {
        for(restb_col in (n_bedder*n_room*n_level+1+n_bedder*n_room*n_level*(bi-1)+4800*(di-1)):(2*n_bedder*n_room*n_level+n_bedder*n_room*n_level*(bi-1)+4800*(di-1)))
        {
          d=d+1
          sim_contact[restb_row,restb_col] <- sample_dormcontact_info[d]
        }
      }
    }
  }
}
# save(sim_contact,file='working/temp_sim_contact_Acacia_C_nointerationAC.Rdata')
# table(sim_contact[4801:(4800+16),4800:(4800+16)])
# table(sim_contact[4801:(4800+192),4801:(4800+192)]) # 1386-12*72
# table(sim_contact[4801:(4800+1152),4801:(4800+1152)]) #13931-1386*6
# table(sim_contact[4801:n,4800:n]) # 38396-8256-13930-11071
# > table(sim_contact[4801:(4800+1152),4801:(4800+1152)]) #13931-1386*6
# 
# 0       1       2 
# 1051850   13930  261324 
# > table(sim_contact[(4801+1152*2):(4800+1152*3),(4801+1152*2):(4800+1152*3)]) #13931-1386*6
# 
# 0       1       2 
# 1094918    8438  223748 
# To fill in the between dorm contact
# sample_between_dormcontact_info <- sample(c(0,1,2),
#                                           choose(n,2)-choose(4800,2)-choose(3456,2),replace=T,
#                                           prob=c(1-0.00005-0.1398,0.00005,0.1398))
sample_between_dormcontact_info <- sample(c(0,1,2),
                                          choose(n,2)-choose(4800,2)-choose(3456,2),replace=T,
                                          prob=c(1-0.0005*4/53-0.1398/20,0.0005*4/53,0.1398/20))
e = 0
for(rest_row in 1:4800)
{
  for(rest_col in 4801:n)
  {
    e=e+1
    sim_contact[rest_row,rest_col] <- sample_between_dormcontact_info[e]
  }
}

t <- 1*as.matrix(sim_contact==1)

a<-c();for(k in 1:n){temp<-sum(t[k,]);a<-c(a,temp)}
mean(a)

t <- 1*as.matrix(sim_contact==2)

a<-c();for(k in 1:n){temp<-sum(t[k,]);a<-c(a,temp)}
mean(a)
save(sim_contact,file='working/simulate contact matrix/simulate_contact.Rdata')
write.csv(sim_contact,file='working/simulate contact matrix/simulate_contact.csv',row.names = F) # for Joel
