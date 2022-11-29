simulation_time = 1#commandArgs(trailingOnly=TRUE)[1]
total_sim_times = 10
a=Sys.time()
for(simulation_time in 101:(total_sim_times+100))
{
 source('code/effectiveness of new dorm.R') 
}
b=Sys.time()