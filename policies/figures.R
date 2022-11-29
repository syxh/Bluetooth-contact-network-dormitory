
library(grid)
xaxes_simulate_time=function(xlm)
{
  
  # Mo-Fr
  xtk1 = 60*seq(xlm[1],xlm[2]/30,1)
  xtk2 = xtk1+30
  
  i=which((xtk1>=xlm[1])&(xtk2<=xlm[2]))
  xtk1touse=xtk1[i]
  xtk2touse=xtk2[i]
  i=which((xtk1<xlm[2])&(xtk2>=xlm[2]))
  if(length(i)>0)
  {
    xtk1touse=c(xtk1touse,xtk1[i])
    xtk2touse=c(xtk2touse,xlm[2])
  }
  # for(k in seq(xtk1touse)) print(paste0(xtk1touse[k],',',xtk2touse[k]))
  
  for(k in seq(xtk1touse))grid.polygon(unit(c(xtk1touse[k],xtk1touse[k],xtk2touse[k],xtk2touse[k]),'native'),
                                       unit(c(0,1,1,0),'npc'),
                                       gp=gpar(col=NA,fill='#f0f7fe'))
  
  
  xtk = seq(0,xlm[2],30)
  xtk = xtk[xtk<=xlm[2]]
  grid.xaxis(at=xtk,label=FALSE,main = TRUE)
  grid.text(xtk,x=unit(xtk,'native'),y=unit(-1,'lines'))
  # grid.text('M',x=unit(0.0+as.Date('2020-05-15')-as.Date('2020-01-01'),'native'),y=unit(-1,'lines'))

  
}


yaxis = function(ylm,log=FALSE)
{
  if(!log)
  {
    ytk = pretty(ylm)
    ylb = as.character(ytk*100)
    for(ij in seq(ylb)){if(ytk[ij]>999)ylb[ij]=paste0(ytk[ij]/1000,'k')}
    ij=which(ytk<ylm[2])
    ytk=ytk[ij]
    ylb=ylb[ij]
    grid.yaxis(ytk,ylb)
  }
  if(log)
  {
    ytk = 0:10#pretty(ylm)
    ylb = as.character(10^ytk)
    for(ij in seq(ylb)){if(ytk[ij]>3)ylb[ij]=paste0(10^ytk[ij]/1000,'k')}
    ylb[ytk<0]='0'
    ij=which(ytk<=ylm[2])
    ytk=ytk[ij]
    ylb=ylb[ij]
    ij=which(ytk>=ylm[1])
    ytk=ytk[ij]
    ylb=ylb[ij]
    
    grid.yaxis(ytk,ylb)
    z=1:9
    Z = log10(c(z,10*z,100*z,1000*z,10000*z,100000*z))
    Z = Z[Z<ylm[2]]
    if(length(Z)>0)grid.yaxis(at=Z,label=FALSE,#c(" ",seq(2,9,1)," ",seq(20,90,10)),
                              gp=gpar(cex=0.5))
    for(i in seq_along(Z)){
      grid.lines(x = unit(c(0,1),'npc'),
                 y = unit(Z[i],'native'),
                 gp=gpar(col='#82bcf7',lty='dashed',lwd=1))
    }
  }
  
}


# trunc_daily_dat <- read.csv('data/dorm_daily_data.csv')
# trunc_daily_dat <- trunc_daily_dat[which(trunc_daily_dat$date>='3/20/2020'),]
# model_data <-  trunc_daily_dat[which(trunc_daily_dat$doy<=113),]

paneller_CI = function(file,dormID,scenario_id)
{
  if(dormID%in%c(1,'Acacia','acacia'))
  {
    output <- read.csv(paste0('output/data/Scenario_',scenario_id,'_acacia.csv'))
    ylabel = "Proportion of Infected under Mild scenario"
  }else{
    output <- read.csv(paste0('output/data/Scenario_',scenario_id,'_cassia.csv'))
    ylabel = "Proportion of Infected under Severe scenario"
  }
  t  = c(1:nrow(output))
  oC = output$central/8256
  oL= output$lower/8256
  oU = output$upper/8256
  png(file,height=10,width=10,units='cm',res=900,pointsize=12)
    
  xlm = c(0,range(t)[2])
  logy=FALSE
  top=TRUE
  ylabel_gap=-2.5
  
  ylm=c(0,max(oU)*1.05)
  grid.newpage()
  pushViewport(plotViewport(c(3,3,1,1),
                            xscale=xlm,yscale=ylm))
  xaxes_simulate_time(xlm)
  
  mastercol = '#ed2939'
  # if(version==2)mastercol= '#2171b5'#'#003171'#
  # if(version==3)mastercol= '#607c3c'
  graycol = '#969696'
  
  q=as.vector(col2rgb(mastercol))/255
  colbg = rgb(q[1],q[2],q[3],0.5)
  
  # q1=as.vector(col2rgb(graycol))/255
  # colbg1 = rgb(q1[1],q1[2],q1[3],0.5)
  grid.polygon(c(c(1:range(t)[2]),rev(c(1:range(t)[2]))),
               c(oL,rev(oU)),
               default.units = 'native',
               gp=gpar(col=NA,fill=colbg,alpha=0.8))
  
  grid.lines(c(1:range(t)[2]),oC,default.units = 'native',
             gp=gpar(col=mastercol,lwd=2))
  

  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  if(!logy)yaxis(ylm)
  if(logy)yaxis(ylm,log = TRUE)
  grid.text(ylabel,x=unit(0,'npc')+unit(ylabel_gap,'lines'),y=unit(0.5,'npc')-unit(0.5,'lines'),rot=90)
  grid.text('Time (day)', x= unit(0.5,'npc'), y=unit(0,'npc')+unit(-2,'lines'))
  popViewport()
  dev.off()
}

paneller = function(file,dormID,scenario_id)
{
  if(dormID%in%c(1,'Acacia','acacia'))
  {
    output <- read.csv(paste0('working/combine_Scenario_',scenario_id,'_theta_acacia.csv'))
    ci <- read.csv(paste0('output/data/Scenario_',scenario_id,'_acacia.csv'))
    ylabel = "Proportion of Infected (%)"
    leg_label = "Mild Scenario"
  }else{
    output <- read.csv(paste0('working/combine_Scenario_',scenario_id,'_theta_cassia.csv'))
    ci <- read.csv(paste0('output/data/Scenario_',scenario_id,'_cassia.csv'))
    ylabel = "Proportion of Infected (%)"
    leg_label = "Severe Scenario"
    
  }
  t  = c(1:nrow(output))
  oC = log10(100*ci$central/8256)
  oC[which(oC<0)]=0
  oL= log10(100*ci$lower/8256)
  oL[which(oL<0)]=0
  oU = log10(100*ci$upper/8256)
  oU[which(oU<0)]=0
  
  png(file,height=10,width=10,units='cm',res=900,pointsize=12)
  
  xlm = c(0,range(t)[2])
  logy=TRUE
  top=TRUE
  ylabel_gap=-3
  xlabel_gap=-2
  # ylabel = "Proportion of Infected (%)"
  xlabel = "Time (Day)"
  
  ylm=c(0,log10(100))#max(oU)*1.25)#c(0,1)
  grid.newpage()
  pushViewport(plotViewport(c(3,3.5,1,1),
                            xscale=xlm*1.1,yscale=ylm))
  xaxes_simulate_time(xlm)
  if(dormID==1)
  {
    mastercol= 'darkorange'#'2171b5'
    light_mastercol = '#ffc988'#'#71c1ff'
  }
  if(dormID==2)
  {
    mastercol = 'steelblue' #'#ed2939'
    # if(version==2)mastercol= '#2171b5'#'#003171'#
    # if(version==3)mastercol= '#607c3c'
    light_mastercol = '#B0C4DE'#'#ffc0b8'
  }

  yaxis(ylm,log = T)
  
  for(i in 2:ncol(output))
  {
    tempy <- log10(100*output[,i]/8256)
    tempy[which(tempy<0)]=0
    grid.lines(c(1:range(t)[2]),tempy,default.units = 'native',
               gp=gpar(col=light_mastercol,lwd=1))
  }
  
  grid.lines(c(1:range(t)[2]),oC,default.units = 'native',
             gp=gpar(col=mastercol,lwd=2))
  
  grid.points(x=unit(120,'native')+unit(0.7,'lines'),y=oC[120],pch=20,
              gp=gpar(col=mastercol))
  grid.lines(x=unit(120,'native')+unit(0.7,'lines'),
             y=unit(c(oL[120],oU[120]),'native'),gp=gpar(col=mastercol,lwd=2))
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  # grid.text(label=paste0("Scenario ",scenario_id), x = unit(1,'npc')-unit(2.5,'lines'),
  #           y = unit(1,'lines'), gp=gpar(fontsize=12))
  # if(!logy)yaxis(ylm)
  # if(logy)yaxis(ylm,log = TRUE)
  # grid.yaxis(at=unit(c(0,pretty(ylm[2])[2]),'native'),label = c(0,pretty(ylm[2])[2]*100))
  # if(scenario_id %in% c(3,4)){
  #   # grid.yaxis(at=seq(0,0.2,0.01),label=seq(0,20,1))
  #   yaxis(ylm)
  #   
  # }else{
    # grid.yaxis(at=seq(0,1,0.1),label=seq(0,100,10))
  # }
  grid.text(ylabel,x=unit(0,'npc')+unit(ylabel_gap,'lines'),y=unit(0.5,'npc')-unit(0.5,'lines'),rot=90)
  grid.text('Time (day)', x= unit(0.5,'npc'), y=unit(0,'npc')+unit(-2,'lines'))
  # grid.text(leg_label, x= unit(1,'npc')+unit(-2.5,'lines'), y=unit(1,'npc')+unit(-1,'lines'))
  
  # grid.text("CI",x=unit(120,'native')+unit(0.7,'lines'),y=min(unit(oU[120],'native')+unit(3,'lines'),unit(0.98,'npc')),gp=gpar(fontface='bold'))
  popViewport()
  dev.off()
}


paneller(file='output/plot/acacia_newdorm_s1_ci.png',dormID = 1,scenario_id = 1)
paneller(file='output/plot/cassia_newdorm_s1_ci.png',dormID = 2,scenario_id = 1)

paneller(file='output/plot/acacia_newdorm_s2_ci.png',dormID = 1,scenario_id = 2)
paneller(file='output/plot/cassia_newdorm_s2_ci.png',dormID = 2,scenario_id = 2)

paneller(file='output/plot/acacia_newdorm_s3_ci.png',dormID = 1,scenario_id = 3)
paneller(file='output/plot/cassia_newdorm_s3_ci.png',dormID = 2,scenario_id = 3)

paneller(file='output/plot/acacia_newdorm_s4_ci.png',dormID = 1,scenario_id = 4)
paneller(file='output/plot/cassia_newdorm_s4_ci.png',dormID = 2,scenario_id = 4)

paneller(file='output/plot/acacia_newdorm_s5_ci.png',dormID = 1,scenario_id = 5)
paneller(file='output/plot/cassia_newdorm_s5_ci.png',dormID = 2,scenario_id = 5)

paneller(file='output/plot/acacia_newdorm_s6_ci.png',dormID = 1,scenario_id = 6)
paneller(file='output/plot/cassia_newdorm_s6_ci.png',dormID = 2,scenario_id = 6)

paneller(file='output/plot/acacia_newdorm_s7_ci.png',dormID = 1,scenario_id = 7)
paneller(file='output/plot/cassia_newdorm_s7_ci.png',dormID = 2,scenario_id = 7)



# check the number in paper
dormID=1
for(scenario_id in 1:7)
{
  if(dormID%in%c(1,'Acacia','acacia'))
  {
    output <- read.csv(paste0('working/combine_Scenario_',scenario_id,'_theta_acacia.csv'))
    ci <- read.csv(paste0('output/data/Scenario_',scenario_id,'_acacia.csv'))
  }else{
    output <- read.csv(paste0('working/combine_Scenario_',scenario_id,'_theta_cassia.csv'))
    ci <- read.csv(paste0('output/data/Scenario_',scenario_id,'_cassia.csv'))
  }
  oC = ci$central/8256
  oL= ci$lower/8256
  oU = ci$upper/8256
  cat("Scenario",scenario_id,'\n')
  cat('Proportion of infected in 4 months:',round(oC[120],3),'(95% CI:',round(oL[120],3),'-',round(oU[120],3),')\n')
}


dormID=2
for(scenario_id in 1:7)
{
  if(dormID%in%c(1,'Acacia','acacia'))
  {
    output <- read.csv(paste0('working/combine_Scenario_',scenario_id,'_theta_acacia.csv'))
    ci <- read.csv(paste0('output/data/Scenario_',scenario_id,'_acacia.csv'))
  }else{
    output <- read.csv(paste0('working/combine_Scenario_',scenario_id,'_theta_cassia.csv'))
    ci <- read.csv(paste0('output/data/Scenario_',scenario_id,'_cassia.csv'))
  }
  oC = ci$central/8256
  oL= ci$lower/8256
  oU = ci$upper/8256
  cat("Scenario",scenario_id,'\n')
  cat('Proportion of infected in 4 months:',round(oC[120],3),'(95% CI:',round(oL[120],3),'-',round(oU[120],3),')\n')
}




