
library(grid)
xaxes_simulate_time=function(xlm)
{
  
  # Mo-Fr
  xtk1 = 30*seq(xlm[1],xlm[2]/30,2)
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
    ylb = as.character(ytk)
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
    output <- read.csv(paste0('new dorm/output/data/Scenario_',scenario_id,'_acacia_daily.csv'))
    ylabel = "New cases"
    mastercol = 'darkorange'
  }else{
    output <- read.csv(paste0('new dorm/output/data/Scenario_',scenario_id,'_cassia_daily.csv'))
    ylabel = "New cases"
    mastercol = 'steelblue'
  }
  t  = c(1:nrow(output))
  oC = output$central
  oC[1] = oC[2]
  # oC[1] = output$central[2]
  oL= output$lower
  oL[1] = oL[2]
  oU = output$upper
  oU[1]  = oU[2]
  png(file,height=10,width=10,units='cm',res=900,pointsize=12)
  
  xlm = c(1,range(t)[2])
  logy=FALSE
  top=TRUE
  ylabel_gap=-3
  
  ylm=c(0,max(oU)*1.05)
  grid.newpage()
  pushViewport(plotViewport(c(3,3.5,1,1),
                            xscale=xlm,yscale=ylm))
  xaxes_simulate_time(xlm)
  
  # mastercol = '#ed2939'
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

paneller_CI(file='new dorm/output/plot/acacia_newdorm_s1_ci_epi.png',dormID = 1,scenario_id = 1)
paneller_CI(file='new dorm/output/plot/cassia_newdorm_s1_ci_epi.png',dormID = 2,scenario_id = 1)
