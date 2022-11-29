

xaxes2=function(xlm)
{
    
    # Mo-Fr
    xtk1 = 48*seq(0,6,1)#c(1:(xlm[2]+1))#rep(1:24,7)
    xtk2 = xtk1+24
    
    i=which((xtk1>=xlm[1])&(xtk2<xlm[2]))
    xtk1touse=xtk1[i]
    xtk2touse=xtk2[i]
    i=which((xtk1<xlm[2])&(xtk2>=xlm[2]))
    if(length(i)>0)
    {
      xtk1touse=c(xtk1touse,xtk1[i])
      xtk2touse=c(xtk2touse,xlm[2])
    }
    
    for(k in seq(xtk1touse))grid.polygon(unit(c(xtk1touse[k],xtk1touse[k],
                                                xtk2touse[k],xtk2touse[k]),'native'),
                                         unit(c(0,1,1,0),'npc'),
                                         gp=gpar(col=NA,fill='#e0dad4'))
    
    
    xtk =  xtk1
    xtk = xtk[xtk<=xlm[2]]
    grid.xaxis(at=xtk,label=FALSE,main = TRUE)
    grid.text('00',x=unit(0.0 +0 +0*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('06',x=unit(0.0 +6 +0*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('12',x=unit(0.0+12+0*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('18',x=unit(0.0+18+0*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('Mon',x=unit(0.0+0*24+12+0*24,'native'),y=unit(-2,'lines'))
    
    grid.text('00',x=unit(0.0 +0 +1*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('06',x=unit(0.0 +6 +1*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('12',x=unit(0.0+12+1*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('18',x=unit(0.0+18+1*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('Tue',x=unit(0.0+1*24+12,'native'),y=unit(-2,'lines'))
    
    grid.text('00',x=unit(0.0 +0 +2*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('06',x=unit(0.0 +6 +2*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('12',x=unit(0.0+12+2*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('18',x=unit(0.0+18+2*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('Wed',x=unit(0.0+2*24+12,'native'),y=unit(-2,'lines'))
    
    grid.text('00',x=unit(0.0 +0 +3*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('06',x=unit(0.0 +6 +3*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('12',x=unit(0.0+12+3*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('18',x=unit(0.0+18+3*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('Thu',x=unit(0.0+3*24+12,'native'),y=unit(-2,'lines'))
    
    grid.text('00',x=unit(0.0 +0 +4*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('06',x=unit(0.0 +6 +4*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('12',x=unit(0.0+12+4*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('18',x=unit(0.0+18+4*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('Fri',x=unit(0.0+4*24+12,'native'),y=unit(-2,'lines'))
    
    grid.text('00',x=unit(0.0 +0 +5*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('06',x=unit(0.0 +6 +5*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('12',x=unit(0.0+12+5*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('18',x=unit(0.0+18+5*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('Sat',x=unit(0.0+5*24+12,'native'),y=unit(-2,'lines'))
    
    grid.text('00',x=unit(0.0 +0 +6*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('06',x=unit(0.0 +6 +6*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('12',x=unit(0.0+12+6*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    # grid.text('18',x=unit(0.0+18+6*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    grid.text('Sun',x=unit(0.0+6*24+12,'native'),y=unit(-2,'lines'))
    grid.text('00',x=unit(0.0 +0 +7*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
    
    grid.text('Day of Week',x=unit(0.5,'npc'),y=unit(-3,'lines'),gp=gpar(fontface='bold'))
    
}

xaxes3=function(xlm)
{
  
  # Mo-Fr
  xtk1 = seq(0,48,2)#c(1:(xlm[2]+1))#rep(1:24,7)
  xtk2 = xtk1+1
  
  i=which((xtk1>=xlm[1])&(xtk2<xlm[2]))
  xtk1touse=xtk1[i]
  xtk2touse=xtk2[i]
  i=which((xtk1<xlm[2])&(xtk2>=xlm[2]))
  if(length(i)>0)
  {
    xtk1touse=c(xtk1touse,xtk1[i])
    xtk2touse=c(xtk2touse,xlm[2])
  }
  
  for(k in seq(xtk1touse))grid.polygon(unit(c(xtk1touse[k],xtk1touse[k],
                                              xtk2touse[k],xtk2touse[k]),'native'),
                                       unit(c(0,1,1,0),'npc'),
                                       gp=gpar(col=NA,fill='#e0dad4'))
  
  
  xtk =  xtk1
  xtk = xtk[xtk<=xlm[2]]
  grid.xaxis(at=xtk,label=FALSE,main = TRUE)
  grid.text('00',x=unit(0.0 +0 +0*24,'native'),y=unit(-1,'lines'))#,gp=gpar(fontsize=12))
  grid.text('06',x=unit(0.0 +6 +0*24,'native'),y=unit(-1,'lines'))#,gp=gpar(fontsize=12))
  grid.text('12',x=unit(0.0+12+0*24,'native'),y=unit(-1,'lines'))#,gp=gpar(fontsize=12))
  grid.text('18',x=unit(0.0+18+0*24,'native'),y=unit(-1,'lines'))#,gp=gpar(fontsize=12))
  grid.text('24',x=unit(0.0+24+0*24,'native'),y=unit(-1,'lines'))#,gp=gpar(fontsize=12))
  
  # grid.text('Mon',x=unit(0.0+0*24+12+0*24,'native'),y=unit(-2,'lines'))
  # 
  # grid.text('00',x=unit(0.0 +0 +1*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('06',x=unit(0.0 +6 +1*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('12',x=unit(0.0+12+1*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('18',x=unit(0.0+18+1*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('Tue',x=unit(0.0+1*24+12,'native'),y=unit(-2,'lines'))
  # 
  # grid.text('00',x=unit(0.0 +0 +2*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('06',x=unit(0.0 +6 +2*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('12',x=unit(0.0+12+2*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('18',x=unit(0.0+18+2*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('Wed',x=unit(0.0+2*24+12,'native'),y=unit(-2,'lines'))
  # 
  # grid.text('00',x=unit(0.0 +0 +3*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('06',x=unit(0.0 +6 +3*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('12',x=unit(0.0+12+3*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('18',x=unit(0.0+18+3*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('Thu',x=unit(0.0+3*24+12,'native'),y=unit(-2,'lines'))
  # 
  # grid.text('00',x=unit(0.0 +0 +4*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('06',x=unit(0.0 +6 +4*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('12',x=unit(0.0+12+4*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('18',x=unit(0.0+18+4*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('Fri',x=unit(0.0+4*24+12,'native'),y=unit(-2,'lines'))
  # 
  # grid.text('00',x=unit(0.0 +0 +5*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('06',x=unit(0.0 +6 +5*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('12',x=unit(0.0+12+5*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('18',x=unit(0.0+18+5*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('Sat',x=unit(0.0+5*24+12,'native'),y=unit(-2,'lines'))
  # 
  # grid.text('00',x=unit(0.0 +0 +6*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('06',x=unit(0.0 +6 +6*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('12',x=unit(0.0+12+6*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # # grid.text('18',x=unit(0.0+18+6*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # grid.text('Sun',x=unit(0.0+6*24+12,'native'),y=unit(-2,'lines'))
  # grid.text('00',x=unit(0.0 +0 +7*24,'native'),y=unit(-1,'lines'),gp=gpar(fontsize=8))
  # 
  grid.text('Hour of Day',x=unit(0.5,'npc'),y=unit(-2,'lines'),gp=gpar(fontface='bold'))
  
}

yaxis = function(ylm,log=FALSE)
{
  if(!log)
  {
    
    ytk = pretty(ylm,n=10)
    if(min(ytk)==-100) ytk=ytk[-1]
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
    for(ij in seq(ylb)){if(ytk[ij]>=3)ylb[ij]=paste0(10^ytk[ij]/1000,'k')}
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
    if(length(Z)>0)grid.yaxis(at=Z,label=FALSE,gp=gpar(cex=0.5))
  }
  
}