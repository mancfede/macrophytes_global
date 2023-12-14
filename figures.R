###FIGURES

library(viridis)
library(akima)
library(fields)
library(ggplot2)

se<-function(x){
  return(sd(x)/sqrt(length(x)))
}


############ GLOBAL SPECIES DIVERSITY

##Fig.1: maps of present diversity, hotpots losses & gains (Python)

##Fig.2: biodiversity loss (SSP3-7.0)
##[Supplementary Fig.6-7 (SSP2-4.5, SSP5-8.5)]

a<-read.csv("lat_trends.csv",header=T)

t_res<-5
scenarios<-unique(a$scenario)
taxa<-unique(a$taxon)

#a,b: trajectories

for(ssp in scenarios){
  
  a_sc<-a[a$scenario==ssp,]

  par(mfrow=c(1,2))
  
  sc_col<-1
    for (ta in taxa){
      a_ta<-subset(a_sc,taxon==ta)
      x_m<-aggregate(a_ta$loss~round(a_ta$year/t_res),FUN="mean")
      x_se<-aggregate(a_ta$loss~round(a_ta$year/t_res),FUN="se")
      x_mye<-aggregate(a_ta$loss~round(a_ta$year),FUN="mean")   #mean local diversity per year
      x_seye<-aggregate(a_ta$loss~round(a_ta$year),FUN="se")    ##s.e. per year
      ll<-x_m[,2]-1.96*x_se[,2]
      ul<-x_m[,2]+1.96*x_se[,2]
      plot(x_m[,1]*t_res,x_m[,2],type="l",col=plasma(3)[sc_col],main=ta,ylab="Diversity change",xlab="Year",ylim=c(min(ll)-0.2,max(ul)+0.2))
      polygon(c(rev(x_m[,1]*t_res),x_m[,1]*t_res),c(rev(ll),ul),col=plasma(3,alpha=0.2)[sc_col],border=NA)  
      sc_col<-sc_col+1
      print(c(ssp,ta))
      print(cbind(x_mye,x_seye))
    }

#c,d: colorplots   

  par(mfrow=c(1,2))
  cm<-colorRampPalette(colors = c("red", "white", "blue"))

    for (ta in taxa){
      a_ta<-subset(a_sc,taxon==ta)
      lat_m<-aggregate(a_ta$loss~a_ta$lat+a_ta$year,FUN='mean')  #average loss per latitude
      max_z<-max(abs(lat_m[,3]))
      min_z<-(-max_z)
      image.plot(interp(a_ta$year,a_ta$lat,a_ta$loss,duplicate='mean',nx=(2100-2015),ny=(max(a_ta$lat)-min(a_ta$lat))),
                col=cm(100),main=ta,xlab='year',ylab='latitude',
                cex.axis=1.2,cex.lab=1.2,las=1,legend.lab='change',zlim=c(min_z,max_z))
    }
}


#% of localities with >X% diversity loss in 2100 (3 scenarios)

thrs<-seq(from=-100,to=-10,by=10)
tab<-c()
for(ssp in scenarios){
  a_sc<-a[a$scenario==ssp,]
  for (ta in taxa){
    a_ta<-subset(a_sc,taxon==ta)
    a_ta2100<-subset(a_ta,year==2100)
    nloc<-nrow(a_ta2100)
    loc_loss<-c()
    for (thr in thrs){
      loc_loss<-c(loc_loss,round(sum(a_ta2100$loss<=thr)/nloc*100,2))
    }
    tab<-rbind(tab,c(ssp,ta,loc_loss))
  }
}

tab<-data.frame(tab)
names(tab)<-c("scenario","kingdom",thrs)


############ GENERIC GLOBAL MACROPHYTE HABITAT

##Fig.3: global macrophyte habitat trajectories
##[Supplementary Fig.9-10 (SSP2-4.5, SSP5-8.5)]

a<-read.csv("./hab_changes.csv",header=T)
a<-a[order(a$year),]
measures<-names(a)[5:15]

#% changes relative to 2015
h_0<-a[a$year==2015,measures]   

rel_change<-c()
for (y in c(2015:2100)){
  h_y<-a[a$year==y,measures]
  h_c<-100*(h_y-h_0)/h_0  
  rel_change<-rbind(rel_change,h_c)
  print (y)
}

rel_change[which(is.na(rel_change),arr.ind=T)]<-0  #0/0 (NaN)
rel_change[sapply(rel_change,simplify='matrix',is.infinite)]<-NA   #n/0 (Inf)
#remove NAs (Inf) from final plots

names<-c(colnames(a),paste0("rel_",measures))
a<-cbind(a,rel_change)
colnames(a)<-names

t_res<-5
scenarios<-unique(a$ssp)
taxa<-unique(a$kingdom)

glob<-subset(a,region=="global")

par(mfrow=c(2,3),mar=c(8,2,2,8),oma=c(2,2,2,2)) 

for (ta in taxa){
  a_ta<-glob[glob$kingdom==ta,]
 for (sc in scenarios){
  a_sc<-a_ta[a_ta$ssp==sc,]
  print(c(ta,sc,round(min(a_sc$rel_X90),1)))  #print max suitability loss (p>0.9)
  plot(0,0,type='n',xlim=c(2015,2100),ylim=c(-100,30),
       las=1,cex.axis=1.2,cex.lab=1.2,xlab='year',ylab='relative change (%)',main=paste(ta,sc)) 
  sc_col<-1
    for (measure in paste0("rel_",measures[7:11])){
      x_m<-aggregate(a_sc[,measure]~round(a_sc$year/t_res),FUN='mean') 
      lines(x_m[,1]*t_res,x_m[,2],lwd=2,col=viridis(5)[sc_col])
      sc_col<-sc_col+1
    }
 }
}
legend("right",title="Threshold",legend=c(0.5,0.6,0.7,0.8,0.9),col=viridis(5),pch=15,cex=1.5,pt.cex=1.5,
       y.intersp=1,inset=-0.5,xpd=T)


##Fig.4: 

#panel a:marine regions macrophyte habitat trajectories (km2)

regions<-unique(a$region)
regions<-subset(regions,regions!="global")

c<-c()

for (sc in scenarios){
  par(mfrow=c(3,5))
  a_sc<-a[a$ssp==sc,]
  
  #use "X90" instead of "X60" to plot Supplementary Fig. 8
    for (re in regions){
      a_reg<-a_sc[a_sc$region==re,]     
      plot(0,0,type='n',xlim=c(2015,2100),ylim=c(0.001*min(a_reg[,"X60"],na.rm=T),0.001*max(a_reg[,"X60"],na.rm=T)),
         las=1,cex.axis=1,cex.lab=1,xlab='year',ylab='Habitat extent (km2*10e-3)',main=re)
      sc_col<-1
        for (ta in taxa){
          a_ta<-a_reg[a_reg$kingdom==ta,]
          x_m<-aggregate(0.001*a_ta$X60~round(a_ta$year/t_res),FUN='mean',na.rm=T)
          lines(x_m[,1]*t_res,x_m[,2],lwd=2,col=plasma(3)[sc_col])
          sc_col<-sc_col+1
        }
      }
  

#panel b:% of global macrophyte habitat hosted in each marine regions (macrophyte habitat trajectories (years 2015,2100)

  cols<-c("#17becfff","#888888ff","#caca50ff","#b0e7a5ff","#d2c1deff","#ffc790ff",
                     "#1f77b4ff","#905c52ff","#ffb0afff","#aec7e8ff","#e995cfff","#9edae5ff")

  par(mfrow=c(2,1))

      for (ta in taxa){
        a_ta<-a_sc[a_sc$kingdom==ta,]
        b<-cbind(rep(sc,length(regions)),
                 rep(ta,length(regions)),
                 regions)
        for (ye in c(2015,2100)){
          a_ye<-a_ta[a_ta$year==ye,]
          gl<-a_ye[a_ye$region=="global",c("X60")]
          re<-a_ye[a_ye$region!="global",c("X60")]
          pe<-100*(re/gl)
          b<-cbind(b,pe)
          tot_ext<-round(sum(re),2)
          print(c(sc,ta,ye,tot_ext))  #total global suitable habitat extension (add to barplots with square brackets)
        }
      c<-rbind(c,b)
 
    barplot(b[,5:4],
            names.arg=c(2100,2015),
            horiz=T,
            main=ta,
            col=cols)
  }
}


c<-data.frame(scenario=c[,1],
              taxon=c[,2],
              regions=c[,3],
              year2015=round(as.numeric(c[,4]),1),
              year2100=round(as.numeric(c[,5]),1))
              

############ MACROPHYTE AREA OF OCCUPANCY

##Fig.5: average macrophyte area-of-occupancy (range) trajectories

a<-read.csv("range_changes.csv",header=T)
a<-subset(a,region=="global")
a<-a[order(a$year),]

lsp_0<-a[a$year==2015,c('loss','suit','prop')]

#% changes relative to 2015
rel_change<-c()
for (y in c(2015:2100)){
  lsp_y<-a[a$year==y,c('loss','suit','prop')]
  lsp_c<-100*(lsp_y-lsp_0)/lsp_0   
  rel_change<-rbind(rel_change,lsp_c)
  print (y)
}

rel_change[which(is.na(rel_change),arr.ind=T)]<-0

names<-c(colnames(a),'rel_loss','rel_suit','rel_prop')
a<-cbind(a,rel_change)
colnames(a)<-names

#remove floating Sargassum spp.
to_del<-which(a$species=="Sargassum_fluitans" | a$species== "Sargassum_natans" | a$species=="Sargassum_pusillum")
a<-a[-to_del,]

t_res<-5
scenarios<-unique(a$ssp)
taxa<-unique(a$kingdom)

par(mfrow=c(2,3))

#km2

for (sc in scenarios){
  a_sc<-a[a$ssp==sc,]
  plot(0,0,type='n',xlim=c(2015,2100),ylim=c(200,700),
       las=1,cex.axis=1.2,cex.lab=1.2,xlab='year',ylab='range (km2*10e-3)',main=sc) 
  sc_col<-1
  for (ta in taxa){
    a_ta<-a_sc[a_sc$kingdom==ta,]
    lty=1
    for (measure in c('loss','prop')){
      x_m<-aggregate(0.001*a_ta[,measure]~round(a_ta$year/t_res),FUN='mean')
      x_se<-aggregate(0.001*a_ta[,measure]~round(a_ta$year/t_res),FUN='se')
      ll<-x_m[,2]-1.96*x_se[,2]
      ul<-x_m[,2]+1.96*x_se[,2]
      if (lty==1){
        polygon(c(rev(x_m[,1]*t_res), x_m[,1]*t_res), c(rev(ll),ul), col = plasma(3,alpha=0.2)[sc_col],border=NA)}  
      lines(x_m[,1]*t_res,x_m[,2],lwd=2,col=plasma(3)[sc_col],lty=lty)
      lty<-lty+1
    }
    sc_col<-sc_col+1
  }    
}

#% variation

for (sc in scenarios){
  a_sc<-a[a$ssp==sc,]
  plot(0,0,type='n',xlim=c(2015,2100),ylim=c(-15,10),
       las=1,cex.axis=1.2,cex.lab=1.2,xlab='year',ylab='range loss (%)',main=sc) 
  sc_col<-1
  for (ta in taxa){
    a_ta<-a_sc[a_sc$kingdom==ta,]
    lty=1
    for (measure in c('rel_loss','rel_prop')){
      x_m<-aggregate(a_ta[,measure]~round(a_ta$year/t_res),FUN='mean')
      x_se<-aggregate(a_ta[,measure]~round(a_ta$year/t_res),FUN='se')
      x_mye<-round(aggregate(a_ta[,measure]~round(a_ta$year),FUN="mean"),1)  #annual relative change (mean)
      x_seye<-round(aggregate(a_ta[,measure]~round(a_ta$year),FUN="se"),1)   #annual relative change (se)
      ll<-x_m[,2]-1.96*x_se[,2]
      ul<-x_m[,2]+1.96*x_se[,2]
      polygon(c(rev(x_m[,1]*t_res), x_m[,1]*t_res), c(rev(ll),ul), col = plasma(3,alpha=0.2)[sc_col], border = NA)  
      lines(x_m[,1]*t_res,x_m[,2],lwd=2,col=plasma(3)[sc_col],lty=lty)
      lty<-lty+1
      print(c(sc,ta,measure))
      print(cbind(x_mye,x_seye))
    }
    sc_col<-sc_col+1
  }    
}


################ SUPPLEMENTARY INFORMATION

##Supplementary Fig. 11: range loss distribution

final<-a[a$year==2100,]

par(mfrow=c(3,2))
for (sc in scenarios){
  final_sc<-final[final$ssp==sc,]
  for (ta in taxa){
    final_ta<-final_sc[final_sc$kingdom==ta,]
    final_ta<-subset(final_ta,rel_loss<0)   #select only range loss (rel_loss<0)
    hist(final_ta$rel_loss,
         col="grey",xlab="range loss (%)",main=paste(sc,ta,'loss'))
  } 
}  


##Supplementary Table 1-2: list of modelled macrophyte species, present and future area-of-occupancy (range)
#extension, environmental variables selected in the models

tab_si<-c()
for (ta in taxa){
  a_ta<-a[a$kingdom==ta,]
  for (sc in scenarios){
    a_sc<-a_ta[a_ta$ssp==sc,]
    for (sp in unique(a_sc$species)){
      a_sp<-a_sc[a_sc$species==sp,]
      spec<-chartr("_"," ",sp)
      r_2015<-round(a_sp[a_sp$year==2015,]$loss,2)
      r_2100<-round(a_sp[a_sp$year==2100,]$loss,2)
      rel<-round(a_sp[a_sp$year==2100,]$rel_loss,2)
      tab_si<-rbind(tab_si,c(sc,ta,spec,r_2015,r_2100,rel)) 
    }
  }
}

tab_si<-data.frame(scenario=as.numeric(tab_si[,1]),
                   kingdom=tab_si[,2],
                   species=tab_si[,3],
                   range_2015=as.numeric(tab_si[,4]),
                   range_2100=as.numeric(tab_si[,5]),
                   range_change=as.numeric(tab_si[,6]))

tab_si[tab_si$scenario=="245",1]<-"SSP2-4.5"
tab_si[tab_si$scenario=="370",1]<-"SSP3-7.0"
tab_si[tab_si$scenario=="585",1]<-"SSP5-8.5"

our_spp<-unique(tab_si$species)

##add information on variable importance

a<-read.csv("RF_var_imp.csv",header=T)
head(a)

a[a$scenario=="ssp245",1]<-"SSP2-4.5"
a[a$scenario=="ssp370",1]<-"SSP3-7.0"
a[a$scenario=="ssp585",1]<-"SSP5-8.5"

#remove species not used with alpha-hulls + floating Sargassum
a<-subset(a,species %in% our_spp)

#aggregate monthly variables
#for each species and scenario order variables based on importance

vars<-c("depth","solar","siconc","sos","tos","tas","intpp")

b<-c()
for (i in 1:nrow(a)){
  aggr<-c()
  for (j in vars){
    aggr<-c(aggr,sum(a[i,grepl(j,names(a))]))
  }
  b<-rbind(b,c(a[i,]$scenario,a[i,]$species,aggr))
}

b<-data.frame(scenario=b[,1],
              species=b[,2],
              depth=as.numeric(b[,3]),
              solar=as.numeric(b[,4]),
              siconc=as.numeric(b[,5]),
              sos=as.numeric(b[,6]),
              tos=as.numeric(b[,7]),
              tas=as.numeric(b[,8]),
              intpp=as.numeric(b[,9]))

#for each row, rank variable based on their importance (removing those with 0)
#make string vector with ordered variables and add it to table S1-S2

#change variable names
vars2<-c("depth","radiation","sea-ice cover","salinity","sea surface temperature","air temperature","primary productivity")

var_imp<-c()
for (i in 1:nrow(b)){
  all_vars<-b[i,3:9]
  names(all_vars)<-vars2
  to_del<-which(all_vars==0)
  if (length(to_del)>0){all_vars<-all_vars[-to_del]}
  ord_vars<-names(all_vars)[order(all_vars,decreasing=T)]
  var_imp<-rbind(var_imp,paste(ord_vars,collapse=", "))
}

c<-cbind(b[,1:2],var_imp)

#merge with table S1-S2
tab_si<-merge(tab_si,c,by=c("scenario","species"))
tab_si<-tab_si[order(tab_si$kingdom),]


##Supplementary Table 3: 10 most impacted macrophyte species
tab_sii<-tab_si[order(tab_si$range_change)[1:10],]


##Supplementary Fig.2: variable importance

#add kingdom
d<-merge(b,tab_si[,1:3],by=c("scenario","species"))
head(d)
scenarios<-unique(d$scenario)
taxa<-unique(d$kingdom)

summ<-c()
for (sc in scenarios){
  d_sc<-d[d$scenario==sc,]
  for (ta in taxa){
    d_ta<-d_sc[d_sc$kingdom==ta,]
    for (j in names(d[3:9])){
      x_m<-mean(d_ta[,j])
      x_se<-se(d_ta[,j])
      summ<-rbind(summ,c(sc,ta,j,x_m,x_se))
    }
  }
  print(summ)
}

summ<-data.frame(scenario=summ[,1],
                 kingdom=as.factor(summ[,2]),
                 var=factor(summ[,3]),
                 mean=as.numeric(summ[,4]),
                 se=as.numeric(summ[,5]))

p<-ggplot(data=summ,
          aes(x=var,
              y=mean,
              fill=kingdom))+
  geom_bar(stat="identity",position="dodge",width=0.6)+
  geom_errorbar(aes(ymin=mean-se,
                    ymax=mean+se),
                colour="black",
                width=.2,
                position=position_dodge(.5))+
  theme_classic()+
  coord_flip()+
  scale_fill_manual(values=c('#0d0887ff','#cc4678ff'))+
  facet_wrap(~scenario)  

p


