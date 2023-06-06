library(raster)
library(parallel)
rm.out<-function(vec,outinf=0.01,outsup=0.99){
  to.rm.out<-which(vec<quantile(vec,outinf) | vec>quantile(vec,outsup))
  return(to.rm.out)
}


source.dir<-"F:/InfoFlora/InfoFlora_EvaluationTool/BP/Structure_generic" #define the directory of the generic structure
my.dir<-paste(source.dir,"/data",sep="")
my.dir.sp<-paste(my.dir,"/sp/sp_list/",sep="")
wd.load<-paste(my.dir,"/env/",sep="")

##### Generating the bias map based on the pool of all the observations
load(paste(my.dir,"/sp/synthesis.Rdata",sep=""))
my.mask<-raster(paste(wd.load,"mask.tif",sep=""))
pts.ras<-rasterize(sp.pool.table[,1:2],my.mask,field=sp.pool.table[,5],fun='sum',background=0)*my.mask #Rasterize the whole extraction
mov.win<-focalWeight(pts.ras,500,type='circle') #Generate a bias matrix for a 1 km radius 
bias<-focal(test, mov.win,fun='sum',na.rm=TRUE) # Generate a bias map
writeRaster(bias,filename=paste(my.dir,"/sp/bias.tif",sep=""),format="GTiff",datatype="INT2U",overwrite=TRUE, 
            options=c("COMPRESS=DEFLATE", "PREDICTOR=2"))

##### Sampling 10000 pseudoabsences following the sample bias
bias<-raster(paste(my.dir,"/sp/bias.tif",sep="")) #load the bias map 
bck<-rasterToPoints(bias)#convert the bias map into points
bck.sample<-sample(1:nrow(bck),size=10000,prob=bck[,3]/sum(bck[,3])) # sample 10000 points following the bias map
my.abs<-bck[bck.sample,1:2]
names(my.abs)<-c("x","y")
save(my.abs,file=paste(my.dir,"/sp/PseudoAbsences_bias.Rdata",sep=""),compress="xz")


##### Generate a stratified sampling

# 1-generate the strata
env.dir<-'P:/InfoFlora_EvaluationTool/BP/V2/data/formatted/'
load(paste0(env.dir,'bck'))
env<-stack(paste0(env.dir,'bio19_pcoldq','_8110_formatted_round.tif'),
               paste0(env.dir,'bio18_pwarmq','_8110_formatted_round.tif'),
               paste0(env.dir,'bio17_pdryq','_8110_formatted_round.tif'),
               paste0(env.dir,'bio16_pwetq','_8110_formatted_round.tif'),
               paste0(env.dir,'bio15_ps','_8110_formatted_round.tif'),
               paste0(env.dir,'bio14_pdry','_8110_formatted_round.tif'),
               paste0(env.dir,'bio13_pwet','_8110_formatted_round.tif'),
               paste0(env.dir,'bio12_p','_8110_formatted_round.tif'),
               paste0(env.dir,'bio11_tcoldq','_8110_formatted_round.tif'),
           paste0(env.dir,'bio10_twarmq','_8110_formatted_round.tif'),
           paste0(env.dir,'bio9_tdryq','_8110_formatted_round.tif'),
           paste0(env.dir,'bio8_twetq','_8110_formatted_round.tif'),
           paste0(env.dir,'bio7_tar','_8110_formatted_round.tif'),
           paste0(env.dir,'bio6_tminc','_8110_formatted_round.tif'),
           paste0(env.dir,'bio5_tmaxw','_8110_formatted_round.tif'),
           paste0(env.dir,'bio4_ts','_8110_formatted_round.tif'),
           paste0(env.dir,'bio3_iso','_8110_formatted_round.tif'),
           paste0(env.dir,'bio2_dr','_8110_formatted_round.tif'),
           paste0(env.dir,'bio1_tmean','_8110_formatted_round.tif'),
           paste0(env.dir,'arridity_formatted_round.tif'),
           paste0(env.dir,'gdd3Y_8110_ngb5_mwconic_','formatted_round.tif'),
           paste0(env.dir,'sdiryy_','formatted_round.tif'),
           paste0(env.dir,'slp25_','formatted_round.tif'), 
           paste0(env.dir,'topos_','formatted_round.tif'),
           paste0(env.dir,'vgh_','formatted_round.tif'),
           paste0(env.dir,'soil_','formatted_round.tif'),
           paste0(env.dir,'alt_','formatted_round.tif'),
           paste0(env.dir,'bias_formatted_round.tif'),
           paste0(env.dir,'mask.tif'))

env<-env[[3:1]]
env.pts<-na.exclude(cbind(extract(env,bck[,1:2]),bck[,1:2]))

alt.class<-c(800, 1600, 2400, 3200,4800)
fac.class<-c(1,2)
topo.class<-0
slope.class<-c(5,20,40,90)
rad.class<-c(6000,12000,18000,22000)

env.strat<-env.pts[,22:27]

env.strat[env.strat[,1]<6000,1]<-1
env.strat[env.strat[,1]>=6000 & env.strat[,1]<12000,1]<-2
env.strat[env.strat[,1]>=12000 & env.strat[,1]<18000,1]<-3
env.strat[env.strat[,1]>=18000 & env.strat[,1]<22000,1]<-4

env.strat[env.strat[,2]<=5,2]<-1
env.strat[env.strat[,2]>5 & env.strat[,2]<=20,2]<-2
env.strat[env.strat[,2]>20 & env.strat[,2]<=40,2]<-3
env.strat[env.strat[,2]>40 & env.strat[,2]<90,2]<-4

env.strat[env.strat[,3]<0,3]<-1
env.strat[env.strat[,3]>=0,3]<-2

env.strat[env.strat[,6]<=80,6]<-1
env.strat[env.strat[,6]>80 & env.strat[,6]<=160,6]<-2
env.strat[env.strat[,6]>160 & env.strat[,6]<=240,6]<-3
env.strat[env.strat[,6]>240 & env.strat[,6]<=320,6]<-4
env.strat[env.strat[,6]>320,6]<-5

my.strat<-apply(env.strat,1,toString)

my.strat.cat<-unique(my.strat)

k<-10000
n.full.strat<-length(my.strat.cat)
dispo<-round(k/(length(my.strat.cat)))

my.pts<-c()
for (i in 1:length(my.strat.cat)){
  a<-length(which(my.strat==my.strat.cat[i]))
  my.ptsi<-env.pts[which(my.strat==my.strat.cat[i]),28:30]
  if (a>=dispo){
    my.ptsi[,1]<-my.ptsi[,1]+1
    my.pts2<-my.ptsi[sample(1:nrow(my.ptsi),size =dispo,prob=my.ptsi[,1]),]
  }else{
    my.pts2<-my.ptsi
  }
  my.pts<-rbind(my.pts,my.pts2)
}


abs.strat<-my.pts
abs.rand.bias<-env.pts[sample(1:nrow(env.pts),size = 10000,prob = env.pts[,28]+1),28:30]
abs.rand<-env.pts[sample(1:nrow(env.pts),size = 10000),28:30] # we used this one

save (abs.strat,abs.rand.bias,abs.rand, file = 'F:/InfoFlora/InfoFlora_EvaluationTool/BP/V2/data/formatted/abs', compress ='xz')




