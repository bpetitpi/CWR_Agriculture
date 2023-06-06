############################################################################
##### initial settings (home directory, functions and library loading) #####
############################################################################

library(rgdal)
library(rgeos)
library(parallel)
library(raster)
library(perm)

source.dir<-"F:/InfoFlora/InfoFlora_EvaluationTool/BP/V2" #define the directory of the generic structure
my.dir<-paste(source.dir,"/data/",sep="")
my.dir.fct<-paste(source.dir,"/Code R/",sep="")
source(paste(my.dir.fct,"functions_predGam_work.R",sep=""))
source(paste(my.dir.fct,"functions_extractSp.R",sep=""))


env<-stack(
  paste(my.dir,"env/env100.tif",sep='')
)

my.sp.list<-list.files(paste(my.dir,'sp/Sp_list/',sep="")) # list of the species 
#sp.eval<-list.files(paste(my.dir,'results/evaluations/',sep=""))
#my.sp.list<-my.sp.list[which(!is.na((charmatch(my.sp.list,sp.eval))))]
#length(my.sp.list)
load(paste(my.dir,'sp/abs.Rdata',sep=""))
load(paste(my.dir,"env/bck",sep=""))
names(env)<-c('bio19_pcoldq',
           'bio18_pwarmq',
           'bio17_pdryq',
           'bio16_pwetq',
           'bio15_ps',
           'bio14_pdry',
           'bio13_pwet',
           'bio12_p',
           'bio11_tcoldq',
           'bio10_twarmq',
           'bio9_tdryq',
           'bio8_twetq',
           'bio7_tar',
           'bio6_tminc',
           'bio5_tmaxw',
           'bio4_ts',
           'bio3_iso',
           'bio2_dr',
           'bio1_tmean',
           'arridity',
           'gdd3Y_8110_ngb5_mwconic_',
           'sdiryy',
           'slp25', 
           'topos',
           'soilPH',
           'ndviSD',
           'ndviQ80',
           'ndviQ50',
           'ndviMIN',
           'ndviMAX',
           'ndviMEAN',
           'ForestQ95',
           'ForestQ25',
           'alt',
           'bias',
           'mask.tif')

  env.bck<-na.exclude(extract(env, abs.rand[,2:3]))

  ##### Download and disaggregation of the desired taxa
  ncores<-detectCores()-1
  if (ncores == 0) {ncores = 1}
  my.files.na<-list.files(paste0(my.dir,'env/varPval_NODATA'))
  my.names.na<-sort(as.numeric(sapply(my.files.na,list.name,code.rank=1,split="_")))
  my.files<-list.files(paste0(my.dir,'env/varPval'))
  my.names<-sort(as.numeric(sapply(my.files,list.name,code.rank=1,split="_")))
  total.list<-sort(as.numeric(sapply(my.sp.list,list.name,code.rank=1,split="_")))
  to.do<-which(is.na(match(total.list,c(my.names,my.names.na))))
  
  cl <- makeCluster(ncores)
  clusterExport(cl,varlist=ls(),envir=.GlobalEnv)
  clusterEvalQ(cl,{
    library(rgdal)
    library(rgeos)
    library(raster)
    library(perm)
  })
  
  parLapply(cl,to.do,var.test,my.sp.list=my.sp.list,
                   my.dir=my.dir,env.bck = env.bck,env =env)
  

stopCluster(cl)

