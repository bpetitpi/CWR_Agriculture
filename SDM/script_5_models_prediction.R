############################################################################
##### initial settings (home directory, functions and library loading) #####
############################################################################

library(rgdal)
library(rgeos)
library(FNN)
library(Hmisc)
library(mgcv)
library(parallel)
library(biomod2)
library(raster)
library(scales)
library(dplyr)

source.dir<-"F:/InfoFlora/InfoFlora_EvaluationTool/BP/V2" #define the directory of the generic structure
my.dir<-paste(source.dir,"/data/",sep="")
my.dir.fct<-paste(source.dir,"/CodesR/",sep="")
source(paste(my.dir.fct,"functions_pred.R",sep=""))
source(paste(my.dir.fct,"functions_evaluation.R",sep=""))

env<-stack(
  paste(my.dir,"env/env100.tif",sep='')
)

my.sp.list<-list.files(paste(my.dir,'sp/Sp_list/',sep="")) # list of the species 
#sp.eval<-list.files(paste(my.dir,'results/evaluations/',sep=""))
#my.sp.list<-my.sp.list[which(!is.na((charmatch(my.sp.list,sp.eval))))]
#length(my.sp.list)
mod.ini<-read.table(paste0(source.dir,'/data/inidat/model_inidat.txt'),h=T,stringsAsF=F)
mod.ini<-mod.ini[mod.ini$Nobs>=10,]
load(paste(my.dir,'sp/abs.Rdata',sep=""))
load(paste(my.dir,"env/bck",sep=""))
var.list<-c('bio19_pcoldq',
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

names(env)<-var.list

Index<- 1 # index of the species to run
my.dir<-paste(source.dir,"/data/",sep="")

my.sp.id<-mod.ini$code_espece[Index]
sp.inidat<-mod.ini[Index,]
print(my.sp.id)

if (sp.inidat$modeling=='EM'){
  sp.th<-read.csv(paste0(my.dir,'results/thresholds/',my.sp.id,'_thresholds.csv'))
  Ew<-read.csv(paste0(my.dir,'results/Eweights/',my.sp.id,'_Eweights.csv'))
  sp.w<-Ew[,-1]
  vari<-which(sp.inidat[2:34]==TRUE)
  
  ks=3 
  min.dist=99
  my.mask.pts<-SpatialPoints(coord=bck[,1:2])
  exclude=500
  geo.exclude=5000
  
  load(paste(my.dir,'Sp/Sp_list/',my.sp.id,sep=""))
  vcol<-19
  to.rm<- which(is.na(match(my.sp.tot[,1],unique(my.sp.dist[,4]))))
  if(length(to.rm)>0){
    my.sp.xy<-my.sp.tot[-to.rm,c(11:12,1,26)]
  }else{
    my.sp.xy<-my.sp.tot[,c(11:12,1,26)]
  }
  
  vld.xy<-my.sp.tot[my.sp.tot[,26]>exclude & my.sp.tot[,26]<geo.exclude & my.sp.tot[,19]>=300
                    & (my.sp.tot[,23]!=3 | my.sp.tot[,23]!=4),c(11,12,1,26)]
  my.sp.xy<-rbind(my.sp.xy,vld.xy) 
  
  if(length(unique(my.sp.dist[,4])) >= 10 | nrow(my.sp.xy)>=10 ){
    
    
    ##### remove aggregated data
    if(min.dist>0){
      my.sp.xy.dup<-zerodist(SpatialPoints(coord=my.sp.xy[,1:2]),zero=min.dist)
      prec<-cbind(my.sp.xy[my.sp.xy.dup[,1],4],my.sp.xy[my.sp.xy.dup[,2],4])
      to.rm<-unique(apply(cbind(my.sp.xy.dup,prec),1,dsg.accurate))
      
      if(length(to.rm)>0){
        my.sp.xy<-my.sp.xy[-to.rm,]
        my.sp.dist<-my.sp.dist[which(!is.na(match(my.sp.dist[,4],my.sp.xy[,3]))),]
      }
      
    }
    
    my.abs<-abs.rand[,2:3]
    colnames(my.abs)<-colnames(my.sp.xy)[1:2]
    knn.tot<-log(round(knnx.dist(my.sp.xy[,1:2],coordinates(my.mask.pts),k=2)/100)[,1]+1)
    Ndist<-rasterize(my.mask.pts,env,field=knn.tot)
    env<-stack(env[[vari]],Ndist)
    names(env)[length(names(env))]<-'Ndist'
    Ndist<-c(log(round(knn.dist(my.sp.xy[,1:2],k=2)/100)[,1]+1),
             log(round(knnx.dist(my.sp.xy[,1:2], my.abs,k=2)/100)[,1]+1))
    sp.geo<-c(rep(1,nrow(my.sp.xy)),rep(0,nrow(my.abs)))
    
    data.env<-data.frame(extract(env,rbind(my.sp.dist[,1:2],setNames(my.abs,colnames(my.sp.dist[,1:2])))))
    names(data.env)[ncol(data.env)]<-"Ndist"
    
    data.resp<-c(rep(1,nrow(my.sp.dist)),rep(0,nrow(my.abs)))
    df.tmp.input<-cbind(data.resp,data.env)
    colnames(df.tmp.input)[1]<-"sp"
    
    w2<-tapply(my.sp.dist[,4], as.factor(my.sp.dist[,4]),re.w)
    wi<-match(my.sp.dist[,4],names(w2))
    w3<-w2[wi]
    w.pres<-sum(w3)/100
    w.abs<-nrow(my.abs)/w.pres
    w<-c(w3*round(w.abs),rep(100,nrow(my.abs)))
    w.geo<-c(rep(length(which(sp.geo==0))/length(which(sp.geo==1)),length(which(sp.geo==1))),rep(1,length(which(sp.geo==0))))
    
    m1.geo<-gam(sp.geo ~ s(Ndist, k = 3),family = binomial(link = "logit"), control = gam.control(maxit = 500),w=round(10*w.geo))
    sp.pred.geo<-round(1000*predict(env,m1.geo,type='response'))
    geo.th<-round(sp.th[3,2])
    sp.dist.geo<-round(1000*m1.geo$fitted[1:sum(sp.geo)])
    sp.dist.geo<-sp.dist.geo[which(sp.dist.geo>=geo.th)]

    # the name of studied species
    myRespName <- paste0('sp_',strsplit(my.sp.id,'_')[[1]][c(1)],'Full')
    #myRespName<- "mySpB"
    
    if (dir.exists(paste0(source.dir,'/data/models/',myRespName))){
      unlink(paste0(source.dir,'/data/models/',myRespName))
    }
    dir.create(paste0(source.dir,'/data/models/',myRespName),recursive = TRUE)
    setwd(paste0(source.dir,'/data/models/',myRespName))
    
    
    # the presence/absences data for our species 
    myResp <- na.exclude(df.tmp.input)$sp
    myExpl <- na.exclude(df.tmp.input)[,2:((ncol(df.tmp.input)-1))]
    myXY<-na.exclude(cbind(rbind(my.sp.dist[,2:3],my.abs),df.tmp.input))[,1:2]
    wz<-na.exclude(cbind(w,df.tmp.input))[,1]
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = myExpl,
                                         resp.name = myRespName,
                                         resp.xy = myXY)
    
    myBiomodOption <- BIOMOD_ModelingOptions(GAM = list(k=ks,control = gam.control(maxit = 500)),
                                             MAXENT = list(Threshold=FALSE,
                                                           path_to_maxent.jar = #'/home/petitpie/tools/maxent.jar',
                                                             'F:/InfoFlora/tools',
                                                           betamultiplier =1.5))
    
    myBiomodModelOut <- BIOMOD_Modeling( 
      myBiomodData, 
      models = c('MAXENT','GAM','RF'), 
      bm.options = myBiomodOption, 
      nb.rep=1, 
      data.split.perc=100, 
      weights  = wz, 
      var.import=3,
      metric.eval = c('ROC'),
      save.output = FALSE,
      scale.models = FALSE,
      do.full.models = TRUE,
      modeling.id = paste(myRespName,sep=""))
    
    var.import<-get_variables_importance(myBiomodModelOut) %>% 
      group_by(expl.var) %>% 
      summarise(var.imp = mean(var.imp)) 
    colnames(var.import)[ncol(var.import)]<-'E'
    if (!dir.exists(paste0(my.dir,"/results/varImport"))){
      dir.create(paste0(my.dir,"/results/varImport"),recursive = TRUE)
    }
    write.csv(var.import,paste(my.dir,"results/varImport/",my.sp.id,"_varImport.csv",sep=""),quote=FALSE,
              row.names=TRUE)
    
    pred.all<- BIOMOD_Projection(
      bm.mod = myBiomodModelOut,
      new.env = env[[-dim(env)[3]]],
      proj.name = 'Full',
      models.chosen = 'all',
      metric.binary = NULL,
      compress = 'xz',
      build.clamping.mask = F )
    
    m_w_name<-substring(Ew[,1],2)
    
    my.proj<-stack(paste0(getwd(),'/',"sp.",strsplit(my.sp.id,'_')[[1]][c(1)],'Full/proj_Full/proj_Full_sp.',strsplit(my.sp.id,'_')[[1]][c(1)],'Full.tif'))
    my.proj.name<-names(my.proj)
    model_sucess<-sapply(m_w_name,grep,my.proj.name)
    my.proj.E<-w.mean(my.proj, sp.w[lengths(model_sucess)>0])
    my.proj.E<-round(my.proj.E/10)
    
    eco.th<-round(sp.th[3,3]/10)
    sp.pred.eco<-round(get_predictions(myBiomodModelOut,model.as.col=TRUE)[1:sum(df.tmp.input$sp),]/10)
    sp.pred.eco<-round(apply(sp.pred.eco,1,weighted.mean,sp.w[lengths(model_sucess)>0]))
    sp.dist.eco<-sp.pred.eco[which(sp.pred.eco>=eco.th)]

    if(!dir.exists(paste0(my.dir,"results/rawPred_geo"))){
      dir.create(paste0(my.dir,"results/rawPred_geo"))
    }
    if(!dir.exists(paste0(my.dir,"results/rawPred_eco"))){
      dir.create(paste0(my.dir,"results/rawPred_eco"))
    }

    writeRaster(sp.pred.geo,filename=paste(my.dir,"results/rawPred_geo/",my.sp.id,"_geo",sep=""),format="GTiff",datatype="INT2U",overwrite=TRUE, 
                options=c("COMPRESS=DEFLATE", "PREDICTOR=2"))
    writeRaster(my.proj.E,filename=paste(my.dir,"results/rawPred_eco/",my.sp.id,"_eco",sep=""),format="GTiff",datatype="INT2U",overwrite=TRUE, 
                options=c("COMPRESS=DEFLATE", "PREDICTOR=2"))
    unlink(paste0(source.dir,'/data/models/',myRespName),rec=TRUE)    
    
  }else{
    nb.eco=length(unique(my.sp.dist[,3]))
    nb.geo=nrow(my.sp.xy)
    
    if(!dir.exists(paste(my.dir,"results/Sp_excluded",sep=""))){
      dir.create(paste(my.dir,"results/Sp_excluded",sep=""),showWarnings = FALSE)
    }
    save(nb.eco,nb.geo,file=paste(my.dir,"results/Sp_excluded/",my.sp.id,sep=""))
  }
  
  
  
}
if (sp.inidat$modeling=='ESM'){
  
  sp.th<-read.csv(paste0(my.dir,'results/thresholds/',my.sp.id,'_thresholds.csv'))
  sp.w<-read.csv(paste0(my.dir,'results/Eweights/',my.sp.id,'_Eweights.csv'))[,-1]
  ini.ESM<- load(paste(my.dir,"results/ESM/",my.sp.id,"_ESM",sep=""))
  sp.w<-sp.w[retained.sm,]
  
  vari<-unique(sort(ESM.list[retained.sm,]))
  var.import.esm<-as.data.frame(matrix(0,nrow=nrow(ESM.list),ncol=length(vari)))
  colnames(var.import.esm)<-names(env)[vari]
  ks=3 
  min.dist=99
  my.mask.pts<-SpatialPoints(coord=bck[,1:2])
  exclude=500
  geo.exclude=5000
  
  load(paste(my.dir,'Sp/Sp_list/',my.sp.id,sep=""))
  vcol<-19
  to.rm<- which(is.na(match(my.sp.tot[,1],unique(my.sp.dist[,4]))))
  if(length(to.rm)>0){
    my.sp.xy<-my.sp.tot[-to.rm,c(11:12,1,26)]
  }else{
    my.sp.xy<-my.sp.tot[,c(11:12,1,26)]
  }
  
  vld.xy<-my.sp.tot[my.sp.tot[,26]>exclude & my.sp.tot[,26]<geo.exclude & my.sp.tot[,19]>=300
                    & (my.sp.tot[,23]!=3 | my.sp.tot[,23]!=4),c(11,12,1,26)]
  my.sp.xy<-rbind(my.sp.xy,vld.xy) 
  
  
  ##### remove aggregated data
  if(min.dist>0){
    my.sp.xy.dup<-zerodist(SpatialPoints(coord=my.sp.xy[,1:2]),zero=min.dist)
    prec<-cbind(my.sp.xy[my.sp.xy.dup[,1],4],my.sp.xy[my.sp.xy.dup[,2],4])
    to.rm<-unique(apply(cbind(my.sp.xy.dup,prec),1,dsg.accurate))
    
    if(length(to.rm)>0){
      my.sp.xy<-my.sp.xy[-to.rm,]
      my.sp.dist<-my.sp.dist[which(!is.na(match(my.sp.dist[,4],my.sp.xy[,3]))),]
    }
    
  }
  
  #if(length(unique(my.sp.dist[,4])) < 8 | nrow(my.sp.xy)<8 ){ ### stop the process in case of too few datapoints
  # nb.eco=length(unique(my.sp.dist[,3]))
  #nb.geo=nrow(my.sp.xy)
  #dir.create(paste(my.dir,"results/Sp_excluded",sep=""),showWarnings = FALSE)
  #save(nb.eco,nb.geo,file=paste(my.dir,"results/Sp_excluded/",my.sp.id,sep=""))
  #} else{
  
  my.abs<-abs.rand[,2:3]
  colnames(my.abs)<-colnames(my.sp.xy)[1:2]
  knn.tot<-log(round(knnx.dist(my.sp.xy[,1:2],coordinates(my.mask.pts),k=2)/100)[,1]+1)
  Ndist<-rasterize(my.mask.pts,env,field=knn.tot)
  env<-stack(env[[vari]],Ndist)
  Ndist<-c(log(round(knn.dist(my.sp.xy[,1:2],k=2)/100)[,1]+1),
           log(round(knnx.dist(my.sp.xy[,1:2], my.abs,k=2)/100)[,1]+1))
  sp.geo<-c(rep(1,nrow(my.sp.xy)),rep(0,nrow(my.abs)))
  
  data.env<-data.frame(extract(env,rbind(my.sp.dist[,1:2],setNames(my.abs,colnames(my.sp.dist[,1:2])))))
  names(data.env)[ncol(data.env)]<-"Ndist"
  names(env[[dim(env)[3]]])<-names(data.env)[length(data.env)]
  
  data.resp<-c(rep(1,nrow(my.sp.dist)),rep(0,nrow(my.abs)))
  df.tmp.input<-cbind(data.resp,data.env)
  colnames(df.tmp.input)[1]<-"sp"
  
  w2<-tapply(my.sp.dist[,4], as.factor(my.sp.dist[,4]),re.w)
  wi<-match(my.sp.dist[,4],names(w2))
  w3<-w2[wi]
  w.pres<-sum(w3)/100
  w.abs<-nrow(my.abs)/w.pres
  w<-c(w3*round(w.abs),rep(100,nrow(my.abs)))
  w.geo<-c(rep(length(which(sp.geo==0))/length(which(sp.geo==1)),length(which(sp.geo==1))),rep(1,length(which(sp.geo==0))))
  
  m1.geo<-gam(sp.geo ~ s(Ndist, k = 3),family = binomial(link = "logit"), control = gam.control(maxit = 500),w=round(10*w.geo))
  sp.pred.geo<-round(1000*predict(env,m1.geo,type='response'))
  geo.th<-round(sp.th[3,2])
  sp.dist.geo<-round(1000*m1.geo$fitted[1:sum(sp.geo)])
  sp.dist.geo<-sp.dist.geo[which(sp.dist.geo>=geo.th)]

  # the name of studied species
  myRespName <- paste0(strsplit(my.sp.id,'_')[[1]][c(2)],strsplit(my.sp.id,'_')[[1]][c(3)],strsplit(my.sp.id,'_')[[1]][c(1)],'Full')
  
  if (dir.exists(paste0(source.dir,'/data/models/',myRespName))){
    unlink(paste0(source.dir,'/data/models/',myRespName))
  }
  
  
  dir.create(paste0(source.dir,'/data/models/',myRespName),recursive = TRUE)
  setwd(paste0(source.dir,'/data/models/',myRespName))
  
  
  # the presence/absences data for our species 
  var.import.esm<-array(0,dim=c(nrow(ESM.list),ncol(ESM.list),3))
  var.import.esm[,,1]<-ESM.list
  
  for (z in retained.sm){
    print(z)
    retained.pred<-ESM.list[z,]
    myResp <- na.exclude(df.tmp.input[,c(1,1+match(retained.pred,vari))])$sp
    myExpl <- na.exclude(df.tmp.input[,1+match(retained.pred,vari)])
    myXY<-na.exclude(cbind(rbind(my.sp.dist[,1:2],my.abs),df.tmp.input[,1+match(retained.pred,vari)]))[,1:2]
    wz<-na.exclude(cbind(w,df.tmp.input[,1+match(retained.pred,vari)]))[,1]
    
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = myExpl,
                                         resp.name = myRespName,
                                         resp.xy = myXY)
    
    myBiomodOption <- BIOMOD_ModelingOptions(GAM = list(k=ks,control = gam.control(maxit = 500)),
                                             MAXENT.Phillips = list(Threshold=FALSE,
                                                                    path_to_maxent.jar = #'/home/petitpie/tools/maxent.jar',
                                                                      'F:/InfoFlora/tools',
                                                                    betamultiplier =1.5))
    
    myBiomodModelOut <- BIOMOD_Modeling( 
      myBiomodData, 
      models = c('MAXENT.Phillips','GAM','RF'), 
      models.options = myBiomodOption, 
      NbRunEval=1, 
      DataSplit=100, 
      Yweights = wz, 
      VarImport=3,
      models.eval.meth = c('ROC'),
      SaveObj = FALSE,
      rescal.all.models = FALSE,
      do.full.models = FALSE,
      modeling.id = as.character(Index) )#myRespName
    
    Mw<-sp.w[match(z,row.names(sp.w)),1:3]
    Mw[which(is.na(Mw))]<-0
    Mw[Mw[]<0.4]<-0
    var.import<-myBiomodModelOut@variables.importances@val[,,1,1]
    var.import<-round(apply(var.import,1,weighted.mean,Mw,na.rm=TRUE),3)
    var.import.esm[z,,2]<-var.import
    
    
    pred.all<- BIOMOD_Projection(
      modeling.output = myBiomodModelOut,
      new.env = env[[match(retained.pred,vari)]],
      proj.name = 'Full',
      selected.models = 'all',
      binary.meth = NULL,
      compress = 'xz',
      build.clamping.mask = TRUE)
    
    my.proj<-stack(paste0(getwd(),'/',myRespName,'/proj_Full/proj_Full_',myRespName,'.grd'))
    
    for (s in 1:dim(my.proj)[3]){
      my.proj[[s]]<-setValues(my.proj[[s]],round(rescale(getValues(my.proj[[s]]),c(0,1000))))
    }
    
    if (z==retained.sm[1]){
      my.proj.E<-w.mean(my.proj,Mw)
      my.proj.E<-round(my.proj.E/10)
      
    }else{
      my.proj.E<-stack(my.proj.E,round(w.mean(my.proj, Mw)/10))
    }
  }
  
  var.import.esm[retained.sm,,3]<-cbind(sp.w[,4],sp.w[,4])
  my.proj.esm<-round(10*w.mean(my.proj.E,sp.w[,4]))
  var.import.synth<-rbind(var.import.esm[,1,],var.import.esm[,2,])
  var.import.synth<-var.import.w(var.import.synth)
  var.import.synth[,1]<-var.list[var.import.synth[,1]]
  
  eco.th<-round(sp.th[3,3])
  sp.pred.eco<-read.csv(paste0(my.dir,'results/innerness_eco/',my.sp.id,'_innerness_eco.csv'),h=T)[,2]
  sp.dist.eco<-sp.pred.eco[which(sp.pred.eco>=eco.th)]

  if (!dir.exists(paste0(my.dir,"/results/varImport"))){
    dir.create(paste0(my.dir,"/results/varImport"),recursive = TRUE)
  }
  write.csv(var.import.synth,paste(my.dir,"results/varImport/",my.sp.id,"_varImport.csv",sep=""),quote=FALSE,
            row.names=FALSE)
  

  if(!dir.exists(paste0(my.dir,"results/rawPred_geo"))){
    dir.create(paste0(my.dir,"results/rawPred_geo"))
  }
  if(!dir.exists(paste0(my.dir,"results/rawPred_eco"))){
    dir.create(paste0(my.dir,"results/rawPred_eco"))
  }

  my.proj.esm<-round(my.proj.esm/10)
  
  writeRaster(sp.pred.geo,filename=paste(my.dir,"results/rawPred_geo/",my.sp.id,"_geo",sep=""),format="GTiff",datatype="INT2U",overwrite=TRUE, 
              options=c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  writeRaster(my.proj.esm,filename=paste(my.dir,"results/rawPred_eco/",my.sp.id,"_eco",sep=""),format="GTiff",datatype="INT2U",overwrite=TRUE, 
              options=c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  unlink(paste0(source.dir,'/data/models/',myRespName),rec=TRUE)    
  
}


