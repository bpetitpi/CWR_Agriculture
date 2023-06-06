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

source.dir<-"F:/InfoFlora/InfoFlora_EvaluationTool/BP/V2" #define the directory of the generic structure
my.dir<-paste(source.dir,"/data/",sep="")
my.dir.fct<-paste(source.dir,"/CodesR/",sep="")
source(paste(my.dir.fct,"functions_pred_models.R",sep=""))

k<- 1 # index of the species to run model

source.dir<-"F:/InfoFlora/InfoFlora_EvaluationTool/BP/V2" #define the directory of the generic structure
my.dir<-paste(source.dir,"/data/",sep="")
my.dir.fct<-paste(source.dir,"/CodesR/",sep="")
source(paste(my.dir.fct,"functions_pred_work.R",sep=""))
source(paste(my.dir.fct,"functions_evaluation.R",sep=""))
source(paste(my.dir.fct,"functions_innerness.r",sep=""))

env<-stack(
  paste(my.dir,"env/env100.tif",sep='')
)

my.sp.list<-list.files(paste(my.dir,'sp/Sp_list/',sep="")) # list of the species 
#sp.eval<-list.files(paste(my.dir,'results/evaluations/',sep=""))
#my.sp.list<-my.sp.list[which(!is.na((charmatch(my.sp.list,sp.eval))))]
#length(my.sp.list)
inidat.model<-read.table(file=paste0(my.dir,'inidat/model_inidat.txt'),header=TRUE,sep='\t')
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

### serialized version : choose this option if you have a low number of cores in your computer

my.sp.id<-my.sp.list[k]
print(my.sp.id)
sp.inidat<-inidat.model[match(my.sp.id,inidat.model[,1]),]


vari<-which(sp.inidat[2:34]==TRUE)
ks=3 # degree of smoothing
min.dist=99
my.mask.pts<-SpatialPoints(coord=bck[,1:2])
exclude=500
geo.exclude=5000

load(paste(my.dir,'sp/Sp_list/',my.sp.id,sep=""))
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

####delimit prediction data (for the whole set of occurences)
my.sp.poly<-try(do.call(bind,lapply(my.sp.tot[,28],readWKT)), silent=T)
if (attributes(my.sp.poly)$class=='try-error'){
  to.rm<-c()
  for (eri in 1 : nrow(my.sp.tot)){
    eri.res<-try(readWKT(my.sp.tot[eri,28]), silent=T)
    if (attributes(eri.res)$class=="try-error"){
      to.rm<-c(to.rm,eri) 
    }
  }
  my.sp.tot<-my.sp.tot[-to.rm,]
  my.sp.poly<-do.call(bind,lapply(my.sp.tot[,28],readWKT))
}
crs(my.sp.poly)<-crs(my.mask.pts)
my.sp.pts<-my.mask.pts %over% my.sp.poly
pts.extract<-my.mask.pts[which(!is.na(my.sp.pts)),]
my.data.pts<-extract(env[[vari]],pts.extract)
my.data.pts<-cbind(my.sp.tot[na.exclude(my.sp.pts),1],coordinates(pts.extract),my.data.pts)
colnames(my.data.pts)[1]<-"obs.id"
id.miss.pts<-which(is.na(match(my.sp.tot[,1],unique(my.data.pts[,1])))) #Add points for smaller polygons
my.data.pts.miss<-extract(env[[vari]],my.sp.tot[id.miss.pts,11:12])
my.data.pts.miss<-cbind(my.sp.tot[id.miss.pts,c(1,11,12)],my.data.pts.miss)
colnames(my.data.pts.miss)<-colnames(my.data.pts)
my.data.pts<-na.exclude(rbind(my.data.pts,my.data.pts.miss))

my.abs.strat<-abs.rand[,2:3]
my.abs.rand<-abs.rand[,2:3]
my.abs.strat<-na.exclude(cbind(rep(0,nrow(my.abs.strat)),my.abs.strat,extract(env[[vari]],my.abs.strat)))
colnames(my.abs.strat)<-colnames(my.data.pts)
my.abs.rand<-na.exclude(cbind(rep(0,nrow(my.abs.rand)),my.abs.rand, extract(env[[vari]],my.abs.rand)))
colnames(my.abs.rand)<-colnames(my.data.pts)

my.sp.dist<-na.exclude(cbind(my.sp.dist[,c(4,1,2)],extract(env[[vari]],my.sp.dist[,1:2]),my.sp.dist[,5]))    
colnames(my.sp.dist)[1:(3+length(vari))]<-colnames(my.data.pts)
colnames(my.sp.dist)[ncol(my.sp.dist)]<-'w'

### EM modeling -> if there are more than 10 obervation per predictor

if (sp.inidat$modeling=='EM'){
  CV = 3
  #####Generate the CV partitions
  pres.id<-my.sp.xy[,3]
  pres.id<-sample(pres.id,length(pres.id))
  length.part.pres<-length(pres.id)%/%CV
  cv.part.pres.geo<-split(pres.id,c(rep(1:CV,each=length.part.pres),rep(CV,length(pres.id)%%CV)))
  
  pres.id<-unique(my.sp.dist[,1])
  pres.id<-sample(pres.id,length(pres.id))
  length.part.pres<-length(pres.id)%/%CV
  cv.part.pres.eco<-split(pres.id,c(rep(1:CV,each=length.part.pres),rep(CV,length(pres.id)%%CV)))
  
  row.abs.strat<-sample(1:nrow(my.abs.strat))
  length.part.abs<-length(row.abs.strat)%/% CV
  cv.part.abs.strat<-split(row.abs.strat,c(rep(1:CV,each=length.part.abs),rep(CV,length(row.abs.strat)%%CV)))
  
  row.abs.rand<-sample(1:nrow((my.abs.rand)))
  length.part.abs<-length(row.abs.rand)%/% CV
  cv.part.abs.rand<-split(row.abs.rand,c(rep(1:CV,each=length.part.abs),rep(CV,length(row.abs.rand)%%CV)))
  
  cv.part<-mapply(FUN = cv.format, cv.part.pres.geo, cv.part.pres.eco,cv.part.abs.strat,cv.part.abs.rand, SIMPLIFY=FALSE)
  cv.part.tot<-cv.part
  
  dir.create(paste0(source.dir,'/data/models/',my.sp.id),recursive = TRUE)
  setwd(paste0(source.dir,'/data/models/',my.sp.id))
  
  #      my.th<-array(NA,dim = c(6,2,CV))
  
  pred.pres.eco.cal.cv<-c()
  pred.abs.eco.cal.cv<-c()
  pred.pres.geo.cal.cv<-c()
  pred.abs.geo.cal.cv<-c()
  pred.pres.eco.eval.cv<-c()
  pred.abs.eco.eval.cv<-c()
  pred.pres.geo.eval.cv<-c()
  pred.abs.geo.eval.cv<-c()
  
  pred.remain.geo.cv<-array(NA,dim = c(nrow(my.sp.tot),2,CV))
  pred.remain.eco.cv<-array(NA,dim = c(nrow(my.sp.tot),4,CV))
  for (k in 1:CV){
    
    print(paste0('CV#',k))
    cv.part<-cv.part.tot[[k]]
    
    
    
    ##### Split into evaluation, calibration and projection dataset
    to.pck.geo<-which(!is.na(match(my.sp.xy[,3],cv.part[[1]])))
    pres.data.eval.geo<-my.sp.xy[to.pck.geo,]
    pres.data.cal.geo<-my.sp.xy[-to.pck.geo,]
    
    to.pck.eco<-which(!is.na(match(my.sp.dist[,1],cv.part[[2]])))
    pres.data.cal.eco<-my.sp.dist[-to.pck.eco,]
    pres.data.eval.eco<-my.sp.dist[to.pck.eco,]
    
    abs.data.cal<-my.abs.strat[-cv.part[[3]],]
    abs.data.eval<-my.abs.rand[cv.part[[4]],]
    
    my.data.cal.eco<-as.data.frame(rbind(pres.data.cal.eco[,-ncol(pres.data.cal.eco)],abs.data.cal))
    my.data.eval.eco<-as.data.frame(rbind(pres.data.eval.eco[,-ncol(pres.data.eval.eco)],abs.data.eval))
    if (length(which(colnames(my.data.cal.eco)=='soil')==1)){
      my.data.cal.eco$soil<-as.factor(my.data.cal.eco$soil)
      my.data.eval.eco$soil<-as.factor(my.data.eval.eco$soil)
    }
    if (length(which(colnames(my.data.cal.eco)=='vgh')==1)){
      my.data.cal.eco$vgh<-as.factor(my.data.cal.eco$vgh)
      my.data.eval.eco$vgh<-as.factor(my.data.eval.eco$vgh)
    }
    
    my.data.cal.eco<-na.exclude(my.data.cal.eco)
    my.data.eval.eco<-na.exclude(my.data.eval.eco)
    ##### Make the remaining data independent from the calibration partition
    to.rm.geo<-which(!is.na(match(my.sp.tot[,1],pres.data.cal.geo[,3]))) #Select occurence not in the calibration data set
    my.data.remain.geo<-my.sp.tot[-to.rm.geo,c(11,12,1,26)]
    
    
    my.sp.xy.eco<-my.sp.tot[which(!is.na(match(my.sp.tot[,1],pres.data.cal.eco[,1]))),c(11,12,1,26)]
    to.rm.eco<-which(is.na(match(my.sp.tot[,1],my.sp.xy.eco[,3])))
    my.data.remain.eco<-my.data.pts[which(!is.na(match(my.data.pts[,1],my.sp.tot[to.rm.eco,1]))),]
    
    if (length(which(colnames(my.data.cal.eco)=='soil')==1)){
      my.data.remain.eco$soil<-as.factor(my.data.remain.eco$soil)
    }
    if (length(which(colnames(my.data.cal.eco)=='vgh')==1)){
      my.data.remain.eco$vgh<-as.factor(my.data.remain.eco$vgh)
    }
    
    
    ##### Build weight for the calibration of the model
    w.pres<-length(unique(pres.data.cal.eco[,1]))
    w.abs<-nrow(abs.data.cal)/w.pres
    w2<-tapply(pres.data.cal.eco[,1], as.factor(pres.data.cal.eco[,1]),re.w)
    wi<-match(pres.data.cal.eco[,1],names(w2))
    w3<-w2[wi]
    w.eco<-round(c(w3*w.abs,rep(100,nrow(abs.data.cal))))
    w.geo<-round(100*c(rep(nrow(abs.data.cal)/nrow(pres.data.cal.geo),nrow(pres.data.cal.geo)),rep(1,nrow(abs.data.cal))))
    
    
    ##### Nearest neighbour, with the calibration dataset as a reference
    knn.cal.pres<-log(round(knn.dist(pres.data.cal.geo[,1:2],k=2)/100)[,1]+1)
    knn.cal.abs<-log(round(knnx.dist(pres.data.cal.geo[,1:2],abs.data.cal[,2:3],k=2)/100)[,1]+1)
    
    knn.cal<-c(knn.cal.pres,knn.cal.abs)
    
    knn.eval.pres<-log(round(knnx.dist(pres.data.cal.geo[,1:2],pres.data.eval.geo[,1:2],k=2)/100)[,1]+1)
    knn.eval.abs<-log(round(knnx.dist(pres.data.cal.geo[,1:2],abs.data.eval[,2:3],k=2)/100)[,1]+1)
    knn.eval<-c(knn.eval.pres,knn.eval.abs)
    knn.remain<-log(round(knnx.dist(pres.data.cal.geo[,1:2], my.data.remain.geo[,1:2],k=2)/100)[,1]+1)
    my.data.remain.geo<-as.data.frame(cbind(my.data.remain.geo,knn.remain))
    colnames(my.data.remain.geo)[length(colnames(my.data.remain.geo))]<-"Ndist"
    
    cal.geo<-as.data.frame(cbind(c(rep(1,nrow(pres.data.cal.geo)),rep(0,nrow(abs.data.cal))),knn.cal))
    eval.geo<-as.data.frame(cbind(c(rep(1,nrow(pres.data.eval.geo)),rep(0,nrow(abs.data.eval))),knn.eval))
    my.data.remain<-as.data.frame(c(knn.remain))
    colnames(my.data.remain)[ncol(my.data.remain)]<-"Ndist"
    names(cal.geo)<-c("sp","Ndist")
    names(eval.geo)<-c("sp","Ndist")
    
    ##### BIOMOD
    
    # the name of studied species
    myRespName <- paste0('sp_',strsplit(my.sp.id,'_')[[1]][c(1)])
    
    # the presence/absences data for our species 
    myResp <- c(rep(1,nrow(pres.data.cal.eco)),rep(0,nrow(abs.data.cal)))
    
    myExpl <- my.data.cal.eco[,-c(1:3)]
    
    w.eco<-na.exclude(cbind(myExpl,w.eco))[,ncol(myExpl)+1]
    
    myXY<-my.data.cal.eco[,2:3]
    
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = myExpl,
                                         resp.name = myRespName,
                                         resp.xy = myXY)
    
    myBiomodOption <- BIOMOD_ModelingOptions(GAM = list(k=ks,control = gam.control(maxit = 500)),
                                             MAXENT = list(Threshold=FALSE,
                                                           path_to_maxent.jar = 'F:/InfoFlora/tools',
                                                           betamultiplier =1.5))
    
    myBiomodModelOut <- BIOMOD_Modeling( 
      myBiomodData, 
      models = c('MAXENT','GAM','RF'), 
      bm.options = myBiomodOption, 
      nb.rep=1, 
      data.split.perc=100, 
      weights  = w.eco, 
      var.import=0,
      metric.eval = c('ROC'),
      save.output = FALSE,
      scale.models = FALSE,
      do.full.models = TRUE,
      modeling.id = paste(myRespName,"FirstModeling",sep=""))
    
    my.pred.cal<-get_predictions(myBiomodModelOut, model.as.col = TRUE)
    pred.pres.cal<-apply(my.pred.cal[1:nrow(pres.data.cal.eco),],2,pred.max,my.data.cal.eco$obs.id[1:nrow(pres.data.cal.eco)])
    pred.abs.cal<-my.pred.cal[which(my.data.cal.eco$obs.id==0),]
    pred.pres.eco.cal.cv<-rbind(pred.pres.eco.cal.cv,pred.pres.cal)
    pred.abs.eco.cal.cv<-rbind(pred.abs.eco.cal.cv,pred.abs.cal)
    #      pred.cal<-rbind(pred.pres.cal,pred.abs.cal) 
    #      pres.abs.cal<-c(rep(1,length(unique(my.data.cal.eco$obs.id))-1),rep(0,length(which(my.data.cal.eco$obs.id==0))))
    
    pred.eco<- BIOMOD_Projection(
      myBiomodModelOut,
      modeling.output = myBiomodModelOut,
      new.env = my.data.eval.eco[,-c(1:3)],
      proj.name = 'eval',
      models.chosen = 'all',
      metric.binary = NULL,
      compress = TRUE,
      build.clamping.mask = F,
      output.format = '.RData')
    
    my.pred.eco <- get_predictions(pred.eco, model.as.col = TRUE)
    ##### Calibrate geo GAMs
    m1.geo<-gam(sp ~ s(Ndist, k = ks),data=cal.geo,family = binomial(link = "logit"), control = gam.control(maxit = 500), weight=w.geo)
    
    ##### Predictions on the evaluation dataset
    pred.geo<-predict.gam(m1.geo, newdata=eval.geo, type="response")
    pred.pres.eco<-apply(my.pred.eco[1:nrow(pres.data.eval.eco),],2,pred.max,my.data.eval.eco$obs.id[1:nrow(pres.data.eval.eco)])
    pred.abs.eco<-my.pred.eco[which(my.data.eval.eco$obs.id==0),]
    pred.pres.eco.eval.cv<-rbind(pred.pres.eco.eval.cv,pred.pres.eco)
    pred.abs.eco.eval.cv<-rbind(pred.abs.eco.eval.cv,pred.abs.eco)
    pred.pres.geo.eval.cv<-c(pred.pres.geo.eval.cv,pred.geo[which(eval.geo[,1]==1)])
    pred.abs.geo.eval.cv<-c(pred.abs.geo.eval.cv,pred.geo[which(eval.geo[,1]==0)])
    
  ##### Evaluations with TSS, AUC and Boyce
  pres.abs.geo<-c(rep(1, length(pred.pres.geo.eval.cv)), rep(0,length(pred.abs.geo.eval.cv)))
  pred.geo<-c(pred.pres.geo.eval.cv,pred.abs.geo.eval.cv)
  pres.abs.eco<-c(rep(1, nrow(pred.pres.eco.eval.cv)), rep(0,nrow(pred.abs.eco.eval.cv)))
  pred.eco<-rbind(pred.pres.eco.eval.cv,pred.abs.eco.eval.cv)
  TSS.th<-0.4
  TSS.geo<-unlist(KappaRepet(pres.abs.geo,pred.geo, TSS=T))
  TSS.eco.m<-apply(pred.eco,2,TSS.Kappa,pres.abs.eco,TSS=T)
  TSS.w<-TSS.eco.m[1,]*(TSS.eco.m[1,]>=TSS.th)
  if (sum(TSS.w)==0){
    TSS.w<-(TSS.eco.m[1,])
    if (length(which(TSS.w<0))>1){
      TSS.w[which(TSS.w<0)]<-0
      if (sum(TSS.w)==0){
        TSS.w<-c(rep(1,length(TSS.w)))
      }
    }
  }
  
  Se.th<-0.4
  Se.geo<-2*((round(TSS.geo[4])/100)-0.5)
  Se.eco.m<-2*((round(TSS.eco.m[4,1:3])/100)-0.5)
  Se.w<-Se.eco.m*(Se.eco.m>=TSS.th)
  if (sum(Se.w)==0){
    Se.w<-(Se.eco.m[1,])
    if (length(which(Se.w<0))>0){
      Se.w[which(Se.w<0)]<-0
      if (sum(Se.w)==0){
        Se.w<-c(rep(1,length(Se.w)))
      }
    }
  }

  AUC.th<-0.4
  AUC.geo<-somers2(pred.geo,pres.abs.geo)[2]
  AUC.eco.m<-apply(pred.eco,2,somers2,pres.abs.eco)[2,]
  AUC.w<-AUC.eco.m*(AUC.eco.m>=AUC.th)
  if (sum(AUC.w)==0){
    AUC.w<-AUC.eco.m
    if (length(which(AUC.w<0))>0){
      AUC.w[which(AUC.w<0)]<-0
      if (sum(AUC.w)==0){
        AUC.w<-c(rep(1,length(AUC.w)))
      }
    }
  }

  B.th<-0.4
  B.geo<-ecospat.boyce(pred.geo, pred.geo[which(pres.abs.geo==1)], PEplot = F,cor.method = "kendall")$correlation
  B.eco.m<-apply(pred.eco,2,Boyce,which(pres.abs.eco==1))
  if(length(which(is.na(B.eco.m)))>0){B.eco.m[which(is.na(B.eco.m))]<-0}
  B.w<-B.eco.m*(B.eco.m>=B.th)
  if (sum(B.w)==0){
    B.w<-B.eco.m
    if (length(which(B.w<0))>0){
      B.w[which(B.w<0)]<-0
      if (sum(B.w)==0){
        B.w<-c(rep(1,length(B.w)))
      }
    }
  }

  EM.w<-round(apply(rbind(TSS.w,Se.w,AUC.w,B.w),2,mean),2)
  EM.pred<-apply(pred.eco,1,weighted.mean,EM.w)
  EM.TSS<-unlist(KappaRepet(pres.abs.eco,EM.pred, TSS=T))[1]
  EM.Se<-2*((round(unlist(KappaRepet(pres.abs.eco,EM.pred, TSS=T))[4])/100)-0.5)
  EM.AUC<-somers2(EM.pred,pres.abs.eco)[2]
  EM.B<-ecospat.boyce(EM.pred, EM.pred[which(pres.abs.eco==1)], PEplot = F,cor.method = "kendall")$correlation
  
  my.evali<-round(rbind(c(TSS.geo[1], TSS.eco.m[1,],EM.TSS),c(Se.geo,Se.eco.m,EM.Se),c(AUC.geo,AUC.eco.m,EM.AUC),
                        c(B.geo,B.eco.m,EM.B)),2)
  colnames(my.evali)<-c("GAM.geo","ME",'GAM','RF','EM')
  row.names(my.evali)<-c("TSS","Se","AUC","B")
  Eval.consensus<-apply(my.evali[,c(2:ncol(my.evali))],2,mean)
  best.model<-which(Eval.consensus==max(Eval.consensus))
  
  pred.eco.final<-cbind(pred.remain.eco.cv[,1,1],apply(pred.remain.eco.cv[,2:4,],c(1,2),mean,na.rm=T))
  pred.eco.final<-round(cbind(pred.eco.final,apply(pred.eco.final[,2:4],1,weighted.mean,EM.w,na.rm=T)))
  colnames(pred.eco.final)<-c('obs_id','ME','GAM','RF','EM')
  pred.geo.final<-round(cbind(pred.remain.geo.cv[,1,1],apply(pred.remain.geo.cv[,2,],1,mean,na.rm=T)))
  colnames(pred.geo.final)<-c('obs_id','geo')
  
  ####Thresholds
  
  th.geo<-round(c(quantile(1000*pred.geo[which(pres.abs.geo==1)],c(0,0.01,0.05,0.1,0.2)),TSS.Kappa(1000*pred.geo,pres.abs.geo,TSS=T)[2]))
  th.eco<-round(c(quantile(EM.pred[which(pres.abs.eco==1)],c(0,0.01,0.05,0.1,0.2)),TSS.Kappa(EM.pred,pres.abs.eco,TSS=T)[2]))
  
  colnames(results.geo)<-c("id_obs","pred","100%","99%","95%","90%","80%","TSS_CutOff")
  
  colnames(results.eco)<-c("id_obs","pred_maxent", 'pred_gam','pred_rf','pred_EM',"100%","99%","95%","90%","80%","TSS_CutOff")
  

  names(EM.w)<-c('wMAXENT','wGAM','wRF')
  
  colnames(my.evali)<-c('GAMgeo','ME','GAM','RF','EM')
  row.names(my.evali)<-c('TSS','Se_std','AUC_std','B')
  
  my.th<-cbind(th.geo,th.eco)
  colnames(my.th)<-c('GEOth','ECOth')
  row.names(my.th)<-c('0%','1%','5%','10%','20%','TSS_Optimised')

  if (!dir.exists(paste0(my.dir,"/results/evaluations"))){
    dir.create(paste0(my.dir,"/results/evaluations"),recursive = TRUE)
  }
  write.csv(my.evali,paste(my.dir,"results/evaluations/",my.sp.id,"_evaluation.csv",sep=""),quote=FALSE,
            row.names=TRUE)
  
  if (!dir.exists(paste0(my.dir,"/results/thresholds"))){
    dir.create(paste0(my.dir,"/results/thresholds"),recursive = TRUE)
  }
  write.csv(my.th,paste(my.dir,"results/thresholds/",my.sp.id,"_thresholds.csv",sep=""),quote=FALSE,
            row.names=TRUE)
  
  if (!dir.exists(paste0(my.dir,"/results/Eweights"))){
    dir.create(paste0(my.dir,"/results/Eweights"),recursive = TRUE)
  }
  write.csv(EM.w,paste(my.dir,"results/Eweights/",my.sp.id,"_Eweights.csv",sep=""),quote=FALSE,
            row.names=TRUE)
  
  # unlink(getwd(),rec=TRUE)
  
  
  #ESM modeling
  
}else{   #ESM modeling 
  CV = 5
  #####Generate the CV partitions
  pres.id<-my.sp.xy[,3]
  pres.id<-sample(pres.id,length(pres.id))
  length.part.pres<-length(pres.id)%/%CV
  cv.part.pres.geo<-split(pres.id,c(rep(1:CV,each=length.part.pres),rep(CV,length(pres.id)%%CV)))
  
  pres.id<-unique(my.sp.dist[,1])
  pres.id<-sample(pres.id,length(pres.id))
  length.part.pres<-length(pres.id)%/%CV
  cv.part.pres.eco<-split(pres.id,c(rep(1:CV,each=length.part.pres),rep(CV,length(pres.id)%%CV)))
  
  row.abs.strat<-sample(1:nrow(my.abs.strat))
  length.part.abs<-length(row.abs.strat)%/% CV
  cv.part.abs.strat<-split(row.abs.strat,c(rep(1:CV,each=length.part.abs),rep(CV,length(row.abs.strat)%%CV)))
  
  row.abs.rand<-sample(1:nrow((my.abs.rand)))
  length.part.abs<-length(row.abs.rand)%/% CV
  cv.part.abs.rand<-split(row.abs.rand,c(rep(1:CV,each=length.part.abs),rep(CV,length(row.abs.rand)%%CV)))
  
  cv.part<-mapply(FUN = cv.format, cv.part.pres.geo, cv.part.pres.eco,cv.part.abs.strat,cv.part.abs.rand, SIMPLIFY=FALSE)
  cv.part.tot<-cv.part
  
  dir.create(paste0(source.dir,'/data/models/',my.sp.id),recursive = TRUE)
  setwd(paste0(source.dir,'/data/models/',my.sp.id))
  
  
  ####model geo
  pred.pres.geo.cal.cv<-c()
  pred.abs.geo.cal.cv<-c()
  pred.pres.geo.eval.cv<-c()
  pred.abs.geo.eval.cv<-c()
  
  pred.remain.geo.cv<-array(NA,dim = c(nrow(my.sp.tot),2,CV))
  
  
  for (k in 1:CV){
    print (paste0('CV #',k))
    cv.part<-cv.part.tot[[k]]
    
    
    ##### Split into evaluation, calibration and projection dataset
    to.pck.geo<-which(!is.na(match(my.sp.xy[,3],cv.part[[1]])))
    pres.data.eval.geo<-my.sp.xy[to.pck.geo,]
    pres.data.cal.geo<-my.sp.xy[-to.pck.geo,]
    
    abs.data.cal<-my.abs.strat[-cv.part[[3]],]
    abs.data.eval<-my.abs.rand[cv.part[[4]],]
    
    to.rm.geo<-which(!is.na(match(my.sp.tot[,1],pres.data.cal.geo[,3]))) #Select occurence not in the calibration data set
    my.data.remain.geo<-my.sp.tot[-to.rm.geo,c(11,12,1,26)]
    
    w.geo<-round(100*c(rep(nrow(abs.data.cal)/nrow(pres.data.cal.geo),nrow(pres.data.cal.geo)),rep(1,nrow(abs.data.cal))))
    
    ##### Nearest neighbour, with the calibration dataset as a reference
    knn.cal.pres<-log(round(knn.dist(pres.data.cal.geo[,1:2],k=2)/100)[,1]+1)
    knn.cal.abs<-log(round(knnx.dist(pres.data.cal.geo[,1:2],abs.data.cal[,2:3],k=2)/100)[,1]+1)
    
    knn.cal<-c(knn.cal.pres,knn.cal.abs)
    
    knn.eval.pres<-log(round(knnx.dist(pres.data.cal.geo[,1:2],pres.data.eval.geo[,1:2],k=2)/100)[,1]+1)
    knn.eval.abs<-log(round(knnx.dist(pres.data.cal.geo[,1:2],abs.data.eval[,2:3],k=2)/100)[,1]+1)
    knn.eval<-c(knn.eval.pres,knn.eval.abs)
    knn.remain<-log(round(knnx.dist(pres.data.cal.geo[,1:2], my.data.remain.geo[,1:2],k=2)/100)[,1]+1)
    my.data.remain.geo<-as.data.frame(cbind(my.data.remain.geo,knn.remain))
    colnames(my.data.remain.geo)[length(colnames(my.data.remain.geo))]<-"Ndist"
    
    cal.geo<-as.data.frame(cbind(c(rep(1,nrow(pres.data.cal.geo)),rep(0,nrow(abs.data.cal))),knn.cal))
    eval.geo<-as.data.frame(cbind(c(rep(1,nrow(pres.data.eval.geo)),rep(0,nrow(abs.data.eval))),knn.eval))
    my.data.remain<-as.data.frame(c(knn.remain))
    colnames(my.data.remain)[ncol(my.data.remain)]<-"Ndist"
    names(cal.geo)<-c("sp","Ndist")
    names(eval.geo)<-c("sp","Ndist")
    
    ##### Calibrate geo GAMs
    m1.geo<-gam(sp ~ s(Ndist, k = ks),data=cal.geo,family = binomial(link = "logit"), control = gam.control(maxit = 500), weight=w.geo)
    
    ##### Predictions on the evaluation dataset
    pred.geo<-predict.gam(m1.geo, newdata=eval.geo, type="response")
    pred.pres.geo.eval.cv<-c(pred.pres.geo.eval.cv,pred.geo[which(eval.geo[,1]==1)])
    pred.abs.geo.eval.cv<-c(pred.abs.geo.eval.cv,pred.geo[which(eval.geo[,1]==0)])
    
    ####Predict on dependent data
    pred.remain.geo<-predict.gam(m1.geo, newdata=my.data.remain.geo, type="response")  
    
    pred.remain.geo.cv[,1,k]<-my.sp.tot$obs_id
    pred.remain.geo.cv[match(my.data.remain.geo$obs_id,my.sp.tot$obs_id),2,k]<-round(1000*pred.remain.geo)
    
  }
  
  pres.abs.geo<-c(rep(1, length(pred.pres.geo.eval.cv)), rep(0,length(pred.abs.geo.eval.cv)))
  pred.geo<-c(pred.pres.geo.eval.cv,pred.abs.geo.eval.cv)
  TSS.geo<-unlist(KappaRepet(pres.abs.geo,pred.geo, TSS=T))
  
  Se.geo<-2*((round(TSS.geo[4])/100)-0.5)
  
  AUC.geo<-somers2(pred.geo,pres.abs.geo)[2]
  
  B.geo<-ecospat.boyce(pred.geo, pred.geo[which(pres.abs.geo==1)], PEplot = F,cor.method = "kendall")$correlation
  
  pred.geo.final<-round(cbind(pred.remain.geo.cv[,1,1],apply(pred.remain.geo.cv[,2,],1,mean,na.rm=T)))
  colnames(pred.geo.final)<-c('obs_id','geo')
  
  ####Thresholds
  
  th.geo<-round(c(quantile(1000*pred.geo[which(pres.abs.geo==1)],c(0,0.01,0.05,0.1,0.2)),TSS.Kappa(1000*pred.geo,pres.abs.geo,TSS=T)[2]))
  
  colnames(results.geo)<-c("id_obs","pred","100%","99%","95%","90%","80%","TSS_CutOff")
  

  
  #### model ECO
  
  ESM.list<-t(combn(vari,2))
  
  ESM.eval<-array(NA,dim = c(4,5,nrow(ESM.list)))
  ESM.pred.eval<-array(NA,dim = c(sum(length(pres.id),nrow(my.abs.rand)),4,nrow(ESM.list)))
  ESM.pred.final<-array(NA,dim = c(nrow(my.sp.tot),5,nrow(ESM.list)))
  ESM.w<-c()
  EM.th<-0.4
  
  for (z in 1:nrow(ESM.list)){
    vari<-ESM.list[z,]
    pred.pres.eco.cal.cv<-c()
    pred.abs.eco.cal.cv<-c()
    pred.pres.eco.eval.cv<-c()
    pred.abs.eco.eval.cv<-c()
    
    pred.remain.eco.cv<-array(NA,dim = c(nrow(my.sp.tot),4,CV))
    
    
    for (k in 1:CV){
      
      print(paste0('SM#',z,' CV#',k))
      cv.part<-cv.part.tot[[k]]
      
      
      to.pck.eco<-which(!is.na(match(my.sp.dist[,1],cv.part[[2]])))
      pres.data.cal.eco<-my.sp.dist[-to.pck.eco,]
      pres.data.eval.eco<-my.sp.dist[to.pck.eco,]
      
      abs.data.cal<-my.abs.strat[-cv.part[[3]],]
      abs.data.eval<-my.abs.rand[cv.part[[4]],]
      
      my.data.cal.eco<-as.data.frame(rbind(pres.data.cal.eco[,-ncol(pres.data.cal.eco)],abs.data.cal))
      my.data.eval.eco<-as.data.frame(rbind(pres.data.eval.eco[,-ncol(pres.data.eval.eco)],abs.data.eval))
      
      my.data.cal.eco<-na.exclude(my.data.cal.eco)
      my.data.eval.eco<-na.exclude(my.data.eval.eco)
      
      ##### Make the remaining data independent from the calibration partition
      
      my.sp.xy.eco<-my.sp.tot[which(!is.na(match(my.sp.tot[,1],pres.data.cal.eco[,1]))),c(11,12,1,26)]
      to.rm.eco<-which(is.na(match(my.sp.tot[,1],my.sp.xy.eco[,3])))
      my.data.remain.eco<-my.data.pts[which(!is.na(match(my.data.pts[,1],my.sp.tot[to.rm.eco,1]))),]
      
      
      ##### Build weight for the calibration of the model
      w.pres<-length(unique(pres.data.cal.eco[,1]))
      w.abs<-nrow(abs.data.cal)/w.pres
      w2<-tapply(pres.data.cal.eco[,1], as.factor(pres.data.cal.eco[,1]),re.w)
      wi<-match(pres.data.cal.eco[,1],names(w2))
      w3<-w2[wi]
      w.eco<-round(c(w3*w.abs,rep(100,nrow(abs.data.cal))))
      
      
      ##### BIOMOD
      
      # the name of studied species
      myRespName <- paste0(strsplit(my.sp.id,'_')[[1]][c(2)],strsplit(my.sp.id,'_')[[1]][c(3)],strsplit(my.sp.id,'_')[[1]][c(1)])
      
      # the presence/absences data for our species 
      myResp <- c(rep(1,nrow(pres.data.cal.eco)),rep(0,nrow(abs.data.cal)))
      
      myExpl <- my.data.cal.eco[,-c(1:3)]
      myExpl<-myExpl[,match(names(env)[vari],colnames(myExpl))]
      
      w.eco<-na.exclude(cbind(myExpl,w.eco))[,ncol(myExpl)+1]
      
      myXY<-my.data.cal.eco[,2:3]
      
      myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                           expl.var = myExpl,
                                           resp.name = myRespName,
                                           resp.xy = myXY)
      
      myBiomodOption <- BIOMOD_ModelingOptions(GAM = list(k=ks,control = gam.control(maxit = 500)),
                                               MAXENT.Phillips = list(Threshold=FALSE,
                                                                      path_to_maxent.jar = 'F:/InfoFlora/tools',
                                                                      betamultiplier =1.5))
      
      myBiomodModelOut <- BIOMOD_Modeling( 
        myBiomodData, 
        models = c('MAXENT.Phillips','GAM','RF'), 
        models.options = myBiomodOption, 
        NbRunEval=1, 
        DataSplit=100, 
        Yweights = w.eco, 
        VarImport=0,
        models.eval.meth = c('ROC'),
        SaveObj = FALSE,
        rescal.all.models = F,
        do.full.models = TRUE,
        modeling.id = paste(myRespName,'a',sep=""))
      
      my.pred.cal<-get_predictions(myBiomodModelOut)[,,1,1]
      pred.pres.cal<-apply(my.pred.cal[1:nrow(pres.data.cal.eco),],2,pred.max,my.data.cal.eco$obs.id[1:nrow(pres.data.cal.eco)])
      pred.abs.cal<-my.pred.cal[which(my.data.cal.eco$obs.id==0),]
      pred.pres.eco.cal.cv<-rbind(pred.pres.eco.cal.cv,pred.pres.cal)
      pred.abs.eco.cal.cv<-rbind(pred.abs.eco.cal.cv,pred.abs.cal)
      
      pred.eco<- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = my.data.eval.eco[,-c(1:3)][,match(names(env)[vari],colnames(my.data.eval.eco[,-c(1:3)]))],
        proj.name = 'eval',
        selected.models = 'all',
        binary.meth = NULL,
        compress = 'xz',
        clamping.mask = F,
        output.format = '.RData')
      
      my.pred.eco<-as.data.frame(matrix(data=NA,nrow=nrow(my.data.eval.eco),ncol=3))
      colnames(my.pred.eco)<-c('MAXENT.Phillips','GAM','RF')
      to.load<- get_predictions(pred.eco)[,,1,1]
      my.pred.eco[,match(colnames(to.load),colnames(my.pred.eco))]<-to.load
      
      ##### Predictions on the evaluation dataset
      pred.pres.eco<-apply(my.pred.eco[1:nrow(pres.data.eval.eco),],2,pred.max,my.data.eval.eco$obs.id[1:nrow(pres.data.eval.eco)])
      pred.abs.eco<-my.pred.eco[which(my.data.eval.eco$obs.id==0),]
      pred.pres.eco.eval.cv<-rbind(pred.pres.eco.eval.cv,pred.pres.eco)
      pred.abs.eco.eval.cv<-rbind(pred.abs.eco.eval.cv,pred.abs.eco)
      
      ####Predict on dependent data
      
      pred.remain<- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = my.data.remain.eco[,-c(1:3)][,match(names(env)[vari],colnames(my.data.remain.eco[,-c(1:3)]))],
        proj.name = 'remain',
        selected.models = 'all',
        binary.meth = NULL,
        compress = 'xz',
        clamping.mask = F,
        output.format = '.RData')
      
      to.load<-get_predictions(pred.remain)[,,1,1]
      my.pred.remain.eco <- as.data.frame(matrix(data=NA,nrow=nrow(my.data.remain.eco),ncol=3))
      colnames(my.pred.remain.eco)<-c('MAXENT.Phillips','GAM','RF')
      my.pred.remain.eco[,match(colnames(to.load),colnames(my.pred.remain.eco))]<-to.load
      sp.dist.remain<-apply(my.pred.remain.eco,2,pred.max,my.data.remain.eco$obs.id)
      
      pred.remain.eco.cv[,1,k]<-my.sp.tot$obs_id
      pred.remain.eco.cv[match(row.names(sp.dist.remain),my.sp.tot$obs_id),2:4,k]<-sp.dist.remain
      
    }
    pres.abs.eco<-c(rep(1, nrow(pred.pres.eco.eval.cv)), rep(0,nrow(pred.abs.eco.eval.cv)))
    pred.eco<-rbind(pred.pres.eco.eval.cv,pred.abs.eco.eval.cv)
    TSS.th<-0.4
    TSS.eco.m<-apply(pred.eco,2,TSS.Kappa,pres.abs.eco,TSS=T)
    
    
    Se.th<-0.4
    Se.eco.m<-2*((round(TSS.eco.m[4,1:3])/100)-0.5)
    
    AUC.th<-0.4
    AUC.eco.m<-apply(pred.eco,2,somers2,pres.abs.eco)[2,]  
    
    B.th<-0.4
    B.eco.m<-apply(pred.eco,2,Boyce,which(pres.abs.eco==1))
    
    EM.w<-round(apply(rbind(TSS.eco.m[1,],Se.eco.m,AUC.eco.m,B.eco.m),2,mean),2)
    if (length(which(is.na(EM.w)))>0){EM.w[which(is.na(EM.w))]<-0}
    EM.w[which(EM.w<EM.th)]<-0
    if (length(which(EM.w<=0))==length(EM.w)){
      EM.w<-round(apply(rbind(TSS.eco.m[1,],Se.eco.m,AUC.eco.m,B.eco.m),2,mean),2)
      EM.w<-EM.w+1.001
      if (length(which(is.na(EM.w)))>0){EM.w[which(is.na(EM.w))]<-0}
      EM.w[EM.w[]<=1]<-0
    }
    EM.pred<-apply(pred.eco,1,weighted.mean,EM.w,na.rm=T)
    
    if(length(na.exclude(EM.pred))!=length(EM.pred)){
      EM.TSS<-0
      EM.Se<-0
      EM.AUC<-0
      EM.B<-0
    }else{
      EM.TSS<-unlist(KappaRepet(pres.abs.eco,EM.pred, TSS=T))[1]
      EM.Se<-2*((round(unlist(KappaRepet(pres.abs.eco,EM.pred, TSS=T))[4])/100)-0.5)
      EM.AUC<-somers2(EM.pred,pres.abs.eco)[2]
      EM.B<-ecospat.boyce(EM.pred, EM.pred[which(pres.abs.eco==1)], PEplot = F,cor.method = "kendall")$correlation
    }
    
    my.evali<-round(rbind(c(TSS.geo[1], TSS.eco.m[1,],EM.TSS),c(Se.geo,Se.eco.m,EM.Se),c(AUC.geo,AUC.eco.m,EM.AUC),
                          c(B.geo,B.eco.m,EM.B)),2)
    colnames(my.evali)<-c("GAM.geo","ME",'GAM','RF','EM')
    row.names(my.evali)<-c("TSS","Se","AUC","B")
    Eval.consensus<-apply(my.evali[,c(2:ncol(my.evali))],2,mean,na.rm =T)
    best.model<-which(Eval.consensus==max(Eval.consensus))
    
    pred.eco.final<-cbind(pred.remain.eco.cv[,1,1],apply(pred.remain.eco.cv[,2:4,],c(1,2),mean,na.rm=T))
    pred.eco.final<-round(cbind(pred.eco.final,apply(pred.eco.final[,2:4],1,weighted.mean,EM.w,na.rm=T)))
    colnames(pred.eco.final)<-c('obs_id','ME','GAM','RF','EM')
    
    ESM.eval[,,z]<-my.evali
    ESM.pred.eval[,,z]<-as.matrix(cbind(pred.eco,round(EM.pred)))
    ESM.pred.final[,,z]<-pred.eco.final
    ESM.w<-rbind(ESM.w,Eval.consensus)
  }
  
  # Selection of the best SM
  th.w<-EM.th
  retained.sm<-which(ESM.w[,4]>=th.w)
  retained.pred<-ESM.list[retained.sm,]
  
  
  # Evaluation of the ESM
  pred.eco<-apply(as.matrix(ESM.pred.eval[,4,retained.sm]),1,weighted.mean,ESM.w[retained.sm,4])
  
  TSS.eco.esm<-TSS.Kappa(pred.eco,pres.abs.eco,TSS=T)
  Se.eco.esm<-2*((round(TSS.eco.esm[4])/100)-0.5)
  AUC.eco.esm<-somers2(pred.eco,pres.abs.eco)[2]
  B.eco.esm<-Boyce(pred.eco,which(pres.abs.eco==1))
  
  
  my.eval.esm<-round(rbind(c(TSS.geo[1], TSS.eco.esm[1]),c(Se.geo,Se.eco.esm),c(AUC.geo,AUC.eco.esm),
                           c(B.geo,B.eco.esm)),2)
  colnames(my.eval.esm)<-c("GAM.geo","ESM")
  row.names(my.eval.esm)<-c("TSS","Se","AUC","B")
  
  pred.eco.final.esm<-apply(as.matrix(ESM.pred.final[,5,retained.sm]),1,weighted.mean,ESM.w[retained.sm,4])
  pred.eco.final.esm<-round(cbind(ESM.pred.final[,1,1],pred.eco.final.esm))
  colnames(pred.eco.final.esm)<-c('obs_id','ESM')
  
  
  ####Thresholds
  
  th.eco<-round(c(quantile(pred.eco[which(pres.abs.eco==1)],c(0,0.01,0.05,0.1,0.2)),TSS.Kappa(pred.eco,pres.abs.eco,TSS=T)[2]))
  
  colnames(results.eco)<-c("id_obs",'pred_ESM',"100%","99%","95%","90%","80%","TSS_CutOff")
  
  
  colnames(ESM.w)<-c('wMAXENT','wGAM','wRF','wESM')
  
  
  my.th<-cbind(th.geo,th.eco)
  colnames(my.th)<-c('GEOth','ECOth')
  row.names(my.th)<-c('0%','1%','5%','10%','20%','TSS_Optimised')
  
  if (!dir.exists(paste0(my.dir,"/results/evaluations"))){
    dir.create(paste0(my.dir,"/results/evaluations"),recursive = TRUE)
  }
  write.csv(my.eval.esm,paste(my.dir,"results/evaluations/",my.sp.id,"_evaluation.csv",sep=""),quote=FALSE,
            row.names=TRUE)
  
  if (!dir.exists(paste0(my.dir,"/results/thresholds"))){
    dir.create(paste0(my.dir,"/results/thresholds"),recursive = TRUE)
  }
  write.csv(my.th,paste(my.dir,"results/thresholds/",my.sp.id,"_thresholds.csv",sep=""),quote=FALSE,
            row.names=TRUE)
  
  if (!dir.exists(paste0(my.dir,"/results/Eweights"))){
    dir.create(paste0(my.dir,"/results/Eweights"),recursive = TRUE)
  }
  write.csv(ESM.w,paste(my.dir,"results/Eweights/",my.sp.id,"_Eweights.csv",sep=""),quote=FALSE,
            row.names=TRUE)
  
  if (!dir.exists(paste0(my.dir,"/results/ESM"))){
    dir.create(paste0(my.dir,"/results/ESM"),recursive = TRUE)
  }
  save(ESM.list,retained.sm,ESM.w,file = paste(my.dir,"results/ESM/",my.sp.id,"_ESM",sep=""), compress = 'xz')
  
  
  
  #unlink(getwd(),rec=TRUE)
}



