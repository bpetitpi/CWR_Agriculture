############################
##### Helper functions #####
############################

# function to weight observations
re.w<-function(x){
  return(c(round(100*(1/length(x)))))
}

# helper to extract elements from lists
list.pckp<-function(x,y){return(x[[y]])}

# helper to apply innerness simultaneously to a matrix
apply.innerness<-function(th,foc.pop,sp.dist,bck.dist){
  sp.disti<-sp.dist[which(sp.dist>=th)]
  if (length(sp.disti)==0){
    sp.dist<-max(sp.dist)
  }else{
    sp.dist<-sp.disti
  }
  
  inner<-innerouterness.nichedens(foc.pop=foc.pop,sp.dist=sp.dist,bck.dist=bck.dist, test=NA, data.test=NA)  
  inandout<-round(100*(inner$innerness+inner$outerness))
  return(inandout)
}

apply.extract<-function(x,y){
  return(extract(y,x,small=T))
}

dsg.accurate<-function(x){
  y<-x[3:4]
  x<-x[1:2]
  return(x[which(y==max(y))][1]) 
}

TSS.Kappa<-function(my.pred,pres.abs,TSS=T){
  if (length(which(is.na(my.pred)))>0){
    if (length(which(is.na(my.pred)))==length(my.pred)){
      return(rep(NA,6))
    }else{
      pres.abs<-pres.abs[-which(is.na(my.pred))]
      my.pred<-my.pred[-which(is.na(my.pred))]
      return(unlist(KappaRepet(pres.abs,my.pred,TSS)))      
    }
  }else{
  return(unlist(KappaRepet(pres.abs,my.pred,TSS)))
  }
}

Boyce<-function(pred.eco,my.pres){
  if (length(which(is.na(pred.eco)))>0){
    if (length(which(is.na(pred.eco)))==length(pred.eco)){
      return(NA)
    }else{
      my.pres<-cbind(my.pres,pred.eco[my.pres])
      my.pres<-my.pres[-which(is.na(my.pres[,2])),1]
      pred.eco<-pred.eco[-which(is.na(pred.eco))]
      return(ecospat.boyce(pred.eco, pred.eco[my.pres], PEplot = F,cor.method = "kendall")$correlation)
      
    }
  }else{
   return(ecospat.boyce(pred.eco, pred.eco[my.pres], PEplot = F,cor.method = "kendall")$correlation)
  }
}

CV.data<-function(data,Nfold=10){
  suppressWarnings(data.split<-split(sample(1:nrow(data),nrow(data)),1:Nfold))
  
}

cv.format<-function(id.pres.geo,id.pres.eco,id.abs.strat ,id.abs.rand){
  return(list(id.pres.geo=id.pres.geo, id.pres.eco=id.pres.eco,id.abs.strat=id.abs.strat, id.abs.rand=id.abs.rand))
}

logit <- function(x){ exp(x)/(exp(x)+1)}

poly2pts<-function(poly,res,ID = NULL){
  if (area(poly)< res){
    poly2ptsi<-as.data.frame(matrix(coordinates(poly),nrow=1))
    poly2ptsi.id<-ID
    poly2pts.df<-cbind(poly2ptsi,poly2ptsi.id)
  }else{
    cellsize= res/10
    poly2ptsi<-spsample(poly,type="regular",cellsize=cellsize)
    poly2ptsi<-remove.duplicates(poly2ptsi,zero = res)
    poly2ptsi.id<-rep(ID,length(poly2ptsi))
    poly2pts.df<-cbind(coordinates(poly2ptsi),poly2ptsi.id)
  }
  colnames(poly2pts.df)<-c("x","y","id.cell")
  return(poly2pts.df)
}

pred.max<-function(val, fac){
  return(tapply(val,fac,max))
}

gam.pred.total<-function(my.sp.id,my.abs,env,ks=3, min.dist=99,my.dir, my.th=0.05,my.mask.pts, exclude = 500, geo.exclude=5000){
  load(paste(my.dir,'sp/Sp_list/',my.sp.id,sep=""))
  to.rm<-which(is.na(match(my.sp.tot[,1],unique(my.sp.dist[,4]))))
  
  if (length(to.rm>0)){
    my.sp.xy<-my.sp.tot[-to.rm,c(11:12,1,26)]
  } else{
    my.sp.xy<-my.sp.tot[,c(11:12,1,26)]
    }
  
  vld.xy<-my.sp.tot[my.sp.tot[,26]>exclude & my.sp.tot[,26]<geo.exclude & my.sp.tot[,19]>=300,c(11,12,1,26)]
  my.sp.xy<-rbind(my.sp.xy,vld.xy)

  if(length(unique(my.sp.dist[,4])) >= 8 | nrow(my.sp.xy)>=8 ){ 

  ##### remove aggregated data
  if (min.dist>0){
    my.sp.xy.dup<-zerodist(SpatialPoints(coord=my.sp.xy[,1:2]),zero=min.dist)
    prec<-cbind(my.sp.xy[my.sp.xy.dup[,1],4],my.sp.xy[my.sp.xy.dup[,2],4])
    to.rm<-unique(apply(cbind(my.sp.xy.dup,prec),1,dsg.accurate))
    
    if(length(to.rm)>0){
      my.sp.xy<-my.sp.xy[-to.rm,]
      my.sp.dist<-my.sp.dist[which(!is.na(match(my.sp.dist[,4],my.sp.xy[,3]))),]
    }
  }
  
  colnames(my.abs)<-colnames(my.sp.xy)[1:2]
  knn.tot<-log(round(knnx.dist(my.sp.xy[,1:2],coordinates(my.mask.pts),k=2)/100)[,1]+1)
  Ndist<-rasterize(my.mask.pts,env,field=knn.tot)
  env<-stack(env,Ndist)
  names(env)<-c("ddeg","prec","bio6","soil","vgh", "Ndist")
  Ndist<-c(log(round(knn.dist(my.sp.xy[,1:2],k=2)/100)[,1]+1),
           log(round(knnx.dist(my.sp.xy[,1:2], my.abs,k=2)/100)[,1]+1))
  sp.geo<-c(rep(1,nrow(my.sp.xy)),rep(0,nrow(my.abs)))
  
  data.env<-data.frame(extract(env,rbind(my.sp.dist[,1:2],setNames(my.abs,colnames(my.sp.dist[,1:2])))))
  names(data.env)<-c("ddeg","prec","bio6","soil","vgh", "Ndist")
  data.env$vgh<-factor(data.env$vgh)
  data.env$soil<-factor(data.env$soil)
  data.env$bio6<-factor(data.env$bio6)
  
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
  m1.eco<-gam(sp~s(ddeg, k = 3) + s(prec, k = 3) + factor(soil) + factor(vgh),data=df.tmp.input,family = binomial(link = "logit"), control = gam.control(maxit = 500), weight=w)
  
  sp.pred.geo<-round(100*predict(env,m1.geo,type='response'))
  sp.pred.eco<-round(100*predict(env,m1.eco,type='response'))
  
  sp.dist.geo<-100*m1.geo$fitted[1:sum(sp.geo)]
  sp.dist.geo<-sp.dist.geo[which(sp.dist.geo>=quantile(sp.dist.geo,my.th))]
  
  sp.dist.eco<-100*m1.eco$fitted[1:nrow(my.sp.dist)]
  sp.dist.eco<-tapply(sp.dist.eco,as.factor(my.sp.dist[,4]),max, na.rm = T)
  sp.dist.eco[which(sp.dist.eco==-Inf)]<-NA
  sp.dist.eco<-sp.dist.eco[which(sp.dist.eco>=quantile(sp.dist.eco,my.th))]
  
  inner.geo<-innerouterness.nichedens(foc.pop=values(sp.pred.geo),sp.dist=sp.dist.geo,bck.dist=values(sp.pred.geo), test=NA, data.test=NA)  
  inner.eco<-innerouterness.nichedens(foc.pop=values(sp.pred.eco),sp.dist=sp.dist.eco,bck.dist=values(sp.pred.eco), test=NA, data.test=NA)  
  
  inner.geo.r<-setValues(sp.pred.geo,round(100*(inner.geo$innerness+inner.geo$outerness)))
  inner.eco.r<-setValues(sp.pred.eco,round(100*(inner.eco$innerness+inner.eco$outerness)))
  
  writeRaster(sp.pred.geo,filename=paste(my.dir,"results/rawPred_geo/",my.sp.id,"_geo",sep=""),format="GTiff",datatype="INT1U",overwrite=TRUE, 
              options=c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  writeRaster(sp.pred.eco,filename=paste(my.dir,"results/rawPred_eco/",my.sp.id,"_eco",sep=""),format="GTiff",datatype="INT1U",overwrite=TRUE, 
              options=c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  writeRaster(inner.geo.r,filename=paste(my.dir,"results/map_innerness_geo/",my.sp.id,"_Innergeo",sep=""),format="GTiff",datatype="INT2S",overwrite=TRUE, 
              options=c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  writeRaster(inner.eco.r,filename=paste(my.dir,"results/map_innerness_eco/",my.sp.id,"_Innereco",sep=""),format="GTiff",datatype="INT2S",overwrite=TRUE, 
              options=c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  }else{
    nb.eco=length(unique(my.sp.dist[,3]))
    nb.geo=nrow(my.sp.xy)
    dir.create(paste(my.dir,"results/Sp_excluded",sep=""),showWarnings = FALSE)
    save(nb.eco,nb.geo,file=paste(my.dir,"results/Sp_excluded/",my.sp.id,sep=""))
  }
}

get.max.pred<-function(x,fac){
  max.pred<-tapply(x,fac,FUN=max)
  return(max.pred)
}

var.import.w<-function(var.import.synth){
  wm<-c()
  for (i in unique(var.import.synth[,1])){
    vari<-var.import.synth[var.import.synth[,1]==i,]
    wmi<-c(i,weighted.mean(vari[,2],vari[,3]))
    wm<-rbind(wm,wmi)
  }
  wm[which(is.nan(wm[,2])),2]<-0
  return(wm)
}


w.mean<-function(r,w){
  check.na<-c()
  for (i in 1:length(w)){
    r[[i]]<-as.numeric(w[i])*r[[i]]
    check.na<-c(check.na,maxValue(r[[i]]))
  }
  if (length(which(is.na(check.na)))>0){
    r<-r[[-which(is.na(check.na))]]
    w<-w[-which(is.na(check.na))]
  }
  return(sum(r,na.rm = T)/sum(w,na.rm =T))
}


