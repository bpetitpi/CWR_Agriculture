source.dir<-"F:/InfoFlora/InfoFlora_EvaluationTool/BP/V2" #define the directory of the generic structure
my.dir<-paste(source.dir,"/data/",sep="")
pred.group<-read.table(paste0(my.dir,'env/predictor_groups.txt'),h=T,sep='\t')
source(paste(my.dir.fct,"functions_extractSp.R",sep=""))


var.names<-c('bio19_pcoldq',
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
                         'ForestQ25')
                         

my.sp.list<-list.files(paste(my.dir,'sp/Sp_list/',sep="")) # list of the species 

sp.var<-c()
for (i in 1:length(my.sp.list)){
  my.sp.name<-my.sp.list[i]
  if (file.exists(paste0(my.dir,'env/varPval/',my.sp.name,'.txt'))){
    my.var<-cbind(pred.group,read.table(paste0(my.dir,'env/varPval/',my.sp.name,'.txt'),h=T,sep="\t"))
    sp.var<-c(sp.var,list(c(my.sp.name,var.select(my.var,core=c('T', 'Season', 'Ar')))))
  }
}

#colonnes : code_nom_espece, pred, expert, nombre d'obs, nombre de prédicteurs, ratio, modeling
inidat0<-as.data.frame(matrix(nrow=length(my.sp.list),ncol=6+length(var.names)))
colnames(inidat0)<-c('code_espece',var.names,'expert_selection','Nobs','Npred','Nobs/Npred','modeling')
inidat0[,1]<-my.sp.list
inidat0[,34]<-rep(F, length(my.sp.list))
for (i in 1:length(sp.var)){
  my.preds<-sp.var[[i]]
  sp.line<-match(my.preds[1],inidat0[,1])
  inidat0[sp.line,2:34]<-rep(F,length(var.names))
  inidat0[sp.line,as.numeric(my.preds[-1])+1]<-T
  load(paste0(my.dir,'Sp/Sp_list/',my.preds[1]))
  inidat0$Nobs[sp.line]<-length(unique(my.sp.dist[,4]))
  inidat0$Npred[sp.line]<-length(my.preds)-1
  inidat0$expert_selection[sp.line]<-FALSE
  inidat0$`Nobs/Npred`[sp.line]<-round(inidat0$Nobs[sp.line]/inidat0$Npred[sp.line],1)
  if (inidat0$`Nobs/Npred`[sp.line]>=10){
    inidat0$modeling[sp.line]<-'EM'
    }else{
    inidat0$modeling[sp.line]<-'ESM'
    }
  }

for (sp.line in which(is.na(inidat0$modeling))){
  inidat0[sp.line,2:34]<-rep(F,length(var.names))
  my.preds<-c(21,14,20)
  inidat0[sp.line,as.numeric(my.preds)+1]<-T
  load(paste0(my.dir,'Sp/Sp_list/',inidat0[sp.line,1]))
  inidat0$Nobs[sp.line]<-length(unique(my.sp.dist[,4]))
  inidat0$Npred[sp.line]<-length(my.preds)
  inidat0$expert_selection[sp.line]<-FALSE
  inidat0$`Nobs/Npred`[sp.line]<-round(inidat0$Nobs[sp.line]/inidat0$Npred[sp.line],1)
  if (inidat0$`Nobs/Npred`[sp.line]>=10){
    inidat0$modeling[sp.line]<-'EM'
  }else{
    inidat0$modeling[sp.line]<-'ESM'
  }
  
}

if (file.exists(paste0(my.dir,'inidat/model_inidat.txt'))){
  inidat<-read.table(paste0(my.dir,'inidat/model_inidat.txt'),h=T,sep='\t')
  to.keep<-inidat[which(inidat$expert_selection==TRUE),]
  inidat0[match(to.keep[,1],inidat0),c(1:35)]<-to.keep[,c(1:35)]
}

inidat<-inidat0
write.table(inidat,file=paste0(my.dir,'inidat/model_inidat.txt'),row.names=FALSE,quote=FALSE,sep='\t')

