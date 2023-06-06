########## SP Richness ##########


# observed richness -------------------------------------------------------


library(preprocessCore)
library(sf)
library(terra)

kt<-st_read('F:/InfoFlora/CWR/sig_data/kt_ch.shp')
kt<-kt[-which(kt$kuerzel==0),]

my.sp.list<-read.table("F:/InfoFlora/CWR/Analyses/Liste_finale_ID2_checklist.txt", h=T, sep="\t")
sp.ids<-my.sp.list$taxon_id_Checklist

exclude<-500 #uncertainty threshold
xcol<-11 # Column with the longitude
ycol<-12 # Column with the latitude
rcol<-26 # Column with the radius
vcol<-19 # Column with validation status of the observation
shpcol<-27 # Column with the geometry (without buffer)
bufcol<-28 # Column with the geometry (with buffer)
date.col<-6 # # Column with the date
date.lim<-'2001-12-31'
dist.pop<-25 # minimal distance between two populations

grid.extent<- st_sf(a = 1:2, geom = st_sfc(st_point(x = c(5.959163471,45.816873303)), # these coordinates are provided by infoflora API
                                           st_point(x = c(10.576523538944393,47.80822259322885))),
                    crs = 4326) %>% 
  st_transform(21781)
r_1km<-st_coordinates(grid.extent) %>% 
  vect %>% 
  rast(res = 1000,crs = "epsg:21781") # create a 100 m raster

Nsp<-setValues(r_1km,0)# empty raster
for (i in 1:nrow(data)){ #load the observation for each species and rasterization
  print(i)
  load(paste0('F:/infoflora/CWR/Analyses/data_sp/',my.sp.list[i,13])) # load observations for the species
  if (NCOL(my.sp)>1){
    my.sp<-my.sp[my.sp[,rcol]<=exclude,] #remove unprecise observations
    if(nrow(my.sp)>0){
      my.sp<-my.sp[my.sp[,date.col]>date.lim,] #remove unprecise observations
      if(nrow(my.sp)>0){
        my.sp<-SpatialPointsDataFrame(coords=my.sp[,11:12], data=my.sp) # transforms the dataframe into spatial points dataframe
        
        crs(my.sp)<-crs(sq)
        my.sp<-remove.duplicates(my.sp, zero = dist.pop) # remove aggregated data
        
        my.sp<-st_as_sf(my.sp)
        my.sp<-my.sp[which(st_intersects(my.sp,st_union(kt),sparse=F)),] #keep swiss observations only
        
        if(nrow(my.sp)>0){
          my.spi<-rep(1,nrow(my.sp)) %>% as.data.frame
          st_geometry(my.spi) <-st_geometry(my.sp)
          r.spi<-vect(my.spi) %>% 
            rasterize(r_1km,background = 0)
          
          Nsp<-Nsp + r.spi
          
        }
      }
    }
  }
}

writeRaster(Nsp, 'D:/infoflora/CWR/Analyses/sp_richness2/nsp_obs.tif',
            datatype = 'INT2U',overwrite=TRUE,
            gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"))


# Distribution for each priority species ----------------------------------

library(terra)
library(RColorBrewer)
#library(rasterVis)
library(sf)
#library(viridis)
#library(PresenceAbsence)
#library(preprocessCore)

exclude<-500 #uncertainty threshold
geo.exclude=5000
xcol<-11 # Column with the longitude
ycol<-12 # Column with the latitude
rcol<-26 # Column with the radius
vcol<-19 # Column with validation status of the observation
shpcol<-27 # Column with the geometry (without buffer)
bufcol<-28 # Column with the geometry (with buffer)
date.col<-6 # # Column with the date
date.lim<-'2001-12-31' #older date
dist.pop<-25 # minimal distance between two populations


ch<-st_read('F:/InfoFlora/CWR/sig_data/ch_aussengrenze.shp',crs = 21781)
lake<-st_read('F:/InfoFlora/CWR/sig_data/see_ch_n.shp',crs = 21781)
kt<-st_read('F:/InfoFlora/CWR/sig_data/kt_ch.shp',crs = 21781)
kt<-kt[-which(kt$kuerzel==0),]
st_crs(kt)<-21781
# kt.sp<-as(kt,'Spatial')
# lake.sp<-as(lake,'Spatial')

data<-read.delim('F:/InfoFlora/CWR/Analyses/Liste_finale_ID2_checklist.txt', h=T, sep='\t', stringsAsFactors = F) #list of the priority CWR
models<-read.delim('F:/InfoFlora/InfoFlora_EvaluationTool/BP/V2/data/results/results_synth.txt', h=T, sep='\t', stringsAsFactors = F) # results of the models
mod.ids<-unlist(lapply(sapply(models[,1],strsplit,"_"), `[[`, 1))
my.models<-models[match(data[,13],mod.ids),]

evaluation.table<-data.frame(sp_id = data$taxon_id_Checklist, 
                             sp_name = data$nom.sans.auteur,
                             Nobs = my.models$Nobs,
                             modeling_type = my.models$modeling,
                             predictors = rep(NA,nrow(my.models)),
                             AUC = my.models$AUC_std_eco/2 + 0.5,
                             TSS = my.models$TSS_eco, 
                             B = my.models$B_eco,
                             SE = my.models$Se_std_eco/2 + 0.5) # models evaluation
my.preds<-c()
for(i in 1:nrow(my.models)){ # create the list of predictors used for SDMs
  my.preds<-c(my.preds,colnames(my.models[,10:42][which(my.models[i,10:42]!= "FALSE")]))
  evaluation.table$predictors[i] = paste(colnames(my.models[,10:42][which(my.models[i,10:42]!= "FALSE")]),collapse = ", ")
}

my.preds<-unique(my.preds)
write.table(evaluation.table,"F:/InfoFlora/CWR/Analyses/Resultats_200428/model_evaluation.txt", sep ="\t",quote =F, row.names =F)
write.table(my.preds,"F:/InfoFlora/CWR/Analyses/Resultats_200428/pred_list.txt", sep ="\t", quote =F, row.names = F)

# load the mask of the species richness map
r.nsp<-rast('F:/infoflora/CWR/Analyses/sp_richness2/nsp_obs.tif')

# specify the directory of the SDMs data
source.dir<-"F:/InfoFlora/InfoFlora_EvaluationTool/BP/V2" #define the directory of the generic structure
load(paste(source.dir,'/data/sp/abs.Rdata',sep="")) # background data used in the SDM

mask100<-rast(paste0(source.dir,'/data/results/rawPred_eco/1000340_Achillea_millefolium_eco.tif')) #100 m projection mask
mask100[!is.na(mask100)]<-1
mask1000<-aggregate(mask100,10,max) 

k<-0
for (i in 1:nrow(my.models)){
  if (!is.na(my.models[i,2]))   { # check if there are species distribution predictions
    x<-as.character(my.models[i,1])
    load(paste0(source.dir,'/data/sp/sp_list/',x))
    
    if (NCOL(my.sp)>1){
      my.sp<-my.sp[my.sp[,rcol]<=exclude,] #remove unprecise observations
      if(nrow(my.sp)>0){
        my.sp<-my.sp[my.sp[,date.col]>date.lim,] #remove unprecise observations
        if(nrow(my.sp)>0){
          my.sp<-SpatialPointsDataFrame(coords=my.sp[,11:12], data=my.sp) # transforms the dataframe into spatial points dataframe
          
          my.sp<-remove.duplicates(my.sp, zero = dist.pop)
          my.sp<-st_as_sf(my.sp) %>% 
            st_set_crs (21781)
          my.sp<-my.sp[which(st_intersects(my.sp,st_union(kt),sparse=F)),] #keep swiss observations only
          
        }
      }
    }
    
    my.var<-my.models[i,10:42]
    my.var<-my.var[which(my.var!=F)]
    names(my.var)<-gsub("gdd3Y_8110_ngb5_mwconic_","gdd3y",names(my.var))# shorten vrariable name
    png(paste0('F:/InfoFlora/CWR/Analyses/distributions/',x,'.png'),
        width = 3940, height = 1920,pointsize =30)
    par(mfrow=c(2,4))
    plot(st_geometry(my.sp),cex=0.8,col='red',main=paste(my.models[i,1],'distribution'),pch=19)
    plot(st_geometry(kt),border='grey45',add=T)
    plot(st_geometry(lake),border=NA,col='lightblue',add=T)
    
    if(!is.na(my.models[i,2])){
      if(my.models[i,2]=='done'){
        k<-k+1
        proj.eco<-rast(paste0(source.dir,'/data/results/rawPred_eco/',x,'_eco.tif'))*(mask100>=0)
        
        multiplier<-c(round((1000-minmax(proj.eco)[2])/1000))*10
        multiplier[multiplier[]==0]<-1
        if (minmax(proj.eco)[2]<=10){
          multiplier = 100
        }
        ex.values<-tapply(terra::extract(proj.eco,my.sp.dist[,c(1:2)])[,1],INDEX = my.sp.dist[,4],FUN = max)
        ex.values<-ex.values * multiplier
        
        abs.values<-terra::extract(proj.eco,abs.rand[,c(2:3)])[,1]
        abs.values<-abs.values *multiplier
        
        DATA<-as.data.frame(na.exclude(cbind(1:(length(abs.values)+length(ex.values)),c(rep(1,length(ex.values)),rep(0,length(abs.values))),
                                             c(ex.values,abs.values))))
        DATA[,3]<-DATA[,3]/1000
        colnames(DATA)<-c('ID','OBS','ECO')
        
        opt.th<-1000*optimal.thresholds(DATA,opt.methods=c(4,9))[,2]
        names(opt.th)<-c('MaxKappa','MinRocDist')
        
        my.func<-function(x){quantile(x,0.95,na.rm =TRUE)}
        proj.eco<-aggregate(proj.eco,10,fun="my.func")*multiplier
        
        YlGn<-colorRampPalette(brewer.pal(9, "YlGn"))
        plot(proj.eco,col=YlGn(1001),axes=F,main=paste(my.models[i,1],'suitability'))
        plot(st_geometry(kt),border='grey45',add=T)
        plot(st_geometry(lake),border=NA,col='lightblue',add=T)
        
        
        plot(proj.eco>=opt.th[1],col=YlGn(2),axes=F,main=paste(my.models[i,1],'MaxKappa'))
        plot(st_geometry(kt),border='grey45',add=T)
        plot(st_geometry(lake),border=NA,col='lightblue',add=T)
        
        
        plot(proj.eco>=opt.th[2],col=YlGn(2),axes=F,main=paste(my.models[i,1],'RocOpt'))
        plot(st_geometry(kt),border='grey45',add=T)
        plot(st_geometry(lake),border=NA,col='lightblue',add=T)
        
        par (mar=c(7,4.1,4.1,2.1))
        barplot(as.numeric(unlist(my.var)),names.arg=names(my.var),las=2,main=paste0('Model evaluation = ',my.models[i,8]))
        
        #plot(1, type="n", axes=F, xlab="", ylab="") #blank plot
        par (mar=c(5.1,4.1,4.1,2.1))
        
        plot(proj.eco>=(quantile(ex.values,0.05,na.rm=T)),col=YlGn(2),axes=F,main=paste(my.models[i,1],'q95'))
        plot(st_geometry(kt),border='grey45',add=T)
        plot(st_geometry(lake),border=NA,col='lightblue',add=T)
        
        plot(proj.eco>=(quantile(ex.values,0.1,na.rm=T)),col=YlGn(2),axes=F,main=paste(my.models[i,1],'q90'))
        plot(st_geometry(kt),border='grey45',add=T)
        plot(st_geometry(lake),border=NA,col='lightblue',add=T)
        
        plot(proj.eco>=(quantile(ex.values,0.2,na.rm=T)),col=YlGn(2),axes=F,main=paste(my.models[i,1],'q85'), )
        plot(st_geometry(kt),border='grey45',add=T)
        plot(st_geometry(lake),border=NA,col='lightblue',add=T)
        
        
        
        if(k==1){
          cum.suit.kappa<-proj.eco>=opt.th[1]
          cum.suit.roc<-proj.eco>=opt.th[2]
          cum.suit.95<-proj.eco>=(quantile(ex.values,0.05,na.rm=T))
          cum.suit.90<-proj.eco>=(quantile(ex.values,0.1,na.rm=T))
          cum.suit.85<-proj.eco>=(quantile(ex.values,0.15,na.rm=T))
          
        }else{
          cum.suit.kappa<-cum.suit.kappa+(proj.eco>=opt.th[1])
          cum.suit.roc<-cum.suit.roc+(proj.eco>=opt.th[2])
          cum.suit.95<-cum.suit.95+(proj.eco>=(quantile(ex.values,0.05,na.rm=T)))
          cum.suit.90<-cum.suit.90+(proj.eco>=(quantile(ex.values,0.1,na.rm=T)))
          cum.suit.85<-cum.suit.85+(proj.eco>=(quantile(ex.values,0.15,na.rm=T)))
          
        }
      }
    } 
    dev.off()
    
  }
}
dir.rich<-'F:/InfoFlora/CWR/Analyses/sp_richness2/'
writeRaster(cum.suit.85,filename = paste0(dir.rich,'cum_suit_85.tif'), overwrite=TRUE)
writeRaster(cum.suit.90,filename = paste0(dir.rich,'cum_suit_90.tif'), overwrite=TRUE)
writeRaster(cum.suit.95,filename = paste0(dir.rich,'cum_suit_95.tif'), overwrite=TRUE)
writeRaster(cum.suit.kappa,filename = paste0(dir.rich,'cum_suit_kappa.tif'), overwrite=TRUE)
writeRaster(cum.suit.roc,filename = paste0(dir.rich,'cum_suit_roc.tif'), overwrite=TRUE)

#generate cum.suit only with good models
k<-0
for (i in which(my.models$eval_eco_mean>=0.6)){
  if (!is.na(my.models[i,2]))   {
    x<-as.character(my.models[i,1])
    load(paste0(source.dir,'/data/sp/sp_list/',x))
    
    if (NCOL(my.sp)>1){
      my.sp<-my.sp[my.sp[,rcol]<=exclude,] #remove unprecise observations
      if(nrow(my.sp)>0){
        my.sp<-my.sp[my.sp[,date.col]>date.lim,] #remove unprecise observations
        if(nrow(my.sp)>0){
          my.sp<-SpatialPointsDataFrame(coords=my.sp[,11:12], data=my.sp) # transforms the dataframe into spatial points dataframe
          
          my.sp<-remove.duplicates(my.sp, zero = dist.pop)
          my.sp<-st_as_sf(my.sp) %>% 
            st_set_crs (21781)
          my.sp<-my.sp[which(st_intersects(my.sp,st_union(kt),sparse=F)),] #keep swiss observations only
          
        }
      }
    }
    
    my.var<-my.models[i,10:42]
    my.var<-my.var[which(my.var!=F)]
    
    if(!is.na(my.models[i,2])){
      if(my.models[i,2]=='done'){
        k<-k+1
        proj.eco<-rast(paste0(source.dir,'/data/results/rawPred_eco/',x,'_eco.tif'))*(mask100>=0)
        
        multiplier<-c(round((1000-minmax(proj.eco)[2])/1000))*10
        multiplier[multiplier[]==0]<-1
        if (minmax(proj.eco)[2]<=10){
          multiplier = 100
        }
        ex.values<-tapply(terra::extract(proj.eco,my.sp.dist[,c(1:2)])[,1],INDEX = my.sp.dist[,4],FUN = max)
        ex.values<-ex.values * multiplier
        
        abs.values<-terra::extract(proj.eco,abs.rand[,c(2:3)])[,1]
        abs.values<-abs.values *multiplier
        
        DATA<-as.data.frame(na.exclude(cbind(1:(length(abs.values)+length(ex.values)),c(rep(1,length(ex.values)),rep(0,length(abs.values))),
                                             c(ex.values,abs.values))))
        DATA[,3]<-DATA[,3]/1000
        colnames(DATA)<-c('ID','OBS','ECO')
        
        opt.th<-1000*optimal.thresholds(DATA,opt.methods=c(4,9))[,2]
        names(opt.th)<-c('MaxKappa','MinRocDist')
        
        my.func<-function(x){quantile(x,0.95,na.rm =TRUE)}
        proj.eco<-aggregate(proj.eco,10,fun="my.func")*multiplier
        
        if(k==1){
          cum.suit.kappa<-proj.eco>=opt.th[1]
          cum.suit.roc<-proj.eco>=opt.th[2]
          cum.suit.95<-proj.eco>=(quantile(ex.values,0.05,na.rm=T))
          cum.suit.90<-proj.eco>=(quantile(ex.values,0.1,na.rm=T))
          cum.suit.85<-proj.eco>=(quantile(ex.values,0.15,na.rm=T))
          
        }else{
          cum.suit.kappa<-cum.suit.kappa+(proj.eco>=opt.th[1])
          cum.suit.roc<-cum.suit.roc+(proj.eco>=opt.th[2])
          cum.suit.95<-cum.suit.95+(proj.eco>=(quantile(ex.values,0.05,na.rm=T)))
          cum.suit.90<-cum.suit.90+(proj.eco>=(quantile(ex.values,0.1,na.rm=T)))
          cum.suit.85<-cum.suit.85+(proj.eco>=(quantile(ex.values,0.15,na.rm=T)))
          
        }
      }
    } 
    
  }
}
dir.rich<-'F:/InfoFlora/CWR/Analyses/sp_richness2/'
writeRaster(cum.suit.85,filename = paste0(dir.rich,'cum_suit_85_goodModels.tif'), overwrite=TRUE)
writeRaster(cum.suit.90,filename = paste0(dir.rich,'cum_suit_90_goodModels.tif'), overwrite=TRUE)
writeRaster(cum.suit.95,filename = paste0(dir.rich,'cum_suit_95_goodModels.tif'), overwrite=TRUE)
writeRaster(cum.suit.kappa,filename = paste0(dir.rich,'cum_suit_kappa_goodModels.tif'), overwrite=TRUE)
writeRaster(cum.suit.roc,filename = paste0(dir.rich,'cum_suit_roc_goodModels.tif'), overwrite=TRUE)

cum.suit.85<-rast(paste0(dir.rich,'cum_suit_85_goodModels.tif'))
cum.suit.90<-rast(paste0(dir.rich,'cum_suit_90_goodModels.tif'))
cum.suit.95<-rast(paste0(dir.rich,'cum_suit_95_goodModels.tif'))
cum.suit.kappa<-rast(paste0(dir.rich,'cum_suit_kappa_goodModels.tif'))
cum.suit.roc<-rast(paste0(dir.rich,'cum_suit_roc_goodModels.tif'))
r.nsp<-rast('F:/infoflora/CWR/Analyses/sp_richness2/sp_richness3_withGoodModels.tif')


# quantile normalization --------------------------------------------------

# # alternate function if you cannot install preprocess core
# quantile_normalisation <- function(df){
#   df_rank <- apply(df,2,rank,ties.method="min")
#   df_sorted <- data.frame(apply(df, 2, sort))
#   df_mean <- apply(df_sorted, 1, mean)
#   
#   index_to_mean <- function(my_index, my_mean){
#     return(my_mean[my_index])
#   }
#   
#   df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
#   rownames(df_final) <- rownames(df)
#   return(df_final)
# }

# extract values from species richness
r.nsp<-round(resample(r.nsp,cum.suit.85)) #align rasters

v.nsp<-subst(r.nsp,NA,0) %>% values
v.th.95<-subst(cum.suit.95,NA,0) %>% values
v.th.90<-subst(cum.suit.90,NA,0) %>% values
v.th.85<-subst(cum.suit.85,NA,0) %>% values
v.th.kappa<-subst(cum.suit.kappa,NA,0) %>% values
v.th.roc<-subst(cum.suit.roc,NA,0) %>% values


mat<-cbind(v.nsp,v.th.95,v.th.90,v.th.85,v.th.kappa, v.th.roc)
#mat<-cbind(v.nsp,v.th.95)

# normalization
#std.val<-quantile_normalisation(mat)
std.val<-preprocessCore::normalize.quantiles(mat)
std.nsp<-setValues(r.nsp,std.val[,1])
std.th.95<-setValues(r.nsp,std.val[,2])
std.th.90<-setValues(r.nsp,std.val[,3])
std.th.85<-setValues(r.nsp,std.val[,4])
std.th.kappa<-setValues(r.nsp,std.val[,5])
std.th.roc<-setValues(r.nsp,std.val[,6])

# save normalized values
writeRaster(std.nsp,filename = paste0(dir.rich,'sp_richness3_goodModels_std.tif'), overwrite=TRUE)
writeRaster(std.th.95,filename = paste0(dir.rich,'cum_suit_95_goodModels_std.tif'), overwrite=TRUE)
writeRaster(std.th.90,filename = paste0(dir.rich,'cum_suit_90_goodModels_std.tif'), overwrite=TRUE)
writeRaster(std.th.85,filename = paste0(dir.rich,'cum_suit_85_goodModels_std.tif'), overwrite=TRUE)
writeRaster(std.th.kappa,filename = paste0(dir.rich,'cum_suit_kappa_goodModels_std.tif'), overwrite=TRUE)
writeRaster(std.th.roc,filename = paste0(dir.rich,'cum_suit_roc_goodModels_std.tif'), overwrite=TRUE)

# difference between observed and modelled richness
std.nsp = rast(paste0(dir.rich,'sp_richness3_goodModels_std.tif'))
std.th.95 = rast( paste0(dir.rich,'cum_suit_95_goodModels_std.tif'))
std.th.90 = rast( paste0(dir.rich,'cum_suit_90_goodModels_std.tif'))
std.th.85 = rast( paste0(dir.rich,'cum_suit_85_goodModels_std.tif'))
std.th.kappa = rast(paste0(dir.rich,'cum_suit_kappa_goodModels_std.tif'))
std.th.roc = rast( paste0(dir.rich,'cum_suit_roc_goodModels_std.tif'))

def.95.std<-std.nsp-std.th.95
def.90.std<-std.nsp-std.th.90
def.85.std<-std.nsp-std.th.85
def.kappa.std<-std.nsp-std.th.kappa
def.roc.std<-std.nsp-std.th.roc

r.nsp[is.na(r.nsp)]<-0

def.95<-r.nsp-cum.suit.95
def.90<-r.nsp-cum.suit.90
def.85<-r.nsp-cum.suit.85
def.kappa<-r.nsp-cum.suit.kappa
def.roc<-r.nsp-cum.suit.roc


writeRaster(def.95,filename = paste0(dir.rich,'def_95.tif'), overwrite=TRUE)
writeRaster(def.90,filename = paste0(dir.rich,'def_90.tif'), overwrite=TRUE)
writeRaster(def.85,filename = paste0(dir.rich,'def_85.tif'), overwrite=TRUE)
writeRaster(def.kappa,filename = paste0(dir.rich,'def_kappa.tif'), overwrite=TRUE)
writeRaster(def.roc,filename = paste0(dir.rich,'def_roc.tif'), overwrite=TRUE)

writeRaster(def.95.std,filename = paste0(dir.rich,'def_95_std.tif'), overwrite=TRUE)
writeRaster(def.90.std,filename = paste0(dir.rich,'def_90_std.tif'), overwrite=TRUE)
writeRaster(def.85.std,filename = paste0(dir.rich,'def_85_std.tif'), overwrite=TRUE)
writeRaster(def.kappa.std,filename = paste0(dir.rich,'def_kappa_std.tif'), overwrite=TRUE)
writeRaster(def.roc.std,filename = paste0(dir.rich,'def_roc_std.tif'), overwrite=TRUE)


# correlation between deficit and sampling effort -------------------------


def.85<-rast(paste0(dir.rich,'def_85.tif'))
def.90<-rast(paste0(dir.rich,'def_90.tif'))
def.95<-rast(paste0(dir.rich,'def_95.tif'))
def.kappa<-rast(paste0(dir.rich,'def_kappa.tif'))
def.kappa<-rast(paste0(dir.rich,'def_kappa.tif'))
r.nsp<-rast('F:/infoflora/CWR/Analyses/sp_richness2/sp_richness3_withGoodModels.tif')
r.nsp.sample<-round(resample(r.nsp,cum.suit.85))
# load("F:/InfoFlora/Analyses/extract_observations/data/DB_extract_210714.Rdata")#Load the full set of observations to build bias map
# DB$date<-as.Date(fasttime::fastPOSIXct(DB$date))
# DB<-DB[as.Date(date) >= "2002-01-01" & as.Date(date) <= "2019-12-31", ]
# DB<-DB[v_xy_radius <= 500,]
# DB$x<-DB$x-2000000
# DB$y<-DB$y-1000000
# sample_effort<-terra::rasterize(as.matrix(DB[,c("x","y")]),y = r.nsp.sample,fun = length)
# writeRaster(sample_effort,filename = paste0(dir.rich,'sample_effort.tif'), overwrite=TRUE)
sample_effort<-rast(paste0(dir.rich,'sample_effort.tif'))

my.files<-c('def_95','def_90','def_85','def_kappa','def_roc')
my.source.dir<-'F:/InfoFlora/CWR/Analyses/sp_richness2/'

ch<-st_read('F:/InfoFlora/CWR/sig_data/ch_aussengrenze.shp')
lake<-st_read('F:/InfoFlora/CWR/sig_data/see_ch_n.shp')
kt<-st_read('F:/InfoFlora/CWR/sig_data/kt_ch.shp')
kt<-kt[-which(kt$kuerzel==0),]

# correlation with sampling effort
a<-r.nsp.sample
a[is.na(a)]<-0
a<-values(a*mask1000)
b<-sample_effort
b[is.na(b)]<-0
b<-values(b*mask1000)
cor.test(a,b, method = "spearman")

# correlation with modelled richness
r.nsp.sample<-round(resample(r.nsp,cum.suit.85))
a<-r.nsp.sample
a[is.na(a)]<-0
a<-a*mask1000
a<-values(a)
my.files<-c("95","90", "85", "kappa", "roc")
cor.mod.obs<-c()
cor.mod.eff<-c()
for (i in 1:5){
  b<-values(get(paste0("cum.suit.",my.files[i])))
  cor.mod.obs<-rbind(cor.mod.obs,
                     cor.test(a,b, method = "spearman")[3:4])
  cor.mod.eff<-rbind(cor.mod.eff,
                     cor.test(values(sample_effort),b,method = "spearman")[3:4])
  
}

apply(cor.mod.obs,2,function(x) mean(unlist(x)))
apply(cor.mod.obs,2,function(x) sd(unlist(x)))

apply(cor.mod.eff,2,function(x) mean(unlist(x)))
apply(cor.mod.eff,2,function(x) sd(unlist(x)))


# figures -----------------------------------------------------------------
png('F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/sampling_effort.png',width = 2048, height = 3000,pointsize =50)
par(mfrow = c(2,1))
sample_effort<-rast(paste0(dir.rich,'sample_effort.tif'))
sample_effort[is.na(sample_effort)]<-0
sample_effort_log<-log(sample_effort,10)
plot(sample_effort_log*mask1000,axes =F, main = "a) Sampling effort (log scale)",mar = c(1, 1, 2, 4))
plot(st_geometry(kt),col = NA,add =T,border = "gray50",lwd = 0.25)
plot(st_geometry(ch),add = T, col = NA)
plot(st_geometry(lake), add = T, col = "lightblue",border = NA)


sample_effort.cat<-sample_effort>=500
plot(sample_effort.cat*mask1000, col = c("lightyellow","pink"),axes =F, mar = c(1, 1, 2, 4),
     main = "b) Sample Effort categorized", levels = c("Low", "High"),type = "classes")
plot(st_geometry(kt),col = NA,add =T,border = "gray70",lwd = 0.25)
plot(st_geometry(ch),add = T, col = NA)
plot(st_geometry(lake), add = T, col = "lightblue",border = NA)
dev.off()

png('F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/fig_sampling_effort.png',width = 2048, height = 1500,pointsize =50)
sample_effort<-rast(paste0(dir.rich,'sample_effort.tif'))
sample_effort[is.na(sample_effort)]<-0
plot(log(sample_effort*mask1000,10),axes =F, main = "Sampling effort (log scale)",mar = c(1, 1, 2, 4),cex.main=1.5)
plot(st_geometry(kt),col = NA,add =T,border = "gray50",lwd = 0.25)
plot(st_geometry(ch),add = T, col = NA)
plot(st_geometry(lake), add = T, col = "lightblue",border = NA)
dev.off()

sample.big<-sample_effort>=500
sample.big[sample.big != 1]<- NA
sample.low<-sample_effort<500
sample.low[sample.low!=1]<-NA


YlOrRd<-colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd")[9:3],bias = 0.5)
Blues<-colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[9:1])

png('F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/deficit_all.png',width = 1024, height = 2048,pointsize =30)
par(mfcol = c(5,2))
my.files<-c('95','90','85','kappa','roc')
my.title<-c('OR5','OR10','OR15','kappa','roc')

for (i in 1:length(my.files)){
  my.ras<-rast(paste0(dir.rich,'def_',my.files[i],'_std.tif'))
  my.ras1<-my.ras*sample.big*(my.ras<0)
  plot(my.ras1,col = YlOrRd(100), axes =F, range=c(minmax(my.ras)[1],0),mar = c(1,1, 1.1, 4), main = paste0("Deficit ",my.title[i]," standardized"))
  my.ras2<-my.ras*sample.low*(my.ras<0)
  plot(my.ras2,col = Blues(100),add =T, plg = list(loc = "bottom"),axes =F, range=c(minmax(my.ras)[1],0),mar = c(1, 1, 1.1, 4))
  plot(st_geometry(kt),col = NA,add =T,border = "gray70",lwd = 0.25)
  plot(st_geometry(ch),add = T, col = NA)
  plot(st_geometry(lake), add = T, col = "lightblue",border = NA)
  
}

for (i in 1:length(my.files)){
  my.ras<-rast(paste0(dir.rich,'def_',my.files[i],'.tif'))
  my.ras1<-my.ras*sample.big*(my.ras<0)
  plot(my.ras1,col = YlOrRd(100), axes =F, range=c(minmax(my.ras)[1],0),mar = c(1,1, 1.1, 4), main = paste0("Deficit ",my.title[i]))
  my.ras2<-my.ras*sample.low*(my.ras<0)
  plot(my.ras2,col = Blues(100),add =T, plg = list(loc = "bottom"),axes =F, range=c(minmax(my.ras)[1],0),mar = c(1, 1, 1.1, 4))
  plot(st_geometry(kt),col = NA,add =T,border = "gray70",lwd = 0.25)
  plot(st_geometry(ch),add = T, col = NA)
  plot(st_geometry(lake), add = T, col = "lightblue",border = NA)
  
}

dev.off()

png('F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/deficit_OR10.png',width = 2048, height = 1500,pointsize =50)
my.ras<-rast(paste0(dir.rich,'def_90_std.tif'))*mask1000
my.ras1<-my.ras*sample.big*(my.ras<0)
plot(abs(my.ras1),
     col = YlOrRd(100)[100:1], 
     axes =F, 
     range=c(0,minmax(abs(my.ras))[2]),
     mar = c(1,1, 2, 4.1), 
     main = "Deficit OR10 standardized",
     cex.main=1.5,
     plg=list(cex=1.25))
my.ras2<-my.ras*sample.low*(my.ras<0)
plot(my.ras2,
     col = Blues(100),
     add =T, 
     plg = list(loc = "bottom", cex = 1.25),
     axes =F, 
     range=c(minmax(my.ras)[1],0),
     mar = c(1, 1, 1.1, 4))
plot(st_geometry(kt),col = NA,add =T,border = "gray40")
plot(st_geometry(ch),add = T, col = NA)
plot(st_geometry(lake), add = T, col = "lightblue",border = NA)
dev.off()

YlGn<-colorRampPalette(RColorBrewer::brewer.pal(9, "YlGn")[-1])
png('F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/Nsp_all.png',width = 2048, height = 2048,pointsize =30)
r.nsp<-rast('F:/infoflora/CWR/Analyses/sp_richness2/sp_richness3_withGoodModels.tif')
mat.layout<-matrix(c(1,2,3,7,8,9,4,5,6,10,11,12),nrow=4,byrow =T)
#par(mfcol = c (6,2))
layout(mat.layout)
plot(r.nsp,col=YlGn(100),axes =FALSE, main = "Nsp_obs",mar = c(1, 1, 1.5, 4))
plot(st_geometry(kt),col = NA,add =T,border = "gray70",lwd = 0.25)
plot(st_geometry(ch),add = T, col = NA)
plot(st_geometry(lake), add = T, col = "lightblue",border = NA)

my.files<-c('cum_suit_95_goodModels','cum_suit_90_goodModels',
            'cum_suit_85_goodModels','cum_suit_kappa_goodModels',
            'cum_suit_roc_goodModels')
my.title<-c('Nsp_OR5','Nsp_OR10','Nsp_OR15','Nsp_kappa','Nspt_roc')

for(i in 1:length(my.files)){
  
  my.def<-rast(paste0(my.source.dir,my.files[i],'.tif'))*mask1000
  plot(my.def,col=YlGn(100),axes =FALSE,main = my.title[i],mar = c(1, 1, 1.5, 4))
  plot(st_geometry(kt),col = NA,add =T,border = "gray70",lwd = 0.25)
  plot(st_geometry(ch),add = T, col = NA)
  plot(st_geometry(lake), add = T, col = "lightblue",border = NA)
  
}

plot(std.nsp*mask1000,col=YlGn(100),axes =FALSE, main = "Nsp_obs_std",mar = c(1, 1, 1.5, 4))
plot(st_geometry(kt),col = NA,add =T,border = "gray70",lwd = 0.25)
plot(st_geometry(ch),add = T, col = NA)
plot(st_geometry(lake), add = T, col = "lightblue",border = NA)

my.files<-c('cum_suit_95_goodModels_std','cum_suit_90_goodModels_std',
            'cum_suit_85_goodModels_std','cum_suit_kappa_goodModels_std',
            'cum_suit_roc_goodModels_std')
my.title<-c('Nsp_OR5_std','Nsp_OR10_std','Nsp_OR15_std','Nsp_kappa_std','Nspt_roc_std')

for(i in 1:length(my.files)){
  
  my.def<-rast(paste0(my.source.dir,my.files[i],'.tif'))*mask1000
  plot(my.def,col=YlGn(100),axes =FALSE,main = my.title[i],mar = c(1, 1, 1.5, 4))
  plot(st_geometry(kt),col = NA,add =T,border = "gray70",lwd = 0.25)
  plot(st_geometry(ch),add = T, col = NA)
  plot(st_geometry(lake), add = T, col = "lightblue",border = NA)
  
}
dev.off()

#only with models
to.keep<-my.models[which(!is.na(my.models$proji)),1] %>%
  strsplit("_") %>% 
  sapply(function(x)x[[1]][1]) %>% 
  as.numeric()

sp_sq<-sp.dat[v_accepted_taxon_id %in% to.keep, .N, by=.(Intial_TaxonID, sq1)]
sq_summary<-sp_sq[, .N, by = .(sq1)]
head(sq_summary)
sp.richness<-grid.1km
sp.richness$Nsp<-0
sp.richness$Nsp[sq_summary$sq1]<-sq_summary$N
sp.richness<-filter(sp.richness,Nsp >0)
sp.rich.class<-mutate(sp.richness,my.class = cut (sp.richness$Nsp,c(1,20,40,60,80,106),include.lowest = T, right = TRUE))

r_1km<-st_coordinates(grid.extent) %>% 
  vect %>% 
  rast(res = 1000,crs = "epsg:21781")


sp.richness.r<-rasterize(vect(sp.rich.class),r_1km,field = "Nsp")
writeRaster(sp.richness.r,"F:/InfoFlora/CWR/Analyses/sp_richness2/sp_richness3_withModels.tif",
            datatype = "INT1U",
            overwrite=TRUE,
            gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"))

#only with GOOD models
which(my.models$eval_eco_mean>=0.6)
to.keep<-my.models[which(my.models$eval_eco_mean>=0.6),1] %>%
  strsplit("_") %>% 
  sapply(function(x)x[[1]][1]) %>% 
  as.numeric()

sp_sq<-sp.dat[v_accepted_taxon_id %in% to.keep, .N, by=.(Intial_TaxonID, sq1)]
sq_summary<-sp_sq[, .N, by = .(sq1)]
head(sq_summary)
sp.richness<-grid.1km
sp.richness$Nsp<-0
sp.richness$Nsp[sq_summary$sq1]<-sq_summary$N
sp.richness<-filter(sp.richness,Nsp >0)
sp.rich.class<-mutate(sp.richness,my.class = cut (sp.richness$Nsp,c(1,20,40,60,80,106),include.lowest = T, right = TRUE))

r_1km<-st_coordinates(grid.extent) %>% 
  vect %>% 
  rast(res = 1000,crs = "epsg:21781")

sp.richness.r<-rasterize(vect(sp.rich.class),r_1km,field = "Nsp")
writeRaster(sp.richness.r,"F:/InfoFlora/CWR/Analyses/sp_richness2/sp_richness3_withGoodModels.tif",
            datatype = "INT1U",
            overwrite=TRUE,
            gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"))



