# BE Analysis -------------------------------------------------------------


########## Bern Analysis ##########
set.seed(7)
library(raster)
library(sp)
library(rgeos)
library(sf)
library(tidyverse)
library(rmapshaper)
library(geojsonio)
library(doParallel)
library(mapview)

# over.test<-function(pts,bias,n.my.sp,my.area,nrepet=1000){
#   xy.i<-pts[sample(nrow(pts),prob = bias,size=n.my.sp),]
#   st_intersects(xy,my.area)
# }

kt<-st_read('F:/InfoFlora/CWR/sig_data/kt_ch.shp')
# sau_spb<-st_read('F:/InfoFlora/CWR/sig_data/LANDKULT/LV95/data/LANDKULT_NUTZFL.shp')
# spb.code<-read.delim('F:/InfoFlora/CWR/analyses/spb_code.txt',header=T,sep='\t')
# spb.code.i<-apply(spb.code[,3:5],1,sum)
# spb.code<-spb.code[which(spb.code.i>0),]
# estivage.code<-c(621,622,623,930,931,935,936)#categories corresponding to estivage
# estivage<-sau_spb[sau_spb$LNF_CODE%in%estivage.code,] 
# estivage<-ms_dissolve(estivage)
# st_write(estivage,dsn='F:/InfoFlora/CWR/sig_data/spb_estivage.shp',append=F)
# estivage<-st_read('F:/InfoFlora/CWR/sig_data/spb_estivage.shp')
# 
# #SAU_large<-st_read('F:/InfoFlora/CWR/sig_data/LANDKULT/LV95/data/LANDKULT_BEWE.shp')
# BE_lu<-st_read('F:/InfoFlora/CWR/sig_data/BE_LU/MOPUBE_BBF.shp')#land use in Bern
# sau_lu_cat<-c(8:10,20) #categories corresponding to SAU (without estivage)
# BE_sau_lu<-BE_lu[BE_lu$ART%in%sau_lu_cat,]
# estivage_lu_cat<-c(18,19)
# BE_estivage_lu<-BE_lu[BE_lu$ART%in%estivage_lu_cat,]
# 
# check.val<-st_is_valid(BE_sau_lu)
# for (i in which(check.val ==F)){
#   BE_sau_lu[i,]<-st_make_valid(BE_sau_lu[i,])
# }
# 
# check.val<-st_is_valid(BE_estivage_lu)
# for (i in which(check.val ==F)){
#   BE_estivage_lu[i,]<-st_make_valid(BE_estivage_lu[i,])
# }
# 
# BE_estivage_lu<-ms_dissolve(st_geometry(BE_estivage_lu))
# 
# estivages_union<-st_union(c(st_geometry(BE_estivage_lu),st_geometry(estivage)))
# estivages_union<-ms_dissolve(estivages_union)
# st_write(estivages_union,dsn='F:/InfoFlora/CWR/sig_data/estivages_union.shp',append=F)
# estivages_union<-st_read('F:/InfoFlora/CWR/sig_data/estivages_union.shp')
# 
# spb<-sau_spb[sau_spb$LNF_CODE%in%spb.code$Code,]
# check.val<-st_is_valid(spb)
# for (i in which(check.val ==F)){
#   spb[i,]<-st_make_valid(spb[i,])
# }
# spb<-ms_dissolve(st_geometry(spb))
# 
# st_write(spb,dsn='F:/InfoFlora/CWR/sig_data/spb.shp',append=F)
# 
# SAU<-ms_dissolve(BE_sau_lu)
# SAU<-st_make_valid(SAU)
# st_write(SAU,dsn='F:/InfoFlora/CWR/sig_data/SAU.shp',append=F)
# 
# SAU<-st_difference(SAU,spb)
# SAU<-st_difference(SAU, estivages_union)  
# 
# spb<-st_difference(spb,estivage_union)
# 
# spb<-spb[sau_spb$LNF_CODE%in%spb.code$Code,]
# spb<-ms_dissolve(st_geometry(spb))
# 
# spb<-st_transform(spb,crs=21781)
# st_crs(kt)<-st_crs(spb)
# be<-st_geometry(st_union(kt[which(kt$kuerzel=='BE'),]))
# 
# spb.dissolve<-   ms_dissolve(spb)
# st_write(spb.dissolve,dsn='F:/InfoFlora/CWR/sig_data/LANDKULT/LV95/spb_union.shp',append=F)
# spb.dissolve<-st_read('F:/InfoFlora/CWR/sig_data/LANDKULT/LV95/spb_union.shp')



# spb.2.dissolve<- ms_dissolve()
# st_write(spb.2.dissolve,dsn='F:/InfoFlora/CWR/sig_data/LANDKULT/LV95/spb2_union.shp',append =F)
# spb.2.dissolve<-st_read('F:/InfoFlora/CWR/sig_data/LANDKULT/LV95/spb2_union.shp')

SAU<-st_read('F:/InfoFlora/CWR/sig_data/sauNoSpb.shp', crs = 2056)
spb<-st_read('F:/InfoFlora/CWR/sig_data/spbNoEstivage.shp', crs = 2056)
estivages<-st_read('F:/InfoFlora/CWR/sig_data/estivages_union.shp', crs = 2056)

SAU_simple<-ms_simplify(SAU, keep = 0.05,keep_shapes = T, explode = TRUE)
spb_simple<-ms_simplify(spb, keep = 0.05,keep_shapes = T, explode = TRUE )
mapview(SAU_simple, col.region = 'red')
mapview(SAU_simple, col.region = 'red') + mapview(spb_simple,col.region = 'green')
m1<-mapview(SAU_simple, col.region = 'red') + mapview(spb_simple,col.region = 'green') + 
  mapview(estivages,col.region = 'blue')
mapshot(m1, url = 'F:/InfoFlora/CWR/sig_data/agricultural_area.html')
wv<-st_read('Y:/00_GIS/Schutzgebiete/Wasser_Zugvogelreservate/LV03/wv.shp') # load the inventory of area for migrating birds
st_crs(wv)<-st_crs(be)
wv<-wv[unlist(st_intersects(be,wv)),]
wv<-st_cast(wv,'MULTIPOLYGON')

fm<-st_read('Y:/00_GIS/Inventare/Flachmoore/Flachmoor_LV03/fm.shp') #load the inventory of fens
st_crs(fm)<-st_crs(be)
fm<-fm[unlist(st_intersects(be,fm)),]

hm<-st_read('Y:/00_GIS/Inventare/Hochmoore/Hochmoor_LV03/hm.shp') #load the inventory of fens
st_crs(hm)<-st_crs(be)
hm<-hm[unlist(st_intersects(be,hm)),]

au<-st_read('Y:/00_GIS/Inventare/Auen/Auen_LV03/au.shp')# load the inventory of alluvial area
st_crs(au)<-st_crs(be)
au<-au[unlist(st_intersects(be,au)),]


tww<-st_read('Y:/00_GIS/Inventare/Trockenwiesen_Weiden/tww.shp') #load the inventor of dry pastures and meadows
st_crs(tww)<-st_crs(be)
tww<-tww[unlist(st_intersects(be,tww)),]

aml<-st_read('Y:/00_GIS/Inventare/Amphibienlaichgebiete/Amphibien_LV03/am_l.shp')
st_crs(aml)<-st_crs(be)
aml<-aml[unlist(st_intersects(be,aml)),]


nsg<-st_read('Y:/00_GIS/Schutzgebiete/Pro_natura_2017/NSG_Perimeter/NSG_PN_2017_lv03.shp') #load the protected area of Pro Natura
st_crs(nsg)<-st_crs(be)
nsg<-nsg[unlist(st_intersects(be,nsg)),]

wald.pronat<-st_read('Y:/00_GIS/Schutzgebiete/Pro_natura_2017/NSG_Perimeter/Wald_PN_2017_lv03.shp')
st_crs(wald.pronat)<-st_crs(be)
wald.pronat<-wald.pronat[unlist(st_intersects(be,wald.pronat)),]

wald<-st_read('Y:/00_GIS/Schutzgebiete/waldreservate_delivery_20170329_ae/data/lv03/shp/waldreservate.shp')
st_crs(wald)<-st_crs(be)
wald<-wald[unlist(st_intersects(be,wald)),]

nsg.be<-st_read('F:/InfoFlora/CWR/sig_data/nsg/lv95/data/NSG_NSGZ.shp')
st_crs(nsg.be)<-st_crs(2056)
nsg.be<-st_transform(nsg.be,21781)

bias.file <- raster('F:/InfoFlora/CWR/sig_data/bias.tif')
crs(bias.file)<-crs(nsg.be)
biasxy<-rasterToPoints(bias.file)
biasxy<-data.frame(biasxy)
coordinates(biasxy)<-biasxy[,c('x','y')]
biasxy<-st_as_sf(biasxy)
st_crs(biasxy)<-st_crs(nsg.be)
biasxy.be<-st_intersects(be,biasxy)
biasxy.be<-biasxy[biasxy.be[[1]],]


# total.nsg<-st_combine(st_union(st_geometry(wv),
#                                st_geometry(fm),
#                                st_geometry(hm),
#                                st_geometry(au),
#                                st_geometry(tww),
#                                st_geometry(nsg),
#                                st_geometry(aml),
#                                st_geometry(natpark),
#                                st_geometry(wald),
#                                st_geometry(wald.pronat),
#                                st_geometry(nsg.be)))
# 

# library(velox)
# my.mask<-raster('D:/InfoFlora/InfoFlora_EvaluationTool/BP/V2/data/env/mask.tif')
# bias<-raster('D:/InfoFlora/InfoFlora_EvaluationTool/BP/V2/data/sp/bias.tif')
# bias<-my.mask*bias
# test<-velox(bias)

#Y:/31_Conservation/CWR/Liste_finale_ID2_checklist.txt  #Cloud
# D:/InfoFlora/CWR/Analyses/Liste_finale_ID2_checklist.txt # Local

data<-read.delim('F:/InfoFlora/CWR/Analyses/Liste_finale_ID2_checklist.txt', h=T, sep='\t', stringsAsFactors = F)
exclude<-100 #uncertainty threshold
xcol<-11 # Column with the longitude
ycol<-12 # Column with the latitude
rcol<-26 # Column with the radius
vcol<-19 # Column with validation status of the observation
shpcol<-27 # Column with the geometry (without buffer)
bufcol<-28 # Column with the geometry (with buffer)
date.col<-6 # # Column with the date
date.lim<-'2001-12-31'
dist.pop<-25 # minimal distance between two populations
sp.results<-c()
#save.image(file = 'F:/InfoFlora/CWR/Analyses/envir_PA.Rdata',compress='xz')
load(file = 'F:/InfoFlora/CWR/Analyses/envir_PA.Rdata')
biasxy.be<-biasxy[biasxy.be[[1]],]
biasxy.be$bias<-biasxy.be$bias+1 # create one observation in all cells of bern
spbi<-c()
for(i in seq(1,nrow(spb),by=10000)){
  spbi<-c(spbi,st_union(spb[i:(i+9999),]))
  print(i)
}
spbi<-st_union(spb[1:10000,],by_feature = F)
for (i in 1:nrow(data)){
  print(i)
  load(paste0('F:/InfoFlora/CWR/Analyses/data_sp/',data[i,13]))
  Nsp.protected<-rep(0,5)
  Nsp.be<-0
  if (NCOL(my.sp)>1){
    my.sp<-my.sp[my.sp[,rcol]<=exclude,] #remove unprecise observations
    if(nrow(my.sp)>0){
      my.sp<-my.sp[my.sp[,date.col]>date.lim,] #remove unprecise observations
      if(nrow(my.sp)>0){
        my.sp<-SpatialPointsDataFrame(coords=my.sp[,11:12], data=my.sp) # transforms the dataframe into spatial points dataframe
        
        my.sp<-remove.duplicates(my.sp, zero = dist.pop)
        
        my.sp<-st_as_sf(my.sp)
        st_crs(my.sp)<-st_crs(be)
        
        my.sp.be<-my.sp[lengths(st_intersects(my.sp,be))>0,] #keep  bern observation only
        
        Nsp.be<-nrow(my.sp.be)
        Nsp.protected<-rep(0,5)
        if(nrow(my.sp.be)>=1){
          my.sp.buf<-st_buffer(my.sp.be,dist=my.sp$v_xy_radius)
          my.sp.pts<-my.sp.be
          my.sp.be<-my.sp.pts
          my.sp.rm<-my.sp.be
          
          
          to.rm<-unlist(st_intersects(st_geometry(wv),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          to.rm<-unlist(st_intersects(st_geometry(hm),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          to.rm<-unlist(st_intersects(st_geometry(fm),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          to.rm<-unlist(st_intersects(st_geometry(au),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          to.rm<-unlist(st_intersects(st_geometry(tww),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          to.rm<-unlist(st_intersects(st_geometry(aml),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          my.sp.fedinv<-nrow(my.sp.be)-nrow(my.sp.rm)
          nsp.rm<-nrow(my.sp.rm)
          
          to.rm<-unlist(st_intersects(st_geometry(wald),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          to.rm<-unlist(st_intersects(st_geometry(wald.pronat),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          to.rm<-unlist(st_intersects(st_geometry(nsg),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          my.sp.reserve<-nsp.rm-nrow(my.sp.rm)
          nsp.rm<-nrow(my.sp.rm)
          
          to.rm<-unlist(st_intersects(st_geometry(nsg.be),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          my.sp.cantinv<-nsp.rm-nrow(my.sp.rm)
          nsp.rm<-nrow(my.sp.rm)
          
          to.rm<-unlist(st_intersects(st_geometry(spb.2),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          my.sp.spb2<-nsp.rm-nrow(my.sp.rm)
          nsp.rm<-nrow(my.sp.rm)
          
          to.rm<-unlist(st_intersects(st_geometry(spb),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
          my.sp.spb<-nsp.rm-nrow(my.sp.rm)
          
          
          Nsp.protected<-c(my.sp.fedinv,my.sp.reserve,my.sp.cantinv,my.sp.spb2,my.sp.spb)
          
          
        }
        
      }
    }
  }
  sp.name<-data[i,c(1:4,13:14)]
  Nobs<-nrow(my.sp)
  sp.results<-rbind(sp.results,c(sp.name,Nobs,Nsp.be,Nsp.protected))
  
}
null.model<-function(ntimes= 435,biasxy.be, 
                     wv, hm, fm, au, tww, aml, 
                     wald, wald.pronat, nsg, 
                     nsg.be, 
                     spb.2, 
                     spb){
  my.sp.null<-biasxy.be[sample(1:nrow(biasxy.be),size = ntimes ,prob = biasxy.be$bias),]
  
  to.rm<-unlist(st_intersects(st_geometry(wv),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(hm),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(fm),my.sp.rm)); if (length(to.rm)>0) my.sp.rm<-my.sp.rm[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(au),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(tww),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(aml),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.fedinv<-ntimes-nrow(my.sp.null)
  nsp.rm<-nrow(my.sp.null)
  
  to.rm<-unlist(st_intersects(st_geometry(wald),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(wald.pronat),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(nsg),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.reserve<-nsp.rm-nrow(my.sp.null)
  nsp.rm<-nrow(my.sp.null)
  
  to.rm<-unlist(st_intersects(st_geometry(nsg.be),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.cantinv<-nsp.rm-nrow(my.sp.null)
  nsp.rm<-nrow(my.sp.null)
  
  to.rm<-unlist(st_intersects(st_geometry(spb.2),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.spb2<-nsp.rm-nrow(my.sp.null)
  nsp.rm<-nrow(my.sp.null)
  
  to.rm<-unlist(st_intersects(st_geometry(spb),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.spb<-nsp.rm-nrow(my.sp.null)
  
  Nsp.protected.null<-c(my.sp.fedinv,my.sp.reserve,my.sp.cantinv,my.sp.spb2,my.sp.spb)
  return(Nsp.protected.null)
}

system.time(asd<-null.model(ntimes= 435,biasxy.be, 
                            wv, hm, fm, au, tww, aml, 
                            wald, wald.pronat, nsg, 
                            nsg.be, 
                            spb.2, 
                            spb)
)
cl <- makeCluster(5)
registerDoParallel(cl)
system.time(Nsp.protected.null<-foreach(i=1:1000, .combine = rbind, .packages = c('sf')) %dopar% 
              null.model(ntimes= 435,biasxy.be, 
                         wv, hm, fm, au, tww, aml, 
                         wald, wald.pronat, nsg, 
                         nsg.be, 
                         spb.2, 
                         spb)
)
stopCluster(cl)

colnames(Nsp.protected.null)<-c('Nobs_fedinv','Nobs_reserves','Nobs_cantinv','Nobs_spb2','Nobs_spb')
write.table(Nsp.protected.null,file= 'F:/InfoFlora/CWR/analyses/null_model.txt',quote=F,sep='\t',row.names = F )

Nsp.protected.null<-c()
for (zi in 1:nrep){
  
  my.sp.null<-biasxy.be[sample(1:nrow(biasxy.be),size = nrow(my.sp.be) ,prob = biasxy.be$bias),]
  
  to.rm<-unlist(st_intersects(st_geometry(wv),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(hm),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(fm),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(au),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(tww),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(aml),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.fedinv<-nrow(my.sp.be)-nrow(my.sp.null)
  nsp.rm<-nrow(my.sp.null)
  
  to.rm<-unlist(st_intersects(st_geometry(wald),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(wald.pronat),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(nsg),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.reserve<-nsp.rm-nrow(my.sp.null)
  nsp.rm<-nrow(my.sp.null)
  
  to.rm<-unlist(st_intersects(st_geometry(nsg.be),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.cantinv<-nsp.rm-nrow(my.sp.null)
  nsp.rm<-nrow(my.sp.null)
  
  to.rm<-unlist(st_intersects(st_geometry(spb.2),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.spb2<-nsp.rm-nrow(my.sp.null)
  nsp.rm<-nrow(my.sp.null)
  
  system.time(to.rm<-unlist(st_intersects(st_geometry(spb.2),my.sp.null))); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.spb2<-nsp.rm-nrow(my.sp.null)
  nsp.rm<-nrow(my.sp.null)
  
  
  system.time(to.rm<-unlist(st_intersects(st_geometry(spb),my.sp.null))); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.spb<-nsp.rm-nrow(my.sp.null)
  
  Nsp.protected.null<-rbind(Nsp.protected.null,c(my.sp.fedinv,my.sp.reserve,my.sp.cantinv,my.sp.spb2,my.sp.spb))
}


colnames(sp.results)[7:13]<-c('Nobs_ch','Nobs_BE','Nobs_FedInv','Nobs_Reserve','Nobs_CantInv','Nobs_SPB','Nobs_SAU')

write.table(sp.results,'D:/InfoFlora/CWR/Analyses/spb_analysis_cumsum_pts.txt',sep='\t',)

# Loading
sp.results <- as.data.frame(readxl::read_excel("F:/InfoFlora/CWR/analyses/Resultats_200428/spb_analysis_cumsum_pts.xlsx"))
sp.null<-read.delim('F:/InfoFlora/CWR/analyses/null_model.txt',sep='\t')

#### A FAIRE -> MATCHER LES LIGNE DE SP NULL AVEC SP.RESULTS
prop.null<-sp.null/435

p.val.FedInv<-sapply(sp.results$Prop_FedInv,
                     function(x){if (!is.na(as.numeric(x)))
                       if (as.numeric(x)>0)
                         1-ecdf(prop.null[,1])(x)})
p.val.reserve<-sapply(sp.results$Prop_Reserve,
                      function(x){if (!is.na(as.numeric(x)))
                        if (as.numeric(x)>0)
                          1-ecdf(prop.null[,2])(x)})
p.val.CantInv<-sapply(sp.results$Prop_CantInv,
                      function(x){if (!is.na(as.numeric(x)))
                        if (as.numeric(x)>0)
                          1-ecdf(prop.null[,3])(x)})
p.val.SPB<-sapply(sp.results$Prop_SPB,
                  function(x){if (!is.na(as.numeric(x)))
                    if (as.numeric(x)>0)
                      1-ecdf(prop.null[,4])(x)})
p.val.SAU<-sapply(sp.results$Prob_SAU,
                  function(x){if (!is.na(as.numeric(x)))
                    if (as.numeric(x)>0)
                      1-ecdf(prop.null[,5])(x)})
p.val.protection<-cbind(p.val.FedInv, p.val.reserve,p.val.CantInv,p.val.SPB,p.val.SAU)
write.table(p.val.protection,'F:/InfoFlora/CWR/Analyses/pval_BE.txt',sep='\t',)

#Y:/31_Conservation/CWR/Analyses/PA_analysis  #Cloud
# D:/InfoFlora/CWR/Analyses/PA_analysis # Local

save(data.prot,file=paste0('D:/InfoFlora/CWR/Analyses/PA_analysis'), compress='xz')

write.table(data.prot, file='D:/InfoFlora/CWR/Analyses/PA_analysis.txt', row.names = F,sep="\t")

load(paste0('F:/InfoFlora/CWR/Analyses/PA_analysis'))
