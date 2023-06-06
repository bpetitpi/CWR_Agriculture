########## PA Analysis ##########


# load libraries and data -------------------------------------------------


#library(sp)
library(sf)
library(raster)
library(rgeos)
library(stars)
library(rgdal)
library(dplyr)
library(doParallel)


# #create SAU (from Szerencsits et al. 2018 and https://data.geo.admin.ch/browser/index.html#/collections/ch.blw.landwirtschaftliche-zonengrenzen)
# load('F:/InfoFlora/projet IE/analyses/grid100_sf_with_enviro.Rdata')
# my.grid<-read_stars('F:/InfoFlora/projet IE/analyses/grid100.tif')
# my.grid[my.grid[]>=0]<-NA
# agrizone<-st_read('F:/InfoFlora/CWR/sig_data/zones_agricoles/LWZ_Shape_20201125.shp')
# SAU<-grid_sf[which(grid_sf$SAU == 1),]
# st_geometry(SAU)<-st_geometry(SAU$centro)
# SAU<-SAU[,1]
# agrizone2<-st_transform(agrizone,21781)
# SAU_zone<-st_join(SAU,agrizone2)
# my.values<-SAU_zone
# xy<-st_coordinates(SAU_zone)
# fieldi<-as.numeric(SAU_zone$LZCODE)
# rSAU<-st_rasterize(SAU_zone["LZCODE"], template = my.grid)
# rSAU<-as(rSAU,'raster')
# plot(rSAU==31|rSAU == 41|rSAU == 51|rSAU == 52|rSAU == 53| rSAU == 54)
# plot(rSAU==61)
# estivage<-rSAU==61 
# estivage[estivage==0]<-NA
# lSAU<-rSAU==31|rSAU == 41|rSAU == 51|rSAU == 52|rSAU == 53| rSAU == 54
# lSAU[lSAU == 0]<-NA
# 
# ch_SAU = st_as_sf(lSAU, as_points = FALSE, merge = TRUE)
# ch_estivage <- st_as_sf(estivage, as_points = FALSE, merge = TRUE)

# st_write(ch_SAU,'F:/InfoFlora/CWR/sig_data/ch_sau.shp')
# st_write(ch_estivage,'F:/InfoFlora/CWR/sig_data/ch_estivage.shp')


ch_SAU<-readOGR('F:/InfoFlora/CWR/sig_data', layer='ch_sau',encoding = 'UTF-8', use_iconv = TRUE) # load agricultural areas
ch_estivage<-readOGR('F:/InfoFlora/CWR/sig_data', layer='ch_estivage',encoding = 'UTF-8', use_iconv = TRUE) # load summer grazing areas
wv<-readOGR('Y:/00_GIS/Schutzgebiete/Wasser_Zugvogelreservate/LV03', layer='wv',encoding = 'UTF-8', use_iconv = TRUE) # load the inventory of area for migrating birds
fm<-readOGR('Y:/00_GIS/Inventare/Flachmoore/Flachmoor_LV03', layer='fm',encoding = 'UTF-8', use_iconv = TRUE) #load the inventory of marshes 1
hm<-readOGR('Y:/00_GIS/Inventare/Hochmoore/Hochmoor_LV03', layer='hm',encoding = 'UTF-8', use_iconv = TRUE) #load the inventory of marshes 2
au<-readOGR('Y:/00_GIS/Inventare/Auen/Auen_LV03', layer='au',encoding = 'UTF-8', use_iconv = TRUE)# load the inventory of alluvial area
tww<-readOGR('Y:/00_GIS/Inventare/Trockenwiesen_Weiden', layer='tww',encoding = 'UTF-8', use_iconv = TRUE) #load the inventor of dry pastures and meadows
nsg<-readOGR('Y:/00_GIS/Schutzgebiete/Pro_natura_2017/NSG_Perimeter', layer='NSG_PN_2017_lv03',encoding = 'UTF-8', use_iconv = TRUE) #load the protected area of Pro Natura
aml<-readOGR('Y:/00_GIS/Inventare/Amphibienlaichgebiete/Amphibien_LV03', layer='am_l',encoding = 'UTF-8', use_iconv = TRUE)  #load the inventory of reproduction sites of Amphibians
natpark<-readOGR('Y:/00_GIS/Schutzgebiete/Schweizerischer_Nationalpark', layer='nat_park',encoding = 'UTF-8', use_iconv = TRUE) #load the perimeter of the national park
wald<-readOGR('Y:/00_GIS/Schutzgebiete/waldreservate_delivery_20170329_ae/data/lv03/shp', layer='waldreservate',encoding = 'UTF-8', use_iconv = TRUE) #load the perimeters of the federal forest reserves
kt<-readOGR(paste0('Y:/00_GIS/Grenzen'), layer='kt_ch',encoding = 'UTF-8', use_iconv = TRUE) #load the cantonal administrative boundaries
wald.pronat<-readOGR('Y:/00_GIS/Schutzgebiete/Pro_natura_2017/NSG_Perimeter', layer='Wald_PN_2017_lv03',encoding = 'UTF-8', use_iconv = TRUE) # load the perimeter of other forest reserves

# Harmonize coordinate reference system
crs(natpark)<-crs(tww)
crs(kt)<-crs(tww)
#crs(my.sp)<-crs(tww)
#crs(sp.prio.ch.centro)<-crs(my.sp)

#Y:/31_Conservation/CWR/Liste_finale_ID2_checklist.txt  #Cloud
# D:/InfoFlora/CWR/Analyses/Liste_finale_ID2_checklist.txt # Local

data<-read.delim('F:/InfoFlora/CWR/Analyses/Liste_finale_ID2_checklist.txt', h=T, sep='\t', stringsAsFactors = F) #load the priority CWR list
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

# convert to sf format
nsg<-st_as_sf(nsg)
tww<-st_as_sf(tww)
au<-st_as_sf(au)
fm<-st_as_sf(fm)
hm<-st_as_sf(hm)
natpark<-st_as_sf(natpark)
wv<-st_as_sf(wv)
aml<-st_as_sf(aml)
wald<-st_as_sf(wald)
wald.pronat<-st_as_sf(wald.pronat)
st_crs(wald.pronat)<-st_crs(wald)
ch_SAU<-st_as_sf(ch_SAU)
ch_estivage<-st_as_sf(ch_estivage)


# overlap analysis --------------------------------------------------------


# overlap between protected areas and species distribution
data.prot<-c()
for (i in 1:nrow(data)){
  print(i)
  load(paste0('Y:/31_Conservation/CWR/Analyses/data_sp/',data[i,13]))
  if (NCOL(my.sp)>1){
    my.sp<-my.sp[my.sp[,rcol]<=exclude,] #remove unprecise observations
    if(nrow(my.sp)>0){
      my.sp<-my.sp[my.sp[,date.col]>date.lim,] #remove unprecise observations
      if(nrow(my.sp)>0){
        my.sp<-SpatialPointsDataFrame(coords=my.sp[,11:12], data=my.sp) # transforms the dataframe into spatial points dataframe
        
        crs(my.sp)<-crs(tww)
        my.sp<-remove.duplicates(my.sp, zero = dist.pop) # remove aggregated observations
        spkt<-as.character(over(my.sp,kt)[,7]) #list the canton where the species is present
        if(length(which(spkt=='0'))>0){
          my.sp<-my.sp[-which(spkt=='0'),]
          spkt<-spkt[-which(spkt=='0')]
        }
        Nkt<-length(unique(spkt))
        my.sp.buf<-gBuffer(my.sp,byid=T,width=my.sp$v_xy_radius)
        my.sp.pts<-my.sp
        my.sp<-st_as_sf(my.sp.buf)
        if(nrow(my.sp)>0){
          Nocc<-nrow(my.sp)
          
          over.nsg<-lengths(st_intersects(my.sp,nsg))
          over.tww<-lengths(st_intersects(my.sp,tww))
          over.au<-lengths(st_intersects(my.sp,au))
          over.fm<-lengths(st_intersects(my.sp,fm))
          over.hm<-lengths(st_intersects(my.sp,hm))
          over.natpark<-lengths(st_intersects(my.sp,natpark))
          over.wv<-lengths(st_intersects(my.sp,wv))
          over.aml<-lengths(st_intersects(my.sp,aml))
          over.wald<-lengths(st_intersects(my.sp,wald))
          over.waldPronat<-lengths(st_intersects(my.sp,wald.pronat))
          
          protection<-cbind(over.nsg,over.tww,over.au,over.fm,over.hm, over.natpark,over.wv,over.aml,over.wald, over.waldPronat)
          class(protection)<-'numeric'
          prosum<-apply(protection, 1,sum)
          nosum<-length(which(prosum==0))
          protected<-which(prosum>0)
          protection<-apply(protection, 2,sum)
          
          over.sau<-lengths(st_intersects(st_centroid(my.sp),ch_SAU))
          over.estivage<-lengths(st_intersects(st_centroid(my.sp),ch_estivage))
          if(length(protected)>0){
            data.proti<-unlist(c(data[i,c(2,13,4)],Nocc,Nkt,protection,nosum,
                                 sum(over.sau[-protected]), sum(over.estivage[-protected])))
          }else{
            data.proti<-unlist(c(data[i,c(2,13,4)],Nocc,Nkt,protection,nosum,
                                 sum(over.sau), sum(over.estivage)))
            
          }
          
        }else{data.proti<-c(data[i,c(2,13,4)],rep(0,15))}
      }else{data.proti<-c(data[i,c(2,13,4)],rep(0,15))}
    }else{data.proti<-c(data[i,c(2,13,4)],rep(0,15))}
  }
  if (NCOL(my.sp)==1){
    data.proti<-c(data[i,c(2,13,4)],rep(0,15))
  }
  
  data.prot<-rbind(data.prot,data.proti)
}
colnames(data.prot)[c(4,5,16,17,18)]<-c('Nocc','Nkt','NoProtection', 'InSAU','InEstivage')

save(data.prot,file=paste0('F:/InfoFlora/CWR/Analyses/PA_analysis'), compress='xz')
write.table(data.prot, file='F:/InfoFlora/CWR/Analyses/PA_analysis.txt', row.names = F,sep="\t")



# randomization test ------------------------------------------------------


load('F:/InfoFlora/CWR/Analyses/PA_analysis')
bias.file <- raster('F:/InfoFlora/CWR/sig_data/bias.tif') # load the number of observations in the InfoFlora database
crs(bias.file)<-crs(wv)
# biasxy<-rasterToPoints(bias.file)
# biasxy<-data.frame(biasxy)
# coordinates(biasxy)<-biasxy[,c('x','y')]
# biasxy<-st_as_sf(biasxy)
# st_crs(biasxy)<-st_crs(wv)
# save(biasxy,file ='F:/InfoFlora/CWR/sig_data/bias.Rdata', compress = 'xz')
load('F:/InfoFlora/CWR/sig_data/bias.Rdata')

# function to generate random distributions. It generate 1991 points according to the sampling bias in the InfoFlora database and estimate the proportion in each PA area
null.model<-function(ntimes= 1991,
                     biasxy=biasxy, 
                     wv=wv, 
                     hm=hm, 
                     fm=fm, 
                     au=au, 
                     tww=tww, 
                     aml=aml, 
                     wald=wald, 
                     wald.pronat=wald.pronat, 
                     nsg=nsg, 
                     natpark=natpark, 
                     ch_SAU=ch_SAU, 
                     ch_estivage=ch_estivage
){
  my.sp.null<-dplyr::sample_n(biasxy,ntimes,weight= biasxy$bias)
  
  to.rm<-unlist(st_intersects(st_geometry(wv),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  to.rm<-unlist(st_intersects(st_geometry(hm),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
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
  
  to.rm<-unlist(st_intersects(st_geometry(ch_SAU),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.SAU<-nsp.rm-nrow(my.sp.null)
  nsp.rm<-nrow(my.sp.null)
  
  to.rm<-unlist(st_intersects(st_geometry(ch_estivage),my.sp.null)); if (length(to.rm)>0) my.sp.null<-my.sp.null[-to.rm,]
  my.sp.estivage<-nsp.rm-nrow(my.sp.null)
  
  Nsp.protected.null<-c(my.sp.fedinv,my.sp.reserve,my.sp.SAU,my.sp.estivage)
  return(Nsp.protected.null)
}

system.time(asd<-null.model(ntimes= 1991,biasxy=biasxy, 
                            wv=wv, hm=hm, fm=fm, au=au, tww=tww, aml=aml, 
                            wald=wald, wald.pronat=wald.pronat, nsg=nsg, 
                            natpark=natpark, ch_SAU=ch_SAU, ch_estivage=ch_estivage
))

# generate the random distributions
cl <- makeCluster(6)
registerDoParallel(cl)
system.time(Nsp.protected.null<-foreach(i=1:1000, .combine = rbind, .packages = c('sf')) %dopar% 
              null.model(ntimes= 1991,biasxy=biasxy, 
                         wv=wv, hm=hm, fm=fm, au=au, tww=tww, aml=aml, 
                         wald=wald, wald.pronat=wald.pronat, nsg=nsg, 
                         natpark=natpark, ch_SAU=ch_SAU, ch_estivage=ch_estivage)
)
stopCluster(cl)

colnames(Nsp.protected.null)<-c('Nobs_fedinv','Nobs_reserves','Nobs_SAU','Nobs_estivage')
write.table(Nsp.protected.null,file= 'F:/InfoFlora/CWR/analyses/CH_null_model.txt',quote=F,sep='\t',row.names = F )
Nsp.protected.null<-read.table(file= 'F:/InfoFlora/CWR/analyses/CH_null_model.txt',sep='\t',h =T)
colnames(Nsp.protected.null)<-c('Nobs_fedinv','Nobs_reserves','Nobs_SAU','Nobs_estivage')

# Histogram of the distributions in the different categories
#load(file=paste0('F:/InfoFlora/CWR/Analyses/PA_analysis'))
# sp.results <- as.data.frame(readxl::read_excel("F:/InfoFlora/CWR/analyses/Resultats_200428/spb_analysis_cumsum_pts.xlsx"))
# sp.null<-read.delim('F:/InfoFlora/CWR/analyses/null_model.txt',sep='\t')
sp.results<-as.data.frame(data.prot) # observed distribtion in PA
prop.null<-Nsp.protected.null/1991 # random distribution in PA
par (mfrow = c(1,3),mar = c(2, 2, 2, 2))
hist((prop.null$Nobs_reserves+prop.null$Nobs_fedinv)*100, main = "protected area [%]")
abline(v = quantile ((prop.null$Nobs_reserves+prop.null$Nobs_fedinv)*100,0.95),col = "red")
mtext("a)",adj = 0)
hist(prop.null$Nobs_SAU*100, main = "AA area [%]")
abline(v = quantile (prop.null$Nobs_SAU*100,0.95),col = "red")
mtext("b)",adj = 0)
hist(prop.null$Nobs_estivage*100, main = "SGA area [%]")
abline(v = quantile (prop.null$Nobs_estivage*100,0.95),col = "red")
mtext("c)",adj = 0)

#Tests of better coverage than random in protected area
my.area<-as.numeric(unlist(sp.results$NoProtection))/
  as.numeric(unlist(sp.results$Nocc))
mean.null.reserve<-c(mean(prop.null[,2]+prop.null[,1]),sd(prop.null[,2]+prop.null[,1]))
mean.reserve<-c(1-mean(my.area,na.rm =T), sd(1-my.area,na.rm =T))
nrow.sp<-sample(1:nrow(prop.null),285)
t.test(1-my.area,prop.null[nrow.sp,1]+prop.null[nrow.sp,2],alt = "greater")
p.val.reserve<-sapply(1-my.area,
                      function(x){if (!is.na(x))
                        if (x>0)
                          1-ecdf(prop.null[,2]+prop.null[,1])(x)})
p.val.reserve[which(lengths(p.val.reserve) ==0)]<-1
unprotected.sp<-which(as.numeric(as.character(p.val.reserve))>0.05)
length(unprotected.sp)
mean.reserve<-c(1-mean(my.area,na.rm =T), sd(1-my.area,na.rm =T))
nrow.sp<-sample(1:nrow(prop.null),285)

#Tests of better coverage than random in agricultural area
my.area<-as.numeric(unlist(sp.results$InSAU))/
  as.numeric(unlist(sp.results$Nocc))
mean.null.SAU<-c(mean(prop.null[,3]),sd(prop.null[,3]))
mean.SAU<-c(mean(my.area,na.rm =T), sd(my.area,na.rm =T))
t.test(my.area,prop.null[nrow.sp,3],alt = "greater")
p.val.SAU<-sapply(my.area,
                  function(x){if (!is.na(x))
                    if (x>0)
                      1-ecdf(prop.null[,3])(x)})
inSAU.sp<-which(as.numeric(as.character(p.val.SAU))<0.05)
length(inSAU.sp)

#Tests of better coverage than random in summer grazing area
my.area<-as.numeric(unlist(sp.results$InEstivage))/
  as.numeric(unlist(sp.results$Nocc))
mean.null.estivage<-c(mean(prop.null[,4]),sd(prop.null[,4]))
mean.estivage<-c(mean(my.area,na.rm =T), sd(my.area,na.rm =T))
t.test(my.area,prop.null[nrow.sp,4],alt = "greater")
p.val.estivage<-sapply(my.area,
                       function(x){if (!is.na(x))
                         if (x>0)
                           1-ecdf(prop.null[,4])(x)})
inestivage.sp<-which(as.numeric(as.character(p.val.estivage))<0.05)
length(inestivage.sp)


p.val.protection<-cbind(sp.results$taxon_id_Checklist,sp.results$nom.sans.auteur, p.val.reserve,p.val.SAU,p.val.estivage)
#p.val.protection<-cbind(sp.results, p.val.protection)
write.table(p.val.protection,file ='F:/InfoFlora/CWR/Analyses/pval_CH_Null.txt',sep='\t',)

