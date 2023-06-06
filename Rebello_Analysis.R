# Rasterize PA and short hotspot analysis ---------------------------------

library(fasterize)
library(sf)
library(raster)
library(viridis)
ch<-st_read('F:/InfoFlora/CWR/sig_data/ch_aussengrenze.shp')
lake<-st_read('F:/InfoFlora/CWR/sig_data/see_ch_n.shp', crs =21781)
kt<-st_read('F:/InfoFlora/CWR/sig_data/kt_ch.shp')
kt<-kt[-which(kt$kuerzel==0),]
kt.ch<-kt[-which(kt$kuerzel==0 | kt$kuerzel == "DE" | kt$kuerzel == "IT"),]
communes<-st_read("F:/InfoFlora/SIG/chBoundary/SHAPEFILE_LV95_LN02/swissBOUNDARIES3D_1_3_TLM_HOHEITSGEBIET.shp") %>% 
  st_transform(crs = 21781) %>% 
  dplyr::filter(OBJEKTART == "Gemeindegebiet")

sp.richness<-raster('F:/infoflora/CWR/Analyses/sp_richness2/nsp_obs.tif')
sp.richness[sp.richness==0]<-NA
kt_num<-as.numeric(as.factor(kt.ch$kuerzel))
kt.ch<-cbind(kt.ch,kt_num)
r.kt<-fasterize(kt.ch,sp.richness,field = "kt_num")
# r.kt<-as(st_rasterize(kt.ch$kt_num, st_as_stars(st_bbox(st_as_stars(sp.richness)),dx=1000,dy=1000),
#                       options="ALL_TOUCHED=TRUE"),"Raster")

SAU<-st_read('F:/InfoFlora/CWR/sig_data/sauNoSpb.shp', crs = 2056)
spb<-st_read('F:/InfoFlora/CWR/sig_data/spbNoEstivage.shp', crs = 2056)
estivages<-st_read('F:/InfoFlora/CWR/sig_data/estivages_union.shp', crs = 2056)

wv<-st_read('Y:/00_GIS/Schutzgebiete/Wasser_Zugvogelreservate/LV03/wv.shp') # load the inventory of area for migrating birds
#wv<-st_cast(wv,'MULTIPOLYGON')
fm<-st_read('Y:/00_GIS/Inventare/Flachmoore/Flachmoor_LV03/fm.shp') #load the inventory of fens
hm<-st_read('Y:/00_GIS/Inventare/Hochmoore/Hochmoor_LV03/hm.shp') #load the inventory of fens
au<-st_read('Y:/00_GIS/Inventare/Auen/Auen_LV03/au.shp')# load the inventory of alluvial area
tww<-st_read('Y:/00_GIS/Inventare/Trockenwiesen_Weiden/tww.shp') #load the inventor of dry pastures and meadows
aml<-st_read('Y:/00_GIS/Inventare/Amphibienlaichgebiete/Amphibien_LV03/am_l.shp')
nsg<-st_read('Y:/00_GIS/Schutzgebiete/Pro_natura_2017/NSG_Perimeter/NSG_PN_2017_lv03.shp') #load the protected area of Pro Natura
wald.pronat<-st_read('Y:/00_GIS/Schutzgebiete/Pro_natura_2017/NSG_Perimeter/Wald_PN_2017_lv03.shp')
wald<-st_read('Y:/00_GIS/Schutzgebiete/waldreservate_delivery_20170329_ae/data/lv03/shp/waldreservate.shp')

protected<-c(st_geometry(wv),
             st_geometry(fm),
             st_geometry(hm),
             st_geometry(au),
             st_geometry(tww),
             st_geometry(aml),
             st_geometry(nsg),
             st_geometry(wald.pronat),
             st_geometry(wald)
)
protected<-st_zm(protected)
protected<-st_cast(protected,'MULTIPOLYGON')
my.field<-as.data.frame(rep(1,length(protected)))
st_geometry(my.field)<-protected

# hotspots (national level)
#r.protected<-fasterize(my.field,sp.richness)
r.protected<-as(st_rasterize(my.field,
                             st_as_stars(st_bbox(st_as_stars(sp.richness)),dx=1000,dy=1000),
                             options="ALL_TOUCHED=TRUE"),"Raster")

quantile(sp.richness,0.83)
plot(sp.richness>=106,col='orange')
sp.hotspot<-(sp.richness>=quantile(sp.richness,0.83)) +
  (sp.richness>=quantile(sp.richness,0.70) )
plot(sp.hotspot,col = c('yellow','orange','purple'),colNA='lightgrey',axes =F,box=F)
plot(kt,add=T,col=NA)
plot(r.protected,col = c('yellow'),colNA='lightgrey',axes =F,box=F)
plot(kt,add=T,col=NA)

plot(r.protected*sp.hotspot)
table(getValues(r.protected*sp.hotspot))
table(getValues(sp.hotspot))

# hotspots (cantonal level)
r.protected<-as(st_rasterize(my.field,
                             st_as_stars(st_bbox(st_as_stars(sp.richness)),dx=1000,dy=1000),
                             options="ALL_TOUCHED=TRUE"),"Raster")

hotspots<-setValues(sp.richness,0)
for (i in na.exclude(unique(getValues(r.kt)))){
  print(i)
  kti<-r.kt==i
  sp.richness.i<-kti*sp.richness
  nodsi<- sp.richness.i>=quantile(sp.richness.i[sp.richness.i[]>0],0.83)
  neti<-sp.richness.i>=quantile(sp.richness.i[sp.richness.i[]>0],0.7)
  hotspots<-hotspots + neti + nodsi
}
plot(hotspots,col = c('yellow','orange','purple'),colNA='lightgrey',axes =F,box=F)
plot(st_geometry(lake),border='grey45',col='lightblue',add=T)
plot(kt,add=T,col=NA)
plot(l)
plot(r.protected,col = c("lightgrey",'yellow'),axes =F,box=F)
plot(st_geometry(lake),border='grey45',col='lightblue',add=T)
plot(kt,add=T,col=NA)

plot(r.protected*hotspots)
table(getValues(r.protected*hotspots))
table(getValues(hotspots))

data<-read.table("F:/InfoFlora/CWR/analyses/Resultats_210319/inSAU.txt",h=F)

### merge the distribution of all CWR
library(sf)
library(sp)
library(terra)
library(data.table)
library(dplyr)


exclude<-500 #uncertainty threshold
geo.exclude=5000
xcol<-11 # Column with the longitude
ycol<-12 # Column with the latitude
rcol<-26 # Column with the radius
vcol<-19 # Column with validation status of the observation
shpcol<-27 # Column with the geometry (without buffer)
bufcol<-28 # Column with the geometry (with buffer)
date.col<-6 # # Column with the date
date.lim<-'2001-12-31'
dist.pop<-25 # minimal distance between two populations

ch<-st_read('F:/InfoFlora/CWR/sig_data/ch_aussengrenze.shp', crs =21781)
lake<-st_read('F:/InfoFlora/CWR/sig_data/see_ch_n.shp', crs =21781)
kt<-st_read('F:/InfoFlora/CWR/sig_data/kt_ch.shp', crs =21781)
kt<-kt[-which(kt$kuerzel==0),]
kt.ch<-kt[-which(kt$kuerzel==0 | kt$kuerzel == "DE" | kt$kuerzel == "IT"),]
communes<-st_read("F:/InfoFlora/SIG/chBoundary/SHAPEFILE_LV95_LN02/swissBOUNDARIES3D_1_4_TLM_HOHEITSGEBIET.shp") %>% 
  st_transform(crs = 21781) %>% 
  dplyr::filter(OBJEKTART == "Gemeindegebiet")

sp.richness<-rast('F:/infoflora/CWR/Analyses/sp_richness2/sp_richness3.tif')
sp.richness[sp.richness==0]<-NA
kt_num<-as.numeric(as.factor(kt.ch$kuerzel))
kt.ch<-cbind(kt.ch,kt_num)
#r.kt<-fasterize(kt.ch,sp.richness,field = "kt_num")
r.kt<-rasterize(vect(kt.ch),sp.richness,field = "kt_num")
# r.kt<-as(st_rasterize(kt.ch$kt_num, st_as_stars(st_bbox(st_as_stars(sp.richness)),dx=1000,dy=1000),
#                       options="ALL_TOUCHED=TRUE"),"Raster")

# SAU<-st_read('F:/InfoFlora/CWR/sig_data/sauNoSpb.shp', crs = 2056)
# spb<-st_read('F:/InfoFlora/CWR/sig_data/spbNoEstivage.shp', crs = 2056)
# estivages<-st_read('F:/InfoFlora/CWR/sig_data/estivages_union.shp', crs = 2056)

wv<-st_read('Y:/00_GIS/Schutzgebiete/Wasser_Zugvogelreservate/LV03/wv.shp') # load the inventory of area for migrating birds
#wv<-st_cast(wv,'MULTIPOLYGON')
fm<-st_read('Y:/00_GIS/Inventare/Flachmoore/Flachmoor_LV03/fm.shp') #load the inventory of fens
hm<-st_read('Y:/00_GIS/Inventare/Hochmoore/Hochmoor_LV03/hm.shp') #load the inventory of fens
au<-st_read('Y:/00_GIS/Inventare/Auen/Auen_LV03/au.shp')# load the inventory of alluvial area
tww<-st_read('Y:/00_GIS/Inventare/Trockenwiesen_Weiden/tww.shp') #load the inventor of dry pastures and meadows
aml<-st_read('Y:/00_GIS/Inventare/Amphibienlaichgebiete/Amphibien_LV03/am_l.shp')
nsg<-st_read('Y:/00_GIS/Schutzgebiete/Pro_natura_2017/NSG_Perimeter/NSG_PN_2017_lv03.shp') #load the protected area of Pro Natura
wald.pronat<-st_read('Y:/00_GIS/Schutzgebiete/Pro_natura_2017/NSG_Perimeter/Wald_PN_2017_lv03.shp')
wald<-st_read('Y:/00_GIS/Schutzgebiete/waldreservate_delivery_20170329_ae/data/lv03/shp/waldreservate.shp')

ch_SAU<-st_read('F:/InfoFlora/CWR/sig_data/ch_sau.shp') #
ch_estivage<-st_read('F:/InfoFlora/CWR/sig_data/ch_estivage.shp') #

protected<-c(st_geometry(wv),
             st_geometry(fm),
             st_geometry(hm),
             st_geometry(au),
             st_geometry(tww),
             st_geometry(aml),
             st_geometry(nsg),
             st_geometry(wald.pronat),
             st_geometry(wald)
)
protected<-st_zm(protected)
protected<-st_cast(protected,'MULTIPOLYGON')
my.field<-as.data.frame(rep(1,length(protected)))
st_geometry(my.field)<-protected

# hotspots (national level)
#r.protected<-fasterize(my.field,sp.richness)
# r.protected<-as(st_rasterize(my.field,
#                              st_as_stars(st_bbox(st_as_stars(sp.richness)),dx=1000,dy=1000),
#                              options="ALL_TOUCHED=TRUE"),"Raster")



grid.extent<- st_sf(a = 1:2, geom = st_sfc(st_point(x = c(5.959163471,45.816873303)), # these coordinates are provided by infoflora API
                                           st_point(x = c(10.576523538944393,47.80822259322885))),
                    crs = 4326) %>% 
  st_transform(21781)

#raster mask at a 10 m resolution
r_10m<-st_coordinates(grid.extent) %>% 
  vect %>% 
  rast(res = 10,crs = "epsg:21781")


grid.1km<-st_make_grid(grid.extent, cellsize = 1000)# Generate a 1x1 km grid (in vector format)
grid.1km<-st_sf(id=1:length(grid.1km),grid.1km) #Add a field with ID

data<-read.delim('F:/InfoFlora/CWR/Analyses/Liste_finale_ID2_checklist.txt', h=T, sep='\t', stringsAsFactors = F)

source.dir<-"F:/InfoFlora/InfoFlora_EvaluationTool/BP/V2" #define the directory of the generic structure
sp.files<-list.files("F:/InfoFlora/InfoFlora_EvaluationTool/BP/V2/data/sp/sp_list")

grid.2km<-st_make_grid(grid.extent, cellsize = 2000)# Generate a 1x1 km grid (in vector format)
grid.2km<-st_sf(id=1:length(grid.2km),grid.2km) #Add a field with ID


grid.5km<-st_make_grid(grid.extent, cellsize = 5000)# Generate a 1x1 km grid (in vector format)
grid.5km<-st_sf(id=1:length(grid.5km),grid.5km) #Add a field with ID

grid.10km<-st_make_grid(grid.extent, cellsize = 10000)# Generate a 1x1 km grid (in vector format)
grid.10km<-st_sf(id=1:length(grid.10km),grid.10km) #Add a field with ID

## Generate a hierarchical square structure
# to2<-st_intersects(st_centroid(grid.1km),grid.2km)
# to5<-st_intersects(st_centroid(grid.1km),grid.5km)
# to10<-st_intersects(st_centroid(grid.1km),grid.10km)
# 
# sq2 <- to2[sp.dat$sq1] %>% unlist
# sq5 <- to5[sp.dat$sq1] %>% unlist
# sq10 <- to10[sp.dat$sq1] %>%  unlist
# 
# sp.dat<-cbind(sp.dat,sq2,sq5,sq10)
# fwrite(sp.dat, "F:/InfoFlora/CWR/Analyses/sp_richness2/dist_cwr.csv",quote = F, )

sp.dat<-fread("F:/InfoFlora/CWR/Analyses/sp_richness2/dist_cwr.csv",quote = F, )

fromi<-seq(1,length(protected),by =100)
toi<-c(seq(100,length(protected),by =100),length(protected))
my.seq<-cbind(fromi,toi)

## estimate the proportion of protected area
# prot_area<-rasterize(vect(protected),r_10m,field =1)
# lake_r<-rasterize(vect(lake),r_10m, field = 0,background = 1)
# lake_r[lake_r == 0]<-NA
# prot_area<-prot_area * lake_r
# prot_area_2<-aggregate(prot_area,fact = 200, fun = "sum", na.rm =T)
# prot_area_5<-aggregate(prot_area,fact = 500, fun = "sum", na.rm =T)
# prot_area_10<-aggregate(prot_area,fact = 1000, fun = "sum", na.rm =T)
# land_2<-aggregate(lake_r,fact = 200, fun = "sum", na.rm =T)
# land_5<-aggregate(lake_r,fact = 500, fun = "sum", na.rm =T)
# land_10<-aggregate(lake_r,fact = 1000, fun = "sum", na.rm =T)
# prop_prot_2<-round(100* prot_area_2 / land_2, 1)
# prop_prot_5 <- round(100* prot_area_5 / land_5, 1)
# prop_prot_10 <- round(100* prot_area_10 / land_10, 1)

# writeRaster(prop_prot_2,"F:/InfoFlora/CWR/sig_data/pro_prot_2.tif")
# writeRaster(prop_prot_5,"F:/InfoFlora/CWR/sig_data/pro_prot_5.tif")
# writeRaster(prop_prot_10,"F:/InfoFlora/CWR/sig_data/pro_prot_10.tif")

prop_prot_2<-rast("F:/InfoFlora/CWR/sig_data/pro_prot_2.tif")
prop_prot_5<-rast("F:/InfoFlora/CWR/sig_data/pro_prot_5.tif")
prop_prot_10<-rast("F:/InfoFlora/CWR/sig_data/pro_prot_10.tif")

prot2v<-extract(prop_prot_2,vect(st_centroid(grid.2km)))$layer
prot2v[is.nan(prot2v)]<-0
prot5v<-extract(prop_prot_5,vect(st_centroid(grid.5km)))$layer
prot5v[is.nan(prot5v)]<-0
prot10v<-extract(prop_prot_10,vect(st_centroid(grid.10km)))$layer
prot10v[is.nan(prot10v)]<-0


# Rebello complementary analyisis -----------------------------------------


### National 10 km
sq_kept<-c()
ntax<-length(unique(sp.dat$Intial_TaxonID))
sp.dati<-sp.dat
while(ntax > 0){ #iterative complementarity analysis 
  sp_sq<-sp.dati[, .N, by=.(Intial_TaxonID, sq10)]
  sq_summary<-sp_sq[, .N, by = .(sq10)]
  sq_kepti<-sq_summary[which(sq_summary$N == max(sq_summary$N))]
  Nsp.tot<-length(unique(sp.dat[sq10 %in% sq_kepti$sq10,]$Intial_TaxonID))
  # Nsp.tot<-sp.dat[sp.dat$sq10 %in% sq_kepti$sq10,] %>% 
  #   .[, .N, by=.(Intial_TaxonID, sq10)] %>% 
  #   .[, .N, by=.(sq10)]
  sq_id<-sq_kepti[order(Nsp.tot,decreasing = TRUE),]$sq10
  sq_kepti<-c(sq_id,
              length(unique(sp.dati[sq10%in%sq_id,]$Intial_TaxonID)),
              length(unique(sp.dat[sq10%in%sq_id,]$Intial_TaxonID)))
  
  tax2rm<-unique(sp.dati[sq10 == sq_id,]$Intial_TaxonID)
  sp.dati<-sp.dati[!(Intial_TaxonID %in% tax2rm),]
  sq_kept<-rbind(sq_kept,sq_kepti)
  ntax<-length(unique(sp.dati$Intial_TaxonID))
}
sq_kept <- as.data.table(sq_kept)
names(sq_kept)<-c('sq10','N_Sp_Rebello','N_Sp_CWR')
best_sq10<-cbind(grid.10km[sq_kept$sq10,],
                 sq_kept[,2:3])
hs<-best_sq10
hs_center<-which(lengths(st_intersects(st_centroid(hs),communes))>0)
hs_name<-rep("hs_name",nrow(hs))
hs_name[hs_center]<-st_intersection(st_centroid(hs),communes)$NAME
hs_name[-hs_center]<-communes$NAME[st_nearest_feature(st_centroid(hs),communes)[-hs_center]]
hs_protection<-extract(prop_prot_10,vect(st_centroid(hs)))$layer
hs_protection[is.nan(hs_protection)]<-0
hs_estivage<-st_intersection(hs,ch_estivage) %>% 
  mutate(area_estivage = as.numeric(st_area(.))) %>% 
  group_by(id) %>% 
  summarise(hs_estivage = sum(area_estivage/1000000)) %>% 
  st_drop_geometry()
hs_sau<-st_intersection(hs,ch_SAU) %>% 
  mutate(area_sau = as.numeric(st_area(.))) %>% 
  group_by(id) %>% 
  summarise(hs_sau = sum(area_sau/1000000)) %>% 
  st_drop_geometry()

best_sq10$name<-hs_name
best_sq10$prot_area<-hs_protection
best_sq10$estivage_area<-0
best_sq10$sau_area<-0
best_sq10$estivage_area[!is.na(match(best_sq10$id,hs_estivage$id))]<-round(hs_estivage$hs_estivage[na.exclude(match(best_sq10$id,hs_estivage$id))],1)
best_sq10$sau_area[!is.na(match(best_sq10$id,hs_sau$id))]<-round(hs_sau$hs_sau[na.exclude(match(best_sq10$id,hs_sau$id))],1)

### National 5 km
sq_kept<-c()
ntax<-length(unique(sp.dat$Intial_TaxonID))
sp.dati<-sp.dat
while(ntax > 0){
  sp_sq<-sp.dati[, .N, by=.(Intial_TaxonID, sq5)]
  sq_summary<-sp_sq[, .N, by = .(sq5)]
  sq_kepti<-sq_summary[which(sq_summary$N == max(sq_summary$N))]
  Nsp.tot<-length(unique(sp.dat[sq5 %in% sq_kepti$sq5,]$Intial_TaxonID))
  # Nsp.tot<-sp.dat[sp.dat$sq5 %in% sq_kepti$sq5,] %>% 
  #   .[, .N, by=.(Intial_TaxonID, sq5)] %>% 
  #   .[, .N, by=.(sq5)]
  sq_id<-sq_kepti[order(Nsp.tot,decreasing = TRUE),]$sq5
  sq_kepti<-c(sq_id,
              length(unique(sp.dati[sq5%in%sq_id,]$Intial_TaxonID)),
              length(unique(sp.dat[sq5%in%sq_id,]$Intial_TaxonID)))
  
  tax2rm<-unique(sp.dati[sq5 == sq_id,]$Intial_TaxonID)
  sp.dati<-sp.dati[!(Intial_TaxonID %in% tax2rm),]
  sq_kept<-rbind(sq_kept,sq_kepti)
  ntax<-length(unique(sp.dati$Intial_TaxonID))
}
sq_kept <- as.data.table(sq_kept)
names(sq_kept)<-c('sq5','N_Sp_Rebello','N_Sp_CWR')
best_sq5<-cbind(grid.5km[sq_kept$sq5,],
                sq_kept[,2:3])
hs<-best_sq5
hs_center<-which(lengths(st_intersects(st_centroid(hs),communes))>0)
hs_name<-rep("hs_name",nrow(hs))
hs_name[hs_center]<-st_intersection(st_centroid(hs),communes)$NAME
hs_name[-hs_center]<-communes$NAME[st_nearest_feature(st_centroid(hs),communes)[-hs_center]]
hs_protection<-extract(prop_prot_5,vect(st_centroid(hs)))$layer
hs_protection[is.nan(hs_protection)]<-0
hs_estivage<-st_intersection(hs,ch_estivage) %>% 
  mutate(area_estivage = as.numeric(st_area(.))) %>% 
  group_by(id) %>% 
  summarise(hs_estivage = sum(area_estivage/250000)) %>% 
  st_drop_geometry()
hs_sau<-st_intersection(hs,ch_SAU) %>% 
  mutate(area_sau = as.numeric(st_area(.))) %>% 
  group_by(id) %>% 
  summarise(hs_sau = sum(area_sau/250000)) %>% 
  st_drop_geometry()

best_sq5$name<-hs_name
best_sq5$prot_area<-hs_protection
best_sq5$estivage_area<-0
best_sq5$sau_area<-0
best_sq5$estivage_area[!is.na(match(best_sq5$id,hs_estivage$id))]<-round(hs_estivage$hs_estivage[na.exclude(match(best_sq5$id,hs_estivage$id))],1)
best_sq5$sau_area[!is.na(match(best_sq5$id,hs_sau$id))]<-round(hs_sau$hs_sau[na.exclude(match(best_sq5$id,hs_sau$id))],1)


### National 2 km
sq_kept<-c()
ntax<-length(unique(sp.dat$Intial_TaxonID))
sp.dati<-sp.dat
while(ntax > 0){
  sp_sq<-sp.dati[, .N, by=.(Intial_TaxonID, sq2)]
  sq_summary<-sp_sq[, .N, by = .(sq2)]
  sq_kepti<-sq_summary[which(sq_summary$N == max(sq_summary$N))]
  Nsp.tot<-length(unique(sp.dat[sq2 %in% sq_kepti$sq2,]$Intial_TaxonID))
  # Nsp.tot<-sp.dat[sp.dat$sq2 %in% sq_kepti$sq2,] %>% 
  #   .[, .N, by=.(Intial_TaxonID, sq2)] %>% 
  #   .[, .N, by=.(sq2)]
  sq_id<-sq_kepti[order(Nsp.tot,decreasing = TRUE),]$sq2
  sq_kepti<-c(sq_id,
              length(unique(sp.dati[sq2%in%sq_id,]$Intial_TaxonID)),
              length(unique(sp.dat[sq2%in%sq_id,]$Intial_TaxonID)))
  
  tax2rm<-unique(sp.dati[sq2 == sq_id,]$Intial_TaxonID)
  sp.dati<-sp.dati[!(Intial_TaxonID %in% tax2rm),]
  sq_kept<-rbind(sq_kept,sq_kepti)
  ntax<-length(unique(sp.dati$Intial_TaxonID))
}
sq_kept <- as.data.table(sq_kept)
names(sq_kept)<-c('sq2','N_Sp_Rebello','N_Sp_CWR')
best_sq2<-cbind(grid.2km[sq_kept$sq2,],
                sq_kept[,2:3])

hs<-best_sq2
hs_center<-which(lengths(st_intersects(st_centroid(hs),communes))>0)
hs_name<-rep("hs_name",nrow(hs))
hs_name[hs_center]<-st_intersection(st_centroid(hs),communes)$NAME
hs_name[-hs_center]<-communes$NAME[st_nearest_feature(st_centroid(hs),communes)[-hs_center]]
hs_protection<-extract(prop_prot_2,vect(st_centroid(hs)))$layer
hs_protection[is.nan(hs_protection)]<-0
hs_estivage<-st_intersection(hs,ch_estivage) %>% 
  mutate(area_estivage = as.numeric(st_area(.))) %>% 
  group_by(id) %>% 
  summarise(hs_estivage = sum(area_estivage/40000)) %>% 
  st_drop_geometry()
hs_sau<-st_intersection(hs,ch_SAU) %>% 
  mutate(area_sau = as.numeric(st_area(.))) %>% 
  group_by(id) %>% 
  summarise(hs_sau = sum(area_sau/40000)) %>% 
  st_drop_geometry()

best_sq2$name<-hs_name
best_sq2$prot_area<-hs_protection
best_sq2$estivage_area<-0
best_sq2$sau_area<-0
best_sq2$estivage_area[!is.na(match(best_sq2$id,hs_estivage$id))]<-round(hs_estivage$hs_estivage[na.exclude(match(best_sq2$id,hs_estivage$id))],1)
best_sq2$sau_area[!is.na(match(best_sq2$id,hs_sau$id))]<-round(hs_sau$hs_sau[na.exclude(match(best_sq2$id,hs_sau$id))],1)


write.table(st_drop_geometry(best_sq5), file = "F:/InfoFlora/CWR/Analyses/Resultats_hotspots/best_sq5.txt",sep ="\t",row.names =F, quote =F)
write.table(st_drop_geometry(best_sq10), file = "F:/InfoFlora/CWR/Analyses/Resultats_hotspots/best_sq10.txt",sep ="\t",row.names =F, quote =F)
write.table(st_drop_geometry(best_sq2), file = "F:/InfoFlora/CWR/Analyses/Resultats_hotspots/best_sq2.txt",sep ="\t",row.names =F, quote =F)

### Create map of richness at 1 km
sp_sq<-sp.dat[, .N, by=.(Intial_TaxonID, sq1)]
sq_summary<-sp_sq[, .N, by = .(sq1)]
head(sq_summary)
sp.richness<-grid.1km
sp.richness$Nsp<-0
sp.richness$Nsp[sq_summary$sq1]<-sq_summary$N
sp.richness<-filter(sp.richness,Nsp >0)
sp.rich.class<-mutate(sp.richness,my.class = cut (sp.richness$Nsp,c(1,20,40,60,80,106),include.lowest = T, right = TRUE))

png('F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/nsp_obs.png',
    width = 15, height = 10, units = "cm",res = 300, pointsize = 20)
# png('F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/nsp_obs.png',
#     width = 2048, height = 1500, pointsize = 100)

p.ch <- ggplot() +geom_sf(data = ch,
                          color = NA,
                          fill ="White",
                          lwd = 1
)+ geom_sf(data = sp.rich.class,
           aes(fill = my.class),
           color = NA) +
  scale_fill_brewer(palette = "YlGn",
                    name = "Nsp")
p.kant<-geom_sf(data = kt.ch,
                fill = NA,
                color = "gray70",
                lwd = 0.25)
p.ch2<-geom_sf(data = ch,
               color = "gray20",
               fill =NA,
               lwd = 0.25)

p.lake<-geom_sf(data = lake, 
                fill = "lightblue",
                color = NA)
p.hotspot<-geom_sf(data = best_sq2,
                   color = "red",
                   fill = NA,
                   lwd = 0.25)

p.hotspot5<-geom_sf(data = best_sq5,
                    color = "red",
                    fill = NA,
                    lwd = 0.25)

p.ch + p.kant + p.ch2 + p.lake + p.hotspot + ggtitle("Observed number of species")+ theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
#p.ch + p.kant + p.ch2 + p.lake + p.hotspot5 + theme_minimal()
dev.off()

png('F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/nsp_obs_raster.png',width = 2048, height = 1500,pointsize =50)
r.nsp<-rast('F:/infoflora/CWR/Analyses/sp_richness2/sp_richness3.tif')
plot(r.nsp,col=YlGn(100), axes =F, mar = c(1,1, 2, 4), main = "Observed number of species",cex.main=1.5,plg=list(cex=1.25))
# plot(r.nsp,col=YlGn(5), axes =F, mar = c(1,1, 2, 4), main = "Observed number of species",
#      type = "interval", breaks=c(1,20,40,60,80,106), plg=list(legend=c("1-20", "21-40","41-60","61-80","81-106")))
plot(st_geometry(kt),col = NA,add =T,border = "gray70")
plot(st_geometry(ch),add = T, col = NA)
plot(st_geometry(lake), add = T, col = "lightblue",border = NA)
plot(st_geometry(best_sq2),col = NA, border = "red",add =T,lwd =2)
sbar(50000,divs = 2, type = "bar", label = c(0,25, 50),below = "kilometers")
north(type = 2, xy = "topleft")
dev.off()

### leaflet
m1<-mapview(best_sq10, zcol = "N_Sp_CWR", layer.name = "CWR National Hotspot 10km")
m2<-mapview(best_sq5, zcol = "N_Sp_CWR", layer.name = "CWR National Hotspot 5km")
m3<-mapview(best_sq2, zcol = "N_Sp_CWR", layer.name = "CWR National Hotspot 2km")

m<-m1+m2+m3
mapshot(m, url = paste0(getwd(), "/Resultats_hotspots/CWR_Hotspots_National.html"))


### cantonal 
my.kt<-sort(unique(sp.dat$sp_kt)[-c(27,28)])
### cantonal 10 km

for (i in 1:length(my.kt)){
  print(my.kt[i])
  sp.dati<-sp.dat[sp_kt == my.kt[i]]
  sq_kept<-c()
  ntax<-length(unique(sp.dat$Intial_TaxonID))
  while(ntax > 0){
    sp_sq<-sp.dati[, .N, by=.(Intial_TaxonID, sq10)]
    sq_summary<-sp_sq[, .N, by = .(sq10)]
    sq_kepti<-sq_summary[which(sq_summary$N == max(sq_summary$N))]
    Nsp.tot<-sp.dat[sp.dat$sq10 %in% sq_kepti$sq10,] %>% 
      .[, .N, by=.(Intial_TaxonID, sq10)] %>% 
      .[, .N, by=.(sq10)]
    sq_kepti<-cbind(sq_kepti[order(Nsp.tot$N,decreasing = TRUE),][1],
                    Nsp.tot[order(Nsp.tot$N, decreasing = TRUE),][1,2])
    
    tax2rm<-unique(sp.dati[sq10 == sq_kepti$sq10,]$Intial_TaxonID)
    sp.dati<-sp.dati[!(Intial_TaxonID %in% tax2rm),]
    sq_kept<-rbind(sq_kept,sq_kepti)
    ntax<-length(unique(sp.dati$Intial_TaxonID))
  }
  names(sq_kept)[c(2,3)]<-c('N_Sp_Rebello','N_Sp_CWR')
  best_sq10<-cbind(grid.10km[sq_kept$sq10,],
                   sq_kept[,2:3])
  if (i == 1){
    mi<-mapview(best_sq10, zcol = "N_Sp_Rebello", 
                layer.name = paste0(my.kt[i], 
                                    " Cantonal CWR Hotspot 10km"))    
  }else{
    mii<-mapview(best_sq10, zcol = "N_Sp_Rebello", 
                 layer.name = paste0(my.kt[i], 
                                     " Cantonal CWR Hotspot 10km"),
                 hide = TRUE)
    mi<-mi + mii
  }
}
mapshot(mi, url = paste0(getwd(), "/Resultats_hotspots/Cantonal_10km.html"))

### cantonal 5 km

for (i in 1:length(my.kt)){
  print(my.kt[i])
  sp.dati<-sp.dat[sp_kt == my.kt[i]]
  sq_kept<-c()
  ntax<-length(unique(sp.dat$Intial_TaxonID))
  while(ntax > 0){
    sp_sq<-sp.dati[, .N, by=.(Intial_TaxonID, sq5)]
    sq_summary<-sp_sq[, .N, by = .(sq5)]
    sq_kepti<-sq_summary[which(sq_summary$N == max(sq_summary$N))]
    Nsp.tot<-sp.dat[sp.dat$sq5 %in% sq_kepti$sq5,] %>% 
      .[, .N, by=.(Intial_TaxonID, sq5)] %>% 
      .[, .N, by=.(sq5)]
    sq_kepti<-cbind(sq_kepti[order(Nsp.tot$N,decreasing = TRUE),][1],
                    Nsp.tot[order(Nsp.tot$N, decreasing = TRUE),][1,2])
    
    tax2rm<-unique(sp.dati[sq5 == sq_kepti$sq5,]$Intial_TaxonID)
    sp.dati<-sp.dati[!(Intial_TaxonID %in% tax2rm),]
    sq_kept<-rbind(sq_kept,sq_kepti)
    ntax<-length(unique(sp.dati$Intial_TaxonID))
  }
  names(sq_kept)[c(2,3)]<-c('N_Sp_Rebello','N_Sp_CWR')
  best_sq5<-cbind(grid.5km[sq_kept$sq5,],
                  sq_kept[,2:3])
  if (i == 1){
    mi<-mapview(best_sq5, zcol = "N_Sp_Rebello", 
                layer.name = paste0(my.kt[i], 
                                    " Cantonal CWR Hotspot 5km"))    
  }else{
    mii<-mapview(best_sq5, zcol = "N_Sp_Rebello", 
                 layer.name = paste0(my.kt[i], 
                                     " Cantonal CWR Hotspot 5km"),
                 hide = TRUE)
    mi<-mi + mii
  }
}
mapshot(mi, url = paste0(getwd(), "/Resultats_hotspots/Cantonal_5km.html"))

### cantonal 2 km

for (i in 1:length(my.kt)){
  print(my.kt[i])
  sp.dati<-sp.dat[sp_kt == my.kt[i]]
  sq_kept<-c()
  ntax<-length(unique(sp.dat$Intial_TaxonID))
  while(ntax > 0){
    sp_sq<-sp.dati[, .N, by=.(Intial_TaxonID, sq2)]
    sq_summary<-sp_sq[, .N, by = .(sq2)]
    sq_kepti<-sq_summary[which(sq_summary$N == max(sq_summary$N))]
    Nsp.tot<-sp.dat[sp.dat$sq2 %in% sq_kepti$sq2,] %>% 
      .[, .N, by=.(Intial_TaxonID, sq2)] %>% 
      .[, .N, by=.(sq2)]
    sq_kepti<-cbind(sq_kepti[order(Nsp.tot$N,decreasing = TRUE),][1],
                    Nsp.tot[order(Nsp.tot$N, decreasing = TRUE),][1,2])
    
    tax2rm<-unique(sp.dati[sq2 == sq_kepti$sq2,]$Intial_TaxonID)
    sp.dati<-sp.dati[!(Intial_TaxonID %in% tax2rm),]
    sq_kept<-rbind(sq_kept,sq_kepti)
    ntax<-length(unique(sp.dati$Intial_TaxonID))
  }
  names(sq_kept)[c(2,3)]<-c('N_Sp_Rebello','N_Sp_CWR')
  best_sq2<-cbind(grid.2km[sq_kept$sq2,],
                  sq_kept[,2:3])
  if (i == 1){
    mi<-mapview(best_sq2, zcol = "N_Sp_Rebello", 
                layer.name = paste0(my.kt[i], 
                                    " Cantonal CWR Hotspot 2km"))    
  }else{
    mii<-mapview(best_sq2, zcol = "N_Sp_Rebello", 
                 layer.name = paste0(my.kt[i], 
                                     " Cantonal CWR Hotspot 2km"),
                 hide = TRUE)
    mi<-mi + mii
  }
}
mapshot(mi, url = paste0(getwd(), "/Resultats_hotspots/Cantonal_2km.html"))

### cantonal 1 km

for (i in 1:length(my.kt)){
  print(my.kt[i])
  sp.dati<-sp.dat[sp_kt == my.kt[i]]
  sq_kept<-c()
  ntax<-length(unique(sp.dat$Intial_TaxonID))
  while(ntax > 0){
    sp_sq<-sp.dati[, .N, by=.(Intial_TaxonID, sq1)]
    sq_summary<-sp_sq[, .N, by = .(sq1)]
    sq_kepti<-sq_summary[which(sq_summary$N == max(sq_summary$N))]
    Nsp.tot<-sp.dat[sp.dat$sq1 %in% sq_kepti$sq1,] %>% 
      .[, .N, by=.(Intial_TaxonID, sq1)] %>% 
      .[, .N, by=.(sq1)]
    sq_kepti<-cbind(sq_kepti[order(Nsp.tot$N,decreasing = TRUE),][1],
                    Nsp.tot[order(Nsp.tot$N, decreasing = TRUE),][1,2])
    
    tax2rm<-unique(sp.dati[sq1 == sq_kepti$sq1,]$Intial_TaxonID)
    sp.dati<-sp.dati[!(Intial_TaxonID %in% tax2rm),]
    sq_kept<-rbind(sq_kept,sq_kepti)
    ntax<-length(unique(sp.dati$Intial_TaxonID))
  }
  names(sq_kept)[c(2,3)]<-c('N_Sp_Rebello','N_Sp_CWR')
  best_sq1<-cbind(grid.1km[sq_kept$sq1,],
                  sq_kept[,2:3])
  if (i == 1){
    mi<-mapview(best_sq1, zcol = "N_Sp_Rebello", 
                layer.name = paste0(my.kt[i], 
                                    " Cantonal CWR Hotspot 1km"))    
  }else{
    mii<-mapview(best_sq1, zcol = "N_Sp_Rebello", 
                 layer.name = paste0(my.kt[i], 
                                     " Cantonal CWR Hotspot 1km"),
                 hide = TRUE)
    mi<-mi + mii
  }
}
mapshot(mi, url = paste0(getwd(), "/Resultats_hotspots/Cantonal_1km.html"))




sq_kept<-c()
ntax<-length(unique(sp.dat$Intial_TaxonID))
sp.dati<-sp.dat[sp_kt == "VD"]
while(ntax > 0){
  sp_sq<-sp.dati[, .N, by=.(Intial_TaxonID, sq10)]
  sq_summary<-sp_sq[, .N, by = .(sq10)]
  sq_kepti<-sq_summary[which(sq_summary$N == max(sq_summary$N))]
  Nsp.tot<-sp.dat[sp.dat$sq10 %in% sq_kepti$sq10,] %>% 
    .[, .N, by=.(Intial_TaxonID, sq10)] %>% 
    .[, .N, by=.(sq10)]
  sq_kepti<-cbind(sq_kepti[order(Nsp.tot$N,decreasing = TRUE),][1],
                  Nsp.tot[order(Nsp.tot$N, decreasing = TRUE),][1,2])
  
  tax2rm<-unique(sp.dati[sq10 == sq_kepti$sq10,]$Intial_TaxonID)
  sp.dati<-sp.dati[!(Intial_TaxonID %in% tax2rm),]
  sq_kept<-rbind(sq_kept,sq_kepti)
  ntax<-length(unique(sp.dati$Intial_TaxonID))
}
names(sq_kept)[c(2,3)]<-c('N_Sp_Rebello','N_Sp_CWR')
best_sq10<-cbind(grid.10km[sq_kept$sq10,],
                 sq_kept[,2:3])


sq_kept<-c()
ntax<-length(unique(sp.dat$v_accepted_taxon_id))
sp.dati<-sp.dat
while(ntax > 0){
  sp_sq<-sp.dati[, .N, by=.(v_accepted_taxon_id, sq5)]
  sq_summary<-sp_sq[, .N, by = .(sq5)]
  sq_kepti<-sq_summary[which(sq_summary$N == max(sq_summary$N))][1]
  
  tax2rm<-unique(sp.dati[sq5 == sq_kepti$sq5,]$v_accepted_taxon_id)
  sp.dati<-sp.dati[!(v_accepted_taxon_id %in% tax2rm),]
  sq_kept<-rbind(sq_kept,sq_kepti)
  ntax<-length(unique(sp.dati$v_accepted_taxon_id))
}

best_sq5<-grid.5km[sq_kept$sq5,]
best_sq5<-cbind(best_sq5,sq_kept$N)
mapview(best_sq5, zcol = "sq_kept.N")
