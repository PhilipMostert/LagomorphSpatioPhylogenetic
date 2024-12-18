library(rworldmap)
library(ggplot2)
library(inlabru)
library(sf)
library(ape)
library(terra)

#proj <- INLA::inla.CRS('globe') #"EPSG:4326"
##CHANGE TO KM
#proj <- "+proj=aea +lat_1=52.5 +lat_2=-10 +lon_0=70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"

#proj <- inlabru::fm_crs('mollweide_globe')
#This proj is new
proj <- '+proj=moll +lon_0=40 +x_0=0 +y_0=0 +ellps=sphere +units=km +no_defs'
#proj <- '+proj=laea +lat_0=30 +lon_0=20 +x_0=0 +y_0=0 +ellps=sphere +datum=WGS84 +units=km +no_defs'
#proj <- '+proj=igh +lat_0=30 +lon_0=20 +x_0=0 +y_0=0 +ellps=sphere +units=km +no_defs'
mam <- read.nexus("data-raw/Bunnies/Phylo.txt")
grep("Lago", mam[[1]]$node.label)
hares <- extract.clade(mam[[1]], mam[[1]]$node.label[c(673)])

Lepus <- sf::st_read('~/Downloads/redlist_species_data_8c436bf9-8db4-4cc8-8657-3662399ac09b/data_0.shp')
Ochotona <- sf::st_read('~/Downloads/redlist_species_data_cdd72903-fff8-4ec9-952e-23935fd863ba/')
Buck <- sf::read_sf('~/Downloads/redlist_species_data_71ad295e-319f-44dc-ac01-012e17531499/data_0.shp')
Lago <- rbind(Lepus, Ochotona, Buck)


world <- getMap(resolution = 'less islands') # rnaturalearth::ne_countries(returnclass = 'sf')
studyArea <- world[world$continent %in% c('Africa', 'Eurasia'),] #c('Africa, 'Eurasia)
studyArea <- as(studyArea, 'sf')
#studyArea <- st_transform(studyArea, proj)
#New from Iceland and including Iceland
studyArea <- studyArea[!studyArea$NAME %in% c('Greenland', 'Madagascar', 'Iceland', 'East Timor', 'Malaysia', 'Brunei','Indonesia', 'Philippines', 'Japan', 'Cyprus'),]

specData <- st_filter(st_transform(Lago,crs = st_crs(studyArea)), st_make_valid(studyArea))#st_filter(st_transform(Lago,crs = proj), studyArea)
specData$species <- gsub(' ', '_', specData$SCI_NAME)
specData <- st_transform(specData, proj)

studyArea <- st_transform(studyArea, proj)
studyArea <- rmapshaper::ms_filter_islands(studyArea, 90000)

hares <- drop.tip(hares, c(hares$tip.label[!hares$tip.label %in% specData$species], 'Lepus_hainanus'))
hares <- drop.tip(mam[[1]], mam[[1]]$tip.label[!mam[[1]]$tip.label%in% c(hares$tip.label, 'Tragelaphus_imberbis')])
specData <- specData[specData$species %in% hares$tip.label,]

prec <- mean(geodata::worldclim_global(var = 'prec', res = 5, path = 'data-raw/'))
prec <- terra::project(prec, proj)
temp <- mean(geodata::worldclim_global(var = 'tavg', res = 5, path = 'data-raw/'))
temp <- terra::project(temp, proj)

worldPoly <- st_make_grid(studyArea, cellsize = 100) #100
world <- st_filter(st_as_sf(worldPoly), studyArea)
st_geometry(world) <- 'geometry'
world <- st_cast(world, 'POINT')
world <- st_filter(st_as_sf(world), studyArea)

world$prec <- inlabru::eval_spatial(data = prec, world)
world$prec <- inlabru::bru_fill_missing(prec, world, values = world$prec)

#world$prec <- inlabru::eval_spatial(data = prec, world)
#world$prec <- inlabru::bru_fill_missing(prec, world, values = world$prec)

world$temp <- inlabru::eval_spatial(data = temp, world)
world$temp <- inlabru::bru_fill_missing(temp, world, values = world$temp)

specPres <- list()
specANS <- list()

for (spec in unique(specData$species)) {
  
  indx <- unlist(sf::st_intersects(specData[specData$species == spec,], world)) #world #worldpoly[world] #st_covers
  specp <- world
  specp$species <- spec
  specp$speciesNum <- which(spec == unique(specData$species))
  specp$speciesNumIndex <-  which(spec == unique(specData$species))
  specp[indx, 'Y'] <- 1
  specp[is.na(specp[['Y']]), 'Y'] <- 0
  specPres[[spec]] <- specp
  
}

numTips <- length(unique(specData$SCI_NAME))

for (anc in rev(hares$node.label[-1])) {
  
  specp <- world
  specp$species <- gsub(' ','_', anc)
  specp$speciesNum <- which(anc == rev(hares$node.label)) + numTips
  specp$speciesNumIndex <-  NA#which(anc == rev(hares$node.label)) + numTips
  specp$Y <- NA
  specANS[[anc]] <- specp
   
}

speciesDataSF <- do.call(rbind, c(specPres, specANS))

mesh <- fmesher::fm_rcdt_2d_inla(loc = world,#coords[,c(1,2)],
                                 boundary = studyArea, #st_combine?
                                 #loc.domain = coords[,1:2],
                                 offset = 2*75,
                                 cutoff = 2*20, #70
                                 globe = 10,
                                 max.edge = c(3,20)*100, #70
                                 #loc.domain = st_as_sfc(speciesDataSF), 
                                 crs = proj)

# save(hares, mesh, speciesDataSF, studyArea,
#       file = 'data/speciesDataNEW.rdata')
predData <- FALSE
if (predData) {

  predData <- fm_pixels(mesh, mask = st_union(studyArea), dims = c(250, 250))
  
  predData$bio12 <- inlabru::eval_spatial(data = prec, predData)
  predData$bio12 <- inlabru::bru_fill_missing(prec, predData, values = predData$bio12)
  
  predData$bio1 <- inlabru::eval_spatial(data = temp, predData)
  predData$bio1 <- inlabru::bru_fill_missing(temp, predData, values = predData$bio1)
  
  predData <- fm_cprod(predData, data.frame(speciesNum = 1:(hares$Nnode + length(hares$tip.label))))
  
  predData$Bio1index <- predData$Bio1index2 <- predData$speciesNum
  predData$Bio12index <- predData$Bio12index2 <- predData$speciesNum
  predData$speciesNumIndex <- predData$speciesNum2 <- predData$speciesNumSpat <- predData$speciesNum
  predData[predData$speciesNumIndex > length(hares$tip.label), ]$speciesNumIndex <- NA
  
  saveRDS(predData, file = 'data/predData.rds')

} 

predANC <- FALSE
if (predANC) {
  
  library(ncdf4)
  Coords <- nc_open(filename ='~/Downloads/6b2. Mean Annual Precipitation (nc)/065Ma_precip.nc')
  Temp <- as.matrix(read.csv('~/Downloads/6a1. Mean Annual Temperatures (csv, list)/065_temp_list.csv'))
  Prec <- as.matrix(read.csv('~/Downloads/6b1. Mean Annual Precipitation (csv)/065_precip.csv'))
  lon <- ncvar_get(Coords,"lon")
  lat <- ncvar_get(Coords,"lat")
  
  #image.plot(x = lon, y  = lat, z = xtabs(elevation ~ longitude + latitude, Temp))
  #image.plot(x = lon, y  = lat, z = xtabs(Value ~  X.Longitude + Latitude, Prec))
  
  
  Coords <- nc_open(filename ='~/Downloads/5a. PaleoDEM (nc, 3601x1801)/000_v21144.nc')
  lon <- ncvar_get(Coords,"longitude")
  lat <- ncvar_get(Coords,"latitude")
  
  #Prec Mean 49.38131 SD 45.96677
  #Temp Mean 13.15843 SD 12.83619
  
  #image.plot(x = lon, y  = lat, z = xtabs(elevation ~ longitude + latitude, Temp))
  library(sf)
  qq <- st_as_sf(as.data.frame(Temp), coords = c('longitude', 'latitude'), crs = '+proj=lonlat')
  qq$Temp <- qq$temperature
  qq$bio1 <- qq$temperature
  qq$temperature <- NULL
  qq$bio12 <- Prec[,3] * 12
  qq$Prec <- qq$bio12
  qq$bio12 <- (qq$bio12 - 49.38131)/45.96677
  qq$bio1 <- (qq$bio1 - 13.15843)/12.83619
  qq$bio1s <- qq$bio1^2
  qq$bio12s <- qq$bio12^2
  qq <- st_transform(qq, crs = '+proj=moll +lon_0=40 +x_0=0 +y_0=0 +ellps=sphere +units=km +no_defs')
  qq$speciesNum <- 90 #91
  qq$Bio1index <- 90 #91
  qq$Bio1index2 <- 90 #91
  qq$Bio12index <- 90 #91
  qq$Bio12index2 <- 90 #91
  qq$speciesNumIndex <- 90 #91
  
  #saveRDS(qq, 'ancData.rds')
  
  library(inlabru)
  library(ggplot2)
  
  mapP <- mapast::getmap(65, model = 'PALEOMAP')# model = 'PALEOMAP'
  mapP <- st_as_sf(mapP)
  mapP <- st_transform(mapP, '+proj=moll +lon_0=40 +x_0=0 +y_0=0 +ellps=sphere +units=km +no_defs')
  qqID <- st_intersects(qq, mapP)
  qq2 <- qq[sapply(qqID, function(x) !identical(x, integer(0))),]
  saveRDS(qq2, 'ancData3.rds')
  
}