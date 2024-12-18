library(rworldmap)
library(ggplot2)
library(inlabru)
library(sf)
library(ape)
library(terra)
library(stringr)

#proj <- INLA::inla.CRS('globe') #"EPSG:4326"
##CHANGE TO KM
#proj <- "+proj=aea +lat_1=52.5 +lat_2=-10 +lon_0=70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"

#proj <- inlabru::fm_crs('mollweide_globe')
#This proj is new
proj <- '+proj=moll +lon_0=40 +x_0=0 +y_0=0 +ellps=sphere +units=km +no_defs'

mam <- read.nexus("data-raw/Bunnies/Phylo.txt")
grep("Lago", mam[[1]]$node.label)
hares <- extract.clade(mam[[1]], mam[[1]]$node.label[c(673)])

Lepus <- sf::st_read('~/Downloads/redlist_species_data_8c436bf9-8db4-4cc8-8657-3662399ac09b/data_0.shp')
Lepus$Genus <- word(Lepus$SCI_NAME)
Ochotona <- sf::st_read('~/Downloads/redlist_species_data_cdd72903-fff8-4ec9-952e-23935fd863ba/')
Ochotona$Genus <- word(Ochotona$SCI_NAME)
Lago <- rbind(Lepus,Ochotona)


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
hares <- drop.tip(mam[[1]], mam[[1]]$tip.label[!mam[[1]]$tip.label%in% c(hares$tip.label)])
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

specDt <- specData[specData$species %in% hares$tip.label,]

##Make the order better such that Lepus and Ochotana are plotted first, and then the rest
ggplot() + 
  geom_sf(data = st_boundary(studyArea)) +
  geom_sf(data = specDt[specDt$Genus == 'Lepus',], aes(fill = Genus), alpha = 0.75) +
  geom_sf(data = specDt[specDt$Genus == 'Ochotona',], aes(fill = Genus), alpha = 0.75) +
  geom_sf(data = specDt[specDt$Genus == 'Pronolagus',], aes(fill = Genus), alpha = 0.75) +
  geom_sf(data = specDt[specDt$Genus == 'Oryctolagus',], aes(fill = Genus), alpha = 0.75) +
  geom_sf(data = specDt[specDt$Genus == 'Caprolagus',], aes(fill = Genus), alpha = 0.75) +
  geom_sf(data = specDt[specDt$Genus == 'Poelagus',], aes(fill = Genus), alpha = 0.75) +
  geom_sf(data = specDt[specDt$Genus == 'Bunolagus',], aes(fill = Genus), alpha = 0.75) +
  scale_fill_manual(values = c('Bunolagus' = '#C2B9CB',
                               'Caprolagus' = '#F9C846',
                               'Lepus' = '#F96E46',
                               'Ochotona' = '#00E8FC', 
                               'Oryctolagus' = '#E086D3',
                               'Poelagus' = '#4A8FE7',
                               'Pronolagus' = '#0FFF95')) +
  theme_void() +
  theme(axis.text=element_text(size=40), 
        legend.text = element_text(size=40),
        legend.title = element_text(size = 60),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

ggsave('Map_of_word.png', scale = 3)
library(ggtree)
hares
haresDat <- tidytree::as_tibble(hares)
haresDat$Genus <- word(gsub('_', ' ', haresDat$label))  

hares$tip.label <- gsub('_', ' ', hares$tip.label)
hares$edge.length <- hares$edge.length/10
ggtree(hares, layout = 'fan', size = 2) %<+% haresDat + geom_tiplab2(aes(fill = Genus),
                                            color = "black",
                                            geom = "label", hjust = -0.01,
                                            align=TRUE, linesize=.5,
                                            label.padding = unit(0.015, "lines"), 
                                            label.size =0, size = 12)  +
  #geom_tiplab(size=2, align=TRUE, linesize=.5) +
  #scale_x_ggtree() +
  xlim(NA, 12) +
  #ylim(NA, 0) +
  scale_fill_manual(values = c('Bunolagus' = '#C2B9CB',
                               'Caprolagus' = '#F9C846',
                               'Lepus' = '#F96E46',
                               'Ochotona' = '#00E8FC', 
                               'Oryctolagus' = '#E086D3',
                               'Poelagus' = '#4A8FE7',
                               'Pronolagus' = '#0FFF95')) +
  guides(fill = guide_legend(override.aes = list(label = "   ", size = 2))) +
  theme(legend.position = c(-100,-100),  #c(0.2,0.75)
        #legend.title = element_blank(), # no title
        legend.key = element_blank()) 
ggsave('Tree.png', scale = 3)
