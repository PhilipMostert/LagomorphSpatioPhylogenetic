library(inlabru)
library(dplyr)
library(reshape2)
library(gridExtra)
library(ggplot2)

#Load Model
model <- readRDS('ModelSpatioPhylogenetic.rds')

#Load predData
predData <- readRDS('Data/ancData3.rds')

##Paleo map
StudyArea <- readRDS('Data/studyArea2.rds')

spatPreds <- predict(model, newdata = predData, formula = ~ INLA::inla.link.logit(Intercept + spEffect + RBio1 + RBio1s +
                                                            RBio12 + RBio12s + bio1 + bio12 + bio1s+
                                                            bio12s + RBio1i +
                                                            RBio1si + RBio12i + RBio12si, inverse = TRUE), 
                     seed = 1,   
                     probs = c(0.1,0.5,0.9), n.samples = 250) 

#saveRDS(spatPreds, 'spatPredsRoot.rds')


proj <- '+proj=longlat +datum=WGS84 +no_defs' 

pp = st_intersects(spatPreds, studyArea)
studyArea <- studyArea[unlist(unique(pp)),]

Rich <- spatPreds[spatPreds$speciesNum %in% 90,]
Rich2 <- Rich %>% group_by(geometry) %>%
  mutate(mid = (mean), low = (q0.1), high = (q0.9)) %>%
  select(low, mid, high) 
Rich3 <- st_as_sf(melt(Rich2, 'geometry'))

lowDist <- ggplot() + 
  gg(Rich2, aes(col = low)) + #predPreds
  theme_bw() +
  gg(st_boundary(studyArea)) +
  ggtitle('10% credibility interval \nof distribution') + #quantile
  scale_colour_gradientn(colours = rev(colorspace::heat_hcl(10))) + 
  coord_sf(crs = proj) +
  labs(col = 'Probability  ') +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = 'bold', colour = 'black'),
        axis.text.y = element_text(size=5, face="italic", colour = "black"),
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),
        axis.text.x = element_text(size=5, colour = "black"),
        strip.text = element_text(face="bold", size=5),
        legend.text = element_text(face="bold", size=5.5,
                                   angle = 30, hjust = 0.7, vjust = 0.5), #hjust= 0.5
        legend.title = element_text(face="bold", size=10), legend.position="bottom")

midDist <- ggplot() + 
  gg(Rich2, aes(col = mid)) + #predPreds
  theme_bw() +
  gg(st_boundary(studyArea)) +
  ggtitle('Mean estimate \nof distribution') +
  scale_colour_gradientn(colours = rev(colorspace::heat_hcl(10))) + 
  coord_sf(crs =  proj) +
  labs(col = 'Probability  ') +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = 'bold', colour = 'black'),
        axis.text.y = element_text(size=5, face="italic", colour = "black"),
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),
        axis.text.x = element_text(size=5, colour = "black"),
        strip.text = element_text(face="bold", size=5),
        legend.text = element_text(face="bold", size=5.5, #5
                                   angle = 30, hjust = 0.7, vjust = 0.5), #hjust= 0.5
        legend.title = element_text(face="bold", size=10),legend.position="bottom") #6

highDist <- ggplot() + 
  gg(Rich2, aes(col = high)) +
  theme_bw() +
  gg(st_boundary(studyArea)) +
  ggtitle('90% credibility interval \nof distribution') +
  scale_colour_gradientn(colours = rev(colorspace::heat_hcl(10))) + 
  coord_sf(crs = proj) +
  labs(col = 'Probability  ') +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = 'bold', colour = 'black'),
        axis.text.y = element_text(size=5, face="italic", colour = "black"),
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),
        axis.text.x = element_text(size=5, colour = "black"),
        strip.text = element_text(face="bold", size=5),
        legend.text = element_text(face="bold", size=5.5,
                                   angle = 30, hjust = 0.7, vjust = 0.5), #0.75 #hjust= 0.5
        legend.title = element_text(face="bold", size=10), legend.position="bottom")
cowplot::plot_grid(lowDist, midDist, highDist, nrow = 1, align = 'hv') #, labels = c('B)', '', '')
ggsave('ProbUrSpatMod.png', dpi = 300)

