library(ggplot2)
library(reshape2)


###########Plot spatial simulations
load('spatialResults/Sigma1/spatSim.Rdata')
Len1 <- length(interceptPhyloSpace)
Len2 <- length(interceptPhyloSpace[[1]])
Len3 <- length(interceptPhyloSpace[[1]][[1]])
sites <- rep(c(10, 50, 100, 150, 200), each = Len1 * Len3)
species <- rep(rep(c(10, 50, 100, 150, 200), each = Len3), times = Len1)

intPS <- unlist(c(dplyr::bind_rows(interceptPhyloSpace)))#unlist(interceptPhyloSpace)
intP <- unlist(c(dplyr::bind_rows(interceptPhylo))) #unlist(interceptPhylo)
intDat <- data.frame(SpatialPhylo = intPS, Phylo = intP, Sites = sites, Species = species)
intDat$dif <- intDat$SpatialPhylo - intDat$Phylo
intDat <- melt(intDat,id.vars = c('Sites', 'Species','dif'))
intDat$Var = 'Intercept'
intDat$Species <- as.character(intDat$Species)


envPS <-  unlist(c(dplyr::bind_rows(envPhyloSpace)))#unlist(envPhyloSpace)
envP <- unlist(c(dplyr::bind_rows(envPhylo)))##nlist(envPhylo)
envDat <- data.frame(SpatialPhylo = envPS, Phylo = envP, Sites = sites, Species = species)
envDat$dif <- envDat$SpatialPhylo - envDat$Phylo
envDat <- melt(envDat,id.vars = c('Sites', 'Species', 'dif'))
envDat$Var <- 'Environmental'
envDat$Species <- as.character(envDat$Species)


ggplot() + 
  geom_boxplot(data = intDat, aes(y = (value), x = variable, fill = reorder(Species, as.numeric(Species))), 
               outlier.size = 0.5, outlier.alpha = 0.5) +
  guides(fill=guide_legend(title="Number \nof tips")) +
  scale_fill_manual(values = c("#B8B8D1", "#06D6A0", '#5B5F97','#EF476F', '#EAD637')) +
  facet_grid(~ Sites, scales = 'free',
             labeller = as_labeller(c(`10`='10 sites',
                                      `50`='50 sites',
                                      `100` ='100 sites',
                                      `150`='150 sites',
                                      `200` = '250 sites'))) +
  scale_x_discrete(labels = c('Spatio-\nphylo', 'Phylo')) +
  ylab(label = expression(sqrt(Mean~squared~error))) +
  xlab('Model') +
  ggtitle('RMSE values for the random intercept') +
  theme_bw() + theme(text = element_text(size = 25),
                     strip.background =element_rect(fill='white'),
                     axis.text.x = element_text(size = 15),
                     plot.title = element_text(hjust = 0.5, face = 'bold'),
                     legend.text=element_text(size=16),
                     legend.title = element_text(size = 20))
ggsave('RMSE.png',dpi = 300)





#####Plot phylogenetic simulations 

load('Simulations/phylogeneticResults/results10.rdata')
signal10 <- data.frame(Intercept = unlist(sInt),
                       Environmental = unlist(sEnv),
                       Lambda = as.character(rep(seq(0,1, by = 0.2), each = 500)),
                       Nspecies = as.character(10))

load('Simulations/phylogeneticResults/results50.rdata')
signal50 <- data.frame(Intercept = unlist(sInt),
                       Environmental = unlist(sEnv),
                       Lambda = as.character(rep(seq(0,1, by = 0.2), each = 500)),
                       Nspecies = as.character(50))


load('Simulations/phylogeneticResults/results100.rdata')
signal100 <- data.frame(Intercept = unlist(sInt),
                        Environmental = unlist(sEnv),
                        Lambda = as.character(rep(seq(0,1, by = 0.2), each = 500)),
                        Nspecies = as.character(100))

load('Simulations/phylogeneticResults/results150.rdata')
signal150 <- data.frame(Intercept = unlist(sInt),
                        Environmental = unlist(sEnv),
                        Lambda = as.character(rep(seq(0,1, by = 0.2), each = 500)),
                        Nspecies = as.character(150))

load('Simulations/phylogeneticResults/results200.rdata')
signal200 <- data.frame(Intercept = unlist(sInt),
                        Environmental = unlist(sEnv),
                        Lambda = as.character(rep(seq(0,1, by = 0.2), each = 500)),
                        Nspecies = as.character(200))

signalAll <- rbind(signal10, signal50, signal100, signal150, signal200)
signalAll[signalAll$Lambda == '0',]$Lambda <- '0.0'
signalAll[signalAll$Lambda == '1',]$Lambda <- '1.0'


ggplot(signalAll) + 
  geom_boxplot(aes(x = Lambda, y = Environmental, fill = reorder(Nspecies, as.numeric(Nspecies))), outlier.size = 0.5, outlier.alpha = 0.5) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0,1)) +
  ggtitle('Phylogenetic signal estimates') +
  #geom_hline(yintercept=0, linetype="dashed", color = "black", alpha = 0.5) +
  #geom_hline(yintercept=0.2, linetype="dashed", color = "black", alpha = 0.5) +
  #geom_hline(yintercept=0.4, linetype="dashed", color = "black", alpha = 0.5) +
  #geom_hline(yintercept=0.6, linetype="dashed", color = "black", alpha = 0.5) +
  #geom_hline(yintercept=0.8, linetype="dashed", color = "black", alpha = 0.5) +
  #geom_hline(yintercept=1, linetype="dashed", color = "black", alpha = 0.5) +
  geom_segment(aes(x = 0.5, y = 0, xend = 1.5, yend = 0), col = '#FF6A00', linetype = 'dashed') +
  geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2), col = '#FF6A00', linetype = 'dashed') +
  geom_segment(aes(x = 2.5, y = 0.4, xend = 3.5, yend = 0.4), col = '#FF6A00', linetype = 'dashed') +
  geom_segment(aes(x = 3.5, y = 0.6, xend = 4.5, yend = 0.6), col = '#FF6A00', linetype = 'dashed') +
  geom_segment(aes(x = 4.5, y = 0.8, xend = 5.5, yend = 0.8), col = '#FF6A00', linetype = 'dashed') +
  geom_segment(aes(x = 5.5, y = 1, xend = 6.5, yend = 1), col = '#FF6A00', linetype = 'dashed') +
  xlab(expression(paste('True value of ', lambda))) +
  ylab(expression(paste('Predicted value of ', lambda))) +
  guides(fill=guide_legend(title="Number \nof tips")) +
  scale_fill_manual(values = c("#B8B8D1", "#06D6A0", '#5B5F97','#EF476F', '#EAD637')) +
  theme_bw() + theme(text = element_text(size = 25),
                     plot.title = element_text(hjust = 0.5, face = 'bold'),
                     legend.text=element_text(size=16),
                     legend.title = element_text(size = 20))
ggsave('phySignalENV.png', dpi = 300)





