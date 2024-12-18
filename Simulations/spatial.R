
##Set sigPar to 0.1, 1, 2
sigPar <- 2

set.seed(100)
library(sf)
library(ape)
library(INLA)
library(terra)
library(geiger)
library(inlabru)
library(MCMCglmm)
library(geostatsp)
library(R.utils)

numUsed <- c(10, 50, 100, 150, 200)
nSpec <- c(10, 50, 100, 150, 200)
namesLists <- paste0('sites:', numUsed)

interceptPhyloSpace <- interceptPhylo <- vector(mode = 'list', length = 5)
envPhyloSpace <- envPhylo <-  vector(mode = 'list', length = 5)

names(interceptPhyloSpace) <- names(interceptPhylo) <- namesLists
names(envPhyloSpace) <- names(envPhylo) <- namesLists

interceptAncestorsPhyloSpace <- interceptAncestorsPhylo <- vector(mode = 'list', length = 5)
envAncestorsPhyloSpace <- envAncestorsPhylo <- vector(mode = 'list', length = 5)

names(interceptAncestorsPhyloSpace) <- names(interceptAncestorsPhylo) <- namesLists
names(envAncestorsPhyloSpace) <- names(envAncestorsPhylo) <- namesLists

precIntPhyloSpace <- precEnvPhyloSpace <- vector(mode = 'list', length = 5)
precIntPhylo <- precEnvPhylo <- vector(mode = 'list', length = 5)

names(precIntPhyloSpace) <- names(precEnvPhyloSpace) <- namesLists
names(precIntPhylo) <- names(precEnvPhylo) <- namesLists

rangeSpatPhylo <- sigmaSpatPhylo <- vector(mode = 'list', length = 5)
names(rangeSpatPhylo) <- names(sigmaSpatPhylo) <- namesLists

dicPhyloSpace <- dicPhylo <- vector(mode = 'list', length = 5)
names(dicPhyloSpace) <- names(dicPhylo) <- namesLists

waicPhyloSpace <- waicPhylo <- vector(mode = 'list', length = 5)
names(waicPhyloSpace) <- names(waicPhylo) <- namesLists



for (Nsites in numUsed) {
  
  for (Nspecies in nSpec) {

for (i in 1:100){
  
  phyloSpaceMod <- try('a' +'a', silent = TRUE)
  phyloMod <- try('b' +'b', silent = TRUE)
  covMat <- try('e' +'e', silent = TRUE)
  
  
while (inherits(phyloSpaceMod, 'try-error') | length(phyloSpaceMod) < 5 | 
       inherits(phyloMod, 'try-error') | length(phyloMod) < 5 | 
       inherits(covMat, 'try-error')) {  

#sigPar <- 2
Nspecies <- Nspecies
Nsites <- Nsites 
delta <- 1
Tree <- geiger::sim.bdtree(n=Nspecies, stop = 'taxa') 
betaEnv <- 1     
betaPhy <- 1 
studyDim <- c(0, 10)
nPoints <- 10
phyloSignal <- 1
sigma.u <- 1
range <- 1
proj <- "epsg:3857"


##Create a continuous environmental map here
randomEffect <- rasterMap <- rast(nrows=50,
                 ncols=50,
                 xmin=studyDim[1],
                 xmax=studyDim[2],
                 ymin=studyDim[1],
                 ymax=studyDim[2],
                 crs = proj)

values(rasterMap) <- 1
y0 <- x0 <- seq(0, 10, length.out = sqrt(nrow(terra::crds(rasterMap))))
values(rasterMap) <- outer(y0,x0,
                           function (x,y) -1*x + 1*y)/5 

names(rasterMap) <- 'Environment'


PAGrid <- terra::as.data.frame(rasterMap, xy = TRUE,
                               cells = FALSE)

MeshGrid <- PAGrid <- st_as_sf(x = PAGrid,
                   coords = c('x', 'y'), crs = proj)

xx <- st_as_sf(st_sfc(st_polygon(list(matrix(c(0,0,10,0,10,10,0,10,0,0),ncol=2, byrow=TRUE)))))
st_crs(xx) <- proj

st_crs(PAGrid) <- proj
Mesh <- fmesher:::fm_mesh_2d(#loc = MeshGrid, 
  boundary = xx, loc = st_coordinates(PAGrid),
  max.edge = c(0.7,0.8) * 1,
  cutoff = 1, crs = proj)

SPDE = inla.spde2.pcmatern(Mesh, prior.range = c(.5, .5), prior.sigma = c(.5, .5))

Qu = inla.spde.precision(SPDE, 
                         theta=c(log(range), 
                                 log(sigma.u)))

PAGrid <- PAGrid[sample.int(n = nrow(PAGrid), size = Nsites, replace = FALSE),]


correctOrder <- c(1:Nspecies, rev((Nspecies+1):(Nspecies + Tree$Nnode)))

spEffect <- (phyloSignal^0.5 * rTraitCont(Tree, model = 'BM', sigma = sigPar, ancestor = TRUE) + (1-phyloSignal)^0.5 * rnorm(n = Nspecies + Tree$Nnode))[correctOrder]

Names <-  scale(rTraitCont(Tree, ancestor = TRUE)[correctOrder]) 

spENV <- (phyloSignal^0.5 * rTraitCont(Tree, model = 'BM', sigma = sigPar, ancestor = TRUE) + (1-phyloSignal)^0.5 * rnorm(n = Nspecies + Tree$Nnode))[correctOrder]

spatEf <- rep(1, length(spEffect))

speciesObs <- vector(mode = 'list', length = length(Tree$tip.label) + Tree$Nnode)
names(speciesObs) <- row.names(Names)

u = inla.qsample(n=length(Tree$tip.label) + Tree$Nnode, Q=Qu)
#u = u[ ,1]

A = inla.spde.make.A(mesh=Mesh, loc=PAGrid)
u = drop(A %*% u)

for(j in 1:(length(Tree$tip.label) + Tree$Nnode)) {

  probs <- INLA::inla.link.logit(rnorm(n = nrow(PAGrid)) + u[,j] + betaPhy*spEffect[j] + betaEnv*spENV[j]*(PAGrid$Environment), inverse = TRUE) #betaEnv*PAGrid$Environment + 
  PAGrid$Probs <- probs
  occs <- rbinom(Nsites, 1, prob = probs)
  
  if (j <= length(Tree$tip.label)) {
  
    PAGrid$Status <- occs
    PAGrid$Truth <- NA
    PAGrid$spIndex2 <- j
  
  }
  else {
    
    PAGrid$Status <- NA
    PAGrid$Truth <- occs
    PAGrid$spIndex2 <- NA
    
  }
  
  PAGrid$Name <- row.names(Names)[j]
  PAGrid$spIndex <- j
  PAGrid$spIndexENV <- j
  
  speciesObs[[j]] <- PAGrid
  
}
speciesObs[[length(speciesObs)]] <- NULL
dataSet <- do.call(rbind, speciesObs)

##Create phylo covariance matrix

covMat <- try(inverseA(Tree, nodes = 'ALL', scale = FALSE)$Ainv, silent = TRUE)
covMat <- try(covMat[na.omit(match(unique(dataSet$Name), 
                     row.names(covMat))),na.omit(match(unique(dataSet$Name), row.names(covMat)))], silent = TRUE)
#Re order this

##Run INLA model
priorSpat <- priorSpat <- list(theta = list(prior = "loggamma", param = c(1, 1)))
priorIID <- list(prec = list(initial = 0, fixed = FALSE, prior = 'pc.prec', param = c(1, 0.1))) 
prior <- list(prior="loggamma", param = c(1, 5e-5))


modOptions <- list(safe = TRUE, num.threads = 2,verbose = FALSE,
                                control.inla = list(int.strategy = 'ccd',
                                                    control.vb= list(enable = FALSE),
                                                    strategy="gaussian"),
                             control.predictor = list(compute = FALSE, link = 1), 
                             control.fixed = list(correlation.matrix=FALSE),
                             control.compute = list(config = TRUE,
                                                    return.marginals.predictor=TRUE))


bruLikeSpace <- like(formula = Status ~ .,
                data = dataSet,
                family = 'binomial',
                include = c('Intercept', 'spIntercept', 'spatEffect','spENV', 'Environment'))

bruLike <- like(formula = Status ~ .,
                data = dataSet,
                family = 'binomial',
                include = c('Intercept', 'spIntercept','spENV', 'Environment'))

phyloSpaceComps <- ~ -1 + Intercept(1) +
               Environment +
               spIntercept(main = spIndex, model = 'generic0', Cmatrix = covMat, constr = TRUE, hyper = list(prec = prior)) +
               spENV(main = spIndexENV, Environment, model = 'generic0',Cmatrix = covMat, constr = TRUE, hyper = list(prec = prior)) +
               spatEffect(main = geometry, model = SPDE, 
               group = spIndex2,
               control.group = list(model = 'iid', hyper = priorIID))

phyloSpaceMod <- try(withTimeout(bru(phyloSpaceComps, bruLikeSpace, options = modOptions), timeout = 10000), silent = TRUE)

phyloComps <- ~ -1 + Intercept(1) + 
  Environment +
  spIntercept(main = spIndex, model = 'generic0', Cmatrix = covMat, constr = TRUE, hyper = list(prec = prior)) +
  spENV(main = spIndexENV, Environment, model = 'generic0',Cmatrix = covMat, constr = TRUE, hyper = list(prec = prior))

phyloMod <- try(withTimeout(bru(phyloComps, bruLike, options = modOptions), timeout = 10000), silent = TRUE)

print(class(phyloSpaceMod))
print(class(phyloMod))

}
#1/bruMod$summary.hyperpar[3,]

interceptPhyloSpace[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- (sum(((phyloSpaceMod$summary.random$spIntercept[1:Nspecies,2] + phyloSpaceMod$summary.fixed[1,1] - spEffect[1:Nspecies])^2))/Nspecies)
interceptPhylo[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- (sum(((phyloMod$summary.random$spIntercept[1:Nspecies,3] + phyloMod$summary.fixed[1,1] - spEffect[1:Nspecies])^2))/Nspecies)

envPhyloSpace[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- (sum(((phyloSpaceMod$summary.random$spENV[1:Nspecies,2] + phyloSpaceMod$summary.fixed[2,1] - spENV[1:Nspecies])^2))/Nspecies)
envPhylo[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- (sum(((phyloMod$summary.random$spENV[1:Nspecies,2] + phyloMod$summary.fixed[2,1] - spENV[1:Nspecies])^2))/Nspecies)

anID <- (Nspecies+1):nrow(phyloSpaceMod$summary.random$spIntercept)

interceptAncestorsPhyloSpace[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- (sum(((phyloSpaceMod$summary.random$spIntercept[anID,2] + phyloSpaceMod$summary.fixed[1,1] - spEffect[anID])^2))/Nspecies) 
interceptAncestorsPhylo[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- (sum(((phyloMod$summary.random$spIntercept[anID,2] + phyloMod$summary.fixed[1,1] - spEffect[anID])^2))/Nspecies) 

envAncestorsPhyloSpace[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- (sum(((phyloSpaceMod$summary.random$spENV[anID,5] + phyloSpaceMod$summary.fixed[2,1] - spENV[anID])^2))/Nspecies) 
envAncestorsPhylo[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- (sum(((phyloMod$summary.random$spENV[anID,5] + phyloMod$summary.fixed[2,1] - spENV[anID])^2))/Nspecies) 

#Check these
precIntPhyloSpace[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- phyloSpaceMod$summary.hyperpar[1,1]
precEnvPhyloSpace[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- phyloSpaceMod$summary.hyperpar[2,1]

precIntPhylo[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- phyloMod$summary.hyperpar[1,1]
precEnvPhylo[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- phyloMod$summary.hyperpar[2,1]

rangeSpatPhylo[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- phyloSpaceMod$summary.hyperpar[3,1]
sigmaSpatPhylo[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- phyloSpaceMod$summary.hyperpar[4,1]

dicPhyloSpace[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <-  phyloSpaceMod$dic$dic
dicPhylo[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- phyloMod$dic$dic

waicPhyloSpace[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <- phyloSpaceMod$waic$waic
waicPhylo[[paste0('sites:',Nsites)]][[paste0('species:', Nspecies)]][i] <-  phyloMod$waic$waic
#print(paste('interceptPhyloSpace',interceptPhyloSpace[i]))
#print(paste('interceptPhylo',interceptPhylo[i]))

print(i)
print(Nsites)
print(Nspecies)
}
    
  }


}

resultsInt <- data.frame(phyloSpace = interceptPhyloSpace, phylo = interceptPhylo, type = 'Intercept')
resultsEnv <- data.frame(phyloSpace = envPhyloSpace, phylo = envPhylo, type = 'Environmental')
resAll <- rbind(resultsInt, resultsEnv)

library(ggplot2)
library(reshape2)

plotData <- melt(resAll, 'type')
names(plotData) <- c('Variable', 'Model', 'MSE')
RMSEplot <- ggplot(plotData, aes(x = Variable, y = sqrt(MSE), fill=Model)) + geom_boxplot()

varRes <- data.frame(Intercept = 1/precIntPhyloSpace, Environmental = 1/precEnvPhyloSpace, Range = rangeSpatPhylo, Sigma = sigmaSpatPhylo)
plotVar <- melt(varRes)
hyperPlot <- ggplot(plotVar, aes(x = variable, y = value)) + geom_boxplot() 

spaceAnc <- data.frame(env = envAncestorsPhyloSpace, int = interceptAncestorsPhyloSpace)
nonAnc <- data.frame(env = envAncestorsPhylo, int = interceptAncestorsPhylo)

print(delta)
print(sigPar)


save(interceptPhyloSpace, interceptPhylo, envPhyloSpace, envPhylo, 
     interceptAncestorsPhyloSpace, interceptAncestorsPhylo, 
     envAncestorsPhyloSpace, envAncestorsPhylo, precIntPhyloSpace, 
     precEnvPhyloSpace, precIntPhylo, precEnvPhylo, rangeSpatPhylo, 
     sigmaSpatPhylo, dicPhyloSpace, dicPhylo, waicPhyloSpace, waicPhylo, file = 'Simulations/spatialResults/spatSim.Rdata')

