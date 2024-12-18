
#Set sigPar
sigPar <- 1
#SetSites
Nsites <- 100


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
namesLists <- paste0('num:', numUsed)

precInt <- vector(mode = 'list', length = 5)
precEnv <- vector(mode = 'list', length = 5)
precInti <- vector(mode = 'list', length = 5)
precEnvi <- vector(mode = 'list', length = 5)
sigInt <- vector(mode = 'list', length = 5)
sigEnv <- vector(mode = 'list', length = 5)

rangeParam <- vector(mode = 'list', length = 5)
sigParam <- vector(mode = 'list', length = 5)

names(precInt) <- namesLists
names(precEnv) <- namesLists
names(precInti) <- namesLists
names(precEnvi) <- namesLists
names(sigInt) <- namesLists
names(sigEnv) <- namesLists
names(rangeParam) <- namesLists
names(sigParam) <- namesLists


for (Nspecies in numUsed) {
  
for (phyloSignal in seq(0, 1, by = 0.2)) { 

for (i in 1:500) { 
  
  phyloSpaceMod <- try('a' +'a', silent = TRUE)
  covMat <- try('b' +'b', silent = TRUE)
  
  precInt[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  precEnv[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  precInti[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  precEnvi[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  
  rangeParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  sigParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  
  
while (inherits(phyloSpaceMod, 'try-error') |
       inherits(covMat, 'try-error') |
       is.na(precInt[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]) |
       is.na(precEnv[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]) |
       is.na(precInti[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]) |
       is.na(precEnvi[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]) |
       is.na(rangeParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]) |
       is.na(sigParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i])) {  

sigPar <- 1
Nspecies <- Nspecies#100 #100
#Nsites <- 100
delta <- 1
Tree <- geiger::sim.bdtree(n=Nspecies, stop = 'taxa') 
Tree <- rescale(Tree, 'delta', delta)
betaEnv <- 1  
betaPhy <- 1  
studyDim <- c(0, 10)
nPoints <- 10
rfSD <- sqrt(1)
phyloSignal <- phyloSignal#0.25
proj <- "epsg:3857"

##Create a continuous environmental map here
randomEffect <- rasterMap <- rast(nrows=50,
                 ncols=50,
                 xmin=studyDim[1],
                 xmax=studyDim[2],
                 ymin=studyDim[1],
                 ymax=studyDim[2],
                 crs = proj)

rfParams = c(mean = 0, variance= rfSD^2, #1 #3 + cellSize
             range=1, shape=2)

values(rasterMap) <- 1
y0 <- x0 <- seq(0, 10, length.out = sqrt(nrow(terra::crds(rasterMap))))
values(rasterMap) <- outer(y0,x0,
                           function (x,y) -1*x + 1*y)/5

names(rasterMap) <- 'Environment'


siteLocs <- data.frame(x = runif(nPoints, 0, 10),
                        y = runif(nPoints, 0, 10))

siteLocs <- st_as_sf(x = siteLocs,
                      coords = c('x', 'y'), crs = proj)

PAGrid <- terra::as.data.frame(rasterMap, xy = TRUE,
                               cells = FALSE)


MeshGrid <- PAGrid <- st_as_sf(x = PAGrid,
                   coords = c('x', 'y'), crs = proj)

xx <- st_as_sf(st_sfc(st_polygon(list(matrix(c(0,0,10,0,10,10,0,10,0,0),ncol=2, byrow=TRUE)))))
st_crs(xx) <- proj
st_crs(PAGrid) <- proj
Mesh <- fmesher:::fm_mesh_2d(#loc = MeshGrid, 
  boundary = xx, loc = st_coordinates(PAGrid),
  max.edge = c(0.7,0.8),
  cutoff = 1, crs = proj)

sigma.u = 1
range = 1
SPDE = inla.spde2.pcmatern(Mesh, prior.range = c(.5, .5), prior.sigma = c(.5, .5))

Qu = inla.spde.precision(SPDE, 
                         theta=c(log(range), 
                                 log(sigma.u)))


PAGrid <- PAGrid[sample.int(n = nrow(PAGrid), size = Nsites, replace = FALSE),]

#Create trait using the geiger package

correctOrder <- c(1:Nspecies, rev((Nspecies+1):(Nspecies + Tree$Nnode)))

spEffect <- (phyloSignal^0.5 * rTraitCont(Tree, model = 'BM', sigma = sigPar, ancestor = TRUE) + (1-phyloSignal)^0.5 * rnorm(n = Nspecies + Tree$Nnode))[correctOrder]

Names <-  scale(rTraitCont(rescale(Tree, "delta", 10), ancestor = TRUE)[correctOrder])

spENV <- (phyloSignal^0.5 * rTraitCont(Tree, model = 'BM', sigma = sigPar, ancestor = TRUE) + (1-phyloSignal)^0.5 * rnorm(n = Nspecies + Tree$Nnode))[correctOrder]

#spatEf <- rnorm(n = length(spENV), mean = 0, sd = sqrt(1/exp(-0)))
spatEf <- rep(1, length(spEffect))

speciesObs <- vector(mode = 'list', length = length(Tree$tip.label) + Tree$Nnode)
names(speciesObs) <- row.names(Names)

u = inla.qsample(n=length(Tree$tip.label) + Tree$Nnode, Q=Qu)
#u = u[ ,1]

A = inla.spde.make.A(mesh=Mesh, loc=PAGrid)
u = drop(A %*% u)

for(j in 1:(length(Tree$tip.label) + Tree$Nnode)) {
  #spatEf[j]*PAGrid$Space
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
#list(prec = list(initial = 0)) 
prior <- list(prior="loggamma", param = c(1, 5e-5))



modOptions <- list(safe = TRUE, num.threads = 6,  verbose = FALSE,
                                control.inla = list(int.strategy = 'ccd',
                                                    stupid.search = FALSE,
                                                    control.vb= list(enable = FALSE),
                                                    strategy="gaussian"),
                             control.predictor = list(compute = FALSE, link = 1), 
                             control.fixed = list(correlation.matrix=FALSE),
                             control.compute = list(config = TRUE,
                                                    return.marginals.predictor=TRUE))


bruLikeSpace <- like(formula = Status ~ .,
                data = dataSet,
                family = 'binomial',
                include = c('Intercept', 
                            'spIntercept',
                            'spatEffect',
                            'spENV', 
                            'Environment',
                            'spInti',
                            'spEnvi'))

phyloSpaceComps <- ~ -1 + Intercept(1) +
               Environment +
               spIntercept(main = spIndex, model = 'generic0', Cmatrix = covMat, constr = TRUE, hyper = list(prec = prior)) +
               spENV(main = spIndexENV, Environment, model = 'generic0',Cmatrix = covMat, constr = TRUE, hyper = list(prec = prior)) +
               spInti(main = spIndex, model = 'iid', constr = TRUE, hyper = list(prec = prior)) + #FALSE
               spEnvi(main = spIndexENV, Environment, model = 'iid', constr = TRUE, hyper = list(prec = prior)) + #FALSE
               spatEffect(main = geometry, model = SPDE, 
               group = spIndex2,
               control.group = list(model = 'iid', hyper = priorIID))


phyloSpaceMod <- try(withTimeout(bru(phyloSpaceComps, bruLikeSpace, options = modOptions), timeout = 100000), silent = TRUE)


if (!inherits(phyloSpaceMod, 'try-error')) {

  #Turn these into median value
precInt[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- try(inla.zmarginal(inla.tmarginal(function(x) 1/x, marginal = phyloSpaceMod$marginals.hyperpar$`Precision for spIntercept`), silent = TRUE)[[5]])
precEnv[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- try(inla.zmarginal(inla.tmarginal(function(x) 1/x, marginal = phyloSpaceMod$marginals.hyperpar$`Precision for spENV`), silent = TRUE)[[5]])

precInti[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- try(inla.zmarginal(inla.tmarginal(function(x) 1/x, marginal = phyloSpaceMod$marginals.hyperpar$`Precision for spInti`), silent = TRUE)[[5]])
precEnvi[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- try(inla.zmarginal(inla.tmarginal(function(x) 1/x, marginal = phyloSpaceMod$marginals.hyperpar$`Precision for spEnvi`), silent = TRUE)[[5]])

rangeParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- try(inla.zmarginal(phyloSpaceMod$marginals.hyperpar$`Range for spatEffect`, silent = TRUE)[[5]])
sigParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- try(inla.zmarginal(phyloSpaceMod$marginals.hyperpar$`Stdev for spatEffect`, silent = TRUE)[[5]])

if (!inherits(precInt[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i], 'try-error') & !inherits(as.numeric(precInt[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]), 'character') &
    !inherits(precEnv[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i], 'try-error') & !inherits(as.numeric(precEnv[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]), 'character') &
    !inherits(precInti[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i], 'try-error') & !inherits(as.numeric(precInti[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]), 'character') &
    !inherits(precEnvi[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i], 'try-error') & !inherits(as.numeric(precEnvi[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]), 'character') &
    !inherits(rangeParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i], 'try-error') & !inherits(as.numeric(rangeParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]), 'character') &
    !inherits(sigParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i], 'try-error') & !inherits(as.numeric(sigParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]), 'character')) {

sigInt[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- as.numeric(precInt[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i])/(as.numeric(precInt[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]) + as.numeric(precInti[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]))
sigEnv[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- as.numeric(precEnv[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i])/(as.numeric(precEnv[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i])  + as.numeric(precEnvi[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i]))

rangeParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- as.numeric(inla.zmarginal(phyloSpaceMod$marginals.hyperpar$`Range for spatEffect`, silent = TRUE)[[5]])
sigParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- as.numeric(inla.zmarginal(phyloSpaceMod$marginals.hyperpar$`Stdev for spatEffect`, silent = TRUE)[[5]])

} else {
  
  precInt[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  precEnv[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  precInti[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  precEnvi[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  rangeParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  sigParam[[paste0('num:',Nspecies)]][[paste0('sig:', phyloSignal)]][i] <- NA
  
  
}



print(i)
print(Nspecies)
print(phyloSignal)

}

}
  
}

}
  
  pInt <- precInt[[paste0('num:',Nspecies)]]
  pEnv <- precEnv[[paste0('num:',Nspecies)]]
  pInti <- precInti[[paste0('num:',Nspecies)]]
  pEnvi <- precEnvi[[paste0('num:',Nspecies)]]
  sInt <- sigInt[[paste0('num:',Nspecies)]] 
  sEnv <- sigEnv[[paste0('num:',Nspecies)]]
  rField <- rangeParam[[paste0('num:',Nspecies)]]
  sigField <- sigParam[[paste0('num:',Nspecies)]]

  save(pInt, pEnv, pInti, pEnvi, sInt, sEnv, rField, sigField,
       file = paste0('Simulations/phylogeneticResults/results', Nspecies, '.rdata'))

}
