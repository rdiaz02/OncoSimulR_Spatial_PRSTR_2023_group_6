# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("OncoSimulR")

# install.packages("scatterplot3d")
# install.packages("viridis")
library(viridis)
library(OncoSimulR)
library(parallel)
library(dplyr)
library(scatterplot3d)

## Definimos genotipos: 
fe <- allFitnessEffects(
  data.frame(parent = c("Root", "Root", "i"),
             child = c("u" , "i" , "v"),
             s = c(0.1 , -0.05 , 0.25),
             sh = -1,
             typeDep = "MN"),
  epistasis = c("u:i" = -1,"u:v" = -1))

evalAllGenotypes (fe , order = FALSE, addwt = TRUE)

## Function for creating the output of SpatialOncoSimul and the SpatialOncoSimul 
## object type.

to_SpatialOncoSimul <- function(x, spatialIterMax, ...){
  UseMethod("to_SpatialOncoSimul")
}

## Take a list and return (if possible) a SpatialOncoSimul object.

to_SpatialOncoSimul.list <- function (finalgrid_list, spatialIterMax, ...){
  JoinGrids <- bind_rows(finalgrid_list)
  TotalPop <- sum(JoinGrids$N)
  TotalDems <- length(finalgrid_list)
  TotalIte <- spatialIterMax - 1
  PopGenotypes <- aggregate(N ~ Genotype, data = JoinGrids, sum)
  ListGenotypes <- PopGenotypes$Genotype
  PopPerGenotypes <- PopGenotypes$N
  DemsGenotype <- c()
  for (Genotype in ListGenotypes){
    DemsPerGenotype <- nrow(JoinGrids[which(JoinGrids$Genotype == Genotype
                                            & JoinGrids$N != 0), ])
    DemsGenotype <- c(DemsGenotype, DemsPerGenotype)
  }
  print(paste0("Final Population: ", TotalPop, " (", TotalDems, " demes ",
        TotalIte, " iterations)"))
  for (g in 1:length(ListGenotypes)){
    if (ListGenotypes[g] != "" && PopPerGenotypes[g] != 0) {
      print(paste0("Genotype: ", ListGenotypes[g], " --- Population: ", 
                   PopPerGenotypes[g], " --- Demes: ", DemsGenotype[g]))
    }
  }
  finalgridspatial <- list(Final_Population = TotalPop, Total_Demes = TotalDems,
                           Total_Iterations = TotalIte, 
                           List_of_Genotypes = PopGenotypes,
                           Grids = finalgrid_list)
  finalgrid_spatial <- structure(finalgridspatial, class = "SpatialOncoSimul")
  return(finalgrid_spatial)
}


## Final function for Spatial Simulation (new arguments have been included
## with respect to OncoSimulIndiv).

SpatialOncoSimul <- function(fp, model = "Exp", 
                             numPassengers = 0, mu = 1e-6, muEF = NULL,
                             detectionSize = 1e8, detectionDrivers = 4,
                             detectionProb = NA,
                             sampleEvery = ifelse(model %in% c("Bozic", "Exp"), 
                                                  1, 0.025),
                             initSize = 500, s = 0.1, sh = -1,
                             K = sum(initSize)/(exp(1) - 1), 
                             keepEvery = sampleEvery,
                             minDetectDrvCloneSz = "auto",
                             extraTime = 0,
                             finalTime = 0.25 * 25 * 365, onlyCancer = FALSE,
                             keepPhylog = FALSE,
                             mutationPropGrowth = ifelse(model == "Bozic",
                                                         FALSE, TRUE),
                             max.memory = 2000, max.wall.time = 200,
                             max.num.tries = 500,
                             errorHitWallTime = TRUE,
                             errorHitMaxTries = TRUE,
                             verbosity = 0,
                             initMutant = NULL,
                             AND_DrvProbExit = FALSE,
                             fixation = NULL,
                             seed = NULL,
                             interventions = NULL,
                             userVars = NULL,
                             rules = NULL,
                             migrationProb = 0.5, 
                             largeDistMigrationProb = 1e-6, 
                             maxMigrationPercentage = 0.2,
                             SpatialModel = "3D", 
                             spatialIterMax = 20,
                             mc.cores = detectCores() - 1,
                             SpatialPlot = TRUE) {
  
    if(.Platform$OS.type == "windows") {
      if(mc.cores != 1)
        message("You are running Windows. Setting mc.cores = 1")
      mc.cores <- 1
    }
  
    # 1. Function for tumor porgression simulation in a grid.
    oncoSimulIndiv_grid <- function(grid){
        # If all the cells from a grid disappear N will be 0, so this grid
        # will be deleted from the list of current grids. The growth simulation 
        # inside a grid will only be performed in grids with N > (exp(1) - 1) to 
        # calculate K initial population equilibrium size. 
      
        if (sum(grid$N) > (exp(1) - 1)){
            grid <- grid[which(grid$N > 0),]
            simul_grid <- invisible(oncoSimulIndiv(fp = fp, model = model, 
                                         numPassengers = numPassengers, mu = mu, 
                                         muEF = muEF,
                                         detectionSize = detectionSize, 
                                         detectionDrivers = detectionDrivers,
                                         detectionProb = detectionProb,
                                         sampleEvery = sampleEvery,
                                         initSize = grid$N, s = s, sh = sh,
                                         K = K, keepEvery = keepEvery,
                                         minDetectDrvCloneSz = minDetectDrvCloneSz,
                                         extraTime = extraTime,
                                         finalTime = finalTime, 
                                         onlyCancer = onlyCancer,
                                         keepPhylog = keepPhylog,
                                         mutationPropGrowth = mutationPropGrowth,
                                         max.memory = max.memory, 
                                         max.wall.time = max.wall.time,
                                         max.num.tries = max.num.tries,
                                         errorHitWallTime = errorHitWallTime,
                                         errorHitMaxTries = errorHitMaxTries,
                                         verbosity = verbosity,
                                         initMutant = grid$Genotype,
                                         AND_DrvProbExit = AND_DrvProbExit,
                                         fixation = fixation,
                                         seed = seed,
                                         interventions = interventions,
                                         userVars = userVars,
                                         rules = rules))
            simul_grid <- data.frame(Genotype = simul_grid$GenotypesLabels,
                                     N = simul_grid$pops.by.time[
                                       nrow(simul_grid$pops.by.time), -1],
                                     Coordinates = I(rep(list(
                                       grid$Coordinates[[1]]), 
                                       length(simul_grid$GenotypesLabels))))
          
            # After running OncoSimulIndiv over a grid, if not wild-type 
            # genotypes are in the results, we add them with a population number 
            # of 0. This will let us to distiguish between old grids and newly 
            # created grids in the migration phase.
            if (! "" %in% simul_grid$Genotype){
              simul_grid <- rbind(simul_grid, 
                                     data.frame(Genotype = "",
                                                N = 0,
                                                Coordinates = I(
                                                  list(grid$Coordinates[[1]]))))
            }
            return (simul_grid)
        }
    }
    
    ## Function for selecting the random coordinates for near and remote 
    ## migration.
    
    SimulMigration <- function(grid, migrationProb = 0.5, 
                               largeDistMigrationProb = 1e-6, 
                               maxMigrationPercentage = 0.2){
      # Create a function to define random coordinates for near and remote 
      # migrating cells.
      randomcoordinates <- function(grid, MigrationType = "Near"){
        init_coordinates <- grid$Coordinates[[1]]
        coord_nearmigration <- init_coordinates
        coord_remotemigration <- init_coordinates
        # While loop to guarantee that coordinates for migrating cells are 
        # different from the ones of the grid they belong to. It also ensures 
        # that only positive coordinates are selected for migrating cells.
        
        while ((identical(coord_nearmigration, init_coordinates) ||  
                any(coord_nearmigration < 0)) & MigrationType == "Near"){ 
          coord_nearmigration <- c()
          for (coord in init_coordinates){
            new_coord <- as.numeric(sample(seq(coord - 1, coord + 1), 1))
            coord_nearmigration <- c(coord_nearmigration, new_coord)
          }
        } 
        while ((identical(coord_remotemigration, init_coordinates) ||  
                any(coord_remotemigration < 0)) & MigrationType == "Remote"){
          coord_remotemigration <- c()
          for (coord in init_coordinates){
            new_coord <- as.numeric(sample(c(coord - sample(seq(10, 35), 1),
                                             coord + 
                                               sample(seq(10, 35), 1)), 1)) 
            coord_remotemigration <- c(coord_remotemigration, new_coord)
          }
        }
        if (MigrationType == "Near"){
          return (list(coord_nearmigration))
        } else if (MigrationType == "Remote"){
          return (list(coord_remotemigration))
        }
      } 
      
      TotalPopSize <- sum(grid$N)
      # If the total population corresponds to wild-type ("") genotype, migration 
      # will not be performed, as only tumor cells are allowed to move.
      if ("" %in% grid$Genotype && 
          TotalPopSize == grid$N[which(grid$Genotype == "")]){
        return (grid)
      } else {
        near_probtest <- sample(c("Migration", "NoMigration"), 1, 
                                prob = c(migrationProb, 1 - migrationProb))
        remote_probtest <- sample(c("Migration", "NoMigration"), 1, 
                                  prob = c(largeDistMigrationProb, 
                                           1 - largeDistMigrationProb))
        ## IF THERE IS ANY KIND KIND OF MIGRATION:
        if ((near_probtest == "Migration") || (remote_probtest == "Migration")) 
        {
          finalpopcomp_mut <- grid[which((grid$Genotype != "") 
                                         & (grid$N > 0)),]
          finalpopcomp_mut$N <- finalpopcomp_mut$N/sum(finalpopcomp_mut$N)
          MigPopulation <- sample(seq(from = 1, to = TotalPopSize * 
                                        maxMigrationPercentage/100), 1)
          if (xor((near_probtest == "Migration"),
                  (remote_probtest == "Migration"))) {
            MigPopulationcomp <- as.data.frame(table(
              sample(finalpopcomp_mut$Genotype, 
                     MigPopulation, 
                     replace = TRUE,
                     prob = finalpopcomp_mut$N)))
            colnames(MigPopulationcomp) <-c("Genotype", "N") 
            if (near_probtest == "Migration"){
              MigPopulationcomp$Coordinates <- randomcoordinates(
                grid = grid, "Near")
            } else if (remote_probtest == "Migration"){
              MigPopulationcomp$Coordinates <- randomcoordinates(
                grid = grid, "Remote")
            }
          } else if ((near_probtest == "Migration" & 
                      remote_probtest == "Migration")) {
            # When both migration probabilities are satisfied for a grid, 
            # the number of cells that migrate to near and remote places are
            # randomly determined (at least one of the migrating cells must 
            # go to each site).
            NearMigPop_number <- sample (from = 1, to = MigPopulation - 1, 
                                         1)
            RemoteMigPop_number <- MigPopulation - NearMigPop_number
            NearMigPopulationcomp <- as.data.frame(table(
              sample(finalpopcomp_mut$Genotype, 
                     NearMigPop_number, 
                     replace = TRUE,
                     prob = finalpopcomp_mut$N)))
            RemoteMigPopulationcomp <- as.data.frame(table(
              sample(finalpopcomp_mut$Genotype, 
                     RemoteMigPop_number, 
                     replace = TRUE,
                     prob = finalpopcomp_mut$N)))
            colnames(NearMigPopulationcomp) <-c("Genotype", "N")
            colnames(RemoteMigPopulationcomp) <-c("Genotype", "N")
            NearMigPopulationcomp$Coordinates <- randomcoordinates(
              grid = grid, "Near")
            RemoteMigPopulationcomp$Coordinates <- randomcoordinates(
              grid = grid, "Remote")
            MigPopulationcomp <- rbind(NearMigPopulationcomp, 
                                       RemoteMigPopulationcomp$Coordinates)
          } 
          # For each genotype of migrating cells, the corresponding number of 
          # migrating cells is substrated from the original grid.
          for (gen in unique(MigPopulationcomp$Genotype)) {
            Nmigration_per_genotype <- sum(MigPopulationcomp$N[which(
              MigPopulationcomp$Genotype == gen)])
            grid$N[which(grid$Genotype == gen)] <- (
              grid$N[which(grid$Genotype == gen)] - Nmigration_per_genotype)
          }
          return (list(grid, MigPopulationcomp))
        } else {
          # The original grid is returned without any change in 
          # genotype populations if there is no migration. 
          return (grid)
        }
      }
    }
    
    ## Function for adding wild-type cells to newly created demes.
    
    AddWT_newgrid <- function(grid, initSize) {
      # When a new grid is created, it initializes with initSize wild-type 
      # cells plus the arriving mutant cells. Already existing demes will have 
      # an empty string in one of the genotypes, indicating the WT (although the
      # WT population could be 0).
      if (! "" %in% grid$Genotype) {
        new_grid <- rbind(grid, data.frame(Genotype = "",
                                           N = initSize,
                                           Coordinates = I(list(
                                             grid$Coordinates[[1]]))))
        return (new_grid)
      } else {
        return (grid)
      }
    }
    
    ## Function for generating the plots corresponding to the spatial simulation.
    
    Simulationplot <- function(SpatialOncoSimGrid, dim){
      grid <- SpatialOncoSimGrid$Grids
      maxgen <- lapply(grid, function(x) x[x$N == max(x$N), ])
      joinmaxgen <- bind_rows(maxgen)
      joinmaxgen <- joinmaxgen[joinmaxgen$Genotype != "", ]
      maxgencoor <- select(joinmaxgen, c("Genotype", "Coordinates"))
      
      if (dim == "3D"){
        coordX <- sapply(maxgencoor$Coordinates, function(x) x[1])
        coordY <- sapply(maxgencoor$Coordinates, function(x) x[2])
        coordZ <- sapply(maxgencoor$Coordinates, function(x) x[3])
        
        Edgencoord <- data.frame(Genotype = maxgencoor$Genotype, 
                                 X = coordX, Y = coordY, Z = coordZ)
        scattercolorsu <- plasma(length(unique(Edgencoord$Genotype)))
        scattercolors <- scattercolorsu[as.numeric(as.factor(Edgencoord$Genotype))]
        s3d <- scatterplot3d(Edgencoord[,2:4], pch = 16, color = scattercolors)
        s3d
        legend(s3d$xyz.convert((max(coordX)*1.2), (max(coordY)/2), 
                               (max(coordZ)/2)), 
               legend = levels(as.factor(Edgencoord$Genotype)),
               col =  scattercolorsu, 
               pch = 16)
      }
      
      if (dim == "2D"){
        coordX <- sapply(maxgencoor$Coordinates, function(x) x[1])
        coordY <- sapply(maxgencoor$Coordinates, function(x) x[2])
        
        Zdgencoord <- data.frame(Genotype = maxgencoor$Genotype, 
                                 X = coordX, Y = coordY)
        scattercolorsu <- plasma(length(unique(Zdgencoord$Genotype)))
        scattercolors <- scattercolorsu[as.factor(Zdgencoord$Genotype)]
        plot(Zdgencoord$X, Zdgencoord$Y, pch = 16, col = scattercolors, 
             xlab = "X", ylab = "Y")
        legend((max(coordX))*0.9, (max(coordY)),
               legend = levels(as.factor(Zdgencoord$Genotype)),
               col =  scattercolorsu, 
               pch = 16)
      }
      if (dim == "1D"){
        coordX <- sapply(maxgencoor$Coordinates, function(x) x[1])
        Idgencoord <- data.frame(Genotype = maxgencoor$Genotype, 
                                 X = coordX)
        scattercolorsu <- plasma(length(unique(Idgencoord$Genotype)))
        scattercolors <- scattercolorsu[as.factor(Idgencoord$Genotype)]
        plot(Idgencoord$X, rep(0, length(Idgencoord$X)), 
             xlab = "", ylab = "", 
             pch = 16, col = scattercolors)
        legend((max(coordX))*0.9, 1,
               legend = levels(as.factor(Idgencoord$Genotype)),
               col =  scattercolorsu, 
               pch = 16)
      }
    }
  
    iter <- 1   
    finalgrid_list <- list()
    spatialIterMax <- spatialIterMax + 1  
    while (iter <= spatialIterMax) {
      # 1. Simulation inside each grid. The simulation for each grid is done 
      # during specific time units.
        if (length(finalgrid_list) == 0) {
            simul_grid <- oncoSimulIndiv(fp = fp, model = model, 
                                         numPassengers = numPassengers, 
                                         mu = mu, 
                                         muEF = muEF,
                                         detectionSize = detectionSize, 
                                         detectionDrivers = detectionDrivers,
                                         detectionProb = detectionProb,
                                         sampleEvery = sampleEvery,
                                         initSize = initSize, s = s, sh = sh,
                                         K = K, keepEvery = keepEvery,
                                         minDetectDrvCloneSz = minDetectDrvCloneSz,
                                         extraTime = extraTime,
                                         finalTime = finalTime, 
                                         onlyCancer = onlyCancer,
                                         keepPhylog = keepPhylog,
                                         mutationPropGrowth = mutationPropGrowth,
                                         max.memory = max.memory, 
                                         max.wall.time = max.wall.time,
                                         max.num.tries = max.num.tries,
                                         errorHitWallTime = errorHitWallTime,
                                         errorHitMaxTries = errorHitMaxTries,
                                         verbosity = verbosity,
                                         initMutant = initMutant,
                                         AND_DrvProbExit = AND_DrvProbExit,
                                         fixation = fixation,
                                         seed = seed,
                                         interventions = interventions,
                                         userVars = userVars,
                                         rules = rules)
            simul_grid <- data.frame(Genotype = simul_grid$GenotypesLabels, 
                                     N = simul_grid$pops.by.time[
                                       nrow(simul_grid$pops.by.time), -1])
            if (SpatialModel == "1D"){
                  simul_grid$Coordinates <- rep(list(0), length(
                    simul_grid$Genotype))
                                                
            } else if (SpatialModel == "2D"){
                  simul_grid$Coordinates <- rep(list(c(0,0)), length(
                    simul_grid$Genotype))
            } else if (SpatialModel == "3D"){
                  simul_grid$Coordinates <- rep(list(c(0,0,0)), length(
                    simul_grid$Genotype))
            }
            simul_grid <- AddWT_newgrid(simul_grid, initSize = 0)
            intragrid_list <- list(simul_grid)
        } else {
            intragrid_list <- mclapply(finalgrid_list, function(x) 
                                     oncoSimulIndiv_grid(grid = x), 
                                     mc.cores = mc.cores)
        }
      # After the individual simulation for each grid is completed, the migration 
      # operations will be performed.
      intragrid_list <- Filter(Negate(is.null), intragrid_list) 
      if (iter < spatialIterMax & length(intragrid_list) > 0){
          migration_out_list <- mclapply(
                                intragrid_list, 
                                function(x) 
                                SimulMigration(
                                   x, migrationProb = migrationProb, 
                                   largeDistMigrationProb = largeDistMigrationProb, 
                                   maxMigrationPercentage = maxMigrationPercentage),
                                mc.cores = mc.cores)
          migration_out_df <- suppressWarnings(bind_rows(migration_out_list)) 
          
          if (SpatialModel == "1D"){
              migration_out_df$Coordinates <- unlist(migration_out_df$Coordinates)
              join_grid <- aggregate(N ~ Genotype + Coordinates, 
                                     data = migration_out_df, FUN = sum)
              join_grid <- join_grid[,c("Genotype", "N", "Coordinates")]
              grid_list <- split(join_grid, f = list(join_grid$Coordinates, 
                                                     drop = TRUE))
              conv_to_list <- function(x) {x$Coordinates <- rep(
                                                        list(x$Coordinates[1]), 
                                                        length(x$Genotype)) 
                                            return (x)}
              grid_list <- mclapply(grid_list, function(x) conv_to_list(x), 
                                    mc.cores = mc.cores)
          } else if (SpatialModel == "2D") {
              coord_matrix <- do.call(rbind, migration_out_df$Coordinates)
              colnames(coord_matrix) <- c("Coordinate_x", "Coordinate_y")
              migration_out_df <- cbind(migration_out_df, coord_matrix)
              join_grid <- aggregate(N ~ Genotype + Coordinate_x + Coordinate_y, 
                                     data = migration_out_df, FUN = sum)
              Coordinates <- matrix(c(join_grid$Coordinate_x, 
                                      join_grid$Coordinate_y), 
                                    ncol = 2)
              join_grid$Coordinates <- split(Coordinates, row(Coordinates))
              grid_list <- split(join_grid, f = list(join_grid$Coordinate_x,
                                                     join_grid$Coordinate_y), 
                                                     drop = TRUE)
              grid_list <- mclapply(grid_list, function(x) x[,-c(2,3)], 
                                    mc.cores = mc.cores)
          } else if (SpatialModel == "3D") {
              coord_matrix <- do.call(rbind, migration_out_df$Coordinates)
              colnames(coord_matrix) <- c("Coordinate_x", "Coordinate_y", 
                                          "Coordinate_z")
              migration_out_df <- cbind(migration_out_df, coord_matrix)
              join_grid <- aggregate(N ~ Genotype + Coordinate_x + Coordinate_y 
                                     + Coordinate_z, data = migration_out_df, 
                                     FUN = sum)
              Coordinates <- matrix(c(join_grid$Coordinate_x, 
                                      join_grid$Coordinate_y, 
                                      join_grid$Coordinate_z), ncol = 3)
              join_grid$Coordinates <- split(Coordinates, row(Coordinates))
              grid_list <- split(join_grid, f = list(join_grid$Coordinate_x,
                                                     join_grid$Coordinate_y,
                                                     join_grid$Coordinate_z), 
                                 drop = TRUE)
              grid_list <- mclapply(grid_list, function(x) x[,-c(2:4)], 
                                    mc.cores = mc.cores)
          }
          finalgrid_list <- mclapply(grid_list, function(x) AddWT_newgrid(
            x, initSize = initSize), 
                                     mc.cores = mc.cores)
      } else if (length(intragrid_list) > 0) {
          finalgrid_list <<- intragrid_list
      }
      names(finalgrid_list) <- paste0("Grid", seq(1, length(finalgrid_list)))
      iter <- iter + 1 
    }
    FinalObject <- to_SpatialOncoSimul(finalgrid_list, spatialIterMax)
    if (SpatialPlot == TRUE) {
      Simulationplot(FinalObject, dim = SpatialModel)}
    return (FinalObject)
}

z <- SpatialOncoSimul(fp = fe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 500,
                      mu = 1e-4,
                      initSize = 1000,
                      keepPhylog = FALSE,
                      seed = NULL,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE, initMutant = c("i"), 
                      spatialIterMax = 10, SpatialModel = "1D")
