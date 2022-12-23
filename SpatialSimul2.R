library(OncoSimulR)
library(parallel)

## Definimos genotipos: 
fe <- allFitnessEffects(
  data.frame(parent = c("Root", "Root", "i"),
             child = c("u" , "i" , "v"),
             s = c(0.1 , -0.05 , 0.25),
             sh = -1,
             typeDep = "MN"),
  epistasis = c("u:i" = -1,"u:v" = -1))

evalAllGenotypes (fe , order = FALSE, addwt = TRUE)

# Lanzamos la simulación para el primer deme.
osi <- oncoSimulIndiv(fe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 500,
                      mu = 1e-4,
                      initSize = 1000,
                      keepPhylog = FALSE,
                      seed = NULL,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)
osi

grid1 <- data.frame(Genotype = osi$GenotypesLabels, 
                    N = osi$pops.by.time[nrow(osi$pops.by.time), -1],
                    Coordinate_x = as.integer(rep(0, length(osi$GenotypesLabels))),
                    Coordinate_y = as.integer(rep(0, length(osi$GenotypesLabels))))

grid2 <- data.frame(Genotype = osi$GenotypesLabels, 
                    N = osi$pops.by.time[nrow(osi$pops.by.time), -1],
                    Coordinate_x = as.integer(rep(0, length(osi$GenotypesLabels))),
                    Coordinate_y = as.integer(rep(0, length(osi$GenotypesLabels))))


## Posible función para simulación en primer deme y resto de demes.
oncoSimulIndiv_grid <- function(grid, iter,...){
  if (iter == 1){
    grid1 <- oncoSimulIndiv(...)
    init_grid <- data.frame(Genotype = grid1$GenotypesLabels, 
                       N = init_grid$pops.by.time[nrow(grid1$pops.by.time), -1],
                       Coordinate_x = as.integer(rep(0, length(grid1$GenotypesLabels))),
                       Coordinate_y = as.integer(rep(0, length(grid1$GenotypesLabels))))
    class(init_grid) <- 'SpatialOncosimul'
  } else{
    Ngenotypes <- grid$N
    Genotypes <- grid$Genotypes
    next_grid <- oncoSimulIndiv(initSize = Ngenotypes, 
                                       initMutant = Genotypes)
    
    next_grid <- data.frame(Genotype = next_grid$GenotypesLabels, 
                       N = next_grid$pops.by.time[nrow(next_grid$pops.by.time), -1],
                       Coordinate_x = grid$Coordinate_x,
                       Coordinate_y = grid$Coordinate_y)
    class(next_grid) <- 'SpatialOncosimul'
  }}


# La función OncoSimulIndiv_grid habría que pasarla con mclapply a una lista que 
# contenga los demes creados hasta el momento, para que la aplique sobre cada deme. 
# Y nos devolvería una nueva lista con los nuevos demes (cada uno 1 dataframe)
# tras la simulación. Los argumentos de la función OncoSimulIndiv_grid son los
# mismos que para la función OncoSimulIndiv, añadiendo el grid y la iteración en 
# la que nos encontramos dentro de la simulación espacial (podemos sumar 1 a iter
# cada vez que terminamos una fase interdeme).


# FUNCIÓN FASE INTERDEME (PARA 2D).

SimulMigration <- function(grid, migrationProb = 0.5, 
                           largeDistMigrationProb = 1e-6, 
                           maxMigrationPercentage = 0.2){
  # Create a new data.frame with the genotypes and population number of mutant
  # cells in the grid.
  TotalPopSize <- sum(grid$N)
  finalpopcomp_mut <- grid[which(grid$Genotype != "" & 
                                   grid$N > 0),]
  # Modify the previous dataframe so instead of mutant population number in the
  # grid, it contains the frequency of each mutant cell.
  finalpopcomp_mut$N <- finalpopcomp_mut$N/sum(finalpopcomp_mut$N)
  
  
  # Number and genotypes of cells that migrate to adjacent spaces (finalpop_near
  # migration dataa.frame).
  pop_nearmigration <- sample(seq(from = 1, to = TotalPopSize * 
                                    migrationProb * 
                                    maxMigrationPercentage), 1)
  finalpop_nearmigration <- as.data.frame(table(sample(finalpopcomp_mut$Genotype, 
                                                       pop_nearmigration, 
                                                       replace = TRUE,
                                                       prob = finalpopcomp_mut$N)))
  colnames(finalpop_nearmigration) <-c("Genotype", "N")
  
  # Number of cells that migrate to remote spaces.
  pop_remotemigration <- sample(seq(from = 1, to = TotalPopSize * 
                                      largeDistMigrationProb * 
                                      maxMigrationPercentage), 1)
  finalpop_remotemigration <- as.data.frame(table(sample(finalpopcomp_mut$Genotype, 
                                                         pop_remotemigration, 
                                                         replace = TRUE,
                                                         prob = finalpopcomp_mut$N)))
  colnames(finalpop_remotemigration) <-c("Genotype", "N")
  init_coordinates <- c(grid$Coordinate_x[1], grid$Coordinate_y[1])
  coord_nearmmigration <- init_coordinates
  coord_remotemigration <- init_coordinates
  
  # This if statement ensures that only positive coordinates were selected for
  # cells that migrate to adjacent spaces.
  if (0 %in% (init_coordinates)){
    # While loop to guarantee that coordinates for migrating cells are different 
    # from the ones of the grid they belong to.
    while (identical(coord_nearmmigration, init_coordinates)){ 
      coord_nearmmigration <- c(sample(seq(init_coordinates[1],
                                           init_coordinates[1] + 1), 1), 
                                sample(seq(init_coordinates[2], 
                                           init_coordinates[2] + 1), 1))
    }
  } else {
    while (identical(coord_nearmmigration, init_coordinates)){
      coord_nearmmigration <- c(sample(seq(init_coordinates[1] - 1,
                                           init_coordinates[1] + 1), 1), 
                                sample(seq(init_coordinates[2] - 1, 
                                           init_coordinates[2] + 1), 1))
    }}
  
  # This while statement ensures that only positive coordinates were selected for
  # cells that migrate to remote spaces.
  while (identical(coord_remotemigration, init_coordinates) ||  
         any(coord_remotemigration < 0)){
    coord_remotemigration <- c(sample(c(init_coordinates[1] - 
                                          sample(seq(10, 35), 1),
                                        init_coordinates[1] + 
                                          sample(seq(10, 35), 1)), 1), 
                               sample(c(init_coordinates[2] - 
                                          sample(seq(10, 35), 1), 
                                        init_coordinates[2] + 
                                          sample(seq(10, 35), 1)), 1))
    
  } 
  finalpop_nearmigration$Coordinate_x <- coord_nearmmigration[1]
  finalpop_nearmigration$Coordinate_y <- coord_nearmmigration[2]
  finalpop_remotemigration$Coordinate_x <- coord_remotemigration[1]
  finalpop_remotemigration$Coordinate_y <- coord_remotemigration[2]
  total_migration <- rbind(finalpop_nearmigration, finalpop_remotemigration)
  for (gen in unique(total_migration$Genotype)) {
    Nmigration_per_genotype <- sum(total_migration$N[which(
      total_migration$Genotype == gen)])
    grid$N[which(grid$Genotype == gen)] <- (grid$N[which(grid$Genotype == gen)] 
                                            - Nmigration_per_genotype)
  }
  return (list(grid, total_migration))
}

SimulMigration(grid1)
grid_list <- list(grid1, grid2) 
y <- mclapply(grid_list, SimulMigration)



#z <- lapply(y, function (x) x[[1]])
#m <- lapply(y, function (x) x[[2]])
#c(z, m)



##1. De la lista que devuelva el mcapply (donde cada elemento es un dataframe 
# de células que migran), unir todos los dataframes en un único dataframe.

## Para unir utilizar alguna función como merge para unir por coordenadas 
# y genotipos en un solo deme.

## Buscar las coordenadas que no tengan genotipo WT o "".



