library(OncoSimulR)

# Cosas generales:
  # las células migran en grupo
  # 1 deme es 1 dimensión
  # correr el indiv (intrademe), mover las células (interdeme), indiv otra vez...

# Definimos los datos
genofit <- data.frame(A = c(0, 1, 0, 1),
                      B = c(0, 0, 1, 1),
                      Fitness = c("3 + 5*f_",
                                  "3 + 5*(f_ + f_1 - f_2)",
                                  "3 + 5*(f_ + f_2 - f_1)",
                                  "5 + 6*(f_1 + f_2 + f_1_2)"))
genofit

# Fitness
afe <- allFitnessEffects(genotFitness = genofit,
                         frequencyDependentFitness = TRUE, 
                         frequencyType = "rel")
afe

# Primera simulación
osi <- oncoSimulIndiv(afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 200,
                      mu = 1e-6,
                      initSize = 5000,
                      keepPhylog = FALSE,
                      seed = NULL,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)
osi
  ## Necesitamos que el output de esta función tenga formato que pueda ser 
  ## pasado de nuevo a la función:

# Extraemos el dataframe con los genotipos y el número de células por genotipo
df <- data.frame(Genotype = osi$GenotypesLabels, 
                 N = osi$pops.by.time[nrow(osi$pops.by.time), -1], 
                 stringsAsFactors = TRUE)

# Separamos los datos y eliminamos los genotipos que tienen 0 células
genotypes <- osi$GenotypesLabels[which(df$N > 0)]
N <- osi$pops.by.time[nrow(osi$pops.by.time), -1][which(df$N > 0)]
  ## habría que meter esto generalizado en la función

# Ahora podemos hacer una nueva simulación con el resultado de la primera
osi2 <- oncoSimulIndiv(afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 200,
                      mu = 1e-6,
                      initSize = N,
                      initMutant = genotypes,
                      keepPhylog = FALSE,
                      seed = NULL,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)
osi2


# Posiciones en el espacio.

# 1. crear objeto deme con coordenadas, genotipos y num de celulas de cada genotipo. 
## pasar de df a objeto

df$Coordinates <- list(c(0,0), c(0,0), c(0,0), c(0,0)) 
#para hacer pruebas añado las coordenadas al dataframe

demeobj <- function(x, ...) {
  UseMethod("deme")
}

deme.data.frame <- function(x) {
  cns <- c("Genotype", "N", "Coordinates")
  if (!(all(colnames(x) %in% cns)))
    stop(paste("Column names are not ", cns))
  deme <- x[, cns]
  class(deme) <- c("deme")
  return(deme)
}

uu <- deme(df)


# 2. función de migración a la que se le pasan los demes.

if (num_grid == 1){
  dim2 <-c(0,0) ## el primer grid está en la posición 0,0 
}

#for loop # para cada deme habría que aplicar la función migración con los parámetros
          # que decidamos.
          # Dos tipos de cada migracion (alejada y adyacente). Para adyacente la probabilidad
          # más alta que para migración alejada.
          # Y contar nuevo número de demes (si son demes con nuevas dimensiones).
# loop o algo rollo apply ¿?¿?


# Fase interdeme
num_grid ##redefinir después de cada fase interdeme (la migración puede crear nuevos
         ##demes (grids) o eliminar otros que estaban antes)





#### FASE INTERDEME.

## Definimos genotipos: 
fe <- allFitnessEffects(
  data.frame(parent = c("Root", "Root", "i"),
             child = c("u" , "i" , "v"),
             s = c(0.1 , -0.05 , 0.25),
             sh = -1,
             typeDep = "MN"),
  epistasis = c("u:i" = -1,"u:v" = -1))

evalAllGenotypes (fe , order = FALSE, addwt = TRUE)

## Creamos primero un deme:
grid <- oncoSimulIndiv(fe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 500,
                      mu = 1e-4,
                      initSize = 1000,
                      keepPhylog = FALSE,
                      seed = NULL,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)
## Añadimos atributo coordenadas a la clase oncoSimul y lo podemos transformar
## en una nueva clase que se llame 'SpatialOncosimul'.

grid$coord <- as.integer(c(0,0)) # Las coordenadas son numeros enteros (importante
# comparar coordenadas con identical y obtener TRUE cuando los objetos son iguales)

# SpatialOncosimul <- osi
# class(SpatialOncosimul) <- 'SpatialOncosimul'
# class(SpatialOncosimul)
# str(SpatialOncosimul)



##FUNCIÓN PARA LA MIGRACIÓN DE CÉLULAS DESDE UN DEME/GRID A REGIONES ADYACENTES
## O ALEJADAS.
SimulMigration<- function(grid, migrationProb = 0.5, 
                          largeDistMigrationProb = 1e-6, 
                          maxMigrationPercentage = 0.2){
    finalpopcomp <- data.frame(Genotype = grid$GenotypesLabels, 
                               N = grid$pops.by.time[nrow(grid$pops.by.time), -1])
    # Create a new data.frame with the genotypes and population number of mutant
    # cells in the grid.
    finalpopcomp_mut <- finalpopcomp[which(finalpopcomp$Genotype != "" & 
                                             finalpopcomp$N > 0),]
    # Modify the previous dataframe so instead of mutant population number in the
    # grid, it contains the frequency of each mutant cell.
    finalpopcomp_mut$N <- finalpopcomp_mut$N/sum(finalpopcomp_mut$N)
    
    
    # Number and genotypes of cells that migrate to adjacent spaces (finalpop_near
    # migration dataa.frame).
    pop_nearmigration <- sample(seq(from = 1, to = grid$TotalPopSize * 
                                      migrationProb * 
                                      maxMigrationPercentage), 1)
    finalpop_nearmigration <- as.data.frame(table(
                                            sample(finalpopcomp_mut$Genotype, 
                                                   pop_nearmigration, 
                                                   replace = TRUE,
                                                   prob = finalpopcomp_mut$N)))
    colnames(finalpop_nearmigration) <-c("Genotype", "N")
    
    # Number of cells that migrate to remote spaces.
    pop_remotemigration <- sample(seq(from = 1, to = grid$TotalPopSize * 
                                        largeDistMigrationProb * 
                                        maxMigrationPercentage), 1)
    finalpop_remotemigration <- as.data.frame(table(
                                              sample(finalpopcomp_mut$Genotype, 
                                                     pop_remotemigration, 
                                                     replace = TRUE,
                                                     prob = finalpopcomp_mut$N)))
    colnames(finalpop_remotemigration) <-c("Genotype", "N")
    coord_nearmmigration <- grid$coord
    coord_remotemigration <- grid$coord
    
    # This if statement ensures that only positive coordinates were selected for
    # cells that migrate to adjacent spaces.
    if (0 %in% (grid$coord)){
      # While loop to guarantee that coordinates for migrating cells are different 
      # from the ones of the grid they belong to.
      while (identical(coord_nearmmigration, grid$coord)){ 
        coord_nearmmigration <- c(sample(seq(grid$coord[1],
                                             grid$coord[1] + 1), 1), 
                                  sample(seq(grid$coord[2], 
                                             grid$coord[2] + 1), 1))
      }
    } else {
      while (identical(coord_nearmmigration, grid$coord)){
        coord_nearmmigration <- c(sample(seq(grid$coord[1]-1,
                                             grid$coord[1]+ 1), 1), 
                                  sample(seq(grid$coord[2]-1, 
                                             grid$coord[2]+ 1), 1))
      }}
    
    # This while statement ensures that only positive coordinates were selected for
    # cells that migrate to remote spaces.
    while (identical(coord_remotemigration, grid$coord) ||  
           any(coord_remotemigration < 0)){
      coord_remotemigration <- c(sample(
                                c(grid$coord[1] - sample(seq(10, 35), 1),
                                grid$coord[1] + sample(seq(10, 35), 1)), 1), 
                                sample(
                                c(grid$coord[2] - sample(seq(10, 35), 1), 
                                grid$coord[2] + sample(seq(10, 35), 1)), 1))
    } 
    
    finalpop_nearmigration$Coordinates <- list(coord_nearmmigration)
    finalpop_remotemigration$Coordinates <- list(coord_remotemigration)
    total_migration <- rbind(finalpop_nearmigration, finalpop_remotemigration)
    return (total_migration)
}






# Esta función habría que aplicarla con un mcapply para cada deme que exista en 
# ese momento.
SimulMigration(osi)





finalpopcomp_mut <- finalpopcomp[which(finalpopcomp$Genotype != "" & 
                                         finalpopcomp$N > 0),]
# Modify the previous dataframe so instead of mutant population number in the
# osi, it contains the frequency of each mutant cell.
finalpopcomp_mut$N <- finalpopcom_pmut$N/sum(finalpopcom_pmut$N)


# Number and genotypes of cells that migrate to adjacent spaces (finalpop_near
# migration dataa.frame).
pop_nearmigration <- sample(seq(from = 1, to = osi$TotalPopSize * 
                                  migrationProb * 
                                  maxMigrationPercentage), 1)
finalpop_nearmigration <- as.data.frame(table(
  sample(finalpopcomp_mut$Genotype, 
         pop_nearmigration, 
         replace = TRUE,
         prob = finalpopcomp_mut$N)))
colnames(finalpop_nearmigration) <-c("Genotype", "N")

# Number of cells that migrate to remote spaces.
pop_remotemigration <- sample(seq(from = 1, to = osi$TotalPopSize * 
                                    largeDistMigrationProb * 
                                    maxMigrationPercentage), 1)
finalpop_remotemigration <- as.data.frame(table(
  sample(finalpopcomp_mut$Genotype, 
         pop_remotemigration, 
         replace = TRUE,
         prob = finalpopcomp_mut$N)))
colnames(finalpop_remotemigration) <-c("Genotype", "N")
coord_nearmmigration <- osi$coord
coord_remotemigration <- osi$coord

# This if statement ensures that only positive coordinates were selected for
# cells that migrate to adjacent spaces.
if (0 %in% (osi$coord)){
  # While loop to guarantee that coordinates for migrating cells are different 
  # from the ones of the osi they belong to.
  while (identical(coord_nearmmigration, osi$coord)){ 
    coord_nearmmigration <- c(sample(seq(osi$coord[1],
                                         osi$coord[1] + 1), 1), 
                              sample(seq(osi$coord[2], 
                                         osi$coord[2] + 1), 1))
  }
} else {
  while (identical(coord_nearmmigration, osi$coord)){
    coord_nearmmigration <- c(sample(seq(osi$coord[1]-1,
                                         osi$coord[1]+ 1), 1), 
                              sample(seq(osi$coord[2]-1, 
                                         osi$coord[2]+ 1), 1))
  }}

# This while statement ensures that only positive coordinates were selected for
# cells that migrate to remote spaces.
while (identical(coord_remotemigration, osi$coord) ||  
       any(coord_remotemigration < 0)){
  coord_remotemigration <- c(sample(
    c(osi$coord[1] - sample(seq(10, 35), 1),
      osi$coord[1] + sample(seq(10, 35), 1)), 1), 
    sample(
      c(osi$coord[2] - sample(seq(10, 35), 1), 
        grid$coord[2] + sample(seq(10, 35), 1)), 1))
  
} 

finalpop_nearmigration$Coordinates <- list(coord_nearmmigration)
finalpop_remotemigration$Coordinates <- list(coord_remotemigration)
total_migration <- rbind(finalpop_nearmigration, finalpop_remotemigration)
return (total_migration)
}
