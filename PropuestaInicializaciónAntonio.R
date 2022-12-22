library(OncoSimulR)

# Aquí definimos los efectos del fitness, afe y bfe.
genofit <- data.frame(A = c(0, 1, 0, 1),
                      B = c(0, 0, 1, 1),
                      Fitness = c("3 + 5*f_",
                                  "3 + 5*(f_ + f_1 - f_2)",
                                  "3 + 5*(f_ + f_2 - f_1)",
                                  "5 + 6*(f_1 + f_2 + f_1_2)"))

afe <- allFitnessEffects(genotFitness = genofit,
                         frequencyDependentFitness = TRUE, 
                         frequencyType = "rel")

bfe <- allFitnessEffects(epistasis = c("A:-B" = 0.1, "B:-A" = 0.4, "A:B" = 0.2,
                                      "C:-A:-B" = 0.5, "C:A" = -0.5, 
                                      "C:B" = 0.4))

# La idea es que GrindInitializer sea la función que cree el primer deme. 
# Aceptaría los argumentos que luego irían en OncoSimulIndiv y devolvería un 
# dataframe.
GridInitializer <- function(...){
  siminit <- oncoSimulIndiv(...)
  sumsinit <- data.frame(Genotype = siminit$GenotypesLabels, 
                         N = siminit$pops.by.time[nrow(siminit$pops.by.time), -1],
                         stringsAsFactors = TRUE)
  N <- as.numeric(siminit$pops.by.time[nrow(siminit$pops.by.time), -1])
  N <- N[which(N>0)]
  genotypes <- siminit$GenotypesLabels[which(N>0)]
  return(data.frame(dem = factor(rep(1, length(genotypes))),
                    Coordinate_X = rep(0, length(genotypes)),
                    Coordinate_y = rep(0, length(genotypes)),
                    Poblational_Size = N, 
                    Genotypes = genotypes))
}

grid <- GridInitializer(fp = afe,
                        model = "McFL",
                        onlyCancer = FALSE,
                        finalTime = 200,
                        mu = 1e-6,
                        initSize = 5000,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = FALSE,
                        errorHitWallTime = FALSE)
# El deme se identificaría con un número en forma de factor, y cada coordenada
# x e y iría en una columna diferente. 
grid

# No se usaría GrindInitializer para crear otro deme. Estos irían surgiendo de 
# de la función de migración, aunque aún no he pensado como hacerla. Utilizo
# aquí un apaño para probar a crear un dataframe con dos demes, lo que 
# representaría el grid.

GridInitializer2 <- function(...){
  siminit <- oncoSimulIndiv(...)
  sumsinit <- data.frame(Genotype = siminit$GenotypesLabels, 
                         N = siminit$pops.by.time[nrow(siminit$pops.by.time), -1],
                         stringsAsFactors = TRUE)
  N <- as.numeric(siminit$pops.by.time[nrow(siminit$pops.by.time), -1])
  N <- N[which(N>0)]
  genotypes <- siminit$GenotypesLabels[which(N>0)]
  return(data.frame(dem = factor(rep(2, length(genotypes))),
                    Coordinate_X = rep(1, length(genotypes)),
                    Coordinate_y = rep(1, length(genotypes)),
                    Poblational_Size = N, 
                    Genotypes = genotypes))
}

grid2 <- rbind(grid, GridInitializer2(bfe,
                                      model = "Exp",
                                      onlyCancer = FALSE,
                                      finalTime = 200,
                                      mu = 1e-6,
                                      initSize = 5000,
                                      keepPhylog = FALSE,
                                      seed = NULL,
                                      errorHitMaxTries = FALSE,
                                      errorHitWallTime = FALSE))
grid2
