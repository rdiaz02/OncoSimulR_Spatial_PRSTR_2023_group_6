genofit <- data.frame(A = c(0, 1, 0, 1),
                      B = c(0, 0, 1, 1),
                      Fitness = c("3 + 5*f_",
                                  "3 + 5*(f_ + f_1 - f_2)",
                                  "3 + 5*(f_ + f_2 - f_1)",
                                  "5 + 6*(f_1 + f_2 + f_1_2)"))

afe <- allFitnessEffects(genotFitness = genofit,
                         frequencyDependentFitness = TRUE, 
                         frequencyType = "rel")

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

x <- data.frame(Genotype = osi$GenotypesLabels, N = osi$pops.by.time[nrow(osi$pops.by.time), -1], 
                stringsAsFactors = TRUE)
genotypes <- osi$GenotypesLabels[which(N>0)]
N <- osi$pops.by.time[nrow(osi$pops.by.time), -1]
N<- N[which(N>0)]
osi <- oncoSimulIndiv(afe,
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

## Posiciones en el espacio.
if (num_grid == 1){
  dim2 <-c(0,0)
}
for loop # para cada deme habría que aplicar la función migración con los parámetros
          # que decidamos.
          # Dos tipos de cada migracion (alejada y adyacente). Para adyacente la probabilidad
          # más alta que para migración alejada.
          # Y contar nuevo número de demes (si son demes con nuevas dimensiones).



##Fase interdeme
num_grid ##redefinir

