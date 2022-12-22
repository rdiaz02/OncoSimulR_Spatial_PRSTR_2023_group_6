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

for loop # para cada deme habría que aplicar la función migración con los parámetros
          # que decidamos.
          # Dos tipos de cada migracion (alejada y adyacente). Para adyacente la probabilidad
          # más alta que para migración alejada.
          # Y contar nuevo número de demes (si son demes con nuevas dimensiones).
# loop o algo rollo apply ¿?¿?


# Fase interdeme
num_grid ##redefinir después de cada fase interdeme (la migración puede crear nuevos
         ##demes (grids) o eliminar otros que estaban antes)

