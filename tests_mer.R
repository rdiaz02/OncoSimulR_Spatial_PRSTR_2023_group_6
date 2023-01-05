library(testthat)

inittime <- Sys.time()
cat(paste("\n Starting test.SpatialOncoSimul at", date(), "\n"))


test_that("SpatialModel parameter is not 1D, 2D or 3D", {
  fe <- allFitnessEffects(
    data.frame(parent = c("Root", "Root", "i"),
               child = c("u" , "i" , "v"),
               s = c(0.1 , -0.05 , 0.25),
               sh = -1,
               typeDep = "MN"),
    epistasis = c("u:i" = -1,"u:v" = -1))
  evalAllGenotypes (fe , order = FALSE, addwt = TRUE)
  expect_error(SpatialOncoSimul(fp = fe,
                                model = "McFL",
                                onlyCancer = FALSE,
                                finalTime = 500,
                                mu = 1e-4,
                                initSize = 1000,
                                keepPhylog = FALSE,
                                seed = NULL,
                                errorHitMaxTries = FALSE,
                                errorHitWallTime = FALSE, initMutant = c("i"), 
                                spatialIterMax = 10, SpatialModel = "4D"),
               "numbers of columns of arguments do not match", ignore.case = TRUE)
})

test_that("errors when initSize and initMutant are not length 1", {
  fe <- allFitnessEffects(
    data.frame(parent = c("Root", "Root", "i"),
               child = c("u" , "i" , "v"),
               s = c(0.1 , -0.05 , 0.25),
               sh = -1,
               typeDep = "MN"),
    epistasis = c("u:i" = -1,"u:v" = -1))
  evalAllGenotypes (fe , order = FALSE, addwt = TRUE)
  expect_error(SpatialOncoSimul(fp = fe,
                                model = "McFL",
                                onlyCancer = FALSE,
                                finalTime = 500,
                                mu = 1e-4,
                                initSize = c(1000,1),
                                keepPhylog = FALSE,
                                seed = NULL,
                                errorHitMaxTries = FALSE,
                                errorHitWallTime = FALSE, initMutant = c("i","u"), 
                                spatialIterMax = 10, SpatialModel = "3D"),
               "(*Argument 3 must have names.)|(second argument must be a list)", ignore.case = TRUE)
})

test_that("SpatialOncoSimul return a SpatialOncosimul object", {
  fe <- allFitnessEffects(
    data.frame(parent = c("Root", "Root", "i"),
               child = c("u" , "i" , "v"),
               s = c(0.1 , -0.05 , 0.25),
               sh = -1,
               typeDep = "MN"),
    epistasis = c("u:i" = -1,"u:v" = -1))
  evalAllGenotypes (fe , order = FALSE, addwt = TRUE)
  p <- SpatialOncoSimul(fp = fe,
                        model = "McFL",
                        onlyCancer = FALSE,
                        finalTime = 500,
                        mu = 1e-4,
                        initSize = 1000,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = FALSE,
                        errorHitWallTime = FALSE, initMutant = c("i"), 
                        spatialIterMax = 10, SpatialModel = "3D")
  expect_that(p, is_a("SpatialOncoSimul"))
})

test_that("Final number of grids equals final number of demes", {
  fe <- allFitnessEffects(
    data.frame(parent = c("Root", "Root", "i"),
               child = c("u" , "i" , "v"),
               s = c(0.1 , -0.05 , 0.25),
               sh = -1,
               typeDep = "MN"),
    epistasis = c("u:i" = -1,"u:v" = -1))
  evalAllGenotypes (fe , order = FALSE, addwt = TRUE)
  p <- SpatialOncoSimul(fp = fe,
                        model = "McFL",
                        onlyCancer = FALSE,
                        finalTime = 500,
                        mu = 1e-4,
                        initSize = 1000,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = FALSE,
                        errorHitWallTime = FALSE, initMutant = c("i"), 
                        spatialIterMax = 10, SpatialModel = "3D")
  expect_that(length(finalgrid_list), is_identical_to(p$Total_Demes))
})

cat(paste("\n Ending test.SpatialOncoSimul at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
