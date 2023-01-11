\name{SpatialOncoSimul}
\alias{SpatialOncoSimul}
\alias{print.oncosimul}
\alias{summary.oncosimul}


\title{
  SpatialOncoSimul with OncoSimulR. Spatial simulation of tumoral cell growth.
}

\description{
  SpatialOncoSimul is a function implemented inside the package OncoSimulR for the
  spatial simulation of tumoral cell growth. This function simulates tumor 
  progression from one starting position in a 1D, 2D or 3D grid specified by the user. 
  After the first population of tumoral cells is developed with OncoSimulIndiv, a 
  small population of tumoral cells will migrate to adjacent or remote spaces based 
  on the fitness of the tumor genotypes. The migrating cells will occupy new space 
  coordinates and will be combined with an initial population of wild-type cells
  present in the new location. Then, the OncoSimulIndiv function will be launched 
  again in each space coordinate to simulate tumor progression from the new tumor 
  composition in each position. This process is repeated the number of times specified 
  by the user to analyse the tumor progression in a spatial context.
}

\usage{
SpatialOncoSimul(fp, model = "McFL", 
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
                 SpatialPlot = TRUE)
}

\arguments{

  \item{fp}{

    Either a poset that specifies the order restrictions (see
    \code{\link{poset}} if you want to use the specification as in
    v.1. Otherwise, a fitnessEffects object (see
    \code{\link{allFitnessEffects}}). You must always use a fitnessEffects object
    when you are in a frequency dependent fitness simulation; of course in this case 
    fp$frequencyDependentFitness must be TRUE.

    Other arguments below (s, sh, numPassengers) make sense only if you
    use a poset, as they are included in the fitnessEffects object.

}

\item{model}{ One of "Bozic", "Exp", "Arb", "McFarlandLog", "McFarlandLogD"
  (the last two can be abbreviated to "McFL" and "McFLD",
  respectively). The default is "McFL". (See OncoSimulR vignette for the difference
  between "McFL" and "McFLD": in the former, death rate = \eqn{\log(1 +
  N/K)} where \eqn{K} is the initial equilibrium population size; when
  using "McFLD", death rate = \eqn{\max(1, \log(1 + N/K))}, so that death
  rate never goes below 1.). If "Arb" (arbitrary) model is specified, death
  must be present in \code{\link{allFitnessEffects}}, and vice versa.
  

}

\item{numPassengers}{

  This has no effect if you use the \code{\link{allFitnessEffects}}
  specification. This happens always when you are in a simulation that use
  frequency dependent fitness.
  
  If you use the specification of v.1., the number of passenger
  genes. Note that using v.1 the total number of genes (drivers plus
  passengers) must be smaller than 64.
  
  All driver genes should be included in the poset (even if they depend
  on no one and no one depends on them), and will be numbered from 1 to
  the total number of driver genes. Thus, passenger genes will be
  numbered from (number of driver genes + 1):(number of drivers + number
  of passengers).

}

\item{mu}{
  Mutation rate. Can be a single value or a named vector. If a single
  value, all genes will have the same mutation rate. If a named vector,
  the entries in the vector specify the gene-specific mutation rate. If
  you pass a vector, it must be named, and it must have entries for all
  the genes in the fitness specification. Passing a vector is only
  available when using fitnessEffects objects for fitness specification.
  Mutation rates \code{<10^-20} are not accepted.
  See also \code{mutationPropGrowth}. 

}

\item{muEF}{

  Mutator effects. A mutatorEffects object as obtained from
  \code{\link{allMutatorEffects}}. This specifies how mutations in
  certain genes change the mutation rate over all the genome. Therefore,
  this allows you to specify mutator phenotypes: models where mutation
  of one (or more) gene(s) leads to an increase in the mutation
  rate. This is only available for version 2 (and above) specifications.
  
  
  All the genes specified in \code{muEF} MUST be included in
  \code{fp}. If you want to have genes that have no direct effect on
  fitness, but that affect mutation rate, you MUST specify them in
  \code{fp}, for instance as \code{noIntGenes} with an effect of 0.
  
  If you use mutator effects you must also use fitnessEffects in
  \code{fp}.

  %% You are not allowed to use mutator effects object in a frequency depedent
  %% fitness simulation.

}

\item{detectionSize}{ What is the minimal number of cells for cancer
to be detected. 

If set to NA, \code{detectionSize} plays no role in stopping the
simulations.

} 

\item{detectionDrivers}{ The minimal number of drivers (not modules,
drivers, whether or not they are from the same module) present in any
clone for cancer to be detected.

If set to NA, \code{detectionDrivers} plays no role in stopping the
simulations.

}

\item{detectionProb}{ Vector of arguments for the mechanism where
  probability of detection depends on size. If \code{NA}, this mechanism
  is not used. If \sQuote{default}, the vector will be populated with
  default values.  Otherwise, a named vector with some of the following
  named elements (see \sQuote{Details}):
  \itemize{

    \item{PDBaseline:}{ Baseline size subtracted to total population size
    to compute the probability of detection. If not given explicitly,
    the default is \code{1.2 * initSize} (or \code{1.2 * sum(initSize)}
    when multiple \code{initMutant}s).
    }

    \item{p2:}{ The probability of detection at population size
    \code{n2}. If you specificy \code{p2} you must also specify
    \code{n2} and you must not specify \code{cPDetect}. The fault is 0.1.}

    \item{n2:}{ The population size at which probability of detection is
      \code{p2}. The default is \code{2 * initSize}.}

    \item{cPDetect:}{ The change in probability of detection with size. If
      you specify it, you should not specify either of \code{p2} or
      \code{n2}. See \sQuote{Details}.  }

    \item{checkSizePEvery:}{ Time between successive checks for the
    probability of exiting as a function of population size. If not
    given explicitly, the default is 20. See \sQuote{Details}.}

}
If you only provide some of the elements (except for the pair \code{p2},
\code{n2}, where you must provide both if you provide any), the rest
will be filled with default values.

This option can not be used with v.1 objects.

}


%% \item{PDBaseline}{For the mechanism where probability of detection
%%   depends on size, the baseline. See \sQuote{Details}.
%% }

%% \item{p2}{For the mechanism where probability of detection
%%   depends on size, the probability of detection at size \code{n2}. See \sQuote{Details}.}

%% \item{n2}{For the mechanism where probability of detection depends on
%%   size, the size at which probability of detection is \code{p2}. See
%%   \sQuote{Details}. }

%% \item{cPDetect}{ For the mechanism where probability of detection
%%   depends on size, the change in probability of detection with size. If
%%   you specify it, you should not specify either of \code{p2} or
%%   \code{n2}. See \sQuote{Details}.  }

%% \item{checkSizePEvery}{ For the mechanism where probability of detection
%%   depends on size, time between successive checks. See \sQuote{Details}.  }

\item{sampleEvery}{
  
  How often the whole population is sampled. This is not the same as the
  interval between successive samples that are kept or stored (for that,
  see \code{keepEvery}).

  For very fast growing clones, you might need to have a small value
  here to minimize possible numerical problems (such as huge increase in
  population size between two successive samples that can then lead to
  problems for random number generators). Likewise, for models with
  density dependence (such as McF) this value should be very small.

}
\item{initSize}{
  Initial population size. If you are passing more than one
  \code{initMutant}, the initial population sizes of each
  clone/species/genotype, given in the same order as in the
  \code{initMutant} vector. \code{initMutant} thus allows to start the
  simulation from arbitrary population compositions. Combined with
  \code{mu} it allows for multispecies simulations (see the OncoSimulR vignette 
  for examples).

}


\item{K}{
  Initial population equilibrium size in the McFarland models.

}

\item{keepEvery}{

  Time interval between successive whole population samples that are
actually stored. This must be larger or equal to \code{sampleEvery}. If
keepEvery is not a multiple integer of sampleEvery, the interval between
successive samples that are stored will be the smallest multiple integer
of \code{sampleEvery} that is larger than or equal to \code{keepEvery}.

  %% keepEvery in use will be the
  %% smallest multiple integer of sampleEvery larger than the specified
  %% keepEvery.

  If you want nice plots, set \code{sampleEvery} and \code{keepEvery} to
  small values (say, 5 or 2). Otherwise, you can use a
  \code{sampleEvery} of 1 but a \code{keepEvery} of 15, so that the
  return objects are not huge and the code runs a lot faster.


  Setting \code{keepEvery = NA} means we only keep the very last
  sample. This is useful if you only care about the final state of the
  simulation, not its complete history.
}

%% \item{endTimeEvery}{
%%   If endTimeEvery is > 0, even if conditions for finishing a simulation
%%   (number of drivers or population size) are met at time \emph{t}, the
%%   simulation will run at least until \emph{t + endTimeEvery} and
%%   conditions will be checked again. Only if conditions for finishing a
%%   simulation are still met, will the simulation end.

%%   The reason for this parameter is to ensure that, say, a clone with a
%%   certain number of drivers that would cause the simulation to end has
%%   not just appeared but then gone extinct shortly after. Beware, though,
%%   that in simulations with very fast growth, setting large endTimeEvery
%%   can result in the simulations taking a long time to finish or hitting
%%   the wall time limit.}

\item{minDetectDrvCloneSz}{

  A value of 0 or larger than 0 (by default equal to
  \code{initSize} in the McFarland model). If larger than 0, when
  checking if we are done with a simulation, we verify that the sum of
  the population sizes of all clones that have a number of mutated
  drivers larger or equal to \code{detectionDrivers} is larger or equal
  to this \code{minDetectDrvCloneSz}.

  The reason for this parameter is to ensure that, say, a clone with a
  certain number of drivers that would cause the simulation to end has
  not just appeared and is present in only one individual that might
  then immediately go extinct. This can be relevant in secenarios such
  as the McFarland model.

  If \code{initSize} is larger than 1 (you are passing multiple
  initMutants), the sum is used.

  See also \code{extraTime}.
}

\item{extraTime}{
  A value larger than zero waits those many additional time periods
  before exiting after having reached the exit condition (population
  size, number of drivers).

  The reason for this setting is to prevent the McFL models from always
  exiting at a time when one clone is increasing its size quickly (see
  \code{minDetectDrvCloneSz}). By setting an \code{extraTime} larger than 0,
  we can sample at points when we are at the plateau.
  

}

\item{finalTime}{
  What is the maximum number of time units that the simulation can
  run. Set to NA to disable this limit.

}
\item{onlyCancer}{
  Return only simulations that reach cancer?

  If set to TRUE, only simulations that satisfy the
  \code{detectionDrivers} or the \code{detectionSize} requirements or
  that are "detected" because of the \code{detectionProb} mechanism will
  be returned: the simulation will be repeated, within the limits set by
  \code{max.num.tries} and \code{max.wall.time} (and, for
  \code{oncoSimulSample} also \code{max.num.tries.total} and
  \code{max.wall.time.total}), until one which meets the
  \code{detectionDrivers} or \code{detectionSize} or one which is
  detected stochastically under \code{detectionProb} is
  obtained.

  
  If \code{onlyCancer = FALSE} the simulation is returned regardless of
  final population size or number of drivers in any clone and this
  includes simulations where the population goes extinct.

  The default used to be \code{onlyCancer = TRUE}; on version 3.99.10
  it was changed to \code{onlyCancer = FALSE} as this is the natural
  setting for simulating general scenarios.   \code{onlyCancer = TRUE},
  by design, leads to selection bias in the   simulations returned: we only see those that "reach cancer".

}

\item{keepPhylog}{
  If TRUE, keep track of when and from which clone each clone is
  created. See also \code{\link{plotClonePhylog}}. 
}



\item{mutationPropGrowth}{ If TRUE, make mutation rate proportional to growth
  rate, so clones that grow faster also mutate faster (laso have a larger
  mutation rate): $mutation_rate = mu
  * birth_rate$. With BNB mutation is actually "mutate after division": p.\ 1232 of Mather
  et al., 2012 explains:
  "(...) mutation is simply defined as the creation and subsequent departure of a single individual from the class". Thus,
  if we want to have individuals of clones/genotypes/populations that divide faster to also produce more
  mutants per unit time (per individual) we have to set \code{mutationPropGrowth = TRUE}.
  %% This is a simple way of approximating that mutation happens
  %% at cell division (it is not strictly making mutation happen at cell division,
  %% since mutation is not strictly coupled with division).
  Of course, this only makes sense in models where birth rate changes.
}

%% \item{seed}{The seed to use for the C++ random number generator. If not
%%   passed explicitly (the default) then chosen from R.}

\item{initMutant}{For v.2: a string with the mutations of the initial
  mutant, if any. This is the same format as for
  \code{\link{evalGenotype}}. The default (if you pass nothing) is to
  start the simulation from the wildtype genotype with nothing mutated.
  For v.1 we no longer accept initMutant: it will be ignored.

  (\code{evalGenotype} also accepts the genotype as a numeric vector;
  \code{initMutant} must be a character string.)
}

%% \item{initMutant}{For v.2, a string with the mutations of the initial
%%   mutant, if any. This is the same format as for
%%   \code{\link{evalGenotype}}. For v.1, the single mutation of the
%%   initial clone for the simulations. The default (if you pass nothing)
%%   is to start the simulation from the wildtype genotype with nothing
%%   mutated.}

%%\item{constantPop}{Indicates if the total population must be constant or not.}

\item{max.num.tries}{Only applies when \code{onlyCancer = TRUE}. What is
  the maximum number of times, for an individual simulation, we can
  repeat the simulation for it to reach cancer? There are certain
  parameter settings where reaching cancer is extremely unlikely and you
  might not want to run forever in those cases.  }

\item{max.wall.time}{ Maximum wall time for the simulation of one
  individual (over all \code{max.num.tries}). If the simulation is
  not done in this time, it is aborted.  %% (If \code{max.num.tries = k},
  %% time is kept track of over the k attempts;)
}

\item{errorHitMaxTries}{ If TRUE (the default) a simulation that reaches
  the maximum number of repetitions allowed is considered not to have
  succesfully finished and, thus, an error, and no output from it will
  be reported. This is often what you want.%%  But if you want to see
  %% what is happening for some set of parameters that seem to never
  %% finish, set this to FALSE, and check what the output looks like.
  See  \code{Details}.
}

\item{errorHitWallTime}{
  If TRUE (the default) a simulation that reaches the maximum wall time
  is considered not to have succesfully finished and, thus, an error,
  and no output from it will be reported. This is often what you
  want. %% But if you want to see what is happening for some set of
  %% parameters that seem to never finish, set this to FALSE, and check
  %% what the output looks like.
  See \code{Details}.
}


\item{max.memory}{
  The largest size (in MB) of the matrix of Populations by Time. If it
  creating it would use more than this amount of memory, it is not
  created. This prevents you from accidentally passing parameters that
  will return an enormous object.

}


\item{verbosity}{ If 0, run silently. Iincreasing values of verbosity
  provide progressively more information about intermediate steps,
  possible numerical notes/warnings from the C++ code, etc. Values less
  than 0 supress some default notes: use with care.

}


\item{mc.cores}{Number of cores to use when simulating more than one
  individual (i.e., when calling oncoSimulPop).}

\item{AND_DrvProbExit}{If TRUE, cancer will be considered to be reached
  if both the \code{detectionProb} mechanism and \code{detectionDrivers}
  are satisfied. This is and AND, not an OR condition. Using this option
  with \code{fixation} is not allowed (as it does not make much sense).}

\item{fixation}{If non-NULL, a list or a vector, where each element of
  is a string with a gene or a gene combination or a genotype (see
  below). Simulations will stop as soon as any of the genes or gene
  combinations or genotypes are fixed (i.e., reach a minimal
  frequency). If you pass gene combinations or genotypes, separate genes
  with commas (not '>'); this means order is not (yet?)  supported. This
  way of specifying gene combinations is the same as the one used for
  \code{initMutant} and \code{\link{evalGenotype}}.

  To differentiate between gene combinations and specific genotypes,
  genotypes are specified by prepending them with a "_,". For instance,
  \code{fixation = c("A", "B, C")} specifies stopping on any genotypes
  with those gene combinations. In contrast, \code{fixation = c("_,A",
  "_,B, C" )} specifies stopping only on gentoypes "A" or "B, C". See
  the vignette for further examples.
  

  In addition to the gene combinations or genotypes themeselves, you can
  add to the list or vector the named elements
  \code{fixation_tolerance}, \code{min_successive_fixation} and
  \code{fixation_min_size}. \code{fixation_tolerance}: fixation is
  considered to have happened if the genotype/gene combinations
  specified as genotypes/gene combinations for fixation have reached a
  frequency \code{> 1 - fixation_tolerance}. (The default is 0, so we
  ask for genotypes/gene combinations with a frequency of 1, which might
  not be what you want with large mutation rates and complex fitness
  landscape with genotypes of similar
  fitness.). \code{min_successive_fixation}: during how many successive
  sampling periods the conditions of fixation need to be fulfilled
  before declaring fixation. These must be successive sampling periods
  without interruptions (i.e., a single period when the condition is not
  fulfilled will set the counter to 0). This can help to exclude short,
  transitional, local maxima that are quickly replaced by other
  genotypes. (The default is 50, but this is probably too small for
  ``real life'' usage). \code{fixation_min_size}: you might only want to
  consider fixation to have happened if a minimal size has been reached
  (this can help weed out local maxima that have fitness that is barely
  above that of the wild-type genotype). (The default is 0).

  Using this option with \code{AND_DrvProbExit} is not allowed (as it
  does not make much sense). This option is not allowed either with the
  old v.1 specification.}

\item{s}{
  Selection coefficient for drivers. 
    Only relevant if using a poset as this is included in the
  fitnessEffects object. This will eventually be deprecated.
}
\item{sh}{
  Selection coefficient for drivers with restrictions not satisfied. A
  value of 0 means there are no penalties for a driver appearing in a
  clone when its restrictions are not satisfied.

  To specify "sh=Inf" (in Diaz-Uriarte, 2015) use sh = -1.

  Only relevant if using a poset as this is included in the
  fitnessEffects object. This will eventually be deprecated.

}

\item{seed}{The seed for the C++ PRNG. You can pass a value. If you set
  it to NULL, then a seed will be generated in R and passed to C++. If
  you set it to "auto", then if you are using v.1, the behavior is the
  same as if you set it to NULL (a seed will be generated in R and
  passed to C++) but if you are using v.2, a random seed will be
  produced in C++. %% using the randutils code from M.E.O'Neill.
  If you   need reproducibility, either pass a value or set it to NULL (setting
  it to NULL will make the C++ seed reproducible if you use the same
  seed in R via \code{set.seed}). However, even using the same value of
  \code{seed} is unlikely to give the exact same results between
  platforms and compilers.  Moreover, note that the defaults for
  \code{seed} are not the same in \code{oncoSimulIndiv},
  \code{oncoSimulPop} and \code{oncoSimulSample}.

  When using oncoSimulPop, if you want reproducibility, you might want
  to, in addition to setting \code{seed = NULL}, also do
  \code{RNGkind("L'Ecuyer-CMRG")} as we use
  \code{\link{mclapply}}; look at the vignette of
  \pkg{parallel}.

}

\item{interventions}{

  This has no effect if you do not specify \code{frequencyDependentFitness = TRUE} in 
  \code{allFitnessEffects} function. Also, you must use \code{createInterventions} 
  function to create the correct type of parameter for the function \code{oncoSimulPop} ,
  \code{oncoSimulIndiv} , \code{oncoSimulSample} to process it correctly.

  Use this argument in case you want to intervene in the simulation. With interventions, 
  you can affect the total population size, or just some genotype-specific
  population. You can complicate it as much as you want, or keep it simple, 
  it is really up to you.

  Formally, interventions must be a list of lists, where each 
  "sub-list" must have the following fields:

  * ID: The identifier of the intervention, must be unique.
  * Trigger: The situation in the simulation that triggers/activates 
             the intervention.

  * What Happens: "What happens" in the simulation. Basically, once the
               trigger is satisfied, this defines how the population is
               affected by  the intervention. Please see the vignette
               for details about the differences between when
               interventions that affect a single genotype and those that affect the complete        population occur.
  
  * Periodicity: Defines the periodicity of the intervention.
  * Repetitions: Defines the maximum repetitions of each intervention in the simulation.
}

\item{userVars}{

  This has no effect if you do not specify \code{frequencyDependentFitness = TRUE} in 
  \code{allFitnessEffects} function. Also, you must use \code{createuserVars} 
  function to create the correct type of parameter for the function \code{oncoSimulPop} ,
  \code{oncoSimulIndiv} , \code{oncoSimulSample} to process it correctly.

  Use this argument in case you want to define arbitrary variables that depend on other simulation values. 
  With the yser Variables you can simulate Adaptive therapy by using this defined variables in the intervention's
  whatHappens definition, or simply get more detailed insight by defining some interesting values you desire as an output.

  Formally, userVars must be a list of lists, where each 
  "sub-list" must have the following fields:

  * Name: The name that identifies the new variable, must be unique.
  * Value: The initial numeric value of the variable.

  You must define the \code{rules} in order to determine how this variables will be modified.
}

\item{rules}{

  This has no effect if you do not specify \code{frequencyDependentFitness = TRUE} in 
  \code{allFitnessEffects} function. Also, you must use \code{createRules} 
  function to create the correct type of parameter for the function \code{oncoSimulPop} ,
  \code{oncoSimulIndiv} , \code{oncoSimulSample} to process it correctly. This also requires you to use
  \code{userVars} as these rules operate on them and will not have any effect if thses do not exist.

  Use this argument in order to determine how the defined user variebles will be modified during the simulation.
  You can use any arbitrarily complex expression depending on other simulation parameters such as T, N, Genotype
  populations or genotype rates.

  Formally, rules must be a list of lists, where each 
  "sub-list" must have the following fields:

  * ID: The identifier of the rule, must be unique.
  * Condition: The situation in the simulation that triggers/activates 
             the user variable modification.
  * Action:  The action that will take place once the condition is true. This defines wich user variables will be 
    modified and the expression that defines the new values for them.
}
                 
\item{migrationProb}{

Indicates, for each deme, the migration probability to adjacent spaces.
The default is 0.5.
}
\item{largeDistMigrationProb}{

Indicates, for each deme, the migration probability to remote spaces
(between 10 and 35 spaces for each axis). The default is 1e-6.
}
\item{maxMigrationPercentage}{

Maximum percentage of cells from each deme that would migrate. Every time a migration is performed (in every iteration, in every deme and either to adjacent or remote spaces) a set of cells will be randomly chosen between 1 (if the migration probability has been fulfilled) and this percentage. The default is 0.2.
}
\item{SpatialModel}{

Select the grid (1-D, 2-D, 3-D) for the spatial model.
The default is 3-D.
}

\item{spatialIterMax}{

Maximum number of iterations that can be reached in a simulation. The default is 20.
}

\item{SpatialPlot}{

Select if it is desired to generate a plot with the result of the spatial 
simulation or not. The default is TRUE.
}

}

\details{

  The algorithm of SpatialOncoSimul is based on the main concepts of different spatial
  models (Bartlomiej Waclaw, et al., 2015, A. Sottoriva, et al., 2015, R. Sun, et al., 
  2017, R. Sun, et al., 2021, Noble, R. et al., 2022) and the ideas set out in the 
  Master's degree final project of Alberto Parramón Castillo (2018).
  
  The \code{summary} method for the class returns lenght, class and mode of each 
  of the components as a multiple row data frame. The \code{print} method for the 
  class simply prints the summary.


}
 
\value{
  
  For \code{SpatialOncoSimul} a list, of class "SpatialOncoSimul", with the 
  following components:
  
    \item{Final_Population }{Total number of cells at the end of the simulation.}
  
  \item{Total_Demes   }{Total number of demes at the end of the simulation.}
  
  \item{Total_Iterations   }{Final number of executed iterations.}
  
  \item{List_of_Genotypes   }{Dataframe containing all genotypes and the number 
  of cells with each of them.}
  
  \item{Grids   }{List of dataframes, each of them being one grid, with its 
  genotypes, number of cells with each genotype and coordinates.}

}
\references{

  Parramón Castillo, A. (2018). Simulación de un modelo espacial de evolución 
  tumoral con R y C++. [TFM]. Universidad Complutense de Madrid.
  
  Diaz-Uriarte, R. (2015). Identifying restrictions in the order of
  accumulation of mutations during tumor progression: effects of
  passengers, evolutionary models, and sampling
  \url{http://www.biomedcentral.com/1471-2105/16/41/abstract}
  
  Bozic, I., et al., (2010). Accumulation of driver and passenger
  mutations during tumor progression. \emph{ Proceedings of the National
  Academy of Sciences of the United States of America\/}, \bold{107},
  18545--18550.
  
  R Diaz-Uriarte. (2017). OncoSimulR: genetic simulation with arbitrary epistasis 
  and mutator genes in asexual populations. Bioinformatics, 33, 1898--1899.
  https://doi.org/10.1093/bioinformatics/btx077.

  R Diaz-Uriarte and C. Vasallo. (2019). Every which way? On predicting tumor 
  evolution using cancer progression models 2019 PLoS Computational Biology
  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007246

  R Diaz-Uriarte. (2017). Cancer progression models and fitness landscapes: a 
  many-to-many relationship 2017 Bioinformatics.
  https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btx663/

  Gerstung et al., 2011. The Temporal Order of Genetic and Pathway
  Alterations in Tumorigenesis. \emph{PLoS ONE}, 6.
  
  McFarland, C.~D. et al. (2013). Impact of deleterious passenger
  mutations on cancer progression.  \emph{Proceedings of the National
  Academy of Sciences of the United States of America\/}, \bold{110}(8),
  2910--5.

  Mather, W.~H., Hasty, J., and Tsimring, L.~S. (2012). Fast stochastic
  algorithm for simulating evolutionary population dynamics.
  \emph{Bioinformatics (Oxford, England)\/}, \bold{28}(9), 1230--1238.

  

}

\note{Please, note that the meaning of the fitness effects in the
  McFarland model is not the same as in the original paper; the fitness
  coefficients are transformed to allow for a simpler fitness function
  as a product of terms. This differs with respect to v.1. See the
  OncoSimulR vignette for details.}

\author{Antonio Giráldez Trujillo, Maria González Bermejo, Mercedes Núñez Bayón.}


\examples{

######### A simple example with spatial simulation.

## Fitness Effects of the population genotypes with restrictionsin the order 
## of accumulation of mutations using a DAG.

s1 <- allFitnessEffects(
  data.frame(parent = c("Root", "Root", "i"),
             child = c("u" , "i" , "v"),
             s = c(0.1 , -0.05 , 0.25),
             sh = -1,
             typeDep = "MN"),
  epistasis = c("u:i" = -1,"u:v" = -1))

evalAllGenotypes (s1 , order = FALSE, addwt = TRUE)

## Simul the spatial tumor cell growth in 3D during 20 iterations.
Spatial_3D <- SpatialOncoSimul(fp = s1,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 500,
                      mu = 1e-4,
                      initSize = 1000,
                      keepPhylog = FALSE,
                      seed = NULL,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE, initMutant = c("i"), 
                      spatialIterMax = 20, SpatialModel = "3D")

######### Using an extended poset for pancreatic cancer from Gerstung et al.
###      (s and sh are made up for the example; only the structure
###       and names come from Gerstung et al.)

pancr <- allFitnessEffects(
    data.frame(parent = c("Root", rep("KRAS", 4), 
                   "SMAD4", "CDNK2A", 
                   "TP53", "TP53", "MLL3"),
               child = c("KRAS","SMAD4", "CDNK2A", 
                   "TP53", "MLL3",
                   rep("PXDN", 3), rep("TGFBR2", 2)),
               s = 0.1,
               sh = -0.9,
               typeDep = "MN"),
    drvNames = c("KRAS", "SMAD4", "CDNK2A", "TP53", 
                 "MLL3", "TGFBR2", "PXDN"))
evalAllGenotypes (pancr , order = FALSE, addwt = TRUE)

## Using exponential growth model
Spatial_3D <- SpatialOncoSimul(fp = pancr,
                      model = "Exp",
                      onlyCancer = FALSE,
                      finalTime = 500,
                      mu = 1e-4,
                      initSize = 1000,
                      keepPhylog = FALSE,
                      seed = NULL,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE, 
                      spatialIterMax = 10, SpatialModel = "3D")


######### Frequency dependent birth example. 

## An example with cooperation. Presence of WT favours all clones and all clones 
## have a positive effect on themselves

genofit <- data.frame(i = c(0, 1, 0, 1),
                      u = c(0, 0, 1, 1),
                      Birth = c("3 + 5*f_",
                                  "3 + 5*(f_ + f_i)",
                                  "3 + 5*(f_ + f_u)",
                                  "5 + 6*(f_ + f_i_u)"))

afe <- allFitnessEffects(genotFitness = genofit,
                         frequencyDependentBirth = TRUE)
                         
Spatial_3D <- SpatialOncoSimul(fp = afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 500,
                      mu = 1e-4,
                      initSize = 1000,
                      keepPhylog = FALSE,
                      seed = NULL,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE, 
                      spatialIterMax = 10, SpatialModel = "3D")                         

}
\keyword{misc}
\keyword{iteration}