####################### variables - remove unnecessary ones
max_loop <- 1 # Define how many sets use to learn, max 10.
max_supra_loop <- 1

ENERGY <- 50
ENERGY_EXCHANGE <- 5
ENERGY_BREEDING <- 80

PROBABILITY_BREEDING <- 0.6

## Optimization parameters
maxit_SANN <- 100
maxit_optimx <- 5000
optim_rel_tol = 1e-20

opti_trace <- TRUE # save drive space if FALSE

## Optimization methods
use_GA <- TRUE
maxit_ga <- 10
minusInfForGaFitness <- -1e16
plusInfForGaFitness <- 1e4
maxFitnessToStopGA <- 1e3 # Has to be bigger than PlusInfForGaFitness

####################### file names

# fileName <- "equationList.txt"
fileName <- "equation.txt"
inputTrain <- "inputTrain.txt"
inputTest <- "inputTest.txt"