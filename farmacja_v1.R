####################### libraries
require(compiler)
require(stringr)
library(GA)
library(foreach)

# ####################### variables - remove unnecessary ones
# max_loop <- 1 # Define how many sets use to learn, max 10.
# max_supra_loop <- 1

# ENERGY <- 50
# ENERGY_EXCHANGE <- 5

# ## Optimization parameters
# maxit_SANN <- 100
# maxit_optimx <- 5000
# optim_rel_tol = 1e-20

# opti_trace <- TRUE # save drive space if FALSE

# ## Optimization methods
# use_GA <- TRUE
# maxit_ga <- 10
# minusInfForGaFitness <- -1e16
# plusInfForGaFitness <- 1e4
# maxFitnessToStopGA <- 1e3 # Has to be bigger than PlusInfForGaFitness

# ####################### file names

# # fileName <- "equationList.txt"
# fileName <- "equation.txt"
# inputTrain <- "inputTrain.txt"
# inputTest <- "inputTest.txt"
source("parameters.R")

####################### functions
## Function used by GA, definiton from source

SELECTION <- function(object, k = 3, ...) {
# (unbiased) Tournament selection 
  sel <- rep(NA, object@popSize)
  for(i in 1:object@popSize)
     { s <- sample(1:object@popSize, size = k)
       sel[i] <- s[which.max(object@fitness[s])]
     }
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

CROSSOVER <- function(object, parents, ...) {
# Blend crossover
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  a <- 0.5
  # a <- exp(-pi*iter/max(iter)) # annealing factor
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  for(i in 1:n)
     { x <- sort(parents[,i])
       xl <- max(x[1] - a*(x[2]-x[1]), object@min[i])
       xu <- min(x[2] + a*(x[2]-x[1]), object@max[i])
       children[,i] <- runif(2, xl, xu) 
     }
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

MUTATION <- function(object, parent, ...) {
# Random mutation around the solution
  mutate <- parent <- as.vector(object@population[parent,])
  dempeningFactor <- 1 - object@iter/object@maxiter
  direction <- sample(c(-1,1),1)
  value <- (object@max - object@min)*0.67
  mutate <- parent + direction*value*dempeningFactor
  outside <- (mutate < object@min | mutate > object@max)
  for(j in which(outside))
     { mutate[j] <- runif(1, object@min[j], object@max[j]) }
  return(mutate)
}

## RMSE function
RMSE1 <- function(matrix, parameters, equat) {
  res <- Inf
  C <- parameters
  try (for (i in 1:(dim(matrix)[2] - 1)) {
    assign(paste("In", i, sep = ""), as.double(matrix[, i]))
    out_RMSE <- as.double(matrix[, dim(matrix)[2]])
  }, TRUE)
  
  try (y <- eval(parse(text = equat)), TRUE)
  try (res <- sqrt(mean((y - out_RMSE) ^ 2)), TRUE)
  return (res)
  
}

##RMSE function
RMSE1<-function(matrix, parameters, equat){
    res<-Inf
    C<-parameters
    try (for (i in 1:(dim(matrix)[2]-1)) {
      assign(paste("In", i, sep=""), as.double(matrix[,i]))
      out_RMSE<-as.double(matrix[,dim(matrix)[2]])
    },TRUE)

    try (y<-eval(parse(text=equat)),TRUE)
    try (res<-sqrt(mean((y-out_RMSE)^2)),TRUE)
    return(res)
}

## Function for minimization
fitness <- function(parameters, equat) {
  C <- parameters
  y <- as.numeric(eval(parse(text = equat)))
  res <- mean((y - out) ^ 2)
  if (is.na(res)) {
    res <- Inf
  }
  return (res)
}

## Negative reflection of function for minimization (fitness)
fitnessReflection <- function(x, y) {
  return(pmin( pmax( -fitness(x, y), minusInfForGaFitness), plusInfForGaFitness))
}

## T - res function
tRes1 <- function(matrix, parameters, equat) {
  C <- parameters
  try (for (i in 1:(dim(matrix)[2] - 1)) {
    assign(paste("In", i, sep = ""), as.double(matrix[, i]))
    observed <- as.double(matrix[, dim(matrix)[2]])
  }, TRUE)
  
  try (predicted <- eval(parse(text = equat)), TRUE)
  try (res <- cbind(observed, predicted), TRUE)
  colnames(res) <- c("Observed", "Predicted")
  return (res)
}

MEETING <- function(population, n_params) {
  fitnessIndex = n_params + 2
  energyIndex = n_params + 1
  out = matrix(ncol = n_params + 2, nrow = 0)
  while (nrow(population) >= 2) {
    indexes = sample(1:nrow(population), 2)
    object1 = population[indexes[1],]
    object2 = population[indexes[2],]
    population = population[-c(indexes[1], indexes[2]),]
    if (object1[fitnessIndex] > object2[fitnessIndex]) {
      
      if (object1[energyIndex] <= ENERGY_EXCHANGE) { #not enough energy 0 < energy < ENERGY_EXCHANGE, take what's left and remove
        object2[energyIndex] = object2[energyIndex] + object1[energyIndex]
        out = rbind(out, object2)
      } else {
        object2[energyIndex] = object2[energyIndex] + ENERGY_EXCHANGE
        object1[energyIndex] = object1[energyIndex] - ENERGY_EXCHANGE
        out = rbind(out, object1, object2)
      }
      
    } else { #if energy is equal we still need to make the exchange
      if (object2[energyIndex] <= ENERGY_EXCHANGE) { #not enough energy 0 < energy < ENERGY_EXCHANGE, take what's left and remove
        object1[energyIndex] = object1[energyIndex] + object2[energyIndex]
        out = rbind(out, object1)
      } else {
        object1[energyIndex] = object1[energyIndex] + ENERGY_EXCHANGE
        object2[energyIndex] = object2[energyIndex] - ENERGY_EXCHANGE
        out = rbind(out, object1, object2)
      }
    }
    #if theres one object left and it didnt met any other object out it in out
    if (is.integer(nrow(population))) {
      out = rbind(out, population)
      break()
    }
  }
  return(out)
}

BREEDING <- function(population, n_params) {
  fitnessIndex = n_params + 2
  energyIndex = n_params + 1
  parents_population = matrix(ncol = n_params + 2, nrow = 0)
  rest_population = matrix(ncol = n_params + 2, nrow = 0)
  while(nrow(population) > 0) {
    index = sample(1:nrow(population), 1)
    parent = population[index,]
    population = population[-index,,drop=FALSE]
    if(parent[energyIndex] > ENERGY_BREEDING) {
      parents_population = rbind(parents_population, parent)
    } else {
      rest_population = rbind(rest_population, parent)
    }
  }

  iter = 1
  children = matrix(nrow = 0, ncol = n_params + 2)
  while(iter < nrow(parents_population)) {
    index1 = iter
    index2 = iter + 1
    iter = iter + 2
    if (runif(1) > PROBABILITY_BREEDING)
      next
    parents = rbind(parents_population[index,], parents_population[index2,])
    child <- matrix(as.double(NA), nrow = 1, ncol = n_params + 2)
    a = 0.1
    for(i in 1:n_params) {
      x <- sort(parents[,i])
      xl <- x[1] - a*(x[2]-x[1])
      xu <- x[2] + a*(x[2]-x[1])
      child[,i] <- runif(1, xl, xu) 
    }
    energy_parent1 = as.integer(parents_population[index1,energyIndex]/2)
    parents_population[index1, energyIndex] = parents_population[index1, energyIndex] - energy_parent1
    energy_parent2 = as.integer(parents_population[index2,energyIndex]/2)
    parents_population[index2, energyIndex] = parents_population[index2, energyIndex] - energy_parent1
    child[1,energyIndex] = energy_parent1 + energy_parent2
    children = rbind(children, child)
  }

  population = rbind(rest_population, parents_population, children)
  population = RECALCULATE_FITNESS(population, n_params)
  return(population)
}

ELIMINATION <- function(population, n_params) {
  return(population)
}

RECALCULATE_FITNESS <- function(population, n_params) {
  for (i in 1:nrow(population)) {
    population[i,n_params+2] <- fitness(parameters = head(population[i,], -2), equat = equation)
  }
  return(population)
}

##Optimize functions
RMSE <- cmpfun(RMSE1)
funct <- cmpfun(fitnessReflection)
tRes <- cmpfun(tRes1)

####################### reading equations

con <- file(fileName, open = "r")
equationsRaw <- readLines(con)
close(con)

####################### preprocessing
skel_plik<-paste("./25in_IKr_iv_padelall_pIC_no")
skel_plik1<-paste("./t-25in_IKr_iv_padelall_pIC_no")

# equationBuf <- gsub(" ", "", equationsRaw[2])
# equation <- gsub("ln", "log", equationBuf)

equaionToOptim<-readLines("equation.txt")
equationBuf<-gsub(" ", "", equaionToOptim)
equation<-gsub("ln", "log", equationBuf)
print("Equation for optimization")
print(equation)

N_params <- str_count(equation, "C")

## Prepare parameters and variables
RMSE_val        <- vector(length = max_loop)
all_params      <- matrix(data = NA, nrow = max_loop, ncol = N_params)
RMSE_ind        <- vector(length = max_loop)
best_all_params <-matrix(data = NA, nrow = max_loop, ncol = N_params)
best_RMSE_val   <- vector(length = max_loop)
best_RMSE_ind   <- vector(length = max_loop)
best_RMSE_total <- 1000000000000

###############################
########## FIRST ONE ##########
############ START ############
###############################
# ####################### reading inputs

# trainMatrix <-
#   read.csv(inputTrain,
#            header = FALSE,
#            sep = "\t",
#            colClasses = "numeric")
# testMatrix <-
#   read.csv(inputTest,
#            header = FALSE,
#            sep = "\t",
#            colClasses = "numeric")

# for (i in 1:(dim(trainMatrix)[2] - 1)) {
#   assign(paste("In", i, sep = ""), as.double(trainMatrix[, i]))
#   out <- as.double(trainMatrix[, dim(trainMatrix)[2]])
# }

# population <- 100
# prmsList <- matrix(data = NA, nrow = population, ncol = N_params)

# for (i in 1:population) {
#   prmsList[i, ] <- rnorm(N_params) / 10
# }

# fitness <- function(prms) {
#   return (funct(prms, equation))
# }

# populationFunction <- function() {
#   return (prmsList)
# }

# #initialFit<-funct(params, equation) #fitness?
# #cat("Initial fitness = ",initialFit,"\n")

# GA <-
#   ga(
#     fitness = fitness,
#     popSize = population,
#     type = "real-valued",
#     min = c(-1,-1),
#     max = c(1, 1)
#   )
# summary(GA)
###############################
########## FIRST ONE ##########
############# END #############
###############################


##Supra optimization loop
for (lk_supra_loop in 1: max_supra_loop) {
  RMSE_total<-0

  ##Main optimization function
  for(lk_loop in 1:max_loop) {

    ##Read learn, test files
    plik <- paste(skel_plik,lk_loop,".txt",sep="")
    plik1 <- paste(skel_plik1,lk_loop,".txt",sep="")
    # outfile<-paste(skel_outfile,"_",lk_loop,".RData",sep="")

    cat("Iteration no = ",lk_loop,"\n")
    cat("Training file = ",plik,"\n")
    cat("Test file = ",plik1,"\n")
    # cat("Output file = ",outfile,"\n")

    matryca <- read.csv(plik,header=FALSE,sep="\t", colClasses="numeric")
    matryca1 <- read.csv(plik1,header=FALSE,sep="\t", colClasses="numeric")
        
    for(i in 1:(dim(matryca)[2]-1)) {
        assign(paste("In", i, sep=""), as.double(matryca[,i]))
        out<-as.double(matryca[,dim(matryca)[2]])
    }                         
    #paramFunct <-vector(length=N_params, "numeric")
    paramFunct<-rnorm(N_params)/10
    print("paramFunct")    
    print(paramFunct)
    best_error<-100000000
    cat("Iteration no = ",lk_loop,"\n")
    cat("best_error INIT = ",best_error,"\n")

    para1.backup<-paramFunct
    X<-matryca

    print("Check init values")

    preliminary_output<-funct(paramFunct, equation)                      
    cat("Preliminary output = ",preliminary_output,"\n")   
              
    for_domain <- matrix(data=NA,nrow=length(paramFunct),ncol=2)

    for(i in 1:length(paramFunct)) {
      for_domain[i,1]<--100*max(abs(paramFunct))
      for_domain[i,2]<-100*max(abs(paramFunct))
    }
    
    N_ROWS <- 50
    
    initialPopulation <- matrix(nrow = N_ROWS, ncol = N_params + 2)
    
    for (i in 1:N_ROWS) {
      initialPopulation[i,] <- c(rnorm(N_params), 50, 0) #energy, fitness
      initialPopulation[i,N_params+2] <- fitness(parameters = head(initialPopulation[i,], -2), equat = equation)
    }
    
    ########################################
    population <- initialPopulation
    
    if (use_GA){
      for (i in 1:1000) {
        population <- MEETING(population, N_params)
        population <- BREEDING(population, N_params)
        population <- ELIMINATION(population, N_params)
      }
    }
    
    print(population)
    paramFunct <- head(population[1,], -2)
    
    ########################################

    ## Optim with optim(BFGS)
    print("params on start BFGS")
    print(paramFunct)

    par_optim_NM<-paramFunct

    try(fit1 <- optim(
        paramFunct,
        equat=equation,
        fn=funct,    
        method="BFGS",
        control=list(trace=opti_trace,maxit=maxit_optimx)
    ),TRUE)

    print("FINAL RESULTS OPTIM(BFGS)")
    try(print(fit1$par),TRUE)

    print("WHOLE object")
    try(print(fit1),TRUE)

    try(par_optim_NM<-fit1$par,TRUE)

    RMSE_ucz_NM<-Inf
    RMSE_test_NM<-Inf
    print("learn_error")
    try(RMSE_ucz_NM<-RMSE(matryca,par_optim_NM, equation))
    print(RMSE_ucz_NM)
    print("test_error")
    try(RMSE_test_NM<-RMSE(matryca1,par_optim_NM, equation))
    print(RMSE_test_NM)

    cat("Iteration no = ",lk_loop,"\n")
    print("Final params")
    try(print(par_optim_NM),TRUE)
    try(all_params[lk_loop,]<-par_optim_NM,TRUE)
    print(" ")

    ##This conditions should help with situation where test RMSE is NA what couse problems and stop current job.
    if(is.na(RMSE_test_NM)){
    RMSE_test_NM<-Inf
    }
    rmserror<-RMSE_test_NM
    RMSE_ind[lk_loop]<-rmserror

    cat("Iteration no = ",lk_loop,"\n")
    cat("RMSE_test = ",rmserror,"\n")

    # try(save(fit1,file=outfile),TRUE)

    print("-------------------------------------")

  }         
  ##End of optimization function
  ##----------------------------------------------------


  ## End of loop loop lk_loop

  print(" ")
  print("SUMMARY")
  print(" ")
  for(lk_loop in 1:max_loop) {
    cat("RMSE no",lk_loop," = ",RMSE_ind[lk_loop],"\n")
    RMSE_total<-RMSE_total+RMSE_ind[lk_loop]
  }
  RMSE_total<-RMSE_total/max_loop
  print("-------------------------------------")
  cat("avRMSE = ",RMSE_total,"\n")
  print(" ")
  print("All parameters :")
  print(all_params)
  print("-------------------------------------")
  print("<<<<<<<<<< END OF LOOP >>>>>>>>>>>>")
  print("-------------------------------------")

  if (RMSE_total<best_RMSE_total){
    best_all_params<-all_params
    best_RMSE_val<-RMSE_val
    best_RMSE_ind<-RMSE_ind
    best_RMSE_total<-RMSE_total  
  }
}
#end of supra_optimization_loop

    all_params<-best_all_params
    RMSE_val<-best_RMSE_val
    RMSE_ind<-best_RMSE_ind
    RMSE_total<-best_RMSE_total

    
print(" ")
print("OVERALL SUMMARY")
print(" ")
for(lk_loop in 1:max_loop){
  cat("RMSE no",lk_loop," = ",RMSE_ind[lk_loop],"\n")
}
print("-------------------------------------")
cat("avRMSE = ",RMSE_total,"\n")
print(" ")
print("All parameters :")
print(all_params)
print("-------------------------------------")
print("<<<<<<<<<< END OF LOOP >>>>>>>>>>>>")
print("-------------------------------------")


##Save equation and params into text files
sink(paste("optimizedEquation.txt", sep=""))
cat("Orig equation\n\n")
cat(equaionToOptim, sep=" ")
cat("\n\n")
cat("Equation with C \n\n")
cat(equation, "\n\n")
cat("Parameters for equation after optimization \n\n")

for(i in 1:dim(all_params)[1]) {
  cat(i, ": ", sep="")
  cat(all_params[i,], sep=";")
  cat("\n")
}
cat("\n")
cat("RMSE for each data sets in k-cross validataion method\n\n")

for(i in 1:max_loop) {
  cat("RMSE_", i, ": ")
  cat(RMSE_ind[i], sep="; ")
  cat("\n")
}

cat("\n")
cat("Mean_RMSE: ", RMSE_total, "\n", sep="")

##Make Observed and Predicted table
obsPredFileName<-"Results.txt"
cat(paste("Observed\tPredicted\n", sep=""), file=obsPredFileName, append=FALSE)  

for(i in 1:max_loop) {
  fileName <- paste(skel_plik1,i,".txt",sep="")
  testData <- read.csv(fileName,header=FALSE,sep="\t", colClasses="numeric")
  predObsMat<-tRes(
      matrix=testData,
      parameters=all_params[i,],
      equat=equation)

  write.table(predObsMat, sep="\t", file=obsPredFileName, col.names=FALSE, row.names=FALSE, append=TRUE)
}