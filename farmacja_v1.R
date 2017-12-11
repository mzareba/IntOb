####################### libraries
require(compiler)
require(stringr)
library(GA)
library(foreach)

####################### variables - remove unnecessary ones
max_loop<-10 # cv no
max_supra_loop<-1

## Optimization parameters
maxit_SANN<-10000
maxit_optimx<-5000
optim_rel_tol = 1e-20

opti_trace<-TRUE # save drive space if FALSE

## Optimization methods
use_nloptr<-TRUE
use_rgenoud<-FALSE
use_gensa<-FALSE
use_SANN<-FALSE
use_NM<-FALSE
max_iter_rgenoud<-500
max_iter_gensa<-5000
maxit_NM<-500000
maxit_nloptr<-10000

####################### functions 

## RMSE function
RMSE1<-function(matrix, parameters, equat) {
  res<-Inf
  C<-parameters
  try (for (i in 1: (dim(matrix)[2] - 1)) {
    assign(paste("In", i, sep = ""), as.double(matrix[, i]))
    out_RMSE<-as.double(matrix[, dim(matrix)[2]])
  }, TRUE)
  
  try (y<-eval(parse(text = equat)), TRUE)
  try (res<-sqrt(mean((y - out_RMSE) ^ 2)), TRUE)
  return (res)
  
}

## Function for minimization
funct1<-function(parameters, equat) {
  C<-parameters
  y<-as.numeric(eval(parse(text = equat)))
  res<-mean((y - out) ^ 2)
  if (is.na(res)) {
    res<-Inf
  }
  return (res)
}

## T - res function
tRes1<-function(matrix, parameters, equat) {
  C<-parameters
  try (for (i in 1: (dim(matrix)[2] - 1)) {
    assign(paste("In", i, sep = ""), as.double(matrix[, i]))
    observed<-as.double(matrix[, dim(matrix)[2]])
  }, TRUE)
  
  try (predicted<-eval(parse(text = equat)), TRUE)
  try (res<-cbind(observed, predicted), TRUE)
  colnames(res)<-c("Observed", "Predicted")
  return (res)
}

####################### reading equations

fileName<-"equationList.txt"
con<-file(fileName, open = "r")
equationsRaw<-readLines(con)
close(con)

####################### preprocessing


equationBuf<-gsub(" ", "", equationsRaw)
equation<-gsub("ln", "log", equationBuf)
N_params<-str_count(equation, "C")

## Prepare parameters and variables
RMSE_val<-vector(length = max_loop)
all_params<-matrix(data = NA, nrow = max_loop, ncol = N_params)
RMSE_ind<-vector(length = max_loop)
best_all_params<-matrix(data = NA, nrow = max_loop, ncol = N_params)
best_RMSE_val<-vector(length = max_loop)
best_RMSE_ind<-vector(length = max_loop)
best_RMSE_total<-1000000000000

print("Works")