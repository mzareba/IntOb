# A script for parameters optimization v0.4
# Authors: 
# Adam Pacławski: adam.paclawski@uj.edu.pl
# Aleksander Mendyk: mfmendyk@cyf-kr.edu.pl; aleksander.mendyk@uj.edu.pl
# Jakub Szlęk: j.szlek@uj.edu.pl
# License: LGPLv3
###########################################################################################################################################
##Libraries load
require(compiler)
require(stringr)
##################################################### For user modification ###############################################################
max_loop<-10                                                   #cv no
max_supra_loop<-1

skel_plik<-paste("./25in_IKr_iv_padelall_pIC_no")
skel_plik1<-paste("./t-25in_IKr_iv_padelall_pIC_no")
# skel_outfile<-paste("output")

##Optimization parameters                                                 
maxit_SANN<-10000
maxit_optimx<-5000
optim_rel_tol=1e-20 

opti_trace<-TRUE #save drive space if FALSE

##Optimization methods
use_nloptr<-FALSE
use_GA<-TRUE
use_rgenoud<-FALSE
use_gensa<-FALSE
use_SANN<-FALSE
use_NM<-FALSE
max_iter_rgenoud<-500
max_iter_gensa<-5000
maxit_NM<-500000
maxit_nloptr<-1000
maxit_ga<-1000

##Load equation in particular formalt: param[i] - parameters for optimization, Ini wektors used in equation whe i is number of column in dataset 
equaionToOptim<-readLines("equation.txt")

############################################################# Functions ###################################################################
##RMSE function
RMSE1<-function(matrix, parameters, equat){
    res<-Inf
    C<-parameters
    try(for(i in 1:(dim(matrix)[2]-1)) {
    assign(paste("In", i, sep=""), as.double(matrix[,i]))
    out_RMSE<-as.double(matrix[,dim(matrix)[2]])
                                  },TRUE)

    try(y<-eval(parse(text=equat)),TRUE)
    try(res<-sqrt(mean((y-out_RMSE)^2)),TRUE)
    return(res)
}

##Function for minimization
funct1<-function(parameters, equat){
    C<-parameters
    y<-as.numeric(eval(parse(text=equat)))
    res<-mean((y-out)^2)
    if (is.na(res)){res<-Inf}
    return(res)
                      }
##Function for minimization using by GA
funct2<-function(x, y){
  return(pmax(-funct1(x, y),-1e16))
                      }
##T-res function
tRes1<-function(matrix, parameters, equat){
    C<-parameters
    try(for(i in 1:(dim(matrix)[2]-1)) {
    assign(paste("In", i, sep=""), as.double(matrix[,i]))
    observed<-as.double(matrix[,dim(matrix)[2]])
                                  },TRUE)

    try(predicted<-eval(parse(text=equat)),TRUE)
    try(res<-cbind(observed, predicted), TRUE)
    colnames(res)<-c("Observed", "Predicted")
    return(res)

                            }
                      

##Optimize functions
RMSE<-cmpfun(RMSE1)
funct<-cmpfun(funct1)
functGA<-cmpfun(funct2)
tRes<-cmpfun(tRes1)
############################################################# Preprocesing ################################################################
##Correct input equation if you use and 'ln' symbol for logarithm 
equationBuf<-gsub(" ", "", equaionToOptim)
equation<-gsub("ln", "log", equationBuf)
N_params<-str_count(equation, "C")

##Prepare parameters and variables
RMSE_val<-vector(length=max_loop)
all_params <- matrix(data=NA,nrow=max_loop,ncol=N_params)
RMSE_ind<-vector(length=max_loop)
best_all_params<-matrix(data=NA,nrow=max_loop,ncol=N_params)
best_RMSE_val<-vector(length=max_loop)
best_RMSE_ind<-vector(length=max_loop)
best_RMSE_total<-1000000000000

##Display equation
print("Equation for optimization")
print(equation)

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

if (use_SANN){
## Optim with optim(SANN)
print(equation)
print("params on start SANN")
print(paramFunct)

try(fit0 <- optim(
	paramFunct,
	equat=equation,
 	fn=funct,
 	method="SANN",
 	control=list(trace=opti_trace,maxit=maxit_SANN))
 	,TRUE)

print("summary")
try(summary(fit0),TRUE)


try(if(fit0$convergence==0){
  par_optim_SANN<-as.vector(as.matrix(fit0$par))
  RMSE_ucz_SANN<-RMSE(matryca, par_optim_SANN, equation)
  best_error<-RMSE_ucz_SANN
  RMSE_test_SANN<-RMSE(matryca1, par_optim_SANN, equation)
  paramFunct<-par_optim_SANN
  RMSE_val[lk_loop]<-RMSE_test_SANN
                            }else{
  par_optim_SANN<-paramFunct
  RMSE_ucz_SANN<-Inf
  RMSE_test_SANN<-Inf
     },TRUE)

}     
          
for_domain <- matrix(data=NA,nrow=length(paramFunct),ncol=2)

for(i in 1:length(paramFunct)) {
  for_domain[i,1]<--100*max(abs(paramFunct))
  for_domain[i,2]<-100*max(abs(paramFunct))
}

if (use_nloptr){

require(nloptr)

print("Running nloptr")
  fit0 <- nloptr(x0=paramFunct,
    	eval_f=funct,    	
    	lb=for_domain[,1],
    	ub=for_domain[,2],
    	equat=equation,
        opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=optim_rel_tol,maxeval=maxit_nloptr,print_level=1,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=optim_rel_tol))
  )
  
  paramFunct<-fit0$solution

  print("FINAL RESULTS nloptr")
  print(paramFunct)

}
########################################
if (use_GA){

require(GA)

print("Running GA")
  fit0 <- ga(
      suggestions=paramFunct,
    	fitness=function(x) functGA(x, equation),
      # fitness=functGA,
      # equat=equation,
      type = "real-valued",
    	min=for_domain[,1],
    	max=for_domain[,2],
      maxiter=maxit_ga
  )
  
  paramFunct<-fit0@solution

  print("FINAL RESULTS GA")
  print(paramFunct)

}
########################################
if (use_gensa){

require(GenSA)

print("Running GenSA")
  fit0 <- GenSA(par=paramFunct,
    	#errorMeasure=errorOpt,  
    	#ind=inpIndex, 
    	#model=modelOpt, 
    	#tar=tarOpt,    
    	fn=funct,
    	equat=equation,
    	lower=for_domain[,1],
    	upper=for_domain[,2],
  	control=list(smooth=FALSE,maxit=max_iter_gensa,verbose=TRUE,nb.stop.improvement=max_iter_gensa)
  )
  
  paramFunct<-fit0$par

  print("FINAL RESULTS GenSA")
  print(paramFunct)

}


if (use_rgenoud){

print("Running rgenoud")

require(rgenoud)

fit2 <- genoud(
	fn=funct,	
	nvars=length(paramFunct),
	starting.values=paramFunct,
	equat=equation,
	max.generations=max_iter_rgenoud,
	Domains=for_domain
)
 
 paramFunct<-fit2$par
 print("Current estimates of the parameters - rgenoud1")
 print(paramFunct)

}


if (use_NM){
## Optim with optim(SANN)
print(equation)
print("params on start SANN")
print(paramFunct)

try(fit0 <- optim(
	paramFunct,
	equat=equation,
 	fn=funct,
 	method="Nelder-Mead",
 	control=list(trace=opti_trace,maxit=maxit_NM))
 	,TRUE)

print("summary")
try(summary(fit0),TRUE)

  paramFunct<-fit0$par

  print("FINAL RESULTS NM")
  print(paramFunct)

}     
     
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
for(lk_loop in 1:max_loop){
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



