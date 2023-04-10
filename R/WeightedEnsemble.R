
#' @title Weighted Ensemble for Hybrid Model
#'
#' @param df Data set (training result) with first column as observed value
#' @param Method Method of optimization
#' @param test_data Test result
#' @param forecast Forecast result
#'
#' @import stats metaheuristicOpt
#' @return
#' \itemize{
#'   \item Weights: Optimized weight
#'   \item Optimized_Result: Optimized result
#'   }
#' @export
#'
#' @examples
#' \donttest{
#' y1<-rnorm(100,mean=100,sd=50)
#' y2<- rnorm(100,mean=100,sd=50)
#' y3<- rnorm(100,mean=100,sd=50)
#' y4<-rnorm(100,mean=100,sd=50)
#' y<-rnorm(100,mean=100,sd=50)
#' data<-cbind(y,y1,y2,y3,y4)
#' OptiSemble<-WeightedEnsemble(df=data)
#' }
#' @references
#' J. S. Armstrong. Combining forecasts: The end of the beginning or the beginning of the end? International Journal of Forecasting, 5(4):585–588, 1989.

WeightedEnsemble<-function(df,Method="PSO", test_data=NULL, forecast=NULL){  train_sel<-df
  Optimization<-Method
  numVar <- ncol(train_sel)-1

  rangeVar <- matrix(c(0,1), nrow=2)
  input_data <- train_sel

  if(ncol(train_sel)==2){
    message("Number of selected model is 1")
  } else if(ncol(train_sel)==3){
    sphere <- function(x){

      error_fn <- sum((input_data[,1] - x[1]*input_data[, 2]
                       -x[2]* input_data[,3])^2)

      return(error_fn)
    }
  } else if(ncol(train_sel)==4){
    sphere <- function(x){
      error_fn <- sum((input_data[,1] - x[1]*input_data[, 2]
                       -x[2]* input_data[,3]-x[3]* input_data[, 4])^2)
      return(error_fn)
    }
  } else if(ncol(train_sel)==5){
    sphere <- function(x){
      error_fn <- sum((input_data[,1] - x[1]*input_data[, 2]
                       -x[2]* input_data[,3]-x[3]* input_data[, 4]-x[4]* input_data[,5])^2)
      return(error_fn)
    }
  } else if(ncol(train_sel)==6){
    sphere <- function(x){
      error_fn <- sum((input_data[,1] - x[1]*input_data[, 2]
                       -x[2]* input_data[,3]-x[3]* input_data[, 4]-x[4]* input_data[,5]-x[5]* input_data[,6])^2)
      return(error_fn)
    }
  } else if(ncol(train_sel)==7){
    sphere <- function(x){
      error_fn <- sum((input_data[,1] - x[1]*input_data[, 2]
                       -x[2]* input_data[,3]-x[3]* input_data[, 4]-x[4]* input_data[,5]-x[5]* input_data[,6]
                       -x[6]* input_data[,7])^2)
      return(error_fn)
    }
  } else if(ncol(train_sel)==8){
    sphere <- function(x){
      error_fn <- sum((input_data[,1] - x[1]*input_data[, 2]
                       -x[2]* input_data[,3]-x[3]* input_data[, 4]-x[4]* input_data[,5]-x[5]* input_data[,6]
                       -x[6]* input_data[,7]-x[7]* input_data[,8])^2)
      return(error_fn)
    }
  } else if(ncol(train_sel)==9){
    sphere <- function(x){
      error_fn <- sum((input_data[,1] - x[1]*input_data[, 2]
                       -x[2]* input_data[,3]-x[3]* input_data[, 4]-x[4]* input_data[,5]-x[5]* input_data[,6]
                       -x[6]* input_data[,7]-x[8]* input_data[,9])^2)
      return(error_fn)
    }
  } else if(ncol(train_sel)==10){
    sphere <- function(x){
      error_fn <- sum((input_data[,1] - x[1]*input_data[, 2]
                       -x[2]* input_data[,3]-x[3]* input_data[, 4]-x[4]* input_data[,5]-x[5]* input_data[,6]
                       -x[6]* input_data[,7]-x[8]* input_data[,9]-x[9]* input_data[,10])^2)
      return(error_fn)
    }
  } else if(ncol(train_sel)==11){
    sphere <- function(x){
      error_fn <- sum((input_data[,1] - x[1]*input_data[, 2]
                       -x[2]* input_data[,3]-x[3]* input_data[, 4]-x[4]* input_data[,5]-x[5]* input_data[,6]
                       -x[6]* input_data[,7]-x[8]* input_data[,9]-x[9]* input_data[,10]
                       -x[10]* input_data[,11])^2)
      return(error_fn)
    }
  }else if(ncol(train_sel)==12){
    sphere <- function(x){
      error_fn <- sum((input_data[,1] - x[1]*input_data[, 2]
                       -x[2]* input_data[,3]-x[3]* input_data[, 4]-x[4]* input_data[,5]-x[5]* input_data[,6]
                       -x[6]* input_data[,7]-x[8]* input_data[,9]-x[9]* input_data[,10]
                       -x[10]* input_data[,11]-x[11]* input_data[,12])^2)
      return(error_fn)
    }
  }else{
    message("Number of selected model is more than 11")
  }

  if(Optimization=="ABC"){
    opt <- ABC(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)

  } else if (Optimization=="ALO"){
    opt <- ALO(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="BA"){
    opt <- BA(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
              maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="BHO"){
    opt <- BHO(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="CLONALG"){
    opt <- CLONALG(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
                   maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="CS"){
    opt <- CS(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
              maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  }else if (Optimization=="CSO"){
    opt <- CSO(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="DA"){
    opt <- DA(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
              maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="DE"){
    opt <- DE(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
              maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="FFA"){
    opt <- FFA(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="GA"){
    opt <- GA(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
              maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="GBS"){
    opt <- GBS(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="GOA"){
    opt <- GOA(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="GWO"){
    opt <- GWO(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="HS"){
    opt <- HS(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
              maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="KH"){
    opt <- KH(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
              maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="MFO"){
    opt <- MFO(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="PSO"){
    opt <- PSO(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="SCA"){
    opt <- SCA(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="SFL"){
    opt <- SFL(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else if (Optimization=="WOA"){
    opt <- WOA(sphere, optimType="MIN", numVar, numPopulation=nrow(train_sel),
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  } else {
    message("Please select a valid optimization technique")
  }

    final_train<-NULL
    for (q in c(1:nrow(train_sel))) {
      row_pred<-sum(train_sel[q,-1]*opt)
      final_train<-rbind(final_train,row_pred)
    }
    colnames(final_train)<-"Ensemble"

  return(list(Weights=opt, Optimized_Result=final_train))
}

