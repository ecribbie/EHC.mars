#' make_B
#' 
#' @description Function used in predict.mars to create initial matrix
#'
#' @param X data matrix from mars
#' @param Bfuncs Bfuncs data from mars object
#'
#' @return Initialized matrix for predict.mars function
#' @export
#'
make_B<-function(X,Bfuncs){
  B=matrix(1,nrow=nrow(X),ncol=length(Bfuncs))
  for (i in 2:ncol(B)){
    for (j in 1:nrow(Bfuncs[[i]]))
      B[,i]<-B[,i]*h(X[,Bfuncs[[i]][j,][2]],Bfuncs[[i]][j,][1],Bfuncs[[i]][j,][3])
  }
  colnames(B)<-paste0("B",(0:(ncol(B)-1)))
  return(B)
}

#' predict.mars
#'
#' @description function used to predict response value using mars object on set of data
#'
#' @param object mars class object created with EHC.MARS::mars function
#' @param newdata dataset containing data to predict on (optional, if none provided training data for mars object will be used)
#'
#' @return predictions using mars object on new data set
#' @export
#'
#' @examples
#' mars_object<-mars(formula, data, control = mars.control()) 
#' predict(mars_object,new_data)
predict.mars <- function(object,newdata) {
  if(missing(newdata) || is.null(newdata)) {
    B <- as.matrix(object$B)
  }
  else {
    tt <- terms(object$formula,data=newdata)
    tt <- delete.response(tt)
    mf <- model.frame(tt,newdata)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)[,-1] # remove intercept
    B <- make_B(X,object$Bfuncs)
  }
  beta <- object$coefficients
  drop(B %*% beta)
}
