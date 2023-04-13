#' Title
#'
#' @param X
#' @param Bfuncs
#'
#' @return
#' @export
#'
#' @examples
make_B<-function(X,Bfuncs){
  B=matrix(1,nrow=nrow(X),ncol=length(Bfuncs))
  for (i in 2:ncol(B)){
    for (j in 1:nrow(Bfuncs[[i]]))
      B[,i]<-B[,i]*h(X[,Bfuncs[[i]][j,][2]],Bfuncs[[i]][j,][1],Bfuncs[[i]][j,][3])
  }
  colnames(B)<-paste0("B",(0:(ncol(B)-1)))
  return(B)
}

#' Title
#'
#' @param object
#' @param newdata
#'
#' @return
#' @export
#'
#' @examples
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
