#' Mars Constructor
#'
#'@description
#'Description here
#'
#' @param control
#'
#' @return
#' @export
#'
#' @examples
new_mars.control <- function(control) {
  structure(control,class="mars.control")}


#' Mars Helper
#'
#'@description
#'Description here
#'
#' @param control Mars control object
#'
#' @return
#' @export
#'
#' @examples
validate_mars.control <- function(control) {
  stopifnot(is.integer(control$Mmax),is.numeric(control$d),
            is.logical(control$trace))
  if(control$Mmax < 2) {
    warning("Mmax must be >= 2; Reset it to 2")
    control$Mmax <- 2}
  if(control$Mmax %% 2 > 0) {
    control$Mmax <- 2*ceiling(control$Mmax/2)
    warning("Mmax should be an even integer. Reset it to ",control$Mmax)}
  control
}


#' Mars Constructor
#'
#'@description
#'Description here
#' @param Mmax
#' @param d
#' @param trace
#'
#' @return
#' @export
#'
#' @examples
#'
mars.control <- function(Mmax=2,d=3,trace=FALSE) {
  Mmax <- as.integer(Mmax)
  control <- list(Mmax=Mmax,d=d,trace=trace)
  control <- validate_mars.control(control)
  new_mars.control(control)
}


#' Forward stepwise selection
#'
#'@description
#'Description here
#'
#' @return
#' @export
#'
#' @examples
fwd_stepwise<-function(y,x,control=mars.control()){
  N<-length(y)
  n<-ncol(x)
  B<-init_B(N,control$Mmax)
  Bfuncs<-vector(mode="list",length=control$Mmax+1)
  for (M in seq(1,control$Mmax,by=2)){
    if (control$trace){
      cat("M",M,"\n")
    }
    lof_best<-Inf
    for (m in 1:M){
      svars<-setdiff(1:n,Bfuncs[[m]][,"v"])
      if (control$trace){
        cat("M, m, svars",M,m,svars,"\n")
      }
      for (v in svars){
        tt<-split_points(x[,v],B[,m])
        for (t in tt){
          Bnew<-data.frame(B[,1:M],Btem1=B[,m]*h(x[,v],1,t),Btem2=B[,m]*h(x[,v],-1,t))
          gdat<-data.frame(y=y,Bnew)
          lof<-LOF(y~.,gdat,control)
          if (lof<lof_best){
            lof_best<-lof
            split_best<-c(m=m,v=v,t=t)
          }
        }
      }
    }
    m<-split_best["m"]
    v<-split_best["v"]
    t<-split_best["t"]
    Bfuncs[[M+1]]<-rbind(Bfuncs[[m]],c(s=-1,v,t))
    Bfuncs[[M+2]]<-rbind(Bfuncs[[m]],c(s=1,v,t))
    B[,M+1:2]<-cbind(B[,m]*h(x[,v],-1,t),B[,m]*h(x[,v],1,t))
  }
  colnames(B)<-paste0("B",(0:(ncol(B)-1)))
  return(list(y=y,B=B,Bfuncs=Bfuncs))
}


#' init_B
#'
#' @description
#' Description here
#' @param N
#' @param Mmax
#'
#' @return
#' @export
#'
#' @examples
init_B<-function(N,Mmax){
  B<-data.frame(matrix(nrow=N,ncol=Mmax+1))
  B[,1]<-1
  names(B)<-paste0("B",0:Mmax)
  return(B)
}


#' LOF
#'
#' @param formula
#' @param data
#' @param control
#'
#' @return
#' @export
#'
#' @examples
LOF<-function(formula,data,control){
  mod<-lm(formula,data)
  RSS<-sum(residuals(mod)^2)
  N<-nrow(data)
  M<-length(coef(mod))-1
  c<-sum(diag(hatvalues(mod)))+control$d*M
  out<-RSS*N/(N-c)^2
  return(out)
}


#' h
#'
#' @param x
#' @param s
#' @param t
#'
#' @return
#' @export
#'
#' @examples
h<-function(x,s,t){
  out<-pmax(0,s*(x-t))
  return(out)
}


split_points<-function(x,B){
  out<-sort(unique(x[B>0]))[-length(sort(unique(x[B>0])))]
  return(out)
}


#' Backward stepwise selection
#'
#' @description
#' description here
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#'
bwd_stepwise<-function(fwd,control){
  Mmax<-ncol(fwd$B)-1
  Jstar<-2:(Mmax+1)
  Kstar<-Jstar
  dat<-data.frame(y=fwd$y,fwd$B)
  lofstar<-LOF(y~.-1,dat,control)
  for (M in (Mmax+1):2){
    b<-Inf
    L<-Kstar
    if (control$trace){
      cat("L:",L,"\n")
    }
    for (m in L){
      K<-setdiff(L,m)
      dat<-data.frame(y=fwd$y,fwd$B[,K])
      lof<-LOF(y~.,dat,control)
      if (control$trace){
        cat("M:K:lof ",M,":",K,":",lof,"\n")
      }
      if(lof<b){
        b<-lof
        Kstar<-K
      }
      if(lof<lofstar){
        lofstar<-lof
        Jstar<-K
      }
    }
    if (control$trace){
      cat("M:Jstar:lofstar",M,":",Jstar,":",lofstar,"\n")
    }
  }
  Jstar<-c(1,Jstar)
  return(list(y=fwd$y,B=fwd$B[,Jstar],Bfuncs=fwd$Bfuncs[Jstar]))
}


#' Mars
#'
#' @description
#' Description of what mars is
#'
#' @param formula a formula
#' @param data dataset
#' @param control mars control object (can be created using mars.control)
#'
#' @return
#' Mars returns ...
#' @export
#'
#' @examples
#'
mars<-function(formula,data,control=mars.control()){
  cc<-match.call()
  mf<-model.frame(formula,data)
  y<-model.response(mf)
  mt<-attr(mf, "terms")
  x<-model.matrix(mt,mf)[,-1,drop=F]
  x_names<-colnames(x)
  control<-validate_mars.control(control)
  fwd<-fwd_stepwise(y,x,control)
  bwd<-bwd_stepwise(fwd,control)
  fit<-lm(y~.-1,data=data.frame(y=y,bwd$B))
  out<-c(list(call=cc,formula=formula,y=y,B=bwd$B,Bfuncs=bwd$Bfuncs,x_names=x_names),fit)
  class(out)<-c("mars",class(fit))
  return(out)
}
