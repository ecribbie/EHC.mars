#' mars: Multivariate Adaptive Regression Splines
#'
#' @description
#' `mars` is used for fitting Multivariate Adaptive Regression Splines
#'  (MARS) models, which are non-parametric regression models that capture complex interactions between predictor variables.
#'  
#'  Download from: https://github.com/ecribbie/EHC.mars.git
#'
#' @param formula   a formula for linear model
#' @param data      dataset with response variable (y) and explanatory variables
#' @param control   mars control object (can be created using EHC.MARS::mars.control, See Details for more information)
#'
#' @return `mars` returns an object of class `"mars"`. The functions `summary` and `anova` are used to obtain and print a summary and analysis of variance table of the results.
#'
#'  an object of class `"mars"` is a list contains the following components:
#'
#' * `call`   the matched call.

#' * `formula`  the matched formula.
#' * `y`  the response variable.
#'
#' * `B`  a matrix represents the basis functions (terms) used in the MARS model and their coefficients
#' * `Bfuncs`   a list of functions that represent the basis functions (terms) used in the MARS model
#' *  `x_names`   names of the predictor variables used
#'
#'
#'
#'
#' * `coefficients`   a named vector of coefficients
#' * `residuals`  the residuals, that is response minus fitted values.
#'
#'
#' * `fitted.values`  the fitted mean values
#' * `df.residual`  the residual degrees of freedom.
#' *  `xlevels`   (only where relevant) a record of the levels of the factors used in fitting.
#'
#'
#' * `model`  the fitted model.
#' * `terms`  the [terms] object used.
#'
#'


#'
#' @export
#' @author
#' Evan Cribbie
#'
#' Huong Thi Mai Nguyen
#'
#'
#' Christine Orcullo
#' @references Friedman, J. H. (1991). Multivariate Adaptive Regression Splines. The Annals of Statistics, 19(1), 1-67.

#' @seealso
#'  `summary.mars` for summaries and `anova.mars` for the ANOVA table
#'
#'  `plot.mars` for visualizing the fitted model, `print.mars` for ouputting call and coefficients used for creating MARS object.
#'
#' `predict.mars` for making predictions using the fitted model.
#' @section Details:
#' The `mars` function uses a formula-based approach to specify MARS models.
#' The formula should follow the standard R formula syntax, with 'response' on the
#' left-hand side and 'terms' on the right-hand side. The 'terms' specification can include
#' main effects, interactions, and higher-order terms.
#'
#' The `mars` function fits a MARS model to the data, which involves automatically selecting
#' basis functions (e.g., hinge functions) and their locations in the predictor space, and
#' estimating coefficients for each basis function using a forward-backward algorithm. The resulting
#' MARS model captures complex interactions between predictor variables, making it a powerful tool for modeling nonlinear relationships in data.
#'
#' `mars.control` is used to create a MARS control object with the following arguments:
#'
#' * `Mmax` max number of splits for mars function (default is 2)
#' * `d` d value used in GCV criterion for mars (default is 3)
#' * `trace` Logical value for if additional output should be printed when run (default is False)
#'
#' @examples
#'#https://www.kaggle.com/datasets/joyshil0599/mlb-hitting-and-pitching-stats-through-the-years
#' dat<-read.csv("data/pitcher_data.csv")
#' dat$AVG<-as.numeric(dat$AVG)
#' dat1<-na.omit(dat)[1:500,]
#' mars_obj<-mars(Win~Games.played+Innings.pitched+Hit.Batsmen+base.on.balls+WHIP+AVG,dat1,mars.control(Mmax=8,d=3,trace=T))
#' summary(mars_obj)
#' plot(mars_obj)
#' pred_dat<-na.omit(dat)[500:600,]
#' predict(mars_obj,pred_dat)
#'
#' #https://www.kaggle.com/datasets/neuromusic/avocado-prices
#' dat2<-read.csv("data/avocado.csv")[1:400,]
#' mars_obj2<-mars(AveragePrice~Total.Volume+X4046+X4225+X4770+Small.Bags+Large.Bags+XLarge.Bags+year,dat2)
#' print(mars_obj2)
#' anova(mars_obj2)
#'
#' #https://www.kaggle.com/datasets/grubenm/austin-weather
#' dat3<-read.csv("data/austin_weather.csv")[1:400,]
#' dat3$PrecipitationSumInches<-as.numeric(dat3$PrecipitationSumInches)
#' dat3<-na.omit(dat3)
#' dat3<-dat3[-305,]
#' mars_obj3<-mars(PrecipitationSumInches~HumidityAvgPercent+TempAvgF+SeaLevelPressureAvgInches+WindAvgMPH,dat3)
#' summary(mars_obj3)
mars<-function(formula,data,control=mars.control()){
  cc<-match.call()
  mf<-model.frame(formula,data)
  y<-model.response(mf)
  mt<-attr(mf, "terms")
  x_names<-colnames(x)
  control<-validate_mars.control(control)
  fwd<-fwd_stepwise(y,x,control)
  bwd<-bwd_stepwise(fwd,control)
  fit<-lm(y~.-1,data=data.frame(y=y,bwd$B))
  out<-c(list(call=cc,formula=formula,y=y,B=bwd$B,Bfuncs=bwd$Bfuncs,x_names=x_names),fit)
  class(out)<-c("mars",class(fit))
  return(out)
}


#' Mars Control
#'
#' @param Mmax max number of splits for mars function
#' @param d d value used in GCV criterion for mars
#' @param trace Logical value for if additional output should be printed when run
#'
#' @return mars control object
mars.control <- function(Mmax=2,d=3,trace=FALSE) {
  Mmax <- as.integer(Mmax)
  control <- list(Mmax=Mmax,d=d,trace=trace)
  control <- validate_mars.control(control)
  new_mars.control(control)
}

#' Mars Validator
#'
#' @param control Mars control object created with EHC.MARS::mars.control function
#'
#' @return mars control object
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
#' @param control a mars control object created with EHC.MARS::mars.control function
#'
#' @return mars control object
new_mars.control <- function(control) {
  structure(control,class="mars.control")}



#' Forward stepwise selection
#'
#' @param y response data
#' @param x explanatory data
#' @param control mars control object created with EHC.MARS::mars.control function
#'
#' @return list of y,B, and Bfuncs objects from forward selection for mars

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
        cat("M:m:svars= ",M,":",m,":",svars,"\n")
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

#' Backward stepwise selection
#'
#' @param fwd output from forward selection
#' @param control mars control object
#'
#' @return list of y, B, Bfuncs outputs from backward selection
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

#' init_B
#'
#' @param N number of rows for data matrix
#' @param Mmax number of columns not including intercept
#'
#' @return initialized matrix
init_B<-function(N,Mmax){
  B<-data.frame(matrix(nrow=N,ncol=Mmax+1))
  B[,1]<-1
  names(B)<-paste0("B",0:Mmax)
  return(B)
}

#' LOF
#'
#' @param formula formula used in linear model
#' @param data dataset
#' @param control mars control object, created with EHC.MARS::mars.control function
#'
#' @return GCV criterion value
LOF<-function(formula,data,control){
  mod<-lm(formula,data)
  RSS<-sum(residuals(mod)^2)
  N<-nrow(data)
  M<-length(coef(mod))-1
  c<-sum(diag(hatvalues(mod)))+control$d*M
  out<-RSS*N/(N-c)^2
  return(out)
}

#' hinge function
#'
#' @param x x value(s)
#' @param s side (+/- 1)
#' @param t split value
#'
#' @return hinge function ouput value
h<-function(x,s,t){
  out<-pmax(0,s*(x-t))
  return(out)
}

#' split_points
#'
#' @param x x data
#' @param B B matrix
#'
#' @return possible points to split on
split_points<-function(x,B){
  out<-sort(unique(x[B>0]))[-length(sort(unique(x[B>0])))]
  return(out)
}

