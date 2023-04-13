#' Title
#'
#' @param mars_obj
#'
#' @return
#' @export
#'
#' @examples
anova.mars<-function(mars_obj){
  class(mars_obj)<-setdiff(class(mars_obj),"mars")
  aov<-anova(mars_obj)
  class(mars_obj)<-c(class(mars_obj),"mars")
  return(aov)
}
