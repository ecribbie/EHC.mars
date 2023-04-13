#' anova.mars
#'
#' @param mars_obj mars class object, created using EHC.MARS::mars function
#'
#' @return anova of mars object regression
#' @export
#'
#' @examples
#' mars_object<-mars(formula, data, control = mars.control()) 
#' anova(mars_object)
anova.mars<-function(mars_obj){
  class(mars_obj)<-setdiff(class(mars_obj),"mars")
  aov<-anova(mars_obj)
  class(mars_obj)<-c(class(mars_obj),"mars")
  return(aov)
}
