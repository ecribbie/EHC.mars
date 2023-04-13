#' plot.mars
#'
#' @param mars_obj A mars class object, created using the EHC.MARS::mars() function 
#'
#' @return Gives set of plots for mars object function
#' @export
#'
#' @examples 
#' mars_object<-mars(formula, data, control = mars.control()) 
#' plot(mars_object)
plot.mars<-function(mars_obj){
  class(mars_obj)<-setdiff(class(mars_obj),"mars")
  plot(mars_obj)
  class(mars_obj)<-c(class(mars_obj),"mars")
}
