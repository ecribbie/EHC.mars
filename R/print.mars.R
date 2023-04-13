#' print.mars
#'
#' @param mars_obj A mars class object, created using the EHC.MARS::mars() function 
#'
#' @return Outputs the call and coefficients used for creating mars object
#' @export
#'
#' @examples 
#' mars_object<-mars(formula, data, control = mars.control()) 
#' print(mars_object)
print.mars<-function(mars_obj){
  print.lm(mars_obj)
}
