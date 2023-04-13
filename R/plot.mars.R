plot.mars<-function(mars_obj){
  class(mars_obj)<-setdiff(class(mars_obj),"mars")
  plot(mars_obj)
  class(mars_obj)<-c(class(mars_obj),"mars")
}
