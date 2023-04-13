summary.mars<-function(mars_obj){
  sm<-summary.lm(mars_obj)
  for (i in 1:length(rownames(sm$coefficients))){
    cat(paste0(c("\n",rownames(sm$coefficients)[i],":\n"),collapse=''))
    if (i==1){
      cat("Intercept")
      next
    }
    for (j in 1:(length(mars_obj$Bfuncs[[i]])/3)){
      r<-mars_obj$Bfuncs[[i]][j,]
      cat(c("Split on variable: ",r[2],", at value: ",signif(r[3],4),", with direction: ",r[1],"\n"),sep='')
    }
  }
  print(sm)
}
