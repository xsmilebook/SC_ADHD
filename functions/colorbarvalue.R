colorbarvalues <- function(inputvalue, midprob){
  middle <- quantile(inputvalue, probs=midprob, na.rm=T)
  maxup<- max(inputvalue, na.rm = T)
  minup<- min(inputvalue, na.rm = T)
  values1<-scales::rescale(seq(0.3,1,length.out=100), to=c(middle+abs(maxup/100),maxup))
  values2<-scales::rescale(seq(0.3,1,length.out=100), to=c(minup,middle))
  values<-c(values2,values1)
  values <- scales::rescale(values, to=c(0,1))
  return(values)
}