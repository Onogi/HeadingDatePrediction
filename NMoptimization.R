NMoptimization <- function(Emergingdates, Headingdates, Tempdash, Photodash, int.a.vec, int.b.vec, int.g.vec){
  
  #Emergingdates is a vector with length = number of environments includes the emerging dates at each environment
  #Headingdates is a vector with length = number of environments includes the heading dates at each environment
  #Tempdash is a matrix with all days x environments includes the adjusted temps.
  #Photodash is a matrix with all days x environments includes the adjusted photoperiods.
  
  stopifnot(length(Emergingdates)==length(Headingdates))
  stopifnot(ncol(Tempdash)==ncol(Photodash))
  stopifnot(length(Emergingdates)==ncol(Tempdash))
  Nenv <- length(Emergingdates)
  
  #for standardize
  unif.a<-5
  unif.b<-5
  unif.g<-55
  
  #missing records in Tempdash and Photodash
  #Missing <- NA

  #define the photoperiod sensitivity length
  dvs1<-function(sigma){0.145+0.005*sigma}
  dvs2<-function(sigma){0.345+0.005*sigma}
  
  fdvs<-function(x, Temp, Photo, emergingdate, headingdate){
    #Temp and Photo are the vectors taken from Tempdash and Photodash respectively. 
    
    alpha <-x[1]*unif.a
    beta  <-x[2]*unif.b
    sigma <-x[3]*unif.g
    dvs<-0
    Maxdates <- max(which(!is.na(Temp)))
    
    i <- emergingdate - 1
    
    if(is.na(headingdate)){
      fdvs.out <- 0
    }else{
      repeat{
        i <- i+1
        t <- Temp[i]
        p <- Photo[i]
        
        ##dvrŽZo
        if((dvs1(sigma) < dvs) && (dvs < dvs2(sigma))){
          dvs<-dvs+(t^alpha)*(p^beta)/sigma
        } else {
          dvs<-dvs+(t^alpha)/sigma
        }
#cat(Maxdates, t,alpha,sigma, dvs,"\n")        
        if(i==Maxdates){break}
        if(dvs>=1){break}
      }#end.repeat
      fdvs.out <- abs(headingdate - i)
    }
    fdvs.out
  }  
  
  fdvs.op<-function(x){
    optim.sum <- 0
    for(env in 1:Nenv){
      optim.sum <- optim.sum + fdvs(x, Tempdash[,env], Photodash[,env], Emergingdates[env], Headingdates[env])
    }
    if( x[1] < 0 ){
      optim.sum <- optim.sum - 10000 * x[1]
    }
    if( x[2] < 0 ){
      optim.sum <- optim.sum - 10000 * x[2]
    }
    optim.sum
  }
  
  fd <- function(a,b,g){
    optim.out<-optim(c(a/unif.a,b/unif.b,g/unif.g),fdvs.op,method="Nelder-Mead")
    instantframe<-data.frame(
      INT.a = a,
      INT.b = b,
      INT.G = g,
      OPT.a = optim.out$par[1]*unif.a,
      OPT.b = optim.out$par[2]*unif.b,
      OPT.G = optim.out$par[3]*unif.g,
      OPT.value = optim.out$value,
      OPT.counts = optim.out$counts[1]
    )
    instantframe
  }  
  
  output<-data.frame()
  for(int.a in int.a.vec){
    for(int.b in int.b.vec){
      for(int.g in int.g.vec){
        output <- rbind(output,	fd(int.a,int.b,int.g))
      } #g
    } #b
  } #a
    
  return(output[which(output$OPT.value==min(output$OPT.value)),])
}