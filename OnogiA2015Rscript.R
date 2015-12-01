#Prediction results are stored in the following objects###########################################################################
#175line and 8env denote the number of lines and environments used for training, 
#i.e. corresponding to LOLO and LOEO, respectively.

  #Genomic prediction using extended Bayesian lasso (EBL)
    #Full data	: Direct.EBL.fitting
    #LOLO             : Direct.EBL.175line
    #LOELO            : Direct.EBL.8env.175line

  #C-Nel
    #Full data	: O.fitting
    #LOEO             : O.8env

  #T-Nel
    #Full data    : O.fitting.EBL
    #LOEO             : O.8env.EBL
    #LOLO             : O.175line.EBL
    #LOELO            : O.8env.175line.EBL

  #C-Bay
    #Full data    : Ilog.fitting
    #LOEO             : Ilog.8env

  #T-Bay
    #Full data    : Ilog.fitting.EBL
    #LOEO             : Ilog.8env.EBL
    #LOLO             : Ilog.175line.EBL
    #LOELO            : Ilog.8env.175line.EBL

  #IM
    #Full data    : EBL.fitting
    #LOEO             : EBL.8env
    #LOLO             : EBL.175line
    #LOELO            : EBL.8env.175line

  #C-Bay and IM were conducted using the C scripts in "FilesForBayesianInference"
  #Then the results were retrieved using this R script.
  #Other calculations were conducted using this R script.

  #The package "VIGoR" is required to perform EBL.
  #VIGoR is available at CRAN or at https://github.com/Onogi/VIGoR

#Input information#################################################################################################################
Environments <- c("Tsukuba2007","Fukuoka2008","HaNoi2008","Ishigaki2008","Ishikawa2008",
                  "ThaiNguyen2008","Tsukuba2008E","Tsukuba2008L","Tsukuba2009")

Heading <- as.matrix(read.table("BIL.headingdate.txt"))
Emergence <- as.matrix(read.table("BIL.emergencedate.txt"))
Heading[Heading==-9] <- NA
Emergence[Emergence==-9] <- NA
Ne<-ncol(Heading) #number of environments
Nline<-nrow(Heading) #number of lines
Temp <- as.matrix(read.table("BIL.dailytemp.txt"))
Photo <- as.matrix(read.table("BIL.photoperiod.txt"))
Temp[which(Temp==-9)]<-NA
Photo[which(Photo==-9)]<-NA
colnames(Heading)<-colnames(Emergence)<-colnames(Temp)<-colnames(Photo)<-Environments

#Fixed parameters of the DVR model (Nakagawa et al. 2005) 
Tc <- 42; To <- 30; Tb <- 8
Pc <- 24; Po <- 10; Pb <- 0

TcTo <- Tc-To; ToTb <- To-Tb; Tratio <- TcTo/ToTb
PcPo <- Pc-Po; PoPb <- Po-Pb; Pratio <- PcPo/PoPb

#Transform daily mean temperature and photoperiod.
Temp2 <- Temp; Photo2 <- Photo
for(env in 1:Ne){
  for(i in 1:length(Temp[,env])){
    if(!is.na(Temp[i,env])){
      if(Temp[i,env]>=Tb&Temp[i,env]<=Tc){
        Temp2[i,env] <- (Temp[i,env]-Tb)/ToTb*((Tc-Temp[i,env])/TcTo)^Tratio
      }else{
        Temp2[i,env] <- 0
      }
      if(Photo[i,env]>=Pb&Photo[i,env]<=Pc){
        Photo2[i,env] <- (Photo[i,env]-Pb)/PoPb*((Pc-Photo[i,env])/PcPo)^Pratio
      }else{
        Photo2[i,env] <- 1
      }
    }
  }
}

#Marker genotypes
Geno<-as.matrix(read.table("KKBIL_geno.txt"))

#Plotting DTH distributions#######################################################################################################
range(Heading-Emergence+1,na.rm=T) #53, 167
Mar<-c(2,2,0.5,0.5)
Cex<-0.5
Width<-4.2
Height<-3
Xlim<-c(50,170)
for(i in 1:9){
  tiff(paste("DTH_",Environments[i],".tif",sep=""),unit="cm",width=Width,height=Height,res=600)
  par(mar=Mar)
  par(cex=Cex)
  hist((Heading-Emergence+1)[,i],main="", xlab="",ylab="",breaks=seq(50,170,2))
  dev.off()    
}

#bi-plots of DTH between environments
Heading.plot<-data.frame(Heading[,c(1,7,8,9,5,2,4,6,3)])
colnames(Heading.plot)<-c("T2007","T2008E","T2008L","T2009","IK","FO","IG","TN","HN")
tiff("DTH.biplot.tif",unit="cm",width=17,height=17,res=400)
par(cex=0.5)
par(mar=c(2,1,1,0.5))
plot(Heading.plot,main="")
dev.off()

#Local outlier factor
library(DMwR)
#impute missing records with line means
ImputeWithLineMean<-(Heading-Emergence+1)
for(line in 1:176){
  if(any(is.na(ImputeWithLineMean[line,]))){
    ImputeWithLineMean[line,is.na(ImputeWithLineMean[line,])]<-mean(ImputeWithLineMean[line,!is.na(ImputeWithLineMean[line,])])
  }
}
for(K in seq(10,170,10)){
  cat(K,max(lofactor(ImputeWithLineMean,k=K)),"\n")
}
#K LOF
#10 1.609725 
#20 1.532554 
#30 1.768315 
#40 1.898824 
#50 1.99278 
#60 2.033349 
#70 2.062421 
#80 2.110563 
#90 2.131876 
#100 2.120339 
#110 2.07462 
#120 2.013855 
#130 1.928487 
#140 1.813758 
#150 1.667141 
#160 1.437872 
#170 1.164944 

#impute missing records with environment means
ImputeWithEnvMean<-(Heading-Emergence+1)
for(env in 1:9){
  if(any(is.na(ImputeWithEnvMean[,env]))){
    ImputeWithEnvMean[is.na(ImputeWithEnvMean[,env]),env]<-mean(ImputeWithEnvMean[!is.na(ImputeWithEnvMean[,env]),env])
  }
}
for(K in seq(10,170,10)){
  cat(K,max(lofactor(ImputeWithEnvMean,k=K)),"\n")
}
#K LOF
#10 1.508084 
#20 1.649134 
#30 1.870834 
#40 2.061766 
#50 2.196042 
#60 2.264309 
#70 2.324046 
#80 2.333093 maximum value
#90 2.317962 
#100 2.286599 
#110 2.211828 
#120 2.122092 
#130 2.018004 
#140 1.885568 
#150 1.718859 
#160 1.443116 
#170 1.168238

#Map information of markers########################################################################################################
MAP<-t(as.matrix(read.table("KKBIL_geno_map.txt")))
temp<-MAP[MAP[,1]==1,2]
for(i in 2:12){
  temp<-c(temp,MAP[MAP[,1]==i,2]+temp[length(temp)]+1)
}
MAP<-cbind(MAP,temp)
colnames(MAP)<-c("Chr","Pos","CumPos")
MAP<-data.frame(MAP)

#Markers nearest to heading date genes
Hdnearest<-c(59,88,92,108,112)#correspond to Hd6, Hd3a, Hd1, Hd2, and Hd5, respectively


#function for prediction of heading date using the DVR model#########################################################################
Predict <- function(Temp2,Photo2,g,a,b,emergence){
  
  if(is.na(emergence)){start<-1}else{start<-emergence}
  
  DVS1 <- 0.145*g+0.005*g*g
  DVS2 <- 0.345*g+0.005*g*g
  
  dvs<-0
  Maxdate<-max(which(!is.na(Temp2)))
  for(i in start:Maxdate){
    if(dvs>=DVS1&&dvs<=DVS2){
      dvs<-dvs+Temp2[i]^a*Photo2[i]^b         
    }else{
      dvs<-dvs+Temp2[i]^a        
    }
    if(dvs>g) break
  }
  i
}


#function to calculate MSE###############################################################################################################
rmse<-function(dif){sqrt(mean(dif^2,na.rm=T))}

#Analysis with C-Nel#####################################################################################################################
source("NMoptimization.R") #conduct Nelder-Mead optimization of the DVR model
#initial values
avec<-c(0.01,5,10)
bvec<-c(0.01,5,10)
gvec<-c(40,55,80)

#Full data
O.fitting<-as.list(numeric(2)) #optimized parameters are stored in O.fitting[[1]] and fitted values are in O.fitting[[2]]
O.fitting[[1]]<-matrix(0,nr=Nline,nc=6)
for(line in 1:Nline){
  cat(line,"\n")
  temp<-NMoptimization(Emergence[line,],Heading[line,],Temp2,Photo2,avec,bvec,gvec)
  O.fitting[[1]][line,]<-c(mean(temp$OPT.G),sd(temp$OPT.G),mean(temp$OPT.a),sd(temp$OPT.a),mean(temp$OPT.b),sd(temp$OPT.b))
}
O.fitting[[1]]<-data.frame(O.fitting[[1]])
colnames(O.fitting[[1]])<-c("g","gsd","a","asd","b","bsd")
O.fitting[[2]]<-matrix(0,nr=Nline,nc=2*Ne)
for(line in 1:Nline){
  for(env in 1:Ne){
    O.fitting[[2]][line,2*env-1]<-Predict(Temp2[,env],Photo2[,env],O.fitting[[1]]$g[line],O.fitting[[1]]$a[line],O.fitting[[1]]$b[line],Emergence[line,env])
  }
}


#LOEO
O.8env <- as.list(numeric(Ne+1))
O.8env[[Ne+1]]<-matrix(0,nc=2*Ne,nr=Nline)
#Optimized parameters and predicted (fitted) values at each fold are stored in O.8env[[1]] to O.8env[[Ne]].
#Only predicted values are extracted and stored in O.8env[[Ne+1]]

for(foldenv in 1:Ne){
  O.8env[[foldenv]]<-as.list(numeric(2))
  O.8env[[foldenv]][[1]]<-matrix(0,nr=Nline,nc=6)
  for(line in 1:Nline){
    cat("foldenv",foldenv,"line",line,"\n")
    temp<-NMoptimization(Emergence[line,-foldenv],Heading[line,-foldenv],Temp2[,-foldenv],Photo2[,-foldenv],avec,bvec,gvec)
    O.8env[[foldenv]][[1]][line,]<-c(mean(temp$OPT.G),sd(temp$OPT.G),mean(temp$OPT.a),sd(temp$OPT.a),mean(temp$OPT.b),sd(temp$OPT.b))
  }
  O.8env[[foldenv]][[1]]<-data.frame(O.8env[[foldenv]][[1]])
  colnames(O.8env[[foldenv]][[1]])<-c("g","gsd","a","asd","b","bsd")
  O.8env[[foldenv]][[2]]<-matrix(0,nr=Nline,nc=2*Ne)
  for(line in 1:Nline){
    for(env in 1:Ne){
      O.8env[[foldenv]][[2]][line,2*env-1]<-
        Predict(Temp2[,env],Photo2[,env],O.8env[[foldenv]][[1]]$g[line],O.8env[[foldenv]][[1]]$a[line],O.8env[[foldenv]][[1]]$b[line],Emergence[line,env])
    }
  }
  O.8env[[Ne+1]][,2*foldenv-1]<-O.8env[[foldenv]][[2]][,2*foldenv-1]
}

#Analysis with T-Nel######################################################################################################################
#Full data
library(VIGoR)
O.fitting.EBL<-as.list(numeric(2))
O.fitting.EBL[[1]]<-matrix(0,nr=176,nc=6)
O.fitting.EBL[[2]]<-matrix(0,nr=176,nc=2*Ne)
O.fitting.coef<-matrix(0,nr=162,nc=6)
for(para in c(1,3,5)){
  v<-vigor(O.fitting[[1]][,para],Geno,"EBL",c(0.1,0.1,1,0.1),Threshold=9,Printinfo=F)
  O.fitting.EBL[[1]][,para]<-colSums(t(Geno)*v$Beta)+v$Alpha
  O.fitting.coef[,c(para,para+1)]<-cbind(v$Beta,v$Sd.beta)
}
O.fitting.coef<-data.frame(O.fitting.coef)
colnames(O.fitting.coef)<-c("g","gsd","a","asd","b","bsd")
O.fitting.EBL[[1]][O.fitting.EBL[[1]]<0]<-0
for(env in 1:Ne){
  for(line in 1:Nline){
    O.fitting.EBL[[2]][line,2*env-1]<-
      Predict(Temp2[,env],Photo2[,env],O.fitting.EBL[[1]][line,1],O.fitting.EBL[[1]][line,3],O.fitting.EBL[[1]][line,5],Emergence[line,env])    
  } 
}

#Computational time
t<-proc.time()
for(line in 1:Nline){
  temp<-NMoptimization(Emergence[line,],Heading[line,],Temp2,Photo2,avec,bvec,gvec)
  v<-c(mean(temp$OPT.G),sd(temp$OPT.G),mean(temp$OPT.a),sd(temp$OPT.a),mean(temp$OPT.b),sd(temp$OPT.b))
}
TNelfull<-proc.time()-t
for(para in c(1,3,5)){
  v<-vigor(O.fitting[[1]][,para],Geno,"EBL",c(0.1,0.1,1,0.1),Threshold=9,Printinfo=F)
  temp<-colSums(t(Geno)*v$Beta)+v$Alpha
}
for(env in 1:Ne){
  for(line in 1:Nline){
    v<-Predict(Temp2[,env],Photo2[,env],O.fitting.EBL[[1]][line,1],O.fitting.EBL[[1]][line,3],O.fitting.EBL[[1]][line,5],Emergence[line,env])    
  } 
}


#plot marker effects
#g
tiff("./Figures/O.fitting.coef.g.tif",unit="cm",height=3,width=17,res=600)
par(mar=c(1,2,1,0.5))
par(cex=0.8)
plot(MAP$CumPos,abs(O.fitting.coef$g),col=c("gray60","gray10")[MAP$Chr%%2+1],pch=8,main="",xlab="",ylab="",xaxt="n")
for(i in 1:length(Hdnearest)) lines(c(MAP$CumPos[Hdnearest[i]],MAP$CumPos[Hdnearest[i]]),c(0,100),lty=2,lwd=1)
dev.off()
#a
tiff("./Figures/O.fitting.coef.a.tif",unit="cm",height=3,width=17,res=600)
par(mar=c(1,2,1,0.5))
par(cex=0.8)
plot(MAP$CumPos,abs(O.fitting.coef$a),col=c("gray60","gray10")[MAP$Chr%%2+1],pch=8,main="",xlab="",ylab="",xaxt="n")
for(i in 1:length(Hdnearest)) lines(c(MAP$CumPos[Hdnearest[i]],MAP$CumPos[Hdnearest[i]]),c(0,100),lty=2,lwd=1)
dev.off()
#b
tiff("./Figures/O.fitting.coef.b.tif",unit="cm",height=3,width=17,res=600)
par(mar=c(1,2,1,0.5))
par(cex=0.8)
plot(MAP$CumPos,abs(O.fitting.coef$b),col=c("gray60","gray10")[MAP$Chr%%2+1],pch=8,main="",xlab="",ylab="",xaxt="n")
for(i in 1:length(Hdnearest)) lines(c(MAP$CumPos[Hdnearest[i]],MAP$CumPos[Hdnearest[i]]),c(0,100),lty=2,lwd=1)
dev.off()

#LOEO
O.8env.EBL<-as.list(numeric(Ne+1))
for(env in 1:Ne){
  O.8env.EBL[[env]]<-matrix(0,nr=176,nc=6)
}
O.8env.EBL[[Ne+1]]<-matrix(0,nr=176,nc=Ne*2)
for(foldenv in 1:Ne){
  for(para in c(1,3,5)){
    cat(foldenv,para,"\n")
    temp<-O.8env[[foldenv]][[1]][,para]
    v<-vigor(temp,Geno,"EBL",c(0.1,0.1,1,0.1),Threshold=9,Printinfo=F)
    O.8env.EBL[[foldenv]][,para]<-colSums(t(Geno)*v$Beta)+v$Alpha
  }
  O.8env.EBL[[foldenv]][O.8env.EBL[[foldenv]]<0]<-0
  for(line in 1:Nline){
    O.8env.EBL[[Ne+1]][line,2*foldenv-1]<-
      Predict(Temp2[,foldenv],Photo2[,foldenv],O.8env.EBL[[foldenv]][line,1],O.8env.EBL[[foldenv]][line,3],O.8env.EBL[[foldenv]][line,5],Emergence[line,foldenv])    
  }
}

#LOEO using markers except for R1684
Geno.woR1684<-Geno
Geno.woR1684<-Geno.woR1684[,-162]
O.8env.EBL.woR1684<-as.list(numeric(Ne+1))
for(env in 1:Ne){
  O.8env.EBL.woR1684[[env]]<-matrix(0,nr=176,nc=6)
}
O.8env.EBL.woR1684[[Ne+1]]<-matrix(0,nr=176,nc=Ne*2)
for(foldenv in 1:Ne){
  for(para in c(1,3,5)){
    cat(foldenv,para,"\n")
    temp<-O.8env[[foldenv]][[1]][,para]
    v<-vigor(temp,Geno.woR1684,"EBL",c(0.1,0.1,1,0.1),Threshold=9,Printinfo=F)
    O.8env.EBL.woR1684[[foldenv]][,para]<-colSums(t(Geno.woR1684)*v$Beta)+v$Alpha
  }
  O.8env.EBL.woR1684[[foldenv]][O.8env.EBL.woR1684[[foldenv]]<0]<-0
  for(line in 1:Nline){
    O.8env.EBL.woR1684[[Ne+1]][line,2*foldenv-1]<-
      Predict(Temp2[,foldenv],Photo2[,foldenv],O.8env.EBL.woR1684[[foldenv]][line,1],O.8env.EBL.woR1684[[foldenv]][line,3],O.8env.EBL.woR1684[[foldenv]][line,5],Emergence[line,foldenv])    
  }
}

#LOLO
O.175line.EBL <- as.list(numeric(2))
O.175line.EBL[[1]]<-matrix(0,nc=6,nr=Nline)
O.175line.EBL[[2]]<-matrix(0,nc=2*Ne,nr=Nline)
names(O.175line.EBL)<-c("abg","DH")
for(para in c(1,3,5)){
  v<-vigor(O.fitting[[1]][,para],Geno,"EBL",c(0.1,0.1,1,0.1),"cv",-1,Threshold=9,Printinfo=F)
  O.175line.EBL[[1]][,para]<-v$Prediction$Yhat
}
for(foldline in 1:Nline){
  g<-O.175line.EBL[[1]][foldline,1];if(g<0)g<-0  
  a<-O.175line.EBL[[1]][foldline,3];if(a<0)a<-0
  b<-O.175line.EBL[[1]][foldline,5];if(b<0)b<-0
  for(env in 1:Ne){
    O.175line.EBL[[2]][foldline,2*env-1]<-Predict(Temp2[,env],Photo2[,env],g,a,b,Emergence[foldline,env])   
  }        
}


#LOELO
O.8env.175line.EBL<-as.list(numeric(Ne+1))
O.8env.175line.EBL[[Ne+1]]<-matrix(0,nc=2*Ne,nr=Nline)
for(foldenv in 1:Ne){
  O.8env.175line.EBL[[foldenv]]<-matrix(0,nr=Nline,nc=6)
  for(para in c(1,3,5)){ 
    v<-vigor(O.8env[[foldenv]][[1]][,para],Geno,"EBL",c(0.1,0.1,1,0.1),"cv",-1,Threshold=9,Printinfo=F)
    O.8env.175line.EBL[[foldenv]][,para]<-v$Prediction$Yhat
  }
}
for(foldenv in 1:Ne){
  for(foldline in 1:Nline){
    cat(foldenv,foldline,"\n")
    g<-O.8env.175line.EBL[[foldenv]][foldline,1];if(g<0)g<-0
    a<-O.8env.175line.EBL[[foldenv]][foldline,3];if(a<0)a<-0
    b<-O.8env.175line.EBL[[foldenv]][foldline,5];if(b<0)b<-0
    O.8env.175line.EBL[[Ne+1]][foldline,2*foldenv-1]<-Predict(Temp2[,foldenv],Photo2[,foldenv],g,a,b,Emergence[foldline,foldenv])
  }
}

#Analysis with GP####################################################################################################################################
DTH<-Heading-Emergence+1
colnames(DTH)<-Environments
DTHsqdif<-matrix(0,Ne,Ne) #mean square differences of DTH between environments 
for(env in 1:Ne){
  for(i in 1:Ne){
    temp<-DTH[,i]-DTH[,env]
    temp<-length(temp[!is.na(temp)])
    DTHsqdif[env,i]<-sum((DTH[,i]-DTH[,env])^2,na.rm=T)/temp
  }
}
Pairedenv<-numeric(Ne) #training environments for LOELO
for(env in 1:Ne){temp<-DTHsqdif; diag(temp)<-1000000; Pairedenv[env]<-which.min(temp[env,])}

round(sqrt(DTHsqdif[c(1,7,8,9,5,2,4,6,3),c(1,7,8,9,5,2,4,6,3)]),1)
#[1,]  0.0 28.6 18.6 10.3 11.4 11.1 24.7 26.7 24.7
#[2,] 28.6  0.0 46.8 18.9 17.9 38.8 52.5 54.4 52.5
#[3,] 18.6 46.8  0.0 28.4 29.3  8.9  8.5 10.6  9.1
#[4,] 10.3 18.9 28.4  0.0  2.7 20.4 34.6 36.6 34.6
#[5,] 11.4 17.9 29.3  2.7  0.0 21.2 35.6 37.5 35.6
#[6,] 11.1 38.8  8.9 20.4 21.2  0.0 16.0 18.1 16.2
#[7,] 24.7 52.5  8.5 34.6 35.6 16.0  0.0  3.9  3.7
#[8,] 26.7 54.4 10.6 36.6 37.5 18.1  3.9  0.0  4.2
#[9,] 24.7 52.5  9.1 34.6 35.6 16.2  3.7  4.2  0.0

round(cor(DTH[,c(1,7,8,9,5,2,4,6,3)],use="pairwise.complete.obs"),2)
#                 Tsukuba2007 Tsukuba2008E Tsukuba2008L Tsukuba2009 Ishikawa2008 Fukuoka2008 Ishigaki2008 ThaiNguyen2008 HaNoi2008
#Tsukuba2007           1.00         0.96         0.87        0.97         0.97        0.91         0.45           0.30      0.33
#Tsukuba2008E          0.96         1.00         0.80        0.98         0.97        0.85         0.32           0.18      0.23
#Tsukuba2008L          0.87         0.80         1.00        0.83         0.85        0.96         0.71           0.58      0.59
#Tsukuba2009           0.97         0.98         0.83        1.00         0.98        0.88         0.35           0.20      0.24
#Ishikawa2008          0.97         0.97         0.85        0.98         1.00        0.91         0.40           0.25      0.29
#Fukuoka2008           0.91         0.85         0.96        0.88         0.91        1.00         0.64           0.48      0.52
#Ishigaki2008          0.45         0.32         0.71        0.35         0.40        0.64         1.00           0.91      0.88
#ThaiNguyen2008        0.30         0.18         0.58        0.20         0.25        0.48         0.91           1.00      0.86
#HaNoi2008             0.33         0.23         0.59        0.24         0.29        0.52         0.88           0.86      1.00


#Full data
Direct.EBL.fitting<-Emergence
for(env in 1:Ne){
  v<-vigor(DTH[,env],Geno,"EBL",c(0.1,0.1,1,0.1),Threshold=9,Printinfo=F)
  Direct.EBL.fitting[,env]<-Direct.EBL.fitting[,env]+colSums(t(Geno)*v$Beta)+v$Alpha-1 
}

#LOLO
#train the model at each environment
Direct.EBL.175line<-Emergence
for(env in 1:Ne){
  v<-vigor(DTH[,env],Geno,"EBL",c(0.1,0.1,1,0.1),"cv",-1,Threshold=9,Printinfo=F)
  Direct.EBL.175line[,env]<-Direct.EBL.175line[,env]+v$Prediction$Yhat-1
}
#Use GBLUP
library(rrBLUP)
Direct.GBLUP.175line<-Emergence
for(env in 1:Ne){
  for(foldline in 1:Nline){
    v<-kinship.BLUP(DTH[-foldline,env],Geno[-foldline,],Geno[foldline,,drop=FALSE])
    Direct.GBLUP.175line[foldline,env]<-Direct.GBLUP.175line[foldline,env]+v$g.pred+v$beta-1
  }
}


#LOELO
#train the model at the environment closed to the target with regards to mean squared difference
Direct.EBL.8env.175line<-Emergence
for(foldenv in 1:Ne){
  v<-vigor(DTH[,Pairedenv[foldenv]],Geno,"EBL",c(0.1,0.1,1,0.1),"cv",-1,Threshold=9,Printinfo=F)
  Direct.EBL.8env.175line[,foldenv]<-Direct.EBL.8env.175line[,foldenv]+v$Prediction$Yhat-1
}
#GBLUP
Direct.GBLUP.8env.175line<-Emergence
for(foldenv in 1:Ne){
  for(foldline in 1:Nline){
    v<-kinship.BLUP(DTH[-foldline,Pairedenv[foldenv]],Geno[-foldline,],Geno[foldline,,drop=FALSE])
    Direct.GBLUP.8env.175line[foldline,foldenv]<-Direct.GBLUP.8env.175line[foldline,foldenv]+v$g.pred+v$beta-1
  }
}

#train the model using all the environments except for the target one
Direct.EBL.8env.175line.all<-Emergence
for(foldenv in 1:Ne){
  w<-apply(DTH[,-foldenv],1,mean,na.rm=T)
  v<-vigor(w,Geno,"EBL",c(0.1,0.1,1,0.1),"cv",-1,Threshold=9,Printinfo=F)
  Direct.EBL.8env.175line.all[,foldenv]<-Direct.EBL.8env.175line.all[,foldenv]+v$Prediction$Yhat-1
}




#Analysis with C-Bay#################################################################################################################################
#Function to calculate posterior means from MCMC samples
Samples<-201:1200
Mean.abg<-function(Filename,Samples){
  
  g<-as.matrix(read.table(paste(Filename,"_sampledg.txt",sep=""))[Samples,])
  a<-as.matrix(read.table(paste(Filename,"_sampleda.txt",sep=""))[Samples,])
  b<-as.matrix(read.table(paste(Filename,"_sampledb.txt",sep=""))[Samples,])
  #Because the program Heading output G, a, and b in original scale, transform again
  g[g==0]<-1e-6; g<-log(g)
  a[a==0]<-1e-6; a<-log(a)
  b[b==0]<-1e-6; b<-log(b)
  
  data.frame(g=apply(g,2,mean),gsd=apply(g,2,sd),
             a=apply(a,2,mean),asd=apply(a,2,sd),
             b=apply(b,2,mean),bsd=apply(b,2,sd)
             )
}

#Full data
setwd("ResultsOfBayesianInference/BIL.fitting.Ilog")
Ilog.fitting <- as.list(numeric(2))
Ilog.fitting [[1]]<-Mean.abg("BIL.fitting.Ilog_fitting",Samples) #inferred parameter values (posterior mean and sd)
Ilog.fitting [[2]]<-as.matrix(read.table("BIL.fitting.Ilog_fitting_DH.txt")) #fitted values (posterior mean and sd)

Ilog.fitting.sample<-as.list(numeric(4)) #MCMC samples
Ilog.fitting.sample[[1]]<-scan("BIL.fitting.Ilog_fitting_Loglike.txt")
Ilog.fitting.sample[[2]]<-log(as.matrix(read.table("BIL.fitting.Ilog_fitting_sampledg.txt")))
Ilog.fitting.sample[[3]]<-log(as.matrix(read.table("BIL.fitting.Ilog_fitting_sampleda.txt")))
Ilog.fitting.sample[[4]]<-log(as.matrix(read.table("BIL.fitting.Ilog_fitting_sampledb.txt")))

#plot MCMC samples of likelihood
setwd("../../Figures")
tiff("BIL.fitting.Ilog_Loglike.tif",unit="cm",width=8,height=4,res=600)
par(mar=c(3,2.5,1,0.5))
par(cex=0.7)
plot(Ilog.fitting.sample[[1]]-(9*176)/2*log(2*pi),type="l",main="", xlab="",ylab="",ylim=c(1700,2100))
dev.off()
setwd("../")

#LOEO
setwd("ResultsOfBayesianInference/BIL.8env.Ilog")
Ilog.8env <- as.list(numeric(Ne+1))
Ilog.8env[[Ne+1]]<-matrix(0,nc=2*Ne,nr=Nline)
for(foldenv in 1:Ne){
  Ilog.8env[[foldenv]]<-as.list(numeric(2))
  Ilog.8env[[foldenv]][[1]]<-Mean.abg(paste("BIL.8env.Ilog_foldenv",foldenv,sep=""),Samples)  
  Ilog.8env[[foldenv]][[2]]<-as.matrix(read.table(paste("BIL.8env.Ilog_foldenv",foldenv,"_DH.txt",sep="")))
  names(Ilog.8env[[foldenv]])<-c("abg","DH")
  Ilog.8env[[Ne+1]][,c(2*foldenv-1,2*foldenv)]<-Ilog.8env[[foldenv]][[2]][,c(2*foldenv-1,2*foldenv)]
}
setwd("../../")

#Analysis with T-Bay########################################################################################################################################
#Full data
Ilog.fitting.EBL<-as.list(numeric(2))
Ilog.fitting.EBL[[1]]<-matrix(0,nr=176,nc=6)
Ilog.fitting.EBL[[2]]<-matrix(0,nr=176,nc=2*Ne)
Ilog.fitting.coef<-matrix(0,nr=162,nc=6)
for(para in c(1,3,5)){
  v<-vigor(Ilog.fitting[[1]][,para],Geno,"EBL",c(0.1,0.1,1,0.1),Threshold=9,Printinfo=F)
  Ilog.fitting.EBL[[1]][,para]<-colSums(t(Geno)*v$Beta)+v$Alpha
  Ilog.fitting.coef[,c(para,para+1)]<-cbind(v$Beta,v$Sd.beta)
}
Ilog.fitting.coef<-data.frame(Ilog.fitting.coef)
colnames(Ilog.fitting.coef)<-c("g","gsd","a","asd","b","bsd")
for(env in 1:Ne){
  for(line in 1:Nline){
    Ilog.fitting.EBL[[2]][line,2*env-1]<-
      Predict(Temp2[,env],Photo2[,env],exp(Ilog.fitting.EBL[[1]][line,1]),exp(Ilog.fitting.EBL[[1]][line,3]),exp(Ilog.fitting.EBL[[1]][line,5]),
              Emergence[line,env])    
  } 
}

#plot marker effects
setwd("./Figures")
#g
tiff("Ilog.fitting.coef.logg.tif",unit="cm",height=3,width=17,res=600)
par(mar=c(1,2,1,0.5))
par(cex=0.8)
plot(MAP$CumPos,abs(Ilog.fitting.coef[,1]),col=c("gray60","gray10")[MAP$Chr%%2+1],pch=8,main="",xlab="",ylab="",xaxt="n")
for(i in 1:length(Hdnearest)) lines(c(MAP$CumPos[Hdnearest[i]],MAP$CumPos[Hdnearest[i]]),c(0,100),lty=2,lwd=1)
dev.off()
#a
tiff("Ilog.fitting.coef.loga.tif",unit="cm",height=3,width=17,res=600)
par(mar=c(1,2,1,0.5))
par(cex=0.8)
plot(MAP$CumPos,abs(Ilog.fitting.coef[,3]),col=c("gray60","gray10")[MAP$Chr%%2+1],pch=8,main="",xlab="",ylab="",xaxt="n")
for(i in 1:length(Hdnearest)) lines(c(MAP$CumPos[Hdnearest[i]],MAP$CumPos[Hdnearest[i]]),c(0,100),lty=2,lwd=1)
dev.off()
#b
tiff("Ilog.fitting.coef.logb.tif",unit="cm",height=3,width=17,res=600)
par(mar=c(1,2,1,0.5))
par(cex=0.8)
plot(MAP$CumPos,abs(Ilog.fitting.coef[,5]),col=c("gray60","gray10")[MAP$Chr%%2+1],pch=8,main="",xlab="",ylab="",xaxt="n")
for(i in 1:length(Hdnearest)) lines(c(MAP$CumPos[Hdnearest[i]],MAP$CumPos[Hdnearest[i]]),c(0,100),lty=2,lwd=1)
dev.off()
setwd("../")

#LOEO
Ilog.8env.EBL<-as.list(numeric(Ne+1))
for(env in 1:Ne){
  Ilog.8env.EBL[[env]]<-matrix(0,nr=176,nc=6)
}
Ilog.8env.EBL[[Ne+1]]<-matrix(0,nr=176,nc=Ne*2)
for(foldenv in 1:Ne){
  for(para in c(1,3,5)){
    cat(foldenv,para,"\n")
    v<-vigor(Ilog.8env[[foldenv]][[1]][,para],Geno,"EBL",c(0.1,0.1,1,0.1),Threshold=9,Printinfo=F)
    Ilog.8env.EBL[[foldenv]][,para]<-colSums(t(Geno)*v$Beta)+v$Alpha
  }
  for(line in 1:Nline){
    Ilog.8env.EBL[[Ne+1]][line,2*foldenv-1]<-
      Predict(Temp2[,foldenv],Photo2[,foldenv],exp(Ilog.8env.EBL[[foldenv]][line,1]),exp(Ilog.8env.EBL[[foldenv]][line,3]),
              exp(Ilog.8env.EBL[[foldenv]][line,5]),Emergence[line,foldenv])    
  }
}

#LOLO
Ilog.175line.EBL <- as.list(numeric(2))
Ilog.175line.EBL[[1]]<-matrix(0,nc=6,nr=Nline)
Ilog.175line.EBL[[2]]<-matrix(0,nc=2*Ne,nr=Nline)
names(Ilog.175line.EBL)<-c("abg","DH")
for(para in c(1,3,5)){
  v<-vigor(Ilog.fitting[[1]][,para],Geno,"EBL",c(0.1,0.1,1,0.1),"cv",-1,Threshold=9,Printinfo=F)
  Ilog.175line.EBL[[1]][,para]<-v$Prediction$Yhat
}
for(foldline in 1:Nline){
  g<-Ilog.175line.EBL[[1]][foldline,1]
  a<-Ilog.175line.EBL[[1]][foldline,3]
  b<-Ilog.175line.EBL[[1]][foldline,5]
  for(env in 1:Ne){
    Ilog.175line.EBL[[2]][foldline,2*env-1]<-Predict(Temp2[,env],Photo2[,env],exp(g),exp(a),exp(b),Emergence[foldline,env])   
  }        
}


#LOELO
Ilog.8env.175line.EBL<-as.list(numeric(Ne+1))
Ilog.8env.175line.EBL[[Ne+1]]<-matrix(0,nc=2*Ne,nr=Nline)
for(foldenv in 1:Ne){
  Ilog.8env.175line.EBL[[foldenv]]<-matrix(0,nr=Nline,nc=6)
  for(para in c(1,3,5)){ 
    v<-vigor(Ilog.8env[[foldenv]][[1]][,para],Geno,"EBL",c(0.1,0.1,1,0.1),"cv",-1,Threshold=9,Printinfo=F)
    Ilog.8env.175line.EBL[[foldenv]][,para]<-v$Prediction$Yhat
  }
}
for(foldenv in 1:Ne){
  for(foldline in 1:Nline){
    cat(foldenv,foldline,"\n")
    g<-Ilog.8env.175line.EBL[[foldenv]][foldline,1]
    a<-Ilog.8env.175line.EBL[[foldenv]][foldline,3]
    b<-Ilog.8env.175line.EBL[[foldenv]][foldline,5]
    Ilog.8env.175line.EBL[[Ne+1]][foldline,2*foldenv-1]<-Predict(Temp2[,foldenv],Photo2[,foldenv],exp(g),exp(a),exp(b),Emergence[foldline,foldenv])
  }
}



#Analysis with IM#############################################################################################################################################
#Full data
setwd("ResultsOfBayesianInference/BIL.fitting.EBL")
EBL.fitting.sample <- as.list(numeric(3))
EBL.fitting.sample[[1]]<-as.matrix(read.table("BIL.fitting.EBL_fitting_sampledg.txt"))
EBL.fitting.sample[[2]]<-as.matrix(read.table("BIL.fitting.EBL_fitting_sampleda.txt"))
EBL.fitting.sample[[3]]<-as.matrix(read.table("BIL.fitting.EBL_fitting_sampledb.txt"))
#Because the program Heading output G, a, and b in original scale, transform again
EBL.fitting.sample[[1]][EBL.fitting.sample[[1]]==0]<-1e-6;EBL.fitting.sample[[1]]<-log(EBL.fitting.sample[[1]])
EBL.fitting.sample[[2]][EBL.fitting.sample[[2]]==0]<-1e-6;EBL.fitting.sample[[2]]<-log(EBL.fitting.sample[[2]])
EBL.fitting.sample[[3]][EBL.fitting.sample[[3]]==0]<-1e-6;EBL.fitting.sample[[3]]<-log(EBL.fitting.sample[[3]])

EBL.fitting <- as.list(numeric(2)) #posterior mean and sd of the DVR model parameters
EBL.fitting [[1]]<-data.frame(g=apply(EBL.fitting.sample[[1]],2,mean),gsd=apply(EBL.fitting.sample[[1]],2,sd),
                              a=apply(EBL.fitting.sample[[2]],2,mean),asd=apply(EBL.fitting.sample[[2]],2,sd),
                              b=apply(EBL.fitting.sample[[3]],2,mean),bsd=apply(EBL.fitting.sample[[3]],2,sd)
                              )
EBL.fitting [[2]]<-as.matrix(read.table("BIL.fitting.EBL_fitting_DH.txt")) #fitted values

EBL.fitting.coef <- read.table("BIL.fitting.EBL_fitting_coefficient.txt",header=T) #marker effects
EBL.fitting.coef <- EBL.fitting.coef[-1,] #the first element is the intercept 

#plot marker effects
setwd("../../Figures")
#g
tiff("EBL.fitting.coef.g.tif",unit="cm",height=3,width=17,res=600)
par(mar=c(1,2,1,0.5))
par(cex=0.8)
plot(MAP$CumPos,abs(EBL.fitting.coef$g),col=c("gray60","gray10")[MAP$Chr%%2+1],pch=8,main="",xlab="",ylab="",xaxt="n")
for(i in 1:length(Hdnearest)) lines(c(MAP$CumPos[Hdnearest[i]],MAP$CumPos[Hdnearest[i]]),c(0,100),lty=2,lwd=1)
dev.off()
#a
tiff("EBL.fitting.coef.a.tif",unit="cm",height=3,width=17,res=600)
par(mar=c(1,2,1,0.5))
par(cex=0.8)
plot(MAP$CumPos,abs(EBL.fitting.coef$a),col=c("gray60","gray10")[MAP$Chr%%2+1],pch=8,main="",xlab="",ylab="",xaxt="n")
for(i in 1:length(Hdnearest)) lines(c(MAP$CumPos[Hdnearest[i]],MAP$CumPos[Hdnearest[i]]),c(0,100),lty=2,lwd=1)
dev.off()
#b
tiff("EBL.fitting.coef.b.tif",unit="cm",height=3,width=17,res=600)
par(mar=c(1,2,1,0.5))
par(cex=0.8)
plot(MAP$CumPos,abs(EBL.fitting.coef$b),col=c("gray60","gray10")[MAP$Chr%%2+1],pch=8,main="",xlab="",ylab="",xaxt="n")
for(i in 1:length(Hdnearest)) lines(c(MAP$CumPos[Hdnearest[i]],MAP$CumPos[Hdnearest[i]]),c(0,100),lty=2,lwd=1)
dev.off()
setwd("../")

#assess the reproducibility of IM
EBL.fitting.coef.rep <- EBL.fitting.rep <- as.list(numeric(4))
for(i in 1:4){ #four analyses were conducted with different initial values
  setwd(paste("ResultsOfBayesianInference/BIL.fitting.EBL",i,sep=""))
  EBL.fitting.coef.rep[[i]] <- read.table("BIL.fitting.EBL_fitting_coefficient.txt",header=T)
  EBL.fitting.coef.rep[[i]] <- EBL.fitting.coef.rep[[i]][-1,]
  g<-as.matrix(read.table("BIL.fitting.EBL_fitting_sampledg.txt"));g[g==0]<-1e-6
  a<-as.matrix(read.table("BIL.fitting.EBL_fitting_sampleda.txt"));a[a==0]<-1e-6
  b<-as.matrix(read.table("BIL.fitting.EBL_fitting_sampledb.txt"));b[b==0]<-1e-6
  EBL.fitting.rep[[i]]<-data.frame(g=apply(log(g),2,mean),
                                   gsd=apply(log(g),2,sd),
                                   a=apply(log(a),2,mean),
                                   asd=apply(log(a),2,sd),
                                   b=apply(log(b),2,mean),
                                   bsd=apply(log(b),2,sd))
  setwd("../../")
}
#correlation between replicates
mean(cor(cbind(EBL.fitting.coef$g,EBL.fitting.coef.rep[[1]]$g,EBL.fitting.coef.rep[[2]]$g,EBL.fitting.coef.rep[[3]]$g,EBL.fitting.coef.rep[[4]]$g))[upper.tri(matrix(0,5,5))])
mean(cor(cbind(EBL.fitting.coef$a,EBL.fitting.coef.rep[[1]]$a,EBL.fitting.coef.rep[[2]]$a,EBL.fitting.coef.rep[[3]]$a,EBL.fitting.coef.rep[[4]]$a))[upper.tri(matrix(0,5,5))])
mean(cor(cbind(EBL.fitting.coef$b,EBL.fitting.coef.rep[[1]]$b,EBL.fitting.coef.rep[[2]]$b,EBL.fitting.coef.rep[[3]]$b,EBL.fitting.coef.rep[[4]]$b))[upper.tri(matrix(0,5,5))])
#0.9880665; 0.998862; 0.9967022
mean(cor(cbind(EBL.fitting[[1]]$g,EBL.fitting.rep[[1]]$g,EBL.fitting.rep[[2]]$g,EBL.fitting.rep[[3]]$g,EBL.fitting.rep[[4]]$g))[upper.tri(matrix(0,5,5))])
mean(cor(cbind(EBL.fitting[[1]]$a,EBL.fitting.rep[[1]]$a,EBL.fitting.rep[[2]]$a,EBL.fitting.rep[[3]]$a,EBL.fitting.rep[[4]]$a))[upper.tri(matrix(0,5,5))])
mean(cor(cbind(EBL.fitting[[1]]$b,EBL.fitting.rep[[1]]$b,EBL.fitting.rep[[2]]$b,EBL.fitting.rep[[3]]$b,EBL.fitting.rep[[4]]$b))[upper.tri(matrix(0,5,5))])
#0.9978931; 0.9970721; 0.9989923

#LOEO
setwd("ResultsOfBayesianInference/BIL.8env.EBL")
EBL.8env <- as.list(numeric(Ne+1))
EBL.8env[[Ne+1]]<-matrix(0,nc=2*Ne,nr=Nline)
for(foldenv in 1:Ne){
  EBL.8env[[foldenv]]<-as.list(numeric(2))
  EBL.8env[[foldenv]][[1]]<-Mean.abg(paste("BIL.8env.EBL_foldenv",foldenv,sep=""),1:300)  
  EBL.8env[[foldenv]][[2]]<-as.matrix(read.table(paste("BIL.8env.EBL_foldenv",foldenv,"_DH.txt",sep="")))
  names(EBL.8env[[foldenv]])<-c("abg","DH")
  EBL.8env[[Ne+1]][,c(2*foldenv-1,2*foldenv)]<-EBL.8env[[foldenv]][[2]][,c(2*foldenv-1,2*foldenv)]  
}

#Posterior uncertainty
LOOE.unc<-matrix(0,nr=Nline,nc=Ne)
EBL.sample<-as.list(numeric(3))
EBL.PosteriorSD.LOOE<-matrix(NA,nr=Nline,nc=Ne)
for(env in 1:Ne){    

  DHsamples.EBL<-numeric(300)
  EBL.sample[[1]]<-as.matrix(read.table(paste("BIL.8env.EBL_foldenv",env,"_sampledg.txt",sep="")))
  EBL.sample[[2]]<-as.matrix(read.table(paste("BIL.8env.EBL_foldenv",env,"_sampleda.txt",sep="")))      
  EBL.sample[[3]]<-as.matrix(read.table(paste("BIL.8env.EBL_foldenv",env,"_sampledb.txt",sep="")))
    
  for(line in 1:Nline){
    if(!is.na(Heading[line,env]))
    {   
        for(i in 1:300){
          DHsamples.EBL[i]<-
            Predict(Temp2[,env],Photo2[,env],EBL.sample[[1]][i,line],EBL.sample[[2]][i,line],EBL.sample[[3]][i,line],Emergence[line,env])  
        }     
        if(median(DHsamples.EBL[1:300])<=Heading[line,env]){
          LOOE.unc[line,env]<-rank(c(Heading[line,env],DHsamples.EBL),ties.method="max")[1]/(300+1)
        }else{
          LOOE.unc[line,env]<-rank(c(Heading[line,env],DHsamples.EBL),ties.method="min")[1]/(300+1)        
        }
        EBL.PosteriorSD.LOOE[line,env]<-sd(DHsamples.EBL)        
    }
  }
}

hist(LOOE.unc[LOOE.unc[[3]]>0],breaks=seq(0,1,0.05))
for(i in c(0.95,0.9,0.8,0.7,0.6)){
  cat("EBL",i,length(LOOE.unc[(LOOE.unc>(1-i)/2)&(LOOE.unc<((1-i)/2+i))])/length(LOOE.unc[LOOE.unc>0]),"\n")
}
#EBL 0.95 0.4197531 
#EBL 0.9 0.345679 
#EBL 0.8 0.2586095 
#EBL 0.7 0.1877843 
#EBL 0.6 0.1403509 
mean(EBL.PosteriorSD.LOOE,na.rm=T);sd(EBL.PosteriorSD.LOOE,na.rm=T);
#1.609866; 0.5192332
setwd("../../")

#LOLO
EBL.175line <- as.list(numeric(Nline+1))
EBL.175line[[Nline+1]]<-matrix(0,nc=2*Ne,nr=Nline)
foldline<-1
for(d in 1:35){
  setwd(paste("ResultsOfBayesianInference/BIL.175line.",d,sep="")) 
  for(k in 1:6){
    if((d<35&k<6)|d==35)
    {
      EBL.175line[[foldline]]<-as.list(numeric(2))   
      EBL.175line[[foldline]][[1]]<-Mean.abg(paste("BIL.175line.",d,".EBL_foldline",k,sep=""),1:300)  
      EBL.175line[[foldline]][[2]]<-as.matrix(read.table(paste("BIL.175line.",d,".EBL_foldline",k,"_DH.txt",sep="")))
      names(EBL.175line[[foldline]])<-c("abg","DH")
      EBL.175line[[Nline+1]][foldline,]<-EBL.175line[[foldline]][[2]][foldline,]
      foldline<-foldline+1
    }
  }
  setwd("../../")
}

#Posterior uncertainty
LOOG.unc<-matrix(0,nr=Nline,nc=Ne)
EBL.sample<-as.list(numeric(3))
EBL.PosteriorSD.LOOG<-matrix(NA,nr=Nline,nc=Ne)
line<-1
for(d in 1:35){
  for(k in 1:6){
    if((d<35&k<6)|d==35){
      cat("line",line,"\n")

      DHsamples.EBL<-numeric(300)
      EBL.sample[[1]]<-as.matrix(read.table(paste("BIL.175line.",d,".EBL_foldline",k,"_sampledg.txt",sep="")))
      EBL.sample[[2]]<-as.matrix(read.table(paste("BIL.175line.",d,".EBL_foldline",k,"_sampleda.txt",sep="")))      
      EBL.sample[[3]]<-as.matrix(read.table(paste("BIL.175line.",d,".EBL_foldline",k,"_sampledb.txt",sep="")))

      for(env in 1:Ne){
        if(!is.na(Heading[line,env]))
        {            
          for(i in 1:300){
            DHsamples.EBL[i]<-
              Predict(Temp2[,env],Photo2[,env],EBL.sample[[1]][i,line],EBL.sample[[2]][i,line],EBL.sample[[3]][i,line],
                      Emergence[line,env])  
          }
        
          if(median(DHsamples.EBL[1:300])<=Heading[line,env]){
            LOOG.unc[line,env]<-rank(c(Heading[line,env],DHsamples.EBL),ties.method="max")[1]/(300+1)
          }else{
            LOOG.unc[line,env]<-rank(c(Heading[line,env],DHsamples.EBL),ties.method="min")[1]/(300+1)        
          }
          EBL.PosteriorSD.LOOG[line,env]<-sd(DHsamples.EBL)
        } 
      }#env        
      line<-line+1
    }
  }
}

hist(LOOG.unc[LOOG.unc>0],breaks=seq(0,1,0.05))
for(i in c(0.95,0.9,0.8,0.7,0.6)){
  cat("EBL",i,length(LOOG.unc[(LOOG.unc>(1-i)/2)&(LOOG.unc<((1-i)/2+i))])/length(LOOG.unc[LOOG.unc>0]),"\n")
}
#EBL 0.95 0.8615984 
#EBL 0.9 0.7842755 
#EBL 0.8 0.6861598 
#EBL 0.7 0.5789474 
#EBL 0.6 0.4879792 
mean(EBL.PosteriorSD.LOOG,na.rm=T);sd(EBL.PosteriorSD.LOOG,na.rm=T);
#6.014466 2.450461
setwd("../../")

#LOELO
EBL.8env.175line <- matrix(0,nc=2*Ne,nr=Nline)
foldline<-1
for(d in 1:35){
  setwd(paste("ResultsOfBayesianInference/BIL.8env.175line.",d,sep=""))   
  for(k in 1:5){
    for(foldenv in 1:9){
      temp<-as.matrix(read.table(paste("BIL.8env.175line.",d,".EBL_foldline",k,"_foldenv",foldenv,"_DH.txt",sep="")))
      EBL.8env.175line[foldline,c(2*foldenv-1,2*foldenv)]<-temp[foldline,c(2*foldenv-1,2*foldenv)]
    }
    foldline<-foldline+1
  }
  setwd("../../")
}
k<-6
for(foldenv in 1:9){
  temp<-as.matrix(read.table(paste("BIL.8env.175line.",d,".EBL_foldline",k,"_foldenv",foldenv,"_DH.txt",sep="")))
  EBL.8env.175line[foldline,c(2*foldenv-1,2*foldenv)]<-temp[foldline,c(2*foldenv-1,2*foldenv)]
}

#Posterior uncertainty
LOOEG.unc<-matrix(0,nr=Nline,nc=Ne)
EBL.sample<-as.list(numeric(3))
EBL.PosteriorSD.LOOEG<-matrix(NA,nr=Nline,nc=Ne)
line<-1
for(d in 1:35){
  for(k in 1:6){
    if((d<35&k<6)|d==35){
      for(env in 1:Ne){
        cat("line",line,"env",env,"\n")
        if(!is.na(Heading[line,env]))
        {
          DHsamples.EBL<-numeric(300)
          EBL.sample[[1]]<-as.matrix(read.table(paste("BIL.8env.175line.",d,".EBL_foldline",k,"_foldenv",env,"_sampledg.txt",sep="")))
          EBL.sample[[2]]<-as.matrix(read.table(paste("BIL.8env.175line.",d,".EBL_foldline",k,"_foldenv",env,"_sampleda.txt",sep="")))      
          EBL.sample[[3]]<-as.matrix(read.table(paste("BIL.8env.175line.",d,".EBL_foldline",k,"_foldenv",env,"_sampledb.txt",sep="")))
          
          for(i in 1:300){
            DHsamples.EBL[i]<-
              Predict(Temp2[,env],Photo2[,env],EBL.sample[[1]][i,line],EBL.sample[[2]][i,line],EBL.sample[[3]][i,line],Emergence[line,env])  
          }
          
          if(median(DHsamples.EBL[1:300])<=Heading[line,env]){
            LOOEG.unc[line,env]<-rank(c(Heading[line,env],DHsamples.EBL),ties.method="max")[1]/(300+1)
          }else{
            LOOEG.unc[line,env]<-rank(c(Heading[line,env],DHsamples.EBL),ties.method="min")[1]/(300+1)        
          }
          EBL.PosteriorSD.LOOEG[line,env]<-sd(DHsamples.EBL)       
        }  
      }#env
      line<-line+1
    }
  }
}

hist(LOOEG.unc[LOOEG.unc>0],breaks=seq(0,1,0.05))
for(i in c(0.95,0.9,0.8,0.7,0.6)){
  cat("EBL",i,length(LOOEG.unc[(LOOEG.unc>(1-i)/2)&(LOOEG.unc<((1-i)/2+i))])/length(LOOEG.unc[LOOEG.unc>0]),"\n")
}
#EBL 0.95 0.8362573 
#EBL 0.9 0.7602339 
#EBL 0.8 0.6504224 
#EBL 0.7 0.5490578 
#EBL 0.6 0.4580897 
mean(EBL.PosteriorSD.LOOEG,na.rm=T);sd(EBL.PosteriorSD.LOOEG,na.rm=T);
#6.0678 2.526105
setwd("../../")

#Compare prediction errors########################################################################################################
#used to specify columns
Odds<-seq(1,2*Ne,2)
Even<-seq(2,2*Ne,2)

#Full data
#RMSE
apply(abs(Direct.EBL.fitting-Heading),2,rmse);rmse(abs(Direct.EBL.fitting-Heading))
apply(abs(O.fitting[[2]][,Odds]-Heading),2,rmse);rmse(abs(O.fitting[[2]][,Odds]-Heading))
apply(abs(Ilog.fitting[[2]][,Odds]-Heading),2,rmse);rmse(abs(Ilog.fitting[[2]][,Odds]-Heading))
apply(abs(O.fitting.EBL[[2]][,Odds]-Heading),2,rmse);rmse(abs(O.fitting.EBL[[2]][,Odds]-Heading))
apply(abs(Ilog.fitting.EBL[[2]][,Odds]-Heading),2,rmse);rmse(abs(Ilog.fitting.EBL[[2]][,Odds]-Heading))
apply(abs(EBL.fitting[[2]][,Odds]-Heading),2,rmse);rmse(abs(EBL.fitting[[2]][,Odds]-Heading))

#           T2007      FO      HN        IG      IK        TN    T2008E   T2008L    T2009     All
#GP(EBL) 5.751842 6.188219 4.041799 4.180288 6.861392 3.759361 6.584877 5.154746 6.383042 5.583012
#C-Nel   4.313035 4.370288 3.080342 2.394773 1.673660 2.289813 2.277608 6.430131 2.809278 3.604109
#C-Bay   4.323090 4.993630 3.413094 2.548988 1.743552 3.560448 1.730771 7.463491 3.314021 4.057269
#T-Nel   6.649932 7.912534 4.662149 4.795203 6.475039 4.184915 6.744526 8.089696 6.857726 6.441621
#T-Bay   7.008112 11.34364 4.970025 5.396340 10.59936 4.372828 10.02950 7.37432612.260895 8.686138
#IM      4.337125 4.796980 3.451707 2.535122 1.858361 3.209109 2.090527 7.284266 3.007733 3.966844

#Correlation
diag(cor(cbind(Direct.EBL.fitting,Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Direct.EBL.fitting),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(O.fitting[[2]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(O.fitting[[2]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(Ilog.fitting[[2]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Ilog.fitting[[2]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(O.fitting.EBL[[2]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(O.fitting.EBL[[2]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(Ilog.fitting.EBL[[2]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Ilog.fitting.EBL[[2]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(EBL.fitting[[2]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(EBL.fitting[[2]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")

#GP(EBL) 0.8517118 0.8639958 0.8213335 0.8560473 0.8728057 0.8658527 0.8884947 0.8566027 0.8865676 0.9641027
#C-Nel   0.9757312 0.9666893 0.9108449 0.9684807 0.9939060 0.9510044 0.9917173 0.9471892 0.9879455 0.9853005
#C-Bay   0.9839665 0.9685490 0.8792556 0.9547781 0.9927067 0.9132043 0.9943990 0.9625058 0.9863893 0.9820799
#T-Nel   0.8649093 0.8401366 0.7436821 0.8314384 0.8980867 0.8355793 0.8952710 0.8339637 0.8964494 0.9527208
#T-Bay   0.8074686 0.7416820 0.7107293 0.8321293 0.8254298 0.8423864 0.8449512 0.7751707 0.8108757 0.9316485
#IM      0.9842408 0.9706550 0.8758607 0.9564183 0.9918301 0.9258700 0.9920725 0.9609912 0.9882448 0.9829071

#slope
coef(lm(Direct.EBL.fitting[1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(O.fitting[[2]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(Ilog.fitting[[2]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(O.fitting.EBL[[2]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(Ilog.fitting.EBL[[2]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(EBL.fitting[[2]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))

          #intercept slope
#GP(EBL) 7.648266  0.912508
#C-Nel   3.4621788 0.9644606
#C-Bay   5.2312039 0.9499395 
#T-Nel   11.758474 0.865609
#T-Bay   18.581343 0.7506323 
#IM      5.0836720 0.9515600


#LOEO
#RMSE
apply(abs(O.8env[[Ne+1]][,Odds]-Heading),2,rmse);rmse(abs(O.8env[[Ne+1]][,Odds]-Heading))
apply(abs(Ilog.8env[[Ne+1]][,Odds]-Heading),2,rmse);rmse(abs(Ilog.8env[[Ne+1]][,Odds]-Heading))
apply(abs(O.8env.EBL[[Ne+1]][,Odds]-Heading),2,rmse);rmse(abs(O.8env.EBL[[Ne+1]][,Odds]-Heading))
apply(abs(Ilog.8env.EBL[[Ne+1]][,Odds]-Heading),2,rmse);rmse(abs(Ilog.8env.EBL[[Ne+1]][,Odds]-Heading))
apply(abs(EBL.8env[[Ne+1]][,Odds]-Heading),2,rmse);rmse(abs(EBL.8env[[Ne+1]][,Odds]-Heading))

#            T2007      FO      HN        IG      IK        TN    T2008E   T2008L    T2009     All
#C-Nel     5.247294 5.765901 4.619424 3.907823 2.921317 3.482582 12.59238 7.348856 3.753786 6.242292
#C-Bay     5.184045 7.499569 4.361587 3.544756 2.530542 5.015790 4.088655 7.870381 4.495350 5.223696
#T-Nel     6.897793 8.626086 5.309934 5.536843 6.525300 4.483453 11.04227 8.547328 7.455169 7.460618
#T-Bay     6.550069 13.63003 5.192833 5.656854 11.04407 5.028970 10.77085 7.626583 13.79640 9.492927
#IM        5.151350 6.809989 4.337664 3.388971 2.366832 4.209842 4.918262 7.867051 4.087416 5.064911

#LOEO without R1684
apply(abs(O.8env.EBL.woR1684[[Ne+1]][,Odds]-Heading),2,rmse);rmse(abs(O.8env.EBL.woR1684[[Ne+1]][,Odds]-Heading))
#T-Nel     7.036431 8.783291 5.319666 5.505200 6.771011 4.487972 11.154392 8.596709 7.695926 7.572049

#Correlation
diag(cor(cbind(O.8env[[Ne+1]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(O.8env[[Ne+1]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(Ilog.8env[[Ne+1]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Ilog.8env[[Ne+1]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(O.8env.EBL[[Ne+1]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(O.8env.EBL[[Ne+1]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(Ilog.8env.EBL[[Ne+1]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Ilog.8env.EBL[[Ne+1]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(EBL.8env[[Ne+1]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(EBL.8env[[Ne+1]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")

#C-Nel     0.9719507 0.9309434 0.8300298 0.9254772 0.9842707 0.8816064 0.8301553 0.9264243 0.9808998 0.9554079
#C-Bay     0.9765378 0.9118053 0.8089838 0.9087498 0.9850775 0.8218225 0.9745260 0.9467602 0.9779346 0.9688601
#T-Nel     0.8705566 0.8129008 0.6860544 0.8192625 0.8947017 0.8153760 0.8297193 0.8194291 0.8859120 0.9376822
#T-Bay     0.8131417 0.5425345 0.6741456 0.8252798 0.8383078 0.7824854 0.8420325 0.7641808 0.7847242 0.9168536
#IM        0.9782243 0.9312717 0.8163727 0.9211442 0.9858751 0.8656501 0.9702318 0.9504344 0.9809385 0.9707253

#slope
coef(lm(O.8env[[Ne+1]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(Ilog.8env[[Ne+1]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(O.8env.EBL[[Ne+1]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(Ilog.8env.EBL[[Ne+1]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(EBL.8env[[Ne+1]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))

          #intercept slope
#C-Nel     8.5441772 0.8947598
#C-Bay     7.3078154 0.9208669
#T-Nel     15.029888 0.8186717
#T-Bay     20.808047 0.7214733 
#IM        6.5428450 0.9301512

#LOLO
#RMSE
apply(abs(Direct.EBL.175line-Heading),2,rmse);rmse(abs(Direct.EBL.175line-Heading))
apply(abs(Direct.GBLUP.175line-Heading),2,rmse);rmse(abs(Direct.GBLUP.175line-Heading))
apply(abs(O.175line.EBL$DH[,Odds]-Heading),2,rmse);rmse(abs(O.175line.EBL$DH[,Odds]-Heading))
apply(abs(Ilog.175line.EBL[[2]][,Odds]-Heading),2,rmse);rmse(abs(Ilog.175line.EBL[[2]][,Odds]-Heading))
apply(abs(EBL.175line[[Nline+1]][,Odds]-Heading),2,rmse);rmse(abs(EBL.175line[[Nline+1]][,Odds]-Heading))

#            T2007      FO      HN        IG      IK        TN    T2008E   T2008L    T2009     All
#GP(EBL)  7.004241 7.898256 6.005174 6.175481 8.549422 5.500601 8.298109 6.715000 7.991629 7.234368
#GP(GBLUP)8.119806 8.834194 5.859601 6.311250 9.819537 5.876498 9.803672 7.282159 9.457581 8.125417
#T-Nel    7.600090 9.067325 5.760548 5.860322 8.147434 5.218859 8.295672 8.879240 8.506348 7.653009
#T-Bay    7.745967 12.04450 5.798831 6.284423 11.33302 5.292141 10.83738 8.034726 12.98381 9.408324
#IM       7.310095 8.401885 5.779242 5.774076 7.556609 5.382575 7.731110 9.364695 7.845141 7.385534

#Correlation
diag(cor(cbind(Direct.EBL.175line,Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Direct.EBL.175line),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(Direct.GBLUP.175line,Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Direct.GBLUP.175line),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(O.175line.EBL$DH[,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(O.175line.EBL$DH[,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(Ilog.175line.EBL[[2]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Ilog.175line.EBL[[2]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(EBL.175line[[Nline+1]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(EBL.175line[[Nline+1]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")

#GP(EBL)  0.7664886 0.7589985 0.4971660 0.6377305 0.7905585 0.6705449 0.8117898 0.7332323 0.8120917 0.9386917
#GP(GBLUP)0.6677940 0.6857005 0.5232147 0.6086272 0.7110097 0.6042871 0.7248922 0.6737922 0.7244150 0.9220058
#T-Nel    0.7859122 0.7289697 0.5905511 0.7025250 0.8177778 0.7081397 0.8251733 0.7206670 0.8117228 0.9316071
#T-Bay    0.7308106 0.6376792 0.5801621 0.6977194 0.7554006 0.7098129 0.7849237 0.6712587 0.7383963 0.9145747
#IM       0.8239809 0.7697239 0.5858065 0.6888753 0.8407786 0.7078637 0.8472258 0.7573428 0.8382437 0.9367182


#slope
coef(lm(Direct.EBL.175line[1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(Direct.GBLUP.175line[1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(O.175line.EBL$DH[,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(Ilog.175line.EBL[[2]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(EBL.175line[[Nline+1]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))

          #intercept slope
#GP(EBL)   9.9248190 0.8870839
#GP(GBLUP) 12.339631 0.8590028 
#T-Nel     14.012414 0.8391637
#T-Bay     19.829384 0.7350249
#IM        12.276159 0.8681118


#LOELO
#RMSE
apply(abs(Direct.EBL.8env.175line-Heading),2,rmse);rmse(abs(Direct.EBL.8env.175line-Heading))
apply(abs(Direct.GBLUP.8env.175line-Heading),2,rmse);rmse(abs(Direct.GBLUP.8env.175line-Heading))
apply(abs(Direct.EBL.8env.175line.all-Heading),2,rmse);rmse(abs(Direct.EBL.8env.175line.all-Heading))
apply(abs(O.8env.175line.EBL[[Ne+1]][,Odds]-Heading),2,rmse);rmse(abs(O.8env.175line.EBL[[Ne+1]][,Odds]-Heading))
apply(abs(Ilog.8env.175line.EBL[[Ne+1]][,Odds]-Heading),2,rmse);rmse(abs(Ilog.8env.175line.EBL[[Ne+1]][,Odds]-Heading))
apply(abs(EBL.8env.175line[,Odds]-Heading),2,rmse);rmse(abs(EBL.8env.175line[,Odds]-Heading))

#            T2007      FO      HN        IG      IK        TN    T2008E   T2008L    T2009     All
#GP(EBL)   12.12028 11.45454 6.248535 6.477606 8.560917 6.101026 19.63549 9.511743 8.294040 10.70441
#GP(GBLUP) 12.81274 12.12216 6.121790 6.554352 9.930074 6.203030 20.23028 9.849496 9.609446 11.28748
#GP(EBLall)9.334922 10.24312 22.18773 21.57753 19.92887 23.33470 38.59920 16.43569 18.82514 21.61225
#T-Nel     7.845091 9.691137 6.106337 6.317408 8.346543 5.373232 12.49272 9.196405 8.888834 8.561136
#T-Bay     7.648901 14.65070 6.186284 6.408767 11.57215 5.858696 11.56945 8.512024 14.33190 10.24673
#IM        7.478246 9.847569 6.142419 5.940836 7.760997 5.730730 8.608967 9.781782 8.635301 7.952053

#Correlation
diag(cor(cbind(Direct.EBL.8env.175line,Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Direct.EBL.8env.175line),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(Direct.GBLUP.8env.175line,Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Direct.GBLUP.8env.175line),as.vector(Heading),use="pairwise.complete")
diag(cor(cbind(Direct.EBL.8env.175line.all,Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Direct.EBL.8env.175line.all),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(O.8env.175line.EBL[[Ne+1]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(O.8env.175line.EBL[[Ne+1]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(Ilog.8env.175line.EBL[[Ne+1]][,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(Ilog.8env.175line.EBL[[Ne+1]][,Odds]),as.vector(Heading),use="pairwise.complete.obs")
diag(cor(cbind(EBL.8env.175line[,Odds],Heading),use="pairwise.complete.obs")[(Ne+1):(2*Ne),1:Ne]);cor(as.vector(EBL.8env.175line[,Odds]),as.vector(Heading),use="pairwise.complete.obs")

#GP(EBL)   0.7722683 0.7387117 0.5285189 0.5761220 0.7931362 0.6157740 0.7895968 0.5178414 0.8004131 0.865995
#GP(GBLUP) 0.6804828 0.6695385 0.5171043 0.5625265 0.7073856 0.5850034 0.7114935 0.4645345 0.7180608 0.8493301
#GP(EBLall)0.7050294 0.7256215 0.1914264 0.3705071 0.7205411 0.1909884 0.6772369 0.7026699 0.7064332 0.1058366
#T-Nel     0.7946188 0.6932118 0.5623162 0.6980513 0.8053437 0.6960961 0.7297869 0.7203679 0.8075674 0.9153225
#T-Bay     0.7148260 0.2735767 0.5101267 0.6994475 0.7860739 0.6260976 0.7831049 0.6217288 0.7212464 0.895491
#IM        0.8336310 0.6996155 0.5540926 0.6767816 0.8310955 0.6722134 0.8285574 0.7297034 0.8136272 0.9255645

#slope
coef(lm(Direct.EBL.8env.175line[1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(Direct.GBLUP.8env.175line[1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(Direct.EBL.8env.175line.all[1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(O.8env.175line.EBL[[Ne+1]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(Ilog.8env.175line.EBL[[Ne+1]][,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))
coef(lm(EBL.8env.175line[,Odds][1:(Ne*Nline)]~Heading[1:(Ne*Nline)]))

          #intercept slope
#GP(EBL)   19.4053962 0.7537069
#GP(GBLUP) 21.6946366 0.7269076
#GP(EBLall)84.0624633 0.03967358
#T-Nel     17.3575873 0.7913753
#T-Bay     22.1037759 0.7062496
#IM        13.1565906 0.8538354


