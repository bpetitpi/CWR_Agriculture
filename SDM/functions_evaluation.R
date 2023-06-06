boycei<-function(interval,obs,fit){
  
  fit.bin<-fit
  obs.bin<-obs
  fit.bin[fit[]>=interval[1]&fit[]<=interval[2]]<-"i";fit.bin[fit.bin!="i"]<-0
  obs.bin[obs[]>=interval[1]&obs[]<=interval[2]]<-"i";obs.bin[obs.bin!="i"]<-0
  
  pi<-length(which(obs.bin=="i"))/length(obs)
  ei<-length(which(fit.bin=="i"))/length(fit.bin)
  fi<-pi/ei
  
  return(fi)
}


ecospat.boyce <- function (fit, obs, nclass = 0, window.w = "default", res = 100,
                           PEplot = T, cor.method="spearman")
{
  if(class(fit)=="RasterLayer"){
    if(class(obs)=="data.frame"){
      obs<-extract(fit, obs)}
    fit<-getValues(fit)
    fit <- fit[!is.na(fit)]
  }
  
  if (window.w == "default") {
    window.w <- (max(fit) - min(fit))/10
  }
  interval <- c(min(fit), max(fit))
  mini <- interval[1]
  maxi <- interval[2]
  if (nclass == 0) {
    vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi -
                                                              mini - window.w)/res)
    vec.mov[res + 1] <- vec.mov[res + 1] + 1 #Trick to avoid error with closed interval in R
    interval <- cbind(vec.mov, vec.mov + window.w)
  }
  else if (length(nclass) > 1) {
    vec.mov <- c(mini, nclass)
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }
  else if (nclass > 0 & length(nclass) < 2) {
    vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini)/nclass)
  }
  f <- apply(interval, 1, boycei, obs, fit)
  to.keep<-which(f!="NaN") # index to keep no NaN data
  f<-f[to.keep]
  if (length(f) < 2) {
    b <- NA #at least two points are necessary to draw a correlation
  }
  else {
    r<-c(1:length(f))[f!=c(f[-1],FALSE)] #index to remove successive duplicates
    b <- cor(f[r], vec.mov[to.keep][r], method = cor.method )# calculation of the spearman correlation (i.e. Boyce index) after removing successive duplicated values
  }
  HS <- apply(interval, 1, sum)/2 # mean habitat suitability in the moving window
  HS[length(HS)]<-HS[length(HS)]-1 #Correction of the "trick" to deal with closed interval
  HS<-HS[to.keep] #exlude the NaN
  if (PEplot == T){
    plot(HS, f, xlab = "Habitat suitability",
         ylab = "Predicted/Expected ratio",col="grey",cex=0.75)
    points(HS[r],f[r],pch=19,cex=0.75)
  }
  results <- list(F.ratio = f, correlation = round(b, 3), HS = HS)
  return(results)
}


KappaRepet <-  function(Obs, Fit, TSS=FALSE)
{
  if(sum(Obs)==0) stop("\n The observed data only contains 0")
  tab <- as.data.frame(matrix(0, nrow=101, ncol=2))  ### il faut prÈciser que: "nrow = 101" sinon le data frame se remplit de "NA" au fur et ??? mesure plutÙt que de "0" (et si il y a des NA cela pose problËme plus loin).
  
  if(length(unique(Fit))==1){
    Misc<-table(as.vector(Fit) >= as.numeric(unique(Fit)), Obs) ### Robin modified here...(avant Áa plantait l???)
    if(TSS!=TRUE) a <- KappaStat(Misc)
    else a <- TSS.Stat(Misc)
    TP <- Misc[4]
    TN <- Misc[1]
    ca0 <- (TN * 100)/sum(Misc[,1])
    ca1 <- (TP * 100)/sum(Misc[,2])
    if(is.na(ca0)) ca0<-0
    if(is.na(ca1)) ca1<-0
    if(TSS!=TRUE) invisible(list(Kappa=0, CutOff=unique(Fit), TP=TP, se=ca1, TN=TN, sp=ca0))
    else invisible(list(TSS=0, CutOff=unique(Fit), TP=TP, se=ca1, TN=TN, sp=ca0))
  }
  else{
    Quant <- quantile(Fit)
    for(j in 0:100){
      Seuil <- Quant[1] + (j*((Quant[5] - Quant[1])/100))
      Misc<-table(Fit >= Seuil, Obs)
      if(TSS!=TRUE) a <- KappaStat(Misc) else a <- TSS.Stat(Misc)
      if(!is.na(a)) if(a > 0) {tab[j+1, 1] <- Seuil; tab[j+1, 2] <- a}
      rm(Misc, Seuil)
    }
    
    t <- max(tab[,2],na.rm=TRUE)
    seuil <- tab[tab[,2]==t,1]   ### Note: Ici il se peut qu'on aie plus de 1 seuil...dans ce cas le plus bas est gardÈ.
    if(t > 0) {
      Misc<-table(Fit >= seuil[1], Obs)
      TP <- Misc[4]
      TN <- Misc[1]
      ca0 <- (TN * 100)/sum(Misc[,1])
      ca1 <- (TP * 100)/sum(Misc[,2])
      if(TSS!=TRUE) invisible(list(Kappa = t, CutOff = seuil[1], TP = TP, se = ca1, TN = TN, sp = ca0))
      else invisible(list(TSS = t, CutOff = seuil[1], TP = TP, se = ca1, TN = TN, sp = ca0))
    }
    else {
      if(TSS!=TRUE) invisible(list(Kappa = 0, CutOff = 0, TP = 0, se = 0, TN = 0, sp = 0))
      else invisible(list(TSS = 0, CutOff = 0, TP = 0, se = 0, TN = 0, sp = 0))
    }
  }
}

TSS.Stat<-function(Misc)
{
  if(dim(Misc)[1]==1){
    if(row.names(Misc)[1]=="FALSE") Misc<-rbind(Misc, c(0,0))
    else {
      a<-Misc
      Misc<-c(0,0)
      Misc<-rbind(Misc, a)
      
    }
  }
  n <- sum(Misc)
  a <- Misc[1,1]
  b <- Misc[1,2]
  c <- Misc[2,1]
  d <- Misc[2,2]
  sens<-a/(a+c)
  spec<-d/(b+d)
  K <- (sens + spec) - 1        #TSS
  return(K)
}

KappaStat <-  function(Misc)
{
  if(dim(Misc)[1]==1){
    if(row.names(Misc)[1]=="FALSE") Misc <- rbind(Misc, c(0,0))
    else{
      a <- Misc
      Misc <- c(0,0)
      Misc <- rbind(Misc, a)
    }
  }
  n <- sum(Misc)
  n.1 <- sum(Misc[,1])
  n.2 <- sum(Misc[,2])
  n1. <- sum(Misc[1,])
  n2. <- sum(Misc[2,])
  Po <- (1/n) * (Misc[1,1] + Misc[2,2])
  Pe <- ((1/n)^2) * ((as.numeric(n1.) * as.numeric(n.1)) + (as.numeric(n2.) * as.numeric(n.2)))
  K <- (Po - Pe)/(1 - Pe)
  #cat("\n Kappa=", K, "\n")
  return(K)
}

