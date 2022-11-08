library(xts)
library(zoo)
library(quantmod)
library(PortfolioAnalytics)
library(ROI)
library('knitr')
library(kableExtra)
library(tidyquant)
library(quadprog)

library("ROI")
library("ROML")
library("ROML.portfolio")
library(tidyquant)
library("DEoptim")
library("shiny")
library("shinydashboard")
library("reshape")

# Carga de datos
fechai <- "2015-08-31"
fechaf <- "2020-09-01"
periodicidad <- "monthly"


# Indice de referencia: SP500
# SP500: ^GSPC


indice <- c(tq_index("SP500")[,1])[[1]]
indice[which(indice=="BRK.B")] <- "BRK-B"

sp=c(tq_index("SP500")[,6])[[1]]
boom=sp

boom=ifelse(boom=="Information Technology","A",boom)
boom=ifelse(boom=="Communication Services","A",boom)

boom=ifelse(boom=="Financials","B",boom)

boom=ifelse(boom=="Consumer Discretionary","C",boom)
boom=ifelse(boom=="Real Estate","C",boom)
boom=ifelse(boom=="Consumer Staples","C",boom)
boom=ifelse(boom=="Utilities","C",boom)

boom=ifelse(boom=="Health Care","D",boom)

boom=ifelse(boom=="Industrials","E",boom)
boom=ifelse(boom=="Materials","E",boom)
boom=ifelse(boom=="Energy","E",boom)

boom=boom[-(which(indice=="MRNA"))]
indice=indice[-(which(indice=="MRNA"))]

boom=boom[-(which(indice=="CARR"))]
indice=indice[-(which(indice=="CARR"))]

boom=boom[-(which(indice=="DOW"))]
indice=indice[-(which(indice=="DOW"))]

boom=boom[-(which(indice=="OTIS"))]
indice=indice[-(which(indice=="OTIS"))]

boom=boom[-(which(indice=="CTVA"))]
indice=indice[-(which(indice=="CTVA"))]

boom=boom[-(which(indice=="FTV"))]
indice=indice[-(which(indice=="FTV"))]

boom=boom[-(which(indice=="IR"))]
indice=indice[-(which(indice=="IR"))]

boom=boom[-(which(indice=="HPE"))]
indice=indice[-(which(indice=="HPE"))]

boom=boom[-(which(indice=="BBWI"))]
indice=indice[-(which(indice=="BBWI"))]

boom=boom[-(which(indice=="CDAY"))]
indice=indice[-(which(indice=="CDAY"))]

boom=boom[-(which(indice=="LUMN"))]
indice=indice[-(which(indice=="LUMN"))]

boom=boom[-(which(indice=="FOXA"))]
indice=indice[-(which(indice=="FOXA"))]

boom=boom[-(which(indice=="HWM"))]
indice=indice[-(which(indice=="HWM"))]

boom=boom[-(which(indice=="BF.B"))]
indice=indice[-(which(indice=="BF.B"))]

boom=boom[-(which(indice=="LW"))]
indice=indice[-(which(indice=="LW"))]

boom=boom[-(which(indice=="OGN"))]
indice=indice[-(which(indice=="OGN"))]

boom=boom[-(which(indice=="FOX"))]
indice=indice[-(which(indice=="FOX"))]

boom=boom[-(which(indice=="UA"))]
indice=indice[-(which(indice=="UA"))]

#HWM,BF.B,LW,OGN,FOX,UA

dataframeindice = data.frame("Acciones"=indice, "Retorno Medio"=matrix(0,nrow = length(indice)),"Desviación Estandar"=matrix(0,nrow = length(indice)),"Coeficiente de Sharpe"=matrix(0,nrow = length(indice)),"Treynor"=matrix(0,nrow = length(indice)), "Medida de Sortino"=matrix(0,nrow = length(indice)), "Omega"=matrix(0,nrow = length(indice)),"Sector"=matrix(0,nrow = length(indice)))
accionindice="^GSPC"
pindice <- Ad(getSymbols(accionindice,from=fechai,to=fechaf,
                         periodicity=periodicidad,auto.assign=FALSE))
rindice <- diff(log(pindice))[-1]

for(i in 1:length(indice)){
  aux <- Ad(getSymbols(indice[i],from=fechai,to=fechaf,
                       periodicity=periodicidad,auto.assign=FALSE))
  aux <- diff(log(aux))[-1]
  dataframeindice[i,2] <- mean(aux)
  dataframeindice[i,3] <- sd(aux)
  dataframeindice[i,4] <- mean(aux)/sd(aux)
  dataframeindice[i,5] <- mean(aux)/(lm(aux~rindice)[["coefficients"]][2])
  dataframeindice[i,6] <- mean(aux)/sd(ifelse(aux<0,aux,0))
  dataframeindice[i,7] <- sum(ifelse(aux>0,aux,0))/sum(ifelse(aux<0,-aux,0))
}

dataframeindice[,8]<-boom

dataframemarkowitz <- dataframeindice[order(dataframeindice$Desviación.Estandar),]
x=0
y=0
z=0
w=0
m=0
i=0

## boom es el vector de clasificacion 

boom0<-dataframemarkowitz[,8]

while (x<1 | y<1 | z<1 | w<1 | m<1 | (((x+y+z+w+m)<8)&(x+y+z+w+m)<20)) {
  i=i+1
  
  if(boom0[i]=="A"){
    x=x+1
  }
  else if(boom0[i]=="B"){
    y=y+1
  }
  else if(boom0[i]=="C"){
    z=z+1
  }
  else if(boom0[i]=="D"){
    w=w+1
  }
  else if(boom0[i]=="E"){
    m=m+1
  }
}
dataframemarkowitz<-dataframemarkowitz[c(1:16,i),c(1,3,8)]
activosmarkowitz <- c(dataframemarkowitz[,1])

options(digits=4, width=70)
actm=kable(dataframemarkowitz,booktabs=TRUE)%>%
  kable_styling(bootstrap_options = "striped",full_width = F, position = "center", font_size = 12)

actm=dataframemarkowitz

preciosmarkowitz <- xts()

for(i in 1:length(activosmarkowitz)){
  aux <- Ad(getSymbols(activosmarkowitz[i],from=fechai, to=fechaf, periodicity=periodicidad, auto.assign = FALSE))
  aux <- na.approx(aux,na.rm =FALSE) #na.omit para no interpolar sino eliminar
  preciosmarkowitz <-cbind(preciosmarkowitz,aux)
}

colnames(preciosmarkowitz)<- activosmarkowitz
tclass(preciosmarkowitz)<-"Date"

retornosmarkowitz <- diff(log(preciosmarkowitz)) [-1]

dataframesharpe <- dataframeindice[order(dataframeindice$Coeficiente.de.Sharpe,decreasing=TRUE),]
x=0
y=0
z=0
w=0
m=0
i=0

## boom es el vector de clasificacion 

boom1<-dataframesharpe[,8]

while (x<1 | y<1 | z<1 | w<1 | m<1 | (((x+y+z+w+m)<8)&(x+y+z+w+m)<20)) {
  i=i+1
  
  if(boom1[i]=="A"){
    x=x+1
  }
  else if(boom1[i]=="B"){
    y=y+1
  }
  else if(boom1[i]=="C"){
    z=z+1
  }
  else if(boom1[i]=="D"){
    w=w+1
  }
  else if(boom1[i]=="E"){
    m=m+1
  }
}
dataframesharpe<-dataframesharpe[1:i,c(1,4,8)]
activossharpe <- c(dataframesharpe[,1])

options(digits=4, width=70)
acts=kable(dataframesharpe,booktabs=TRUE)%>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "center", font_size = 12)

acts=dataframesharpe

preciossharpe <- xts()

for(i in 1:length(activossharpe)){
  aux <- Ad(getSymbols(activossharpe[i],from=fechai, to=fechaf, periodicity=periodicidad, auto.assign = FALSE))
  aux <- na.approx(aux,na.rm =FALSE) #na.omit para no interpolar sino eliminar
  preciossharpe <-cbind(preciossharpe,aux)
}

colnames(preciossharpe)<- activossharpe
tclass(preciossharpe)<-"Date"

retornossharpe <- diff(log(preciossharpe)) [-1]

dataframetreynor <-dataframeindice[order(dataframeindice$Treynor,decreasing=TRUE),]
x=0
y=0
z=0
w=0
m=0
i=0

## boom es el vector de clasificacion 

boom2<-dataframetreynor[,8]

while (x<1 | y<1 | z<1 | w<1 | m<1 | (((x+y+z+w+m)<8)&(x+y+z+w+m)<20)) {
  i=i+1
  
  if(boom2[i]=="A"){
    x=x+1
  }
  else if(boom2[i]=="B"){
    y=y+1
  }
  else if(boom2[i]=="C"){
    z=z+1
  }
  else if(boom2[i]=="D"){
    w=w+1
  }
  else if(boom2[i]=="E"){
    m=m+1
  }
}
dataframetreynor<-dataframetreynor[1:i,c(1,5,8)]
activostreynor <- c(dataframetreynor[,1])

actt=options(digits=4, width=70)
kable(dataframetreynor,booktabs=TRUE)%>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "center", font_size = 12)

actt=dataframetreynor

preciostreynor <- xts()

for(i in 1:length(activostreynor)){
  aux <- Ad(getSymbols(activostreynor[i],from=fechai, to=fechaf, periodicity=periodicidad, auto.assign = FALSE))
  aux <- na.approx(aux,na.rm =FALSE) #na.omit para no interpolar sino eliminar
  preciostreynor <-cbind(preciostreynor,aux)
}

colnames(preciostreynor)<- activostreynor
tclass(preciostreynor)<-"Date"

retornostreynor <- diff(log(preciostreynor)) [-1]

dataframesortino <-dataframeindice[order(dataframeindice$Medida.de.Sortino,decreasing=TRUE),]
x=0
y=0
z=0
w=0
m=0
i=0

## boom es el vector de clasificacion 

boom3<-dataframesortino[,8]

while (x<1 | y<1 | z<1 | w<1 | m<1 | (((x+y+z+w+m)<8)&(x+y+z+w+m)<20)) {
  i=i+1
  
  if(boom3[i]=="A"){
    x=x+1
  }
  else if(boom3[i]=="B"){
    y=y+1
  }
  else if(boom3[i]=="C"){
    z=z+1
  }
  else if(boom3[i]=="D"){
    w=w+1
  }
  else if(boom3[i]=="E"){
    m=m+1
  }
}
dataframesortino<-dataframesortino[1:i,c(1,6,8)]

activossortino <- c(dataframesortino[,1])

options(digits=4, width=70)
actso=kable(dataframesortino,booktabs=TRUE)%>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "center", font_size = 12)

actso=dataframesortino

preciossortino <- xts()

for(i in 1:length(activossortino)){
  aux <- Ad(getSymbols(activossortino[i],from=fechai, to=fechaf, periodicity=periodicidad, auto.assign = FALSE))
  aux <- na.approx(aux,na.rm =FALSE) #na.omit para no interpolar sino eliminar
  preciossortino <-cbind(preciossortino,aux)
}

colnames(preciossortino)<- activossortino
tclass(preciossortino)<-"Date"

retornossortino <- diff(log(preciossortino)) [-1]

dataframeomega <- dataframeindice[order(dataframeindice$Omega,decreasing=TRUE),]
x=0
y=0
z=0
w=0
m=0
i=0

## boom es el vector de clasificacion 

boom4<-dataframeomega[,8]

while (x<1 | y<1 | z<1 | w<1 | m<1 | (((x+y+z+w+m)<8)&(x+y+z+w+m)<20)) {
  i=i+1
  
  if(boom4[i]=="A"){
    x=x+1
  }
  else if(boom4[i]=="B"){
    y=y+1
  }
  else if(boom4[i]=="C"){
    z=z+1
  }
  else if(boom4[i]=="D"){
    w=w+1
  }
  else if(boom4[i]=="E"){
    m=m+1
  }
}
dataframeomega<-dataframeomega[1:i,c(1,7,8)]
activosomega <- c(dataframeomega[,1])

options(digits=4, width=70)
acto=kable(dataframeomega,booktabs=TRUE)%>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "center", font_size = 12)

acto=dataframeomega

preciosomega <- xts()

for(i in 1:length(activosomega)){
  aux <- Ad(getSymbols(activosomega[i],from=fechai, to=fechaf, periodicity=periodicidad, auto.assign = FALSE))
  aux <- na.approx(aux,na.rm =FALSE) #na.omit para no interpolar sino eliminar
  preciosomega <-cbind(preciosomega,aux)
}

colnames(preciosomega)<- activosomega
tclass(preciosomega)<-"Date"

retornosomega <- diff(log(preciosomega)) [-1]

modeloMV <- function(ret){
  # Inputs
  rf <- rf
  mu <- colMeans(ret)
  cov <- cov(ret)
  activos <- names(ret)
  
  # Optimizacion sin restricciones en cortos
  if(short == 1){
    ones <- rep(1,n)
    x <- t(mu)%*%solve(cov)%*%mu
    y <- t(mu)%*%solve(cov)%*%ones
    z <- t(ones)%*%solve(cov)%*%ones
    d <- x*z - y*y
    g <- (solve(cov,ones)%*%x-solve(cov,mu)%*%y)%*%solve(d)
    h <- (solve(cov,mu)%*%z-solve(cov,ones)%*%y)%*%solve(d)
    rpmin <- min(mu)
    rpmax <- max(mu)*1.5
    nport <- 1000
    j <- seq(rpmin,rpmax, length=nport) 
    wpo <- matrix(c(0), ncol=n, nrow=nport) 
    rpo <- matrix(c(0), nrow=nport)
    sigmapo <- matrix(c(0), nrow=nport)
    wj <- 0
    cont <- 1
    for(i in 1:nport){
      wj <- g + h*j[i] 
      wpo[cont,] <- t(wj)
      rpo[cont,] <- t(wj)%*%mu
      sigmapo[cont,] <- sqrt(t(wj)%*%cov%*%wj)
      cont <- cont+1
    }
    # PMVG
    cov_inv_1 <- solve(cov, ones) 
    wpmvg <- (1/as.numeric(ones %*% cov_inv_1)) * cov_inv_1
    rpmvg <- mu%*%wpmvg
    sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)
    # Sharpe
    Er <- mu-rf 
    Z <- solve(cov,Er)  
    sumZ <- sum(Z) 
    wpt <- Z/sumZ 
    rpt <- t(wpt)%*%mu
    sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)
    wpmvg <- t(wpmvg)
    wpt <- t(wpt)
    
    MV <- list()
    MV[[1]] <- wpo
    MV[[2]] <- rpo
    MV[[3]] <- sigmapo
    MV[[4]] <- t(wpmvg)
    MV[[5]] <- rpmvg
    MV[[6]] <- sigmapmvg
    MV[[7]] <- t(wpt)
    MV[[8]] <- rpt 
    MV[[9]] <- sigmapt
    return(MV)
  }
  # Con restricciones en corto
  else {
    # FE    
    library(quadprog)
    if(min(mu) > 0){rpmin = min(mu)*1.001}
    else{rpmin = 0.00}
    rpmax <- max(mu)*0.999
    n <- length(mu)
    nport <- 1000
    j <- seq(rpmin,rpmax,length=nport)
    sigmapo <- matrix(0,nrow=nport)
    wpo <- matrix(0,nrow=nport, ncol=n)
    Amat <- t(rbind(rep(1,n),mu,diag(1,nrow=n)))
    dvec <- rep(0,n) 
    Dmat <- 2*cov
    for(i in 1:nport){
      bvec <- c(1,j[i],rep(0,n))
      result <- solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=2)
      wpo[i,] <- result$solution
      sigmapo[i,] <- sqrt(result$value)
    }
    rpo <- j
    colnames(wpo) <- c(activos)
    # PMVG
    pmvg <- cbind(sigmapo,wpo)
    pmvg.sort <- pmvg[order(pmvg[,1]),]
    pmvg.sel <- cbind(pmvg.sort[1,])
    wpmvg <- cbind(round(pmvg.sel[2:length(pmvg.sel)],6))
    rownames(wpmvg) <- c(activos)
    rpmvg <- mu%*%wpmvg
    sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)
    # Sharpe    
    sharpe_port <- (rpo-rf)/sigmapo
    sharpe <- cbind(sharpe_port,wpo)
    sharpe.sort <- sharpe[order(-sharpe[,1]),]
    sharpe.sel <- cbind(sharpe.sort[1,])
    wpt <- round(cbind(sharpe.sel[2:length(sharpe.sel)]),6)
    rownames(wpt) <- c(activos)
    rpt <- mu%*%wpt
    sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)
    
    MV <- list()
    MV[[1]] <- wpo
    MV[[2]] <- rpo
    MV[[3]] <- sigmapo
    MV[[4]] <- wpmvg
    MV[[5]] <- rpmvg
    MV[[6]] <- sigmapmvg
    MV[[7]] <- wpt
    MV[[8]] <- rpt 
    MV[[9]] <- sigmapt
    return(MV)
  }
}


##--------------------------------------------------
## Programa optimizacion de Sortino (Downside Risk)
## Teoria de portafolios - 2021-2
##-------------------------------------------------

m.sortino <- function(retornos,h){
  # Inputs
  rf <- rf
  mu <- colMeans(retornos)
  semiret <- pmin(retornos,h)
  semicov <- cov(semiret) # semi-covarianzas
  #cov <- semicov
  activos <- names(retornos)
  
  # Optimizacion sin restricciones en cortos
  if(short == 1){
    ones <- rep(1,n)
    x <- t(mu)%*%solve(semicov)%*%mu
    y <- t(mu)%*%solve(semicov)%*%ones
    z <- t(ones)%*%solve(semicov)%*%ones
    d <- x*z - y*y
    g <- (solve(cov,ones)%*%x-solve(semicov,mu)%*%y)%*%solve(d)
    h <- (solve(cov,mu)%*%z-solve(semicov,ones)%*%y)%*%solve(d)
    rpmin <- min(mu)
    rpmax <- max(mu)*1.5
    nport <- 1000
    j <- seq(rpmin,rpmax, length=nport) 
    wpo <- matrix(c(0), ncol=n, nrow=nport) 
    rpo <- matrix(c(0), nrow=nport)
    sigmapo <- matrix(c(0), nrow=nport)
    wj <- 0
    cont <- 1
    for(i in 1:nport){
      wj <- g + h*j[i] 
      wpo[cont,] <- t(wj)
      rpo[cont,] <- t(wj)%*%mu
      sigmapo[cont,] <- sqrt(t(wj)%*%semicov%*%wj)
      cont <- cont+1
    }
    # PMVG
    cov_inv_1 <- solve(semicov, ones) 
    wpmvg <- (1/as.numeric(ones %*% cov_inv_1)) * cov_inv_1
    rpmvg <- mu%*%wpmvg
    sigmapmvg <- sqrt(t(wpmvg)%*%semicov%*%wpmvg)
    # Sharpe
    Er <- mu-rf 
    Z <- solve(semicov,Er)  
    sumZ <- sum(Z) 
    wpt <- Z/sumZ 
    rpt <- t(wpt)%*%mu
    sigmapt <- sqrt(t(wpt)%*%semicov%*%wpt)
    wpmvg <- t(wpmvg)
    wpt <- t(wpt)
    
    SMV <- list()
    SMV[[1]] <- wpo
    SMV[[2]] <- rpo
    SMV[[3]] <- sigmapo
    SMV[[4]] <- t(wpmvg)
    SMV[[5]] <- rpmvg
    SMV[[6]] <- sigmapmvg
    SMV[[7]] <- t(wpt)
    SMV[[8]] <- rpt 
    SMV[[9]] <- sigmapt
    return(SMV)
  }
  # Con restricciones en corto
  else {
    # FE    
    library(quadprog)
    if(min(mu) > 0){rpmin = min(mu)*1.001}
    else{rpmin = 0.00}
    rpmax <- max(mu)*0.999
    n <- length(mu)
    nport <- 1000
    j <- seq(rpmin,rpmax,length=nport)
    sigmapo <- matrix(0,nrow=nport)
    wpo <- matrix(0,nrow=nport, ncol=n)
    Amat <- t(rbind(rep(1,n),mu,diag(1,nrow=n)))
    dvec <- rep(0,n) 
    Dmat <- 2*semicov
    for(i in 1:nport){
      bvec <- c(1,j[i],rep(0,n))
      result <- solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=2)
      wpo[i,] <- result$solution
      sigmapo[i,] <- sqrt(result$value)
    }
    rpo <- j
    colnames(wpo) <- c(activos)
    # PMVG
    pmvg <- cbind(sigmapo,wpo)
    pmvg.sort <- pmvg[order(pmvg[,1]),]
    pmvg.sel <- cbind(pmvg.sort[1,])
    wpmvg <- cbind(round(pmvg.sel[2:length(pmvg.sel)],6))
    rownames(wpmvg) <- c(activos)
    rpmvg <- mu%*%wpmvg
    sigmapmvg <- sqrt(t(wpmvg)%*%semicov%*%wpmvg)
    # Sharpe    
    sharpe_port <- (rpo-rf)/sigmapo
    sharpe <- cbind(sharpe_port,wpo)
    sharpe.sort <- sharpe[order(-sharpe[,1]),]
    sharpe.sel <- cbind(sharpe.sort[1,])
    wpt <- round(cbind(sharpe.sel[2:length(sharpe.sel)]),6)
    rownames(wpt) <- c(activos)
    rpt <- mu%*%wpt
    sigmapt <- sqrt(t(wpt)%*%semicov%*%wpt)
    
    SMV <- list()
    SMV[[1]] <- wpo
    SMV[[2]] <- rpo
    SMV[[3]] <- sigmapo
    SMV[[4]] <- t(wpmvg)
    SMV[[5]] <- rpmvg
    SMV[[6]] <- sigmapmvg
    SMV[[7]] <- wpt
    SMV[[8]] <- rpt 
    SMV[[9]] <- sigmapt
    return(SMV)
  }
}





## ---------------------------------
## Programa evaluacion de desempeno
## Teoria de portafolios - 2021-2
## ---------------------------------

performance <- function(ret,indice){
  t <- nrow(ret)
  rport <- matrix(0,nrow=t,ncol=3)
  colnames(rport) <- c("PMVG","Sharpe","Benchmark")
  vport <- matrix(0,nrow=t,ncol=3)
  colnames(vport) <- c("PMVG","Sharpe","Benchmark")
  
  # Retornos
  # PMVG
  rpmv <- ret%*%wpmv
  rport[,1] <- rpmv
  
  #Sharpe
  rpsharpe <- ret%*%wpt
  rport[,2] <- rpsharpe 
  
  # Benchmark
  r.benchmark <- indice
  rport[,3] <- r.benchmark
  
  # Valor del portafolio
  # PMV
  port.mv <- matrix(0, nrow=t)
  port.mv[1] <- valor
  for(i in 2:t){
    port.mv[i] <- port.mv[i-1]*exp(rpmv[i-1])
  }
  vport[,1] <- port.mv
  
  # Sharpe
  port.sharpe <- matrix(0, nrow=t)
  port.sharpe[1] <- valor
  for(i in 2:t){
    port.sharpe[i] <- port.sharpe[i-1]*exp(rpsharpe[i-1])
  }
  vport[,2] <- port.sharpe
  
  # Benchmark
  v.benchmark <- matrix(0, nrow=t)
  v.benchmark[1] <- valor
  
  for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.benchmark[i-1])
  }
  vport[,3] <- v.benchmark
  
  DH <- list()
  DH[[1]] <- vport
  DH[[2]] <- rport
  return(DH)
}



## ----------------------------------
## Programa para importar los precios
## Teoria de portafolios - 2021-2
## ----------------------------------

f.precios <- function(activos,fechai,fechaf,periodicidad){
  precios <- xts()
  for(i in 1:length(activos)){
    aux <- Ad(getSymbols(activos[i],from=fechai,to=fechaf,
                         periodicity=periodicidad,auto.assign=FALSE))
    aux <- na.approx(aux,na.rm=FALSE) # Interpolación de datos con NA
    precios <- cbind(precios,aux)
  }
  colnames(precios) <- activos
  tclass(precios) <- "Date"
  return(precios)
}

##--------------------------------------------------
## Programa optimizacion de Treynor
## Teoria de portafolios - 2021-2
##-------------------------------------------------

m.treynor <- function(retornos,r.indice){
  
  MT=list()
  mu <- colMeans(retornos)
  sigma<-apply(retornos,2,sd)
  cov <- cov(retornos)
  ## Modelo de Treynor
  # Recomendaci???n: al menos 60 datos de historia
  # Para el calculo de las betas
  
  # Insumos: matriz de retornos de los activos y del indice
  
  # Calculos del modelo: beta y varianza del error
  n <- length(mu)
  betas <- matrix(0,ncol=n)
  varerror <- matrix(0,ncol=n)
  
  # Regresi???n iterativa para los n activos
  for(i in 1:n){
    modelo <- lm(retornos[,i]~r.indice)
    betas[i] <- modelo[["coefficients"]][2]
    varerror[i] <- var(modelo[["residuals"]])
  }
  
  # Calculo del coef. de Treynor (pasos 1 y 2)
  
  treynori <- (mu-rf)/betas
  matriz <- t(rbind(treynori,betas,varerror,mu,sigma))
  matriz.ord <- matriz[order(-matriz[,1]),]
  colnames(matriz.ord) <- c("Treynor","Betas","VarError","Mu","Sigma")
  
  # Paso 3: calculo de los ratios y las sumas acumuladas
  sigmam <- sd(r.indice)
  
  ratio1 <- ((matriz.ord[,4]-rf)*matriz.ord[,2])/matriz.ord[,3]
  ratio2 <- matriz.ord[,2]^2/matriz.ord[,3]
  suma1 <- cumsum(ratio1)
  suma2 <- cumsum(ratio2)
  
  tasac <- (sigmam^2*suma1)/(1+sigmam^2*suma2)
  
  diff <- matriz.ord[,1] - tasac
  cond.diff <- diff[!is.na(diff) & diff>0 ]
  n.optimo <- length(cond.diff)
  cmax <- tasac[n.optimo]
  
  zi <- (matriz.ord[,2]/matriz.ord[,3])*(matriz.ord[,1]-cmax)
  zi <- pmax(zi,0)
  
  wpot <- zi/sum(zi)
  rpot <- t(wpot)%*%mu
  sigmapot <- sqrt(t(wpot)%*%cov%*%wpot)  
  
  
  MT[[1]] <- wpot
  MT[[2]] <- rpot
  MT[[3]] <- sigmapot
  MT[[4]] <- n.optimo
  
  
  return(MT)
}





m.omega2 <- function(retornos,h){
  
  # Funcion optimizacion omega
  
  # Formulacion 2: DEoptim
  
  #install.packages("DEoptim")
  library(DEoptim)
  
  # Crear portafolio inicial
  n<-length(colMeans(retornos))
  x <- rep(1/n,n)
  h <- 0
  cov<-cov(retornos)
  
  lower <- rep(0,n) # Permitir cortos: -1
  upper <- rep(1,n)
  
  # Optimizacion
  mu<-colMeans(retornos)
  resultado <- DEoptim(f.omega,lower,upper,retornos=coredata(retornos),
                       h=h, control=c(itermax=2000,strategy=6))
  
  wpomega2 <- round(resultado[["optim"]][["bestmem"]],6)
  rpomega2 <- mu%*%wpomega2
  sigmapomega2 <- sqrt(t(wpomega2)%*%cov%*%wpomega2)
  
  PO <- list()
  PO[[1]] <- cbind(wpomega2)
  PO[[2]] <- rpomega2
  PO[[3]] <- sigmapomega2
  return(PO)
}

# Funcion de optimizacion Omega
f.omega <- function(x,retornos,h){
  
  rhpo <- retornos%*%(x/sum(x))
  omegap <- -(sum(pmax(rhpo-h,0)))/(sum(pmax(h-rhpo,0)))
  penalty <- ((1-sum(x))^2)*100
  return(omegap+penalty)
}





m.omega <- function(retornos,h){
  
  library("ROI")
  library("ROML")
  library("ROML.portfolio")
  
  short <- short
  if(short == 1){lb=-1}
  else{lb=0}
  m <- model()
  m$variable(portfolio, lb = 0) # the portfolio choice vector; 
  m$maximize( omega(portfolio) )
  opt <- optimize(m, solver="glpk", 
                  data=list(returns = coredata(retornos))) 
  wpomega <- round(opt$solution[grep("portfolio", names(opt$solution))]/
                     opt$solution[grep("z", names(opt$solution))], 6)
  activos<-colnames(retornos)
  mu<-colMeans(retornos)
  names(wpomega) <- activos
  rpomega <- mu%*%wpomega
  sigmapomega <- sqrt(t(wpomega)%*%cov%*%wpomega)
  
  PO <- list()
  PO[[1]] <- cbind(wpomega)
  PO[[2]] <- rpomega
  PO[[3]] <- sigmapomega
  return(PO)
}



# PERMORMANCE REMIX


performanceremix <- function(ret1,ret2,ret3,ret4,indice){
  
  # Inputs
  
  
  t <- nrow(ret1)
  rport <- matrix(0,nrow=t,ncol=5)
  colnames(rport) <- c("Sharpe","Treynor","Sortino","Omega","Benchmark")
  vport <- matrix(0,nrow=t,ncol=5)
  colnames(vport) <- c("Sharpe","Treynor","Sortino","Omega","Benchmark")
  
  # Retornos
  
  #Sharpe
  rpsharpe <- ret1%*%wpt
  rport[,1] <- rpsharpe 
  
  # Treynor
  rptreynor <- ret2%*%wpot
  rport[,2] <- rptreynor
  
  # Sortino
  rpsortino <- ret3%*%wpts
  rport[,3] <- rpsortino
  
  # Omega
  rpomega <- ret4%*%wpomega
  rport[,4] <- rpomega
  
  # Benchmark
  r.benchmark <- indice
  rport[,5] <- r.benchmark
  
  # Valor del portafolio
  
  # Sharpe
  port.sharpe <- matrix(0, nrow=t)
  port.sharpe[1] <- valor
  for(i in 2:t){
    port.sharpe[i] <- port.sharpe[i-1]*exp(rpsharpe[i-1])
  }
  vport[,1] <- port.sharpe
  
  # Treynor
  port.treynor <- matrix(0, nrow=t)
  port.treynor[1] <- valor
  for(i in 2:t){
    port.treynor[i] <- port.treynor[i-1]*exp(rptreynor[i-1])
  }
  vport[,2] <- port.treynor
  
  # Sortino
  port.sortino <- matrix(0, nrow=t)
  port.sortino[1] <- valor
  for(i in 2:t){
    port.sortino[i] <- port.sortino[i-1]*exp(rpsortino[i-1])
  }
  vport[,3] <- port.sortino
  
  # Omega
  port.omega <- matrix(0, nrow=t)
  port.omega[1] <- valor
  for(i in 2:t){
    port.omega[i] <- port.omega[i-1]*exp(rpomega[i-1])
  }
  vport[,4] <- port.omega
  
  # Benchmark
  v.benchmark <- matrix(0, nrow=t)
  v.benchmark[1] <- valor
  
  for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.benchmark[i-1])
  }
  vport[,5] <- v.benchmark
  
  DH <- list()
  DH[[1]] <- vport
  DH[[2]] <- rport
  return(DH)
}



# FUNCION DE PRECIOS


f.precios <- function(activos,fechai,fechaf,periodicidad){
  precios <- xts()
  for(i in 1:length(activos)){
    aux <- Ad(getSymbols(activos[i],from=fechai,to=fechaf,
                         periodicity=periodicidad,auto.assign=FALSE))
    aux <- na.approx(aux,na.rm=FALSE) # Interpolación de datos con NA
    precios <- cbind(precios,aux)
  }
  colnames(precios) <- activos
  tclass(precios) <- "Date"
  return(precios)
}


performanceremix <- function(ret1,ret2,ret3,ret4,ret5,indice){
  
  # Inputs
  
  
  t <- nrow(ret1)
  rport <- matrix(0,nrow=t,ncol=6)
  colnames(rport) <- c("Sharpe","Treynor","Sortino","Omega","PMV","Benchmark")
  vport <- matrix(0,nrow=t,ncol=6)
  colnames(vport) <- c("Sharpe","Treynor","Sortino","Omega","PMV","Benchmark")
  
  # Retornos
  
  #Sharpe
  rpsharpe <- ret1%*%wpt
  rport[,1] <- rpsharpe 
  
  # Treynor
  rptreynor <- ret2%*%wpot
  rport[,2] <- rptreynor
  
  # Sortino
  rpsortino <- ret3%*%wpts
  rport[,3] <- rpsortino
  
  # Omega
  rpomega <- ret4%*%wpomega
  rport[,4] <- rpomega
  
  #PMV
  # PMVG
  rpmv <- ret5%*%wpmv
  rport[,5] <- rpmv
  
  # Benchmark
  r.benchmark <- indice
  rport[,6] <- r.benchmark
  
  # Valor del portafolio
  
  # Sharpe
  port.sharpe <- matrix(0, nrow=t)
  port.sharpe[1] <- valor
  for(i in 2:t){
    port.sharpe[i] <- port.sharpe[i-1]*exp(rpsharpe[i-1])
  }
  vport[,1] <- port.sharpe
  
  # Treynor
  port.treynor <- matrix(0, nrow=t)
  port.treynor[1] <- valor
  for(i in 2:t){
    port.treynor[i] <- port.treynor[i-1]*exp(rptreynor[i-1])
  }
  vport[,2] <- port.treynor
  
  # Sortino
  port.sortino <- matrix(0, nrow=t)
  port.sortino[1] <- valor
  for(i in 2:t){
    port.sortino[i] <- port.sortino[i-1]*exp(rpsortino[i-1])
  }
  vport[,3] <- port.sortino
  
  # Omega
  port.omega <- matrix(0, nrow=t)
  port.omega[1] <- valor
  for(i in 2:t){
    port.omega[i] <- port.omega[i-1]*exp(rpomega[i-1])
  }
  vport[,4] <- port.omega
  
  # PMV
  port.mv <- matrix(0, nrow=t)
  port.mv[1] <- valor
  for(i in 2:t){
    port.mv[i] <- port.mv[i-1]*exp(rpmv[i-1])
  }
  vport[,5] <- port.mv
  
  # Benchmark
  v.benchmark <- matrix(0, nrow=t)
  v.benchmark[1] <- valor
  
  for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.benchmark[i-1])
  }
  vport[,6] <- v.benchmark
  
  DH <- list()
  DH[[1]] <- vport
  DH[[2]] <- rport
  return(DH)
}



short=0
rf=0
h=0

#Funciones
MV<-modeloMV(retornossharpe)
MT <- m.treynor(retornostreynor,rindice)
SMV=m.sortino(retornossortino,h)
PO<-m.omega2(retornosomega,h)
PMV<-modeloMV(retornosmarkowitz)

## Sharpe
wpt <- MV[[7]]
rpt <- MV[[8]] 
sigmapt <- MV[[9]]

## Treynor

wpot<-MT[[1]] 
rpot<-MT[[2]] 
sigmapot<-MT[[3]] 
n.optimo <- MT[[4]]

## Portafolio Tangente de Sortino
wpts <- SMV[[7]]
rpts <- SMV[[8]] 
sigmapts <- SMV[[9]]

# Omega
wpomega<-PO[[1]]
row.names(wpomega)<-activosomega
rpomega <-PO[[2]] 
sigmapomega <-PO[[3]]

# Markowitz
wpmv<-PMV[[4]] 
rpmv<-PMV[[5]]
sigmapmv<-MV[[6]] 

library(ggplot2)

df <- data.frame("activos"=c(row.names(wpt)),"pesos"=c(wpt))
act_sharpe=row.names(wpt)
pesos_sharpe=wpt
fua <- ggplot(df, aes(act_sharpe, pesos_sharpe)) 
fua <- fua + geom_bar(stat = "identity")


df <- data.frame("activos"=c(names(wpot)),"pesos"=c(wpot))
act_treynor=names(wpot)
pesos_treynor=wpot
fua2 <- ggplot(df, aes(act_treynor, pesos_treynor)) 
fua2 <- fua2 + geom_bar(stat = "identity")

df <- data.frame("activos"=c(row.names(wpts)),"pesos"=c(wpts))
act_sortino=row.names(wpts)
pesos_sortino=wpts
fua3 <- ggplot(df, aes(act_sortino, pesos_sortino)) 
fua3 <- fua3 + geom_bar(stat = "identity")

df <- data.frame("activos"=c(row.names(wpomega)),"pesos"=c(wpomega))
act_omega=row.names(wpomega)
pesos_omega=wpomega
fua4 <- ggplot(df, aes(act_omega, pesos_omega)) 
fua4 <- fua4 + geom_bar(stat = "identity")

df <- data.frame("activos"=c(row.names(wpmv)),"pesos"=c(wpmv))
act_mv=row.names(wpmv)
pesos_mv=wpmv
fua5 <- ggplot(df, aes(act_mv, pesos_mv)) 
fua5 <- fua5 + geom_bar(stat = "identity")

# Evaluacion de desempeno 
# In-sample
valor <- 100 
DH <- performanceremix(retornossharpe,retornostreynor,retornossortino,retornosomega,retornosmarkowitz, rindice)
Performance <-ts(DH[[1]],start=2015.9, frequency=12)
performanceboom=data.frame(x=seq_along(Performance[,1]),Performance)
performanceboom=melt(performanceboom, id.vars="x")
fua6=ggplot(performanceboom,aes(x = x, y = value, color = variable))+geom_line()


# --------------------------------------------------------------
# Tabla resumen:
rp.hist <- DH[[2]]
retorno <- round(rbind(mean(rp.hist[,1]),mean(rp.hist[,2]),mean(rp.hist[,3]),mean(rp.hist[,4]),mean(rp.hist[,5]),mean(rp.hist[,6])),4)
riesgo <- round(rbind(sd(rp.hist[,1]),sd(rp.hist[,2]),sd(rp.hist[,3]),sd(rp.hist[,4]),sd(rp.hist[,5]),sd(rp.hist[,6])),4)
sharpe <- round(rbind((retorno[1]-rf)/riesgo[1],
                      (retorno[2]-rf)/riesgo[2],
                      (retorno[3]-rf)/riesgo[3],
                      (retorno[4]-rf)/riesgo[4],
                      (retorno[5]-rf)/riesgo[5],
                      (retorno[6]-rf)/riesgo[6]),4)
treynor <-round(rbind(mean(rp.hist[,1])/(lm(rp.hist[,1]~rp.hist[,6])[["coefficients"]][2]),
                      mean(rp.hist[,2])/(lm(rp.hist[,2]~rp.hist[,6])[["coefficients"]][2]),
                      mean(rp.hist[,3])/(lm(rp.hist[,3]~rp.hist[,6])[["coefficients"]][2]),
                      mean(rp.hist[,4])/(lm(rp.hist[,4]~rp.hist[,6])[["coefficients"]][2]),
                      mean(rp.hist[,5])/(lm(rp.hist[,5]~rp.hist[,6])[["coefficients"]][2]),NA),4)
sortino<-round(rbind(mean(rp.hist[,1])/sd(ifelse(rp.hist[,1]<0,rp.hist[,1],0)),
                     mean(rp.hist[,2])/sd(ifelse(rp.hist[,2]<0,rp.hist[,2],0)),
                     mean(rp.hist[,3])/sd(ifelse(rp.hist[,3]<0,rp.hist[,3],0)),
                     mean(rp.hist[,4])/sd(ifelse(rp.hist[,4]<0,rp.hist[,4],0)),
                     mean(rp.hist[,5])/sd(ifelse(rp.hist[,5]<0,rp.hist[,5],0)),
                     mean(rp.hist[,6])/sd(ifelse(rp.hist[,6]<0,rp.hist[,6],0))),4)
omega<-round(rbind(sum(ifelse(rp.hist[,1]>0,rp.hist[,1],0))/sum(ifelse(rp.hist[,1]<0,-rp.hist[,1],0)),
                   sum(ifelse(rp.hist[,2]>0,rp.hist[,2],0))/sum(ifelse(rp.hist[,2]<0,-rp.hist[,2],0)),
                   sum(ifelse(rp.hist[,3]>0,rp.hist[,3],0))/sum(ifelse(rp.hist[,3]<0,-rp.hist[,3],0)),
                   sum(ifelse(rp.hist[,4]>0,rp.hist[,4],0))/sum(ifelse(rp.hist[,4]<0,-rp.hist[,4],0)),
                   sum(ifelse(rp.hist[,5]>0,rp.hist[,5],0))/sum(ifelse(rp.hist[,5]<0,-rp.hist[,5],0)),
                   sum(ifelse(rp.hist[,6]>0,rp.hist[,6],0))/sum(ifelse(rp.hist[,6]<0,-rp.hist[,6],0))),4)
Resumen <- cbind(sharpe,treynor,sortino,omega,retorno,riesgo)
colnames(Resumen) <- c("Sharpe","Treynor","Sortino","Omega","Retorno","Riesgo")
rownames(Resumen) <- c("Sharpe", "Treynor","Sortino","Omega","PMV","Indice")
resumeninsample<-Resumen
options(digits=4, width=70)
tabin=kable(resumeninsample,booktabs=TRUE)%>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "center", font_size = 12)
tabin=resumeninsample
tabin=data.frame("Portafolios"=c(row.names(tabin)),tabin)
##outsample
fechaifm <- '2020-09-01'
fechaffm <- '2021-09-02'

preciossharpefm <- f.precios(activossharpe,fechaifm,fechaffm,periodicidad)
retornossharpefm <- diff(log(preciossharpefm))[-1,]
preciostreynorfm <- f.precios(activostreynor,fechaifm,fechaffm,periodicidad)
retornostreynorfm <- diff(log(preciostreynorfm))[-1,]
preciossortinofm <- f.precios(activossortino,fechaifm,fechaffm,periodicidad)
retornossortinofm <- diff(log(preciossortinofm))[-1,]
preciosomegafm <- f.precios(activosomega,fechaifm,fechaffm,periodicidad)
retornosomegafm <- diff(log(preciosomegafm))[-1,]
preciosmarkowitzfm <- f.precios(activosmarkowitz,fechaifm,fechaffm,periodicidad)
retornosmarkowitzfm <- diff(log(preciosmarkowitzfm))[-1,]

indicefm <- f.precios(accionindice,fechaifm,fechaffm,periodicidad)
r.indicefm <- diff(log(indicefm))[-1,]

DFM <- performanceremix(retornossharpefm,retornostreynorfm,retornossortinofm,retornosomegafm,retornosmarkowitzfm,r.indicefm)
performancefm <-ts(DFM[[1]],start=2020.9, frequency=12)

performanceboomfm=data.frame(x=seq_along(performancefm[,1]),performancefm)
performanceboomfm=melt(performanceboomfm, id.vars="x")
fua7=ggplot(performanceboomfm,aes(x = x, y = value, color = variable))+geom_line()


# Tabla resumen:
rp.fm <- DFM[[2]]
retornofm <- round(rbind(mean(rp.fm[,1]),mean(rp.fm[,2]),mean(rp.fm[,3]),mean(rp.fm[,4]),mean(rp.fm[,5]),mean(rp.fm[,6])),4)
riesgofm <- round(rbind(sd(rp.fm[,1]),sd(rp.fm[,2]),sd(rp.fm[,3]),sd(rp.fm[,4]),sd(rp.fm[,5]),sd(rp.fm[,6])),4)
sharpefm <- round(rbind((retornofm[1]-rf)/riesgofm[1],
                        (retornofm[2]-rf)/riesgofm[2],
                        (retornofm[3]-rf)/riesgofm[3],
                        (retornofm[4]-rf)/riesgofm[4],
                        (retornofm[5]-rf)/riesgofm[5],
                        (retornofm[6]-rf)/riesgofm[6]),4)
treynorfm <-round(rbind(mean(rp.fm[,1])/(lm(rp.fm[,1]~rp.fm[,6])[["coefficients"]][2]),
                        mean(rp.fm[,2])/(lm(rp.fm[,2]~rp.fm[,6])[["coefficients"]][2]),
                        mean(rp.fm[,3])/(lm(rp.fm[,3]~rp.fm[,6])[["coefficients"]][2]),
                        mean(rp.fm[,4])/(lm(rp.fm[,4]~rp.fm[,6])[["coefficients"]][2]),
                        mean(rp.fm[,5])/(lm(rp.fm[,5]~rp.fm[,6])[["coefficients"]][2]),NA),4)
sortinofm<-round(rbind(mean(rp.fm[,1])/sd(ifelse(rp.fm[,1]<0,rp.fm[,1],0)),
                       mean(rp.fm[,2])/sd(ifelse(rp.fm[,2]<0,rp.fm[,2],0)),
                       mean(rp.fm[,3])/sd(ifelse(rp.fm[,3]<0,rp.fm[,3],0)),
                       mean(rp.fm[,4])/sd(ifelse(rp.fm[,4]<0,rp.fm[,4],0)),
                       mean(rp.fm[,5])/sd(ifelse(rp.fm[,5]<0,rp.fm[,5],0)),
                       mean(rp.fm[,6])/sd(ifelse(rp.fm[,6]<0,rp.fm[,6],0))),4)
omegafm<-round(rbind(sum(ifelse(rp.fm[,1]>0,rp.fm[,1],0))/sum(ifelse(rp.fm[,1]<0,-rp.fm[,1],0)),
                     sum(ifelse(rp.fm[,2]>0,rp.fm[,2],0))/sum(ifelse(rp.fm[,2]<0,-rp.fm[,2],0)),
                     sum(ifelse(rp.fm[,3]>0,rp.fm[,3],0))/sum(ifelse(rp.fm[,3]<0,-rp.fm[,3],0)),
                     sum(ifelse(rp.fm[,4]>0,rp.fm[,4],0))/sum(ifelse(rp.fm[,4]<0,-rp.fm[,4],0)),
                     sum(ifelse(rp.fm[,5]>0,rp.fm[,5],0))/sum(ifelse(rp.fm[,5]<0,-rp.fm[,5],0)),
                     sum(ifelse(rp.fm[,6]>0,rp.fm[,6],0))/sum(ifelse(rp.fm[,6]<0,-rp.fm[,6],0))),4)
Resumenfm <- cbind(sharpefm,treynorfm,sortinofm,omegafm,retornofm,riesgofm)
colnames(Resumenfm) <- c("Sharpe","Treynor","Sortino","Omega","Retorno","Riesgo")
rownames(Resumenfm) <- c("Sharpe", "Treynor","Sortino","Omega","PMV","Indice")
resumenoutsample<-Resumenfm
options(digits=4, width=70)
tabout=kable(resumenoutsample,booktabs=TRUE)%>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "center", font_size = 12)
tabout=resumenoutsample
tabout=data.frame("Portafolios"=c(row.names(tabout)),tabout)

riesgo_in=resumeninsample[,6]
retorno_in=resumeninsample[,5]
wow=data.frame(riesgo_in,retorno_in)

fua8=ggplot(wow,aes(x=riesgo_in,y=retorno_in))+geom_point()+geom_text(aes(label=row.names(resumeninsample)),size=3,nudge_x = 0.0005)


riesgo_out=resumenoutsample[,6]
retorno_out=resumenoutsample[,5]
wow2=data.frame(riesgo_out,retorno_out)

fua9=ggplot(wow2,aes(x=riesgo_out,y=retorno_out))+geom_point()+geom_text(aes(label=row.names(resumenoutsample)),size=3,nudge_x = 0.0005)


library(shiny)
library(shinydashboard)



header <- dashboardHeader()
sidebar <- dashboardSidebar(
  
  sidebarMenu(
    menuItem("Información y análisis",
             tabName = "ana",menuSubItem("Introducción y metodología",tabName = "im"),
             menuSubItem("Análisis de la selección de activos",tabName = "asa"),
             menuSubItem("Análisis de desempeño",tabName = "ads"),
             menuSubItem("Recomendación de inversión",tabName = "reci")
             ),
    
    menuItem("Sharpe",
             tabName = "sharpe"),
    menuItem("Treynor",
             tabName = "treynor"),
    menuItem("Omega",
             tabName = "omega"),
    menuItem("Sortino",
             tabName = "sortino"),
    menuItem("PMV",
             tabName = "PMV"),
    menuItem("Performance in-sample",
             tabName = "Performanceinsample"),
    menuItem("Performance out-sample",
             tabName = "Performanceoutsample")
    
  )
  
  
  
)

body <- dashboardBody(
  
  tabItems(   
    tabItem(tabName = "im",
            
            h1("Introducción",align="center",font="bold"),
            
            p("El siguiente trabajo tiene como finalidad generar una recomendación de inversión mediante el uso de las diferentes medidas de desempeño y los distintos modelos para la construcción de portafolios óptimos que se han visto a lo largo del curso. Se recomendará invertir en un portafolio que estará conformado por algunos activos pertenecientes al índice SP500;
            
            
            este se compone por las 505 empresas más grandes que cotizan en la bolsas de NYSE o NASDAQ."),
            
            
            h2("Metodología",align="center",font="bold"),
            
            
            p("1. Información Histórica: La construcción del portafolio a recomendar estará dada por la evaluación (in-sample) de los datos históricos de los activos entre el 1 de septiembre de 2015 y el 1 de septiembre de 2020. Los resultados obtenidos en este periodo serán comparados con una evaluación out-sample del comportamiento de esos mismos activos entre el 1 de septiembre de 2020 y el 1 de septiembre de 2021. Cabe aclarar que toda la información será tratada con periodicidad mensual."),

            
            p("2. Sectores: Con el fin de garantizar diversificación, se tendrá en cuenta el sector económico al que pertenece cada una de las empresas consideradas. A cada sector se le asignará una letra para representarlo; esta aparecerá al lado de cada uno de los activos que compongan los portafolios. Los sectores considerados, y sus respectivas convenciones, son los siguientes:"),

              p("     - Tecnologia y comunicaciones = A"),
            
            
              p("     - Financiero = B"),
            
            
              p("     - Consumo Basico y Discrecional = C"),
            
            
              p("     - Salud y Cuidado = D"),
            
            
              p("     - Industria y Materiales = E "),
            
            

           p("3. Restricciones: Para la construcción de los portafolios se van a manejar 2 restricciones:"),
           
            
              p("     a. Para los portafolios óptimos sólo se puede tomar posiciones en largo. Esto quiere decir que no se admiten pesos negativos en ninguno de los activos que los componen"),
           
           
              p("     b. Cada portafolio debe estar compuesto por al menos una empresa de cada sector"),
           
           
              p("     c. Cada portafolio tendrá un mínimo de 8 activos y un máximo de 20"),
           
            
            p("Además, para la recomendación de inversión se manejará la restricción de que el portafolio a recomendar debe superar el desempeño del Benchmark"),
            
           
            p("4. Portafolios: En total se construiran 5 portafolios óptimos, cada uno con base en los resultados de una medida de desempeño asociada a un modelo. Los portafolios a construir serán los siguientes:"),
            
           
              p("     a. Portafolio de Sharpe: Se construirá con los activos que tengan mayor ratio de Sharpe"),
           
           
              p("     b. Portafolio de Treynor: Se construirá con los activos que tengan mayor coeficiente de Treynor"),
           
           
              p("     c. Portafolio de Sortino; Se construirá con los activos que tengan mayor medida de Sortino"),
           
           
              p("     d. Portafolio de Omega: Se construirá con los activos que tengan mayor Omega"),
           
           
              p("     e. Portafolio Mínima Varianza: Se construirá el PMV según la formulación de Markowitz, usando los activos con menor desviación estándar como filtro inicial."),
            
           
            p("5. Selección de Activos: Se parte de una base inicial de 487 activos (18 activos presentes en el SP500 fueron descartados por no tener valores para todo el periodo a estudiar) a los que se les calcularán diferentes medidas de desempeño.Para cada portafolio se organizarán los datos de acuerdo a sus ratios (de mayor a menor) y se seleccionará el número de activos a partir del cual se cumplen con la restricciones de los sectores y de la cantidad de acciones. Los activos más abajo de ese corte no serán tenidos en cuenta para la construcción del portafolio."),
            
           
            p("6. Portafolios Óptimos: Una vez seleccionados los activos para cada portafolio, se aplicará el modelo correspondiente para encontrar los pesos que tiene cada activo para conformar un portafolio óptimo. Estos portafolios resultantes son a los que se les evaluará el proceso de evaluación y comparación para determinar la recomendación de inversión"),
            
           
            p("7. Tasa Libre de Riesgo: Para el siguiente ejercicio, se asumirá que la rf del mercado es igual a cero"),
            
           
            p("8. Benchmark: Los portafolios óptimos serán comparados a su vez con un benchmark. Este se asumirá como el desempeño de la acción del índice escogido")
    ),
    tabItem(tabName = "asa",
            h2("Análisis de la selección de activos"),
            p("Lo primero que llama la atención tras realizar la selección de activos es la escasa presencia de acciones pertenecientes al top 10 del SP500; en el global de todas las medidas de desempeño, sólo MSFT hace presencia reptetitiva en los activos seleccionados como miembro de este grupo (NVDA y BRK-B aparecen sólo en un portafolio). Esto, en últimas, refleja que no siempre las acciones de las empresas más grandes son las que tienen mejores ratios, y, por ende, no siempre son las que mejor relación tienen entre sus retornos y su riesgo (o sus pérdidas)."),
            p("Pasando a un análisis un poco más individual, se observa que los niveles de diversificación (entendiéndola como la pluralidad en los sectores económicos) varían considerablemente dependiendo de la medida de desempeño.Mientras que, por ejemplo, en la clasificación de Omega hay una selección diversa, en la de Treynor y la de Markowitz hay una concentración importante de activos pertenecientes al consumo básico y discrecional, o en la de Sortino hay una fuerte presencia del sector de tecnología y comunicaciones. Estos niveles de diversificación cambiantes permiten intuir que los sectores económicos pueden llegar a influir en el comportamiento de las empresas ante los ratios, siendo este el motivo por el que se presenta esta variación de acuerdo a la medida de desempeño que se tome como referencia."),
            p("Finalmente, cabe destacar que algunos activos fueron una constante en algunas de las clasificaciones. Activos como MSFT, MSCI, GNRC, PYPL, POOL o ADBE están presentes simultaneamente en las clasificaciones de Sharpe, Sortino y Omega, permitiendo observar que son acciones con una importante solidez en materia de ratios y, por ende, en cuanto a su relación riesgo/retorno. Curiosamente, es en las medidas de Treynor y Markowitz en las que se observan activos inéditos y no los que se repiten en las otras."),
            p("Con todo lo anterior presente, se puede dar a la optimización de los activos seleccionados de acuerdo a los modelos planteados para cada una de las medidas de desempeño. Serán estos portafolios óptimos los que se compararán, evualuarán y analizarán para posteriormente dar la recomendación final.")
            
            ),
    
    
    tabItem(tabName = "ads",
            
            h2("Análisis de desempeño"),
            
            p("El primer aspecto a tener en cuenta son los pesos de los activos para los portafolios óptimos. Dada la restricción de tomar sólamente posiciones en largo, no se observa ningún peso negativo; sin embargo, para los portafolios de Sharpe, Sortino y PMV si se tienen algunos activos con peso 0. Esto refleja que a pesar de que los activos fueron seleccionados por tener altos ratios, ello no implica que al interactuar con otros activos vayan a aportar positivamente al objetivo conjunto del portafolio (minimizar la varianza, minimizar las pérdidas, etc)."),
            
            p("Pasando al análisis netamente de los portafolios, la comparación del desempeño in-sample con el out-sample arroja cosas interesantes. En primer lugar, se observa que los portafolios óptimos de de Treynor y de Mínima Varianza superan el desempeño del Benchmark para el periodo in-sample, pero no son capaces de hacerlo para el periodo out-sample, lo que resulta es su descarte automático de cara al análisis y a la recomendación."), 
            
            p("Los demás portafolios si lograron superar el desempeño del Benchmark para los 2 periodos, en general obteniendo muy buenos resultados. Los 3 resultan ser opciones atractivas de inversión, con desempeños superiores al 400% para el periodo in-sample y superiores al 130% para el out-sample; sin embargo, 
es el portafolio óptimo de Sortino el que ha mostrado mejor desempeño para ambos periodos. Eso sí, todo lo relativo al desempeño debe matizarse con los resultados de los ratios, que nos muestran que los 3 portafolios resaltan por su buena relación riesgo-retorno y retornos(+)/retornos(-), siendo Sortino claramenete superior para el periodo out-sample y estando muy parejo con Omega para el in-sample."),
            
            p("Pero más allá de la evaluación de los ratios y desempeños de portafolio, una de las grandes conclusiones que arroja el análisis de los portafolios óptimos es que la capitalización bursátil no es indicador determinante de rentabilidad. Como se evidenció, con una selección de activos que no necesariamente eran lo más grandes y con diversos métodos de optimización, se logró construir portafolios con mejor desempeño al Benchmark (o sea, al SP500) para distintos periodos de prueba. Tiende a asociarse que los activos que lideran este tipo de índices son los que tienen mejor desempeño y menores niveles de riesgo, sin embargo, asumir esto y optar por invertir sólo por su capitalización bursátil puede terminar limitando la rentabilidad al excluir del portafolio a activos más pequeños pero con mucho mejores ratios (lo que a la larga traería un mejor desempeño si se optimizan correctamente).")
            
            
            ),
    
    tabItem(tabName = "reci",
            
            h2("Recomendación de inversión"),
            
            p("Tras la comparación de los diversos portafolios óptimos generados a partir de diferentes medidas de desempeño, se recomienda invertir en el Portafolio de Sortino que se construyó previamente, partiendo de los activos del SP500. Las acciones que componen este portafolio son PYPL, GNRC, MSFT, MSCI, ADBE, SNPS, CDNS, NEE y DHR; como se comentó anteriormente, estos activos fueron elegidos al ser los que tenían una mejor medida de Sortino para el periodo de evaluación.La conclusión de recomendar este portafolio se alcanza teniendo en cuenta diferentes factores. Tras la comparación de los diversos portafolios óptimos generados a partir de diferentes medidas de desempeño, se recomienda invertir en el Portafolio de Sortino que se construyó previamente, partiendo de los activos del SP500. Las acciones que componen este portafolio son PYPL, GNRC, MSFT, MSCI, ADBE, SNPS, CDNS, NEE y DHR; como se comentó anteriormente, estos activos fueron elegidos al ser los que tenían una mejor medida de Sortino para el periodo de evaluación.La conclusión de recomendar este portafolio se alcanza teniendo en cuenta diferentes factores. "),
            p("En primer lugar, para la evaluación de desempeño in-sample fue el portafolio más rentable, llegando a alcanzar un 500% del valor inicial del portafolio y superando ampliamente el desempeño del Benchmark. Por si estos 5 años de muestra no son lo suficiente convincentes, en la evaluación out-sample (que además contempla un importante periodo de incertidumbre en los mercados) confirmó este patrón al ser nuevamente el portafolio más rentable (alrededor de un 147% del valor inicial). En este orden de ideas, el Portafolio de Sortino no sólo se muestra como una alternativa mucho más rentable a invertir en el SP500, sino que además supera de forma consistente la rentabilidad de los otros portafolios óptimos considerados."),
            p("En adición, sus ratios  lo consolidan como una opción atractiva para invertir. En primer lugar, sus ratios de Sharpe muestran la existencia de una relación aceptable entre el riesgo y el retorno, siendo el  portafolio con mejor relación entre ambos para el periodo out-sample y el tercero para el in-sample (aun así es muy superior el Benchmark). Sin embargo, es la relación entre sus retornos positivos y los negativos la que más lo impulsa como una inversión atractiva. Para el periodo in-sample, fue el portafolio con mayor medida de Sortino y tuvo una buena medida Omega. Para el periodo out-sample tiene el mejor Omega y la mejor medida de Sortino. El análisis de estos ratios nos da a entender que el Portafolio de Sortino tiene menos riesgo de caer en pérdidas, pero que además compensa muy bien ese riesgo que se asume gracias a una relación favorable con los retornos."),
            p("Es cierto que el portafolio de Sortino es el más riesgoso de todos los portafolios óptimos observados (puede verse en los planos-riesgo retorno). Sin embargo, el análisis de sus ratios muestra que este mayor riesgo no está necesariamente ligado a pérdidas (se evidencia con sus buenos ratios de Omega y Sortino) y que en todo caso está bien compensado por sus retornos (se evidencia con su buen ratio de Sharpe). Si a lo anterior se le suma que para dos periodos distintos fue capaz de superar el desempeño de los demás, que es el portafolio con mayores retornos medios para ambos periodos, y que en general la formulación de Sortino corresponde a desarrollos más recientes y actualizados en materia de teoría portafolios, se considera a este portafolio como la mejor alternativa de inversión.")
            
            
            
            ),
    
    tabItem(tabName = "sharpe", 
            
            h2("Pesos activos portafolio Sharpe"),
            
            
            plotOutput("wpt"),
            
            h2("Características activos portafolio Sharpe"),
            
            tableOutput("acts")
            #renderPlot("plt")
    ),
    tabItem(tabName = "treynor", 
            
            h2("Pesos activos portafolio Treynor"),
            
            
            plotOutput("wpot"),
            
            h2("Características activos portafolio Treynor"),
            
            tableOutput("actt")
            #renderPlot("plt")
    ),
    
    tabItem(tabName = "omega", 
            
            h2("Pesos activos portafolio Omega"),
            
            
            plotOutput("wpomega"),
            
            h2("Características activos portafolio Omega"),
            
            tableOutput("acto")
            
    ),
    tabItem(tabName = "sortino", 
            h2("Pesos activos portafolio Sortino"),
            
            
            plotOutput("wpts"),
            
            h2("Características activos portafolio Sortino"),
            
            tableOutput("actso")
            #renderPlot("plt")
    ),
    tabItem(tabName = "PMV", 
            h2("Pesos activos portafolio PMV"),
            
            
            plotOutput("wpmv"),
            
            h2("Características activos portafolio PMV"),
            
            tableOutput("actm")
            #renderPlot("plt")
    ),
    tabItem(tabName = "Performanceinsample", 
            h2("Performance dentro de muestra"),
            
            
            plotOutput("pefin"),
            
            h2("Tabla de resumen in-sample"),
            
            tableOutput("tabin"),
            
            h2("Plano Riesgo Retorno In-sample"),
            
            
            plotOutput("rretin")
            #renderPlot("plt")
    ),
    tabItem(tabName = "Performanceoutsample", 
            h2("Performance fuera de muestra"),
            
            
            plotOutput("perfout"),
            
            h2("Tabla de resumen out-sample"),
            
            tableOutput("tabout"),
            
            h2("Plano Riesgo Retorno Out-sample"),
            
            plotOutput("rretout")
            #renderPlot("plt")
    )
    
  )  
)









ui <- dashboardPage(header, sidebar, body, skin = "black")




# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$wpt <- renderPlot({fua})
  output$wpot <- renderPlot({fua2})
  output$wpts <- renderPlot({fua3})
  output$wpomega <- renderPlot({fua4})
  output$wpmv <- renderPlot({fua5})
  output$pefin <- renderPlot({fua6})
  output$perfout <- renderPlot({fua7})
  output$actm <- renderTable({actm})
  output$actt <- renderTable({actt})
  output$acts <- renderTable({acts})
  output$acto <- renderTable({acto})
  output$actso <- renderTable({actso})
  output$tabin <- renderTable({tabin})
  output$tabout <- renderTable({tabout})
  output$rretin <- renderPlot({fua8})
  output$rretout <- renderPlot({fua9})
}



# Run the application 
shinyApp(ui = ui, server = server)





