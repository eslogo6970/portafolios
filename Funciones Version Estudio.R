
# Funciones


##-------------------------------------------
## Programa optimizacion Media-Varianza (MV)
## Teoria de portafolios - 2021-2
## Con cortos (permite pesos negativos): 1
## Sin cortos (no pesos negativos): 0
##-------------------------------------------

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


# PERMORMANCE REMIX 2


performanceremix2 <- function(ret1,ret2,ret3,ret4,indice){
  
  # Inputs
  
  
  t <- nrow(ret1)
  rport <- matrix(0,nrow=t,ncol=5)
  colnames(rport) <- c("Sharpe","Black-Litterman","Sortino","Omega","Benchmark")
  vport <- matrix(0,nrow=t,ncol=5)
  colnames(vport) <- c("Sharpe","Black-Litterman","Sortino","Omega","Benchmark")
  
  # Retornos
  
  #Sharpe
  rpsharpe <- ret1%*%wpt
  rport[,1] <- rpsharpe 
  
  # Black-Litterman
  rpbl <- ret2%*%wpbl
  rport[,2] <- rpbl
  
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
  
  # Black-Litterman
  port.bl <- matrix(0, nrow=t)
  port.bl[1] <- valor
  for(i in 2:t){
    port.bl[i] <- port.bl[i-1]*exp(rpbl[i-1])
  }
  vport[,2] <- port.bl
  
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



performanceremix3 <- function(ret1,ret2,indice){
  
  # Inputs
  
  
  t <- nrow(ret1)
  rport <- matrix(0,nrow=t,ncol=3)
  colnames(rport) <- c("Sharpe","Black-Litterman","Benchmark")
  vport <- matrix(0,nrow=t,ncol=3)
  colnames(vport) <- c("Sharpe","Black-Litterman","Benchmark")
  
  # Retornos
  
  #Sharpe
  rpsharpe <- ret1%*%wpt
  rport[,1] <- rpsharpe 
  
  # Black-Litterman
  rpbl <- ret2%*%wpbl
  rport[,2] <- rpbl
  
  # Benchmark
  r.benchmark <- indice
  rport[,3] <- r.benchmark
  
  # Valor del portafolio
  
  # Sharpe
  port.sharpe <- matrix(0, nrow=t)
  port.sharpe[1] <- valor
  for(i in 2:t){
    port.sharpe[i] <- port.sharpe[i-1]*exp(rpsharpe[i-1])
  }
  vport[,1] <- port.sharpe
  
  # Black-Litterman
  port.bl <- matrix(0, nrow=t)
  port.bl[1] <- valor
  for(i in 2:t){
    port.bl[i] <- port.bl[i-1]*exp(rpbl[i-1])
  }
  vport[,2] <- port.bl
  
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


performanceremix4 <- function(ret1,ret2,ret3,indice){
  
  # Inputs
  
  
  t <- nrow(ret1)
  rport <- matrix(0,nrow=t,ncol=4)
  colnames(rport) <- c("Sharpe","Black-Litterman","Crystal Ball","Benchmark")
  vport <- matrix(0,nrow=t,ncol=4)
  colnames(vport) <- c("Sharpe","Black-Litterman","Crystal Ball","Benchmark")
  
  # Retornos
  
  # Sortino
  rpsortino <- ret3%*%wpts
  rport[,1] <- rpsortino
  
  # Black-Litterman
  rpbl <- ret2%*%wpbl
  rport[,2] <- rpbl
  
  #Crystal Ball
  rcb<-ret3%*%wcb
  rport[,3]<-rcb
  
  # Benchmark
  r.benchmark <- indice
  rport[,4] <- r.benchmark
  
  # Valor del portafolio
  
  # Sortino
  port.sortino <- matrix(0, nrow=t)
  port.sortino[1] <- valor
  for(i in 2:t){
    port.sortino[i] <- port.sortino[i-1]*exp(rpsortino[i-1])
  }
  vport[,1] <- port.sortino
  
  # Black-Litterman
  port.bl <- matrix(0, nrow=t)
  port.bl[1] <- valor
  for(i in 2:t){
    port.bl[i] <- port.bl[i-1]*exp(rpbl[i-1])
  }
  vport[,2] <- port.bl
  
  # Crystal Ball
  
  port.cb <- matrix(0, nrow=t)
  port.cb[1] <- valor
  for(i in 2:t){
    port.cb[i] <- port.cb[i-1]*exp(rcb[i-1])
  }
  vport[,3] <- port.cb
  
  # Benchmark
  v.benchmark <- matrix(0, nrow=t)
  v.benchmark[1] <- valor
  
  for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.benchmark[i-1])
  }
  vport[,4] <- v.benchmark
  
  DH <- list()
  DH[[1]] <- vport
  DH[[2]] <- rport
  return(DH)
}


performancefeatdonomar <- function(retos){
  
  # Inputs
  t <- nrow(retos)
  vport <- matrix(0,nrow=t,ncol=ncol(retos))
  colnames(vport) <- colnames(retos)
  vport[1,] <- valor
  
  for(i in 2:t){
    vport[i,] <- vport[i-1,]*exp(retos[i-1,])
  }
  
  DH <- list()
  DH[[1]] <- vport
  DH[[2]] <- retos
  return(DH)
}


tablaresumen <- function(omg){
  
  rf=rf
  er=t(colMeans(omg))
  rr=t(apply(omg,2,sd))
  sharp=(er-rf)/rr
  sorti=t(apply(omg,2,function(x) mean(x)/sd(ifelse(x<0,x,0))))
  ome=t(apply(omg,2,function(x) sum(ifelse(x>0,x,0))/-sum(ifelse(x<0,x,0))))
  
  te=matrix(0,ncol=ncol(ome),nrow = nrow(ome))
  for (colu in (1:ncol(omg))) {
    te[colu]=sd(omg[,colu]-omg[,ncol(omg)])*sqrt(12)}
  alpha=er-er[length(er)]
  ir=alpha/te
    
  resu=t(rbind(er,rr,sharp,sorti,ome,te,alpha,ir))  
  rownames(resu)<-colnames(omg)
  colnames(resu)<-c("Retorno","Riesgo","Sharpe","Sortino","Omega","Tracking Error","Alpha","Information Ratio")
  
  return(resu)
}


# BLACK-LITTERMAN

c.blackl <- function(retornos,rf,tau,q,P){
  
  
  mu <- (colMeans(retornos))
  cov <- cov(retornos)
  #rpt <- modeloMV(retornos)[[8]]
  #wpt <- modeloMV(retornos)[[7]]
  #sigmapt <- modeloMV(retornos)[[9]]
  
  wpt <- m.sortino(retornos,h)[[7]]
  rpt <- m.sortino(retornos,h)[[8]] 
  sigmapt <- m.sortino(retornos,h)[[9]]
  
  lambda <- (rpt-rf)/sigmapt^2
  pi <- (mu-rf)*12
  cov <- cov*12
  k <- dim(P)[1] # No. de views
  omega <- matrix(0,k,k)
  for(i in 1:k){
    omega[i,i] <- tau*(t(P[i,])%*%cov%*%P[i,])  # tau * P' V P
  }
  r.BL <- solve(solve(c(tau)*cov)+t(P)%*%solve(omega)%*%P)%*%(solve(c(tau)*cov)%*%pi+t(P)%*%solve(omega)%*%q)
  w.BL <- solve(c(lambda)*cov)%*%r.BL/sum(solve(c(lambda)*cov)%*%r.BL)
  BL <- list()
  BL[[1]] <- r.BL
  BL[[2]] <- w.BL
  return(BL)
}


m.blackl <- function(retornos,rf,tau,q,P,short){
  
  activos <- names(retornos)
  
  ## Con cortos
  r.BL <- c.blackl(retornos,rf,tau,q,P)[[1]] 
  w.BL <- c.blackl(retornos,rf,tau,q,P)[[2]]
  ##Sin cortos
  mu <- (colMeans(retornos))
  cov <- cov(retornos)
  
  cov <- cov*12
  n <- length(mu)
  
  Dmat <- cov*2
  dvec <- rep(0,n)
  Amat <- matrix(c(r.BL,-r.BL,rep(1,n),-rep(1,n),diag(length(r.BL))),n,n+4)
  nport <- 10000
  j <- seq(min(r.BL)+0.01,max(r.BL)-0.01,length=nport)
  sigmapBL <- matrix(0,nrow=nport)
  wpoBL <- matrix(0,nrow=nport,ncol=n)
  sortinopBL <- matrix(0,nrow=nport)
  
  for(i in 1:nport){
    rp <- j[i]
    bvec <- c(rp,-rp,1,-1,rep(0,n))
    poptimo <- solve.QP(Dmat,dvec,Amat,bvec,meq=2)
    wpoBL[i,] <- poptimo$solution
    sigmapBL[i,] <- sqrt(poptimo$value)
    retornosBL <- retornos%*%wpoBL[i,]
    sortinopBL[i,] <- mean(retornosBL)/sd(ifelse(retornosBL<0,retornosBL,0))
  }
  
  # Portafolio tangente del modelo BL
  sharpe_port <- sortinopBL
  sharpe <- cbind(sharpe_port,wpoBL)
  sharpe <- cbind(j,sharpe)
  sharpe.sort <- sharpe[order(-sharpe[,2]),]
  sharpe.sel <- cbind(sharpe.sort[1,])
  rptBL <- round(sharpe.sel[1],6)
  wptBL <- round(cbind(sharpe.sel[3:length(sharpe.sel)]),6)
  rownames(wptBL) <- c(activos)
  
  BL <- list()
  if(short==1){
    BL[[1]] <- w.BL}
  else{
    BL[[1]] <- wptBL
  }
  
  BL[[2]]<- t(BL[[1]])%*%mu
  BL[[3]]<- sd(retornos%*%BL[[1]])
  
  return(BL)
  
}





