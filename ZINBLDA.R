################
#ZINBLDA
#estimated parameter based on the ZINBLDA model
fun <- function(x,mu,disperhatmean){
  
  
  signx <- sign(x==0)
  zinb <- function(p) {
    
    res <- sum(sign(x==0)*log(p[1]+(1-p[1])*(1/(1+p[3]*p[2]))^(1/p[2]))+
                 (1-sign(x==0))*(log(1-p[1])+lgamma(x+1/p[2])-lgamma(x+1)-
                                   lgamma(1/p[2])+x*log(p[3]*p[2])-x*log(1+p[3]*p[2])-(1/p[2])*log(1+p[3]*p[2])))
    
    return(-res)
    
  }
  #nlminb(c(0,1,mu),zinb,lower=c(0.01,0.001,0.01),upper=c(0.999,disperhatmean,9999),control = list(step.max=0.2))$par
  
  nlminb(c(0,1,mu),zinb,lower=c(0.01,0.001,0.01),upper=c(0.999,disperhatmean+0.5,9999),control = list(step.max=0.2))$par
}



#ZINB估计数据中0的概率
estimatepZINB <- function(x,y,xte=NULL,yte,phihat,beta=1,type=c("mle","deseq","quantile"), prior=NULL){
  if(is.null(xte)){
    xte <- x
    warning("Since no xte was provided, testing was performed on training data set.")
  }
  
  type <- match.arg(type)
  if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
  
  null.out <- NullModel(x, type=type)
  ns <- null.out$n
  nste <- NullModelTest(null.out,x,xte,type=type)$nste
  uniq <- sort(unique(y))
  ds <-  GetDn(ns,x,y,beta)
  
  md <- matrix(NA, nrow=nrow(x), ncol=length(x[1,]))
  mdte <- matrix(NA, nrow=nrow(xte), ncol=length(xte[1,]))
  for(i in 1:nrow(x)){
    dstar = ds[y[i],]
    md[i,] <- ns[i,]*dstar
  }
  
  for(i in 1:nrow(xte)){
    dstar = ds[yte[i],]
    mdte[i,] <- nste[i,]*dstar
  }
  
  
  G <- length(x[1,])
  lib <- rowSums(x)
  x1 <- t(x)
  
  md1 <- as.vector(t(md))
  librep <- rep(lib,rep(G,length(lib)))/(lib[1])
  x2 <- as.vector(x1)
  
  y <- x2
  y[y!=0] <- 1
  xreg <- cbind(y,librep,md1)
  glm.out <- glm(y~ librep+ md1,family=binomial("logit"),data=data.frame(xreg))
  summary(glm.out)
  
  coef <- as.matrix(glm.out$coefficients)
  inter <- rep(1,G)
  mdd <- t(mdte)
  libte <- rowSums(xte)
  
  xte1 = t(xte)
  p <- xte1
  for(i in 1:length(xte1[1,])){
    
    libsize <- rep(sum(xte1[,i]),G)/libte[1]
    estx1 <- cbind(inter,libsize,mdd[,i])
    dd <- estx1%*% coef
    dd[dd > 50] <- 50
    dd[dd < (-50)] <- -50
    p1 <- (1/(1+mdd[,i]*phihat))^(1/phihat)
    p2 <- exp(dd)
    p[,i] <- (p2-(1+p2)*p1)/((1+p2)*(1-p1))
  }
  p[p<0] <- 0
  pp <- p
  return(pp)   
}




ZINB.cv <-
  function(x,y,rhos=NULL,beta=1,nfolds=5,phihat=0,prob0=NULL,type=c("mle","deseq","quantile"),folds=NULL, prior=NULL){
    type <- match.arg(type)
    if(is.null(rhos)){
      ns <- NullModel(x,type=type)$n
      uniq <- sort(unique(y))
      maxrho <- rep(NA, length(uniq))
      for(k in 1:length(uniq)){
        a <- colSums(x[y==uniq[k],])+beta
        b <- colSums(ns[y==uniq[k],])+beta
        maxrho[k] <- max(abs(a/b-1)*sqrt(b),na.rm=TRUE)
      }
      rhos <- seq(0, max(maxrho,na.rm=TRUE)*(2/3), len=30)
    }
    if(is.null(folds)) folds <- balanced.folds(y,nfolds=nfolds)
    nfolds <- length(folds)
    errs <- nnonzero <- matrix(NA, nrow=nfolds, ncol=length(rhos))
    for(i in 1:nfolds){
      cat(i,fill=FALSE)
      tr <- -folds[[i]]
      te <- folds[[i]]
      out <-ZINBLDA(x[tr,],y[tr],x[te,],rhos=rhos,phihat=phihat,beta=beta,prob0=prob0,type="mle", prior=prior) # Have already power-transformed x, so don't need to do it again!!!
      for(j in 1:length(rhos)){      
        errs[i,j] <- sum(out[[j]]$ytehat!=y[te])
        nnonzero[i,j] <- sum(colSums(out[[j]]$ds!=1)!=0)
      }
    }
    cat(fill=TRUE)
    save <- list(errs=errs, bestrho=rhos[max(which(colMeans(errs)==min(colMeans(errs))))], rhos=rhos, nnonzero=nnonzero,folds=folds,type=type)
    return(save)
  }


ZINBLDA<-
  function(x,y,xte=NULL,rho=0,beta=1,rhos=NULL,phihat=0,prob0=NULL,type=c("mle","deseq","quantile"), prior=NULL){
    if(is.null(xte)){
      xte <- x
      warning("Since no xte was provided, testing was performed on training data set.")
    }
    if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
    if(is.null(rho)&&is.null(rhos)) stop("Must enter rho or rhos.")
    null.out <- NullModel(x, type=type)
    ns <- null.out$n
    nste <- NullModelTest(null.out,x,xte,type=type)$nste
    uniq <- sort(unique(y))
    signx3<-sign(xte==0)
    if(is.null(rhos)){
      ds <- GetD(ns,x,y,rho,beta)
      discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
      for(k in 1:length(uniq)){
        for(i in 1:nrow(xte)){
          dstar = ds[k,]
          part2=nste[i,]*dstar 
          part3<-(1/(1+part2*phihat))^(1/phihat)
          discriminant[i,k] <-sum(signx3[i,]*log(prob0[i,]+(1-prob0[i,])*part3))+sum(xte[i,]*(1-signx3[i,])*(log(dstar)-log(1+part2*phihat)))-sum((1-signx3[i,])*(1/phihat)*log(1+part2*phihat))+log(prior[k])
          
        }
      }
      save <- list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max),rho=rho,x=x,y=y,xte=xte,type=type)
      return(save)
    } else {
      save <- list()
      ds.list <- GetD(ns,x,y,rho=NULL, rhos=rhos,beta)
      for(rho in rhos){
        ds <- ds.list[[which(rhos==rho)]]
        discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
        for(k in 1:length(uniq)){
          for(i in 1:nrow(xte))   {
            
            dstar = ds[k,]
            part2=nste[i,]*dstar 
            part3<-(1/(1+part2*phihat))^(1/phihat)
            discriminant[i,k] <-sum(signx3[i,]*log(prob0[i,]+(1-prob0[i,])*part3))+sum(xte[i,]*(1-signx3[i,])*(log(dstar)-log(1+part2*phihat)))-sum((1-signx3[i,])*(1/phihat)*log(1+part2*phihat))+log(prior[k])
            
          }
        }
        save[[which(rhos==rho)]] <- (list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max),rho=rho,x=x,y=y,xte=xte,type=type))
      }
      return(save)
    }
  }


