#NBLDA
Soft <- function(x,a){
  return(sign(x)*pmax(abs(x)-a,0))
}









GetD <- function(ns, x, y, rho,beta,rhos=NULL){
  if(!is.null(rho) && !is.null(rhos)) stop("do you want to use rho or rhos in GetD function???")
  if(is.null(rhos)){
    uniq <- sort(unique(y))
    ds <- matrix(1, nrow=length(uniq), ncol=ncol(x))
    for(k in 1:length(uniq)){
      a <- colSums(x[y==uniq[k],])+beta
      b <- colSums(ns[y==uniq[k],])+beta
      ds[k,] <- 1+Soft(a/b-1,rho/sqrt(b))
    }
    return(ds)
  } else {
    uniq <- sort(unique(y))
    ds.list <- list()
    for(rho in rhos){
      ds <- matrix(1, nrow=length(uniq), ncol=ncol(x))
      for(k in 1:length(uniq)){
        a <- colSums(x[y==uniq[k],])+beta
        b <- colSums(ns[y==uniq[k],])+beta
        ds[k,] <- 1+Soft(a/b-1,rho/sqrt(b))
      }
      ds.list[[which(rhos==rho)]] <- ds
    }
    return(ds.list)
  }
}





permute.rows <- function(x){
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = T)
}









balanced.folds <- function(y, nfolds = min(min(table(y)), 10)){
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)
  # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)
  # nice way to get the ids in a list, split by class
  ###Make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) {
    bigmat[seq(totals[i]), i] <- sample(yids[[i]])
  }
  smallmat <- matrix(bigmat, nrow = nfolds) # reshape the matrix
  ### Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat)) ### Now a clever unlisting
  x <- apply(smallmat, 2, function(x) x[!is.na(x)])
  if(is.matrix(x)){
    xlist <- list()
    for(i in 1:ncol(x)){
      xlist[[i]] <- x[,i]
    }
    return(xlist)
  }
  return(x)
}



NBLDA.cv <-
  function(x,y,rhos=NULL,beta=1,nfolds=5,phihat=0,type=c("mle","deseq","quantile"),folds=NULL, prior=NULL){
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
      out <- NBLDA(x[tr,],y[tr],x[te,],rhos=rhos,phihat=phihat,beta=beta,
                   type="mle", prior=prior) # Have already power-transformed x, so don't need to do it again!!!
      for(j in 1:length(rhos)){      
        errs[i,j] <- sum(out[[j]]$ytehat!=y[te])
        nnonzero[i,j] <- sum(colSums(out[[j]]$ds!=1)!=0)
      }
    }
    cat(fill=TRUE)
    save <- list(errs=errs, bestrho=rhos[max(which(colMeans(errs)==min(colMeans(errs))))], 
                 rhos=rhos, nnonzero=nnonzero,folds=folds,type=type)
    return(save)
  }







NBLDA <-
  function(x,y,xte=NULL,rho=0,beta=1,rhos=NULL,phihat=0,type=c("mle","deseq","quantile"), prior=NULL){
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
    if(is.null(rhos)){
      ds <- GetD(ns,x,y,rho,beta)
      discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
      for (k in 1:length(uniq)) {
        
        for(l in 1:nrow(xte))   {
          
          dstar = ds[k,]
          part2=1+nste[l,]*dstar*phihat 
          part1=dstar/part2 
          
          
          discriminant[l, k]<- sum(xte[l,]*log(part1))-sum((1/phihat)*log(part2))+log(prior[k])
          
        }
      }
      save <- list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max), rho=rho,x=x,y=y,xte=xte,type=type)
      return(save)
    } else {
      save <- list()
      ds.list <- GetD(ns,x,y,rho=NULL, rhos=rhos,beta)
      for(rho in rhos){
        ds <- ds.list[[which(rhos==rho)]]
        discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
        for (k in 1:length(uniq)) {
          
          for(l in 1:nrow(xte))   {
            
            dstar = ds[k,]
            part2=1+nste[l,]*dstar*phihat 
            part1=dstar/part2 
            
            
            discriminant[l, k]<- sum(xte[l,]*log(part1))-sum((1/phihat)*log(part2))+log(prior[k])
            
          }
        }
        save[[which(rhos==rho)]] <- (list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max), rho=rho,x=x,y=y,xte=xte,type=type))
      }
      return(save)
    }
  }
