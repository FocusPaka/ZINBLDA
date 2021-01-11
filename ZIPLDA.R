########################
#ZIPLDA分类器
#参数估计函数
# 求d_{kg}的值
GetDn <- function(ns, x, y, beta){
  uniq <- sort(unique(y))
  ds <- matrix(1, nrow=length(uniq), ncol=ncol(x))
  for(k in 1:length(uniq)){
    a <- colSums(x[y==uniq[k],])+beta
    b <- colSums(ns[y==uniq[k],])+beta
    ds[k,] <- a/b
  }
  return(ds)
}

#ZIPLDA估计数据中0的概率
estimatepZIP<-function(x,y,xte=NULL,yte,beta=1,type=c("mle","deseq","quantile"), prior=NULL)
{
  if(is.null(xte)){
    xte <- x
    warning("Since no xte was provided, testing was performed on training data set.")
  }
  type <- match.arg(type)      #匹配参数,允许模糊匹配
  if(is.null(prior)) 
    prior <- rep(1/length(unique(y)), length(unique(y)))     #这个参数后面没有用上，不知道含义是什么
  null.out <- NullModel(x, type=type)     #输出两个数值，n和sizes, n表示训练集size factor*lambda_g的值
  ns <- null.out$n
  nste <- NullModelTest(null.out,x,xte,type=type)$nste     #测试样本size facror*lambda_g的值
  uniq <- sort(unique(y))                 #将类别从小到大排序
  ds <-  GetDn(ns,x,y,beta)               #求d_{kg}的值
  
  mu <- matrix(NA, nrow=nrow(x), ncol=length(x[1,]))
  mute <- matrix(NA, nrow=nrow(xte), ncol=length(xte[1,]))
  for(i in 1:nrow(x)){                  #求出训练集的mu值
    dstar = ds[y[i],]
    mu[i,] <- ns[i,]*dstar
  }
  
  for(i in 1:nrow(xte)){                #求出测试集的mu值
    dstar = ds[yte[i],]
    mute[i,] <- nste[i,]*dstar
  }
  
  x1 <- t(x)              #求出log(P0/(1-p0))的值，即logistic regression中的响应变量，先将片段数不为0的设为1，然后将所有的化为一个向量
  x2 <- as.vector(x1)     #向量是按列排列，即按照每个样本每个样本排列
  y <- x2
  y[y!=0] <- 1
  G <- length(x[1,])     #再求N_{ki_k}/N_{1i_1}
  lib <- rowSums(x)      #求出N_{ki_k}
  librep <- rep(lib,rep(G,length(lib)))/(lib[1])      
  mu1 <- as.vector(t(mu))   #最后求mu_{ki_kg}
  xreg <- cbind(y,librep,mu1)      #将上面三项数据合成一个矩阵，然后利用glm()函数来拟合logistic regression得到系数
  glm.out <- glm(y ~ librep+ mu1,family=binomial("logit"),data=data.frame(xreg))
  summary(glm.out)
  coef <- as.matrix(glm.out$coefficients)
  inter <- rep(1,G)
  muu <- t(mute)
  libte <- rowSums(xte)
  xte1 = t(xte)
  p <- xte1
  for(i in 1:length(xte1[1,])){
    libsize <- rep(sum(xte1[,i]),G)/libte[1]
    estx1 <- cbind(inter,libsize,muu[,i])
    dd <- estx1%*% coef
    dd[dd>50] <- 50
    dd[dd<(-50)] <- -50
    p1<-exp(dd)
    
    #p[,i] <- ((1-exp(-muu[,i])*(1+p1))/((1+p1)*(1-exp(-muu[,i]))))    #此处公式有问题
    p[,i] <- ((p1-exp(-muu[,i])*(1+p1))/((1+p1)*(1-exp(-muu[,i]))))
  }
  p[p<0] <- 0
  pp <- p
  return(pp)   
}










ZIPDA.cv <-
  function(x,y,rhos=NULL,beta=1,nfolds=5,prob0=NULL,type=c("mle","deseq","quantile"),
           folds=NULL,transform=TRUE, alpha=NULL, prior=NULL){
    
    type <- match.arg(type)
    if(!transform && !is.null(alpha)) stop("You have asked for NO 
                                           transformation but have entered 
                                           alpha.")
    if(transform && is.null(alpha)) alpha <- FindBestTransform(x)      
    #Find the power transformation that makes a data set approximately Poisson.
    
    if(transform){
      if(alpha<=0 || alpha>1) stop("alpha must be between 0 and 1")
      x <- x^alpha     #将数据进行变换使得数据集近似poisson分布
    }
    if(is.null(rhos)){
      ns <- NullModel(x,type=type)$n   #训练数据的size factor*lambda_g
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
      out <- ZIPLDA(x[tr,],y[tr],x[te,],rhos=rhos,beta=beta,prob0=prob0,
                    type="mle", prior=prior, transform=FALSE) 
      # Have already power-transformed x, so don't need to do it again!!!
      for(j in 1:length(rhos)){      
        errs[i,j] <- sum(out[[j]]$ytehat!=y[te])
        nnonzero[i,j] <- sum(colSums(out[[j]]$ds!=1)!=0)
      }
    }
    cat(fill=TRUE)
    save <- list(errs=errs, bestrho=rhos[max(which(colMeans(errs)==min(colMeans(errs))))], 
                 rhos=rhos, nnonzero=nnonzero,folds=folds, alpha=alpha,type=type)
    return(save)
  }







ZIPLDA<-function(x,y,xte=NULL,rho=0,beta=1,rhos=NULL,prob0=NULL,
                 type=c("mle","deseq","quantile"), prior=NULL, 
                 transform=TRUE, alpha=NULL){
  if(is.null(xte)){
    xte <- x
    warning("Since no xte was provided, testing was 
              performed on training data set.")
  }
  if(!is.null(rho) && length(rho)>1) stop("Can only enter 1 value of rho. 
                                            If you would like to enter multiple
                                            values, use rhos argument.")
  type <- match.arg(type)
  if(!transform && !is.null(alpha)) stop("You have asked for NO 
                                           transformation but have entered 
                                           alpha.")
  if(transform && is.null(alpha)) alpha <- FindBestTransform(x)
  if(transform){
    if(alpha<=0 || alpha>1) stop("alpha must be between 0 and 1")
    x <- x^alpha
    xte <- xte^alpha
  }  
  if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
  if(is.null(rho)&&is.null(rhos)) stop("Must enter rho or rhos.")
  null.out <- NullModel(x, type=type)
  ns <- null.out$n
  nste <- NullModelTest(null.out,x,xte,type=type)$nste
  uniq <- sort(unique(y))
  signx1<-sign(xte==0)
  if(is.null(rhos)){
    ds <- GetD(ns,x,y,rho,beta)
    discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
    for(k in 1:length(uniq)){
      for(i in 1:nrow(xte)){
        dstar = ds[k,]
        part2=nste[i,]*dstar 
        part1=prob0[i,]+(1-prob0[i,])*exp(-part2) 
        part1[part1==0]=1
        discriminant[i,k] <-sum(signx1[i,]*log(part1))+sum(xte[i,]*(1-signx1[i,])*log(dstar))-sum((1-signx1[i,])*part2)+log(prior[k])
      }
    }
    save <- list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,
                 ytehat=apply(discriminant,1,which.max), alpha=alpha, 
                 rho=rho,x=x,y=y,xte=xte,alpha=alpha,type=type)
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
          part1=prob0[i]+(1-prob0[i])*exp(-part2) 
          part1[part1==0]=1
          discriminant[i,k] <-sum(signx1[i,]*log(part1))+sum(xte[i,]*(1-signx1[i,])*log(dstar))-sum((1-signx1[i,])*part2)+log(prior[k])
          
        }
        
      }
      save[[which(rhos==rho)]] <- (list(ns=ns,nste=nste,ds=ds,
                                        discriminant=discriminant
                                        ,ytehat=apply(discriminant,1,which.max),
                                        alpha=alpha, rho=rho,x=x,y=y,xte=xte,
                                        alpha=alpha,type=type))
    }
    return(save)
  }
}


