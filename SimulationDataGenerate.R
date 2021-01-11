#用负二项分布生成模拟数据
DataSetzeroNB <- function (n, p, K, param, sdsignal, DE,pzero=NA) 
{
  if (n < 4 * K) 
    stop("We require n to be at least 4*K.")
  
  q0 <- rexp(p, rate = 1/25)    #生成p个期望为25的值，即为lambda_g的值
  isDE <-  runif(p)< DE         #生成p个均匀分布的随机变量，标出差异表达基因
  #逻辑值，如果是DE genes，则为TRUE
  classk <- matrix(NA, nrow = K, ncol = p)      #生成一个K行p列的空矩阵
  for (k in 1:K) {
    lfc <- rnorm(p, mean=0, sd = sdsignal)      #生成p个均值为0，方差为sd的正态分布随机变量
    classk[k, ] <- ifelse(isDE, q0 * exp(lfc), q0)  #如果是DE genes，得到lambda_g*d_{kg}的值，
    #否则，得到lambda_g的值（即d_{kg}为1，说明是non-DE gene）
  }
  
  truesf <- runif(n)*20+20      #训练数据的size factor
  truesfte <- runif(n)*20+20    #测试数据的size factor
  
  conds <- sample(c(rep(1:K, 4), sample(1:K, n - 4 * K, replace = TRUE)))    #训练集，随机分配样本的类别，对不整除类别的情况随机抽样
  condste <- sample(c(rep(1:K, 4), sample(1:K, n - 4 * K, replace = TRUE)))  #测试集，随机分配样本的类别，对不整除类别的情况随机抽样
  x <- xte <- matrix(NA, nrow = n, ncol = p)     #生成两个（训练、测试）n行p列的矩阵
  
  if (is.null(pzero)) {            #判断是否设置了样本为0的比例
    prob<- runif(n,0.05,0.5)       #如果没有设置，则从[0.05,0.5]的均匀分布中随机抽n个
  }else{
    prob<-rep(pzero,n)             #否则直接使用设置的比例
  }
  
  for (i in 1:n) {       
    for (k in 1:K) {
      if (conds[i] == k){                 #开始生成训练数据的样本,先是k=1的所有样本，再是k=2的所有样本
        x[i, ] <- rnbinom(p, mu = truesf[i] * classk[k, ], size = param)   #生成k=1的第一个样本，利用负二项分布，均值为mu, dispersion参数为size
        index<-which(x[i,]>0)
        if(sum(x[i,]>0)>(p*prob[i])){      #生成为0的基因，如果大于0的基因的个数大于设定的比例
          ind0<-sample(index,round(p*prob[i]), replace = FALSE)    #随机从大于0的基因中抽出设定比例的基因，令其等于0
          x[i,ind0]<-0
        }else{                    #如果大于0的基因个数小于设定值，则随机抽一个gene的值设为1，其他的都设为0
          x[i,]<-0
          x[i,sample(1:p,1)]<-1
        }
      }
      
      if (condste[i] == k){               #同样的步骤生成测试样本
        xte[i, ] <- rnbinom(p, mu = truesfte[i] * classk[k, ], size = param)
        index<-which(xte[i,]>0)
        if(sum(xte[i,]>0)>(p*prob[i])){
          ind0<-sample(index,round(p*prob[i]), replace = FALSE)
          xte[i,ind0]<-0
        }else{
          xte[i,]<-0
          xte[i,sample(1:p,1)]<-1
          
        }
      }
    }
  }
  
  rm <- apply(x, 2, sum) == 0      #对训练数据，如果列求和为0，说明该gene不表达，则在训练集和测试集中去掉
  return(list(x = x[, !rm], xte = xte[, !rm], y = conds, yte = condste, 
              truesf = truesf, truesfte = truesfte, dkg=classk))
}
