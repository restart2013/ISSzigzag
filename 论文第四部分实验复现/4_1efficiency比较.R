### importance sub-sampling zigzag process算法
### 与unif sub-sampling的zigzag process 的 efficiency比较
### 该文件中的内容是对论文4.1节的复现，比较两种方法的efficiency
###
### 其中均匀分布sub-sampling的zigzag process只需对算法实现.R中
### 第74，78，84行分别修改即可得到
###
### 论文figure1左图
### 由于全部复现论文中figure1左图内容耗时较长，本代码中，只对alpha=0.1和0.3和1
### 三种情况进行复现
### 对于三种不同的rho分布，在第159行代码分别替换对应的rho分布即可
###
### 论文figure1右图
### 按论文方式复现。observation number=10时由于数据x过少，可能会导致xi方差过大越界的情况，
### 即 if (a < mi_(theta, y, x, xi_seq[, step], i, B)/Mi_(x, i)) 处 会抛出错误
### 尝试多次更换随机数种子重新生成新的一组数据x即可
###
### 
### 由于代码运行时间较长，实验结果直接记录在了24行-74行，运行后便能查看
### 正式代码过程在75行之后
###
###
######################以下为代码部分#######################
###################4.1 figure1左图部分#####################

par(mfrow=c(1,1))
x_lab <- c(0.1, 0.3, 1)
# 实验结果
T_imp_norm <- c(3385.607, 1117.641, 311.4739)
T_unif_norm <- c(108498.9, 12757.6, 1462.174)
efficiency_norm <- T_unif_norm/T_imp_norm
T_imp_unif <- c(4181.735, 1562.377, 1178.217)
T_unif_unif <- c(83450.71, 10395.96, 2330.876)
efficiency_unif <- T_unif_unif/T_imp_unif
T_imp_laplace <- c(3672.112, 1154.691, 326.3651)
T_unif_laplace <- c(107027, 16404.04, 2475.31)
efficiency_laplace <- T_unif_laplace/T_imp_laplace
# 画图
plot(log(x_lab,10),log(efficiency_norm,10), type='o', lwd=2, pch=16, cex=1.2, col='blue',
     xaxt='n', yaxt='n',xlim=c(-2,0), ylim=c(0,2), xlab='alpha', ylab='Efficiency gain')
lines(log(x_lab,10),log(efficiency_unif,10), type='o', lwd=2, pch=16, cex=1.2, col='green')
lines(log(x_lab,10),log(efficiency_laplace,10), type='o', lwd=2, pch=16, cex=1.2, col='red')
abline(0, -1, col='gray',lwd=2)
legend('topright',legend=c('laplace','gaussian','uniform'),col=c('red','blue','green'),pch=16,lwd=2)
axis(1,at=c(-2,-1,0),labels=c('10^-2','10^-1','10^0'))
axis(2,at=c(0,1,2),labels=c('10^0','10^1','10^2'))
grid()

###################4.1 figure1右图部分#####################

par(mfrow=c(1,1))
x_lab <- c(10, 30, 100, 300, 1000, 3000)
# 实验结果
T_imp_norm <- c(27972.73, 8877.328, 1903.532, 595.879, 127.4075, 27.92898)
T_unif_norm <- c(60715.58, 24101.68, 7069.308, 2389.372, 773.788, 238.3263)
efficiency_norm <- T_unif_norm/T_imp_norm
T_imp_unif <- c(513163.7, 36158.21, 6071.431, 1803.89, 612.4633, 217.8121)
T_unif_unif <- c(506887.5, 68465.84, 11940.94, 3577.381, 1188.807, 419.7562)
efficiency_unif <- T_unif_unif/T_imp_unif
T_imp_laplace <- c(18917.14, 5276.082, 2583.839, 740.6308, 129.3354, 21.9232)
T_unif_laplace <- c(46210.65,22940.74, 13281.48, 4952.098, 1383.95, 449.8852)
efficiency_laplace <- T_unif_laplace/T_imp_laplace
# 画图
plot(log(x_lab,10), log(efficiency_norm), type='o', lwd=2, pch=16, cex=1.2, col='blue',
     xaxt='n', yaxt='n',xlim=c(1,4), ylim=c(log(1.7),log(10)), xlab='observation number', ylab='Efficiency gain')
lines(log(x_lab,10), log(efficiency_unif), type='o', lwd=2, pch=16, cex=1.2, col='green')
lines(log(x_lab,10), log(efficiency_laplace), type='o', lwd=2, pch=16, cex=1.2, col='red')
lines(log(seq(10,10000,2),10),log(log(seq(10,10000,2))),lwd=2,col='gray')
lines(log(seq(10,10000,2),10),log(log(seq(10,10000,2)/log(seq(10,10000,2)))*(1+1/log(seq(10,10000,2)))),lwd=2,col='gray')
legend('topright',legend=c('laplace','gaussian','uniform'),col=c('red','blue','green'),pch=16,lwd=2)
axis(1,at=c(1,2,3,4),labels=c('10^1','10^2','10^3','10^4'))
axis(2,at=log(c(1,2,4,6,8,10)),labels=c(1,2,4,6,8,10))
grid()

###################4.1 figure1共同部分#####################
# x维度5，样本量500
p = 5
n = 500

# xi真值
xi_true <- rep(-1,p)

# xi先验初值
xi <- rep(20,p)
# xi先验分布为正态N(0,sigma^2)，sigma=10^5
sigma <- 10^5

# 函数定义
# 抽取反弹时间t，mi_o部分
generate_t_mi_o <- function(unif,xi){
  unif <- -log(unif)
  t <- (sqrt(unif+abs(xi))-sqrt(abs(xi)))*sqrt(2*sigma^2)
  return(t)
}
# 抽取反弹时间t，mi_*部分
generate_t_mi_ <- function(unif,x){
  unif <- -log(unif)
  t <- unif/rowSums(abs(x))
  return(t)
}
# rate函数mi_o
mi_o <- function(theta,xi,i){
  mi <- theta[i]*xi[i]/(sigma^2)
  return(max(0,mi))
}
# mi_o上界
Mi_o <- function(xi,i,t){
  Mi <- (abs(xi[i])+t)/sigma^2
  return(Mi)
}
# rate函数mi_*
mi_ <- function(theta,y,x,xi,i,B){
  mi <- 0
  for (b in 1:B){
    j <- sample(1:n,1,prob=abs(x[i,]))
    y_ <- y[j]
    x_ <- x[,j]
    mi_ <- x_[i]*exp(x_%*%xi)/(1+exp(x_%*%xi))-y_*x_[i]
    mi <- mi+sum(abs(x[i,]))/abs(x_[i])*mi_
  }
  return(max(0,theta[i]*mi/B))
}
# rate函数mi_*普通版
mi_n <- function(theta,y,x,xi,i,B){
  mi <- 0
  for (b in 1:B){
    j <- sample(1:n,1)
    y_ <- y[j]
    x_ <- x[,j]
    mi_ <- x_[i]*exp(x_%*%xi)/(1+exp(x_%*%xi))-y_*x_[i]
    mi <- mi+n*mi_
  }
  return(max(0,theta[i]*mi/B))
}
# mi_*上界
Mi_ <- function(x,i){
  return(sum(abs(x[i,])))
}
# mi_*上界普通版
Mi_n <- function(x,i){
  return(n*max(abs(x[i,])))
}

###################4.1 figure1左图部分#####################
# 修改随机数种子50次，结果取平均即可
set.seed(114514)

# efficiency比较
T_imp <- vector()
T_unif <- vector()

# alpha不同值
alpha_list <- c(0.1,0.3,1)
for (alpha_i in 1:3){
  # x和y
  # 混合分布(1-alpha)0 + alpha标准正态 ，alpha=0.1
  alpha <- alpha_list[alpha_i]
  # 从混合分布中抽取500个x
  x <- array( rnorm(p*n,0,1) ,dim=c(p,n)) # rnorm 分别替换均匀分布runif和拉普拉斯分布rlaplace即可，rlaplace需要library(VGAM)
  a <- array( as.integer(runif(p*n,0,1) < alpha) ,dim=c(p,n))
  x <- x*a
  # 根据xi真值产生逻辑回归y
  y <- as.integer( runif(n,0,1) < 1/(1+exp(-xi_true%*%x)) )
  
  # improved sub-sampling zigzag过程
  # 跑k步
  k = 3*10^7
  # 初始值
  xi_seq <- array(0,dim=c(p,k))
  xi_seq[,1] <- xi
  t_seq <- vector()
  t_seq[1] <- 0
  
  # theta初始值
  theta <- rep(1,p)
  
  # 算法参数mini_batch
  B = 1
  
  # 记录theta正负号反转次数
  bouncing_time <- 0
  
  # improved zigzag
  for (step in 2:k){
    t_prior <- generate_t_mi_o(runif(p,0,1),xi_seq[,step-1])
    t_condition <- generate_t_mi_(runif(p,0,1),x)
    i_p <- which.min(t_prior)
    i_c <- which.min(t_condition)
    if(t_prior[i_p]<t_condition[i_c]){
      i <- i_p
      t <- t_prior[i]
      # 更新序列
      t_seq[step] <- t_seq[step-1]+t
      xi_seq[,step] <- xi_seq[,step-1]+t*theta
      # 更新theta
      a <- runif(1,0,1)
      if(a < mi_o(theta,xi_seq[,step],i)/Mi_o(xi_seq[,step-1],i,t) ){
        theta[i] <- theta[i]*-1
        bouncing_time <- bouncing_time+1
      }
    }else{
      i <- i_c
      t <- t_condition[i]
      # 更新序列
      t_seq[step] <- t_seq[step-1]+t
      xi_seq[,step] <- xi_seq[,step-1]+t*theta
      # 更新theta
      a <- runif(1,0,1)
      if(a < mi_(theta,y,x,xi_seq[,step],i,B)/Mi_(x,i) ){
        theta[i] <- theta[i]*-1
        bouncing_time <- bouncing_time+1
      }
    }
    if (bouncing_time == 10^5){
      # theta正负号反转10^5次后停止
      break
    }
  }
  T_imp[alpha_i] <- max(t_seq)
  print('improved zigzag total simulation time:')
  print(max(t_seq))
  print('bouncing time:')
  print(bouncing_time)
  
  # unif sub-sampling zigzag过程
  # 跑k步
  k = 3*10^7
  # 初始值
  xi_seq <- array(0,dim=c(p,k))
  xi_seq[,1] <- xi
  t_seq <- vector()
  t_seq[1] <- 0
  
  # theta初始值
  theta <- rep(1,p)
  
  # 算法参数mini_batch
  B = 1
  
  # 记录theta正负号反转次数
  bouncing_time <- 0
  # normal zigzag
  for (step in 2:k){
    t_prior <- generate_t_mi_o(runif(p,0,1),xi_seq[,step-1])
    t_condition <- generate_t_mi_(runif(p,0,1),x)
    i_p <- which.min(t_prior)
    i_c <- which.min(t_condition)
    if(t_prior[i_p]<t_condition[i_c]){
      i <- i_p
      t <- t_prior[i]
      # 更新序列
      t_seq[step] <- t_seq[step-1]+t
      xi_seq[,step] <- xi_seq[,step-1]+t*theta
      # 更新theta
      a <- runif(1,0,1)
      if(a < mi_o(theta,xi_seq[,step],i)/Mi_o(xi_seq[,step-1],i,t) ){
        theta[i] <- theta[i]*-1
        bouncing_time <- bouncing_time+1
      }
    }else{
      i <- i_c
      t <- t_condition[i]
      # 更新序列
      t_seq[step] <- t_seq[step-1]+t
      xi_seq[,step] <- xi_seq[,step-1]+t*theta
      # 更新theta
      a <- runif(1,0,1)
      if(a < mi_n(theta,y,x,xi_seq[,step],i,B)/Mi_n(x,i) ){
        theta[i] <- theta[i]*-1
        bouncing_time <- bouncing_time+1
      }
    }
    if (bouncing_time == 10^5){
      # theta正负号反转10^5次后停止
      break
    }
  }
  T_unif[alpha_i] <- max(t_seq)
  print('unif zigzag total simulation time:')
  print(max(t_seq))
  print('bouncing time:')
  print(bouncing_time)
}

# 画图
x_lab <- c(0.1, 0.3, 1)
plot(log(x_lab,10),log(T_unif/T_imp,10), type='o', lwd=2, pch=16, cex=1.2, col='blue',
     xaxt='n', yaxt='n',xlim=c(-2,0), ylim=c(0,2), xlab='alpha', ylab='Efficiency gain')
abline(0, -1, col='gray',lwd=2)
legend('topright',legend=c('laplace','gaussian','uniform'),col=c('red','blue','green'),pch=16,lwd=2)
axis(1,at=c(-2,-1,0),labels=c('10^-2','10^-1','10^0'))
axis(2,at=c(0,1,2),labels=c('10^0','10^1','10^2'))

###################4.1 figure1右图部分#####################
# 修改随机数种子50次，结果取平均即可
set.seed(8)

# efficiency比较
T_imp <- vector()
T_unif <- vector()

# n不同值
n_list <- c(10,30,100,300,1000,3000)
for (n_i in 1:6){
  # 不同observation数
  n <- n_list[n_i]
  # 重写函数
  # rate函数mi_*
  mi_ <- function(theta,y,x,xi,i,B){
    mi <- 0
    for (b in 1:B){
      j <- sample(1:n,1,prob=abs(x[i,]))
      y_ <- y[j]
      x_ <- x[,j]
      mi_ <- x_[i]*exp(x_%*%xi)/(1+exp(x_%*%xi))-y_*x_[i]
      mi <- mi+sum(abs(x[i,]))/abs(x_[i])*mi_
    }
    return(max(0,theta[i]*mi/B))
  }
  # rate函数mi_*普通版
  mi_n <- function(theta,y,x,xi,i,B){
    mi <- 0
    for (b in 1:B){
      j <- sample(1:n,1)
      y_ <- y[j]
      x_ <- x[,j]
      mi_ <- x_[i]*exp(x_%*%xi)/(1+exp(x_%*%xi))-y_*x_[i]
      mi <- mi+n*mi_
    }
    return(max(0,theta[i]*mi/B))
  }
  
  # x和y
  # 混合分布(1-alpha)0 + alpha标准正态 ，alpha=0.1
  alpha <- 1
  # 从混合分布中抽取500个x
  x <- array( rnorm(p*n,0,1) ,dim=c(p,n)) # rnorm 分别替换均匀分布runif和拉普拉斯分布rlaplace即可，rlaplace需要library(VGAM)
  a <- array( as.integer(runif(p*n,0,1) < alpha) ,dim=c(p,n))
  x <- x*a
  # 根据xi真值产生逻辑回归y
  y <- as.integer( runif(n,0,1) < 1/(1+exp(-xi_true%*%x)) )
  
  # improved sub-sampling zigzag过程
  # 跑k步
  k = 3*10^7
  # 初始值
  xi_seq <- array(0,dim=c(p,k))
  xi_seq[,1] <- xi
  t_seq <- vector()
  t_seq[1] <- 0
  
  # theta初始值
  theta <- rep(1,p)
  
  # 算法参数mini_batch
  B = 1
  
  # 记录theta正负号反转次数
  bouncing_time <- 0
  
  # improved zigzag
  for (step in 2:k){
    t_prior <- generate_t_mi_o(runif(p,0,1),xi_seq[,step-1])
    t_condition <- generate_t_mi_(runif(p,0,1),x)
    i_p <- which.min(t_prior)
    i_c <- which.min(t_condition)
    if(t_prior[i_p]<t_condition[i_c]){
      i <- i_p
      t <- t_prior[i]
      # 更新序列
      t_seq[step] <- t_seq[step-1]+t
      xi_seq[,step] <- xi_seq[,step-1]+t*theta
      # 更新theta
      a <- runif(1,0,1)
      if(a < mi_o(theta,xi_seq[,step],i)/Mi_o(xi_seq[,step-1],i,t) ){
        theta[i] <- theta[i]*-1
        bouncing_time <- bouncing_time+1
      }
    }else{
      i <- i_c
      t <- t_condition[i]
      # 更新序列
      t_seq[step] <- t_seq[step-1]+t
      xi_seq[,step] <- xi_seq[,step-1]+t*theta
      # 更新theta
      a <- runif(1,0,1)
      if(a < mi_(theta,y,x,xi_seq[,step],i,B)/Mi_(x,i) ){
        theta[i] <- theta[i]*-1
        bouncing_time <- bouncing_time+1
      }
    }
    if (bouncing_time == 10^5){
      # theta正负号反转10^5次后停止
      break
    }
  }
  T_imp[n_i] <- max(t_seq)
  print('improved zigzag total simulation time:')
  print(max(t_seq))
  print('bouncing time:')
  print(bouncing_time)
  
  # unif sub-sampling zigzag过程
  # 跑k步
  k = 3*10^7
  # 初始值
  xi_seq <- array(0,dim=c(p,k))
  xi_seq[,1] <- xi
  t_seq <- vector()
  t_seq[1] <- 0
  
  # theta初始值
  theta <- rep(1,p)
  
  # 算法参数mini_batch
  B = 1
  
  # 记录theta正负号反转次数
  bouncing_time <- 0
  # normal zigzag
  for (step in 2:k){
    t_prior <- generate_t_mi_o(runif(p,0,1),xi_seq[,step-1])
    t_condition <- generate_t_mi_(runif(p,0,1),x)
    i_p <- which.min(t_prior)
    i_c <- which.min(t_condition)
    if(t_prior[i_p]<t_condition[i_c]){
      i <- i_p
      t <- t_prior[i]
      # 更新序列
      t_seq[step] <- t_seq[step-1]+t
      xi_seq[,step] <- xi_seq[,step-1]+t*theta
      # 更新theta
      a <- runif(1,0,1)
      if(a < mi_o(theta,xi_seq[,step],i)/Mi_o(xi_seq[,step-1],i,t) ){
        theta[i] <- theta[i]*-1
        bouncing_time <- bouncing_time+1
      }
    }else{
      i <- i_c
      t <- t_condition[i]
      # 更新序列
      t_seq[step] <- t_seq[step-1]+t
      xi_seq[,step] <- xi_seq[,step-1]+t*theta
      # 更新theta
      a <- runif(1,0,1)
      if(a < mi_n(theta,y,x,xi_seq[,step],i,B)/Mi_n(x,i) ){
        theta[i] <- theta[i]*-1
        bouncing_time <- bouncing_time+1
      }
    }
    if (bouncing_time == 10^5){
      # theta正负号反转10^5次后停止
      break
    }
  }
  T_unif[n_i] <- max(t_seq)
  print('unif zigzag total simulation time:')
  print(max(t_seq))
  print('bouncing time:')
  print(bouncing_time)
}

# 画图
plot(log(x_lab,10), log(T_unif/T_imp), type='o', lwd=2, pch=16, cex=1.2, col='blue',
     xaxt='n', yaxt='n',xlim=c(1,4), ylim=c(log(1.8),log(10)), xlab='observation number', ylab='Efficiency gain')
lines(log(seq(10,10000,2),10),log(log(seq(10,10000,2))),lwd=2,col='gray')
lines(log(seq(10,10000,2),10),log(log(seq(10,10000,2)/log(seq(10,10000,2)))*(1+1/log(seq(10,10000,2)))),lwd=2,col='gray')
legend('topright',legend=c('laplace','gaussian','uniform'),col=c('red','blue','green'),pch=16,lwd=2)
axis(1,at=c(1,2,3,4),labels=c('10^1','10^2','10^3','10^4'))
axis(2,at=log(c(1,2,4,6,8,10)),labels=c(1,2,4,6,8,10))
grid()