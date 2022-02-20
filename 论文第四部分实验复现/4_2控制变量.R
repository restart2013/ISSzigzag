### control variates 的 importance sub-sampling算法
### 与 uncontrol variates 进行比较
### 该文件中的内容是对论文4.2节的复现
###
### 本文件提供了control variates 的 importance sub-sampling zigzag process代码
### 并比较了与uncontrol variates 情况下的efficiency
### 涉及论文补充材料S3的部分
### 
### 由于本文件是对control和uncontrol两种方法比较的复现，
### control variates 的 importance sub-sampling zigzag process代码部分将会单独放在
### 4.2控制变量.txt中,方便使用
### 若仅想查看control variates 的 importance sub-sampling算法请直接前往txt
### 
### 由于代码运行时间较长，实验结果直接记录在了xx行-xx行，运行后便能查看
### 正式代码过程在xx行之后
###
###
######################以下为代码部分#######################
library(stats4)
library(splines)
library(VGAM)
library(matrixStats)

set.seed(114514)
# x维度10，样本量5000
p = 10
n = 5000

# x和y
# 混合分布(1-alpha)0 + alpha标准正态 ，alpha=0.1
alpha = 0.1
# 从混合分布中抽取500个x
x <- array( rlaplace(p*n,0,1) ,dim=c(p,n))
a <- array( as.integer(runif(p*n,0,1) < alpha) ,dim=c(p,n))
x <- x*a
# 根据xi真值产生逻辑回归y
y <- rep(0,n)


# xi先验初值
xi <- rnorm(p,0,1)
xi_star <- rnorm(p,0,1)
# xi先验分布为正态N(0,sigma^2)，sigma=10^5
sigma <- 1

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
# 抽取反弹时间t，mi_*部分控制变量
generate_t_mi_c <- function(unif,x,y,xi,theta){
  unif <- -log(unif)
  A <- 1/4*colSums(t(abs(x))*sqrt(colSums(x^2)))*sqrt(p)
  B <- theta*colSums(t(x)*as.vector(exp(t(x)%*%xi_star)/(1+exp(t(x)%*%xi_star)))-y*t(x))
  B <- rowMaxs(cbind(rep(0,p),B))+1/4*colSums(t(abs(x))*sqrt(colSums(x^2)))*sqrt(sum((xi-xi_star)^2))
  C <- B/(2*sqrt(A))
  t <- (sqrt(unif+C^2)-C)/A
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
# rate函数mi_*控制变量
mi_c <- function(theta,y,x,xi,i,B){
  mi <- 0
  for (b in 1:B){
    j <- sample(1:n,1,prob=abs(x[i,])*sqrt(colSums(x^2)))
    y_ <- y[j]
    x_ <- x[,j]
    mi_ <- x_[i]*exp(x_%*%xi)/(1+exp(x_%*%xi))-y_*x_[i]
    mi_star <- x_[i]*exp(x_%*%xi_star)/(1+exp(x_%*%xi_star))-y_*x_[i]
    mi <- mi+sum(abs(x[i,])*sqrt(colSums(x^2)))/(abs(x[i,])*sqrt(colSums(x^2)))[j]*(mi_-mi_star)
  }
  mi <- mi/B + sum(x[i,]*exp(t(x)%*%xi_star)/(1+exp(t(x)%*%xi_star))-y*x[i,])
  return(max(0,theta[i]*mi))
}
# mi_*上界
Mi_ <- function(x,i){
  return(sum(abs(x[i,])))
}
# mi_*上界控制变量
Mi_c <- function(x,y,theta,i,xi,t){
  Mi <- theta[i]*sum(x[i,]*exp(t(x)%*%xi_star)/(1+exp(t(x)%*%xi_star))-y*x[i,])
  Mi <- max(0,Mi)+1/4*sum(abs(x[i,])*sqrt(colSums(x^2)))*(sqrt(sum((xi-xi_star)^2))+t*sqrt(p))
  return(Mi)
}


# efficiency比较
T_uncon <- vector()
T_con <- vector()

k_number <- c(1, 4, 18, 110, 700, 2000)


for (k_i in 6:6){
  for (y_i in 1:k_number[k_i]){  # k分别等于1，4，18，110，700，2000的情况
    y[y_i] <- 1
  }
  print(sum(y))

  # 非控制变量sub-sampling zigzag过程
  # 跑k步
  k = 10^7
  # 初始值
  xi_seq <- array(NA,dim=c(p,k))
  xi_seq[,1] <- xi
  t_seq <- vector()
  t_seq[1] <- 0
  acf_xi_seq <- xi
  acf_t_seq <- t_seq
  
  # theta初始值
  theta <- rep(1,p)
  
  # 算法参数mini_batch
  B = 1
  
  # 记录theta正负号反转次数
  bouncing_time <- 0
  
  # zigzag
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
        acf_xi_seq <- cbind(acf_xi_seq,xi_seq[,step])
        acf_t_seq[bouncing_time] <- t_seq[step]
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
        acf_xi_seq <- cbind(acf_xi_seq,xi_seq[,step])
        acf_t_seq[bouncing_time] <- t_seq[step]
      }
    }
    if (bouncing_time == 4*10^4){
      # theta正负号反转4*10^4次后停止
      break
    }
  }
  t_candidate <- vector()
  for (p_i in 1:p){
    t_candidate[p_i] <- acf_t_seq[which(abs(acf(acf_xi_seq[p_i,],lag.max=40000,plot=F)$acf)<0.02)[1]]
  }
  T_uncon[k_i] <- max(t_candidate)
  print('uncontrol variates slowest mixing time:')
  print(T_uncon[k_i])
  print('bouncing time:')
  print(bouncing_time)
  
  
  # 控制变量sub-sampling zigzag过程
  # 跑k步
  k = 10^7
  # 初始值
  xi_seq <- array(NA,dim=c(p,k))
  xi_seq[,1] <- xi
  t_seq <- vector()
  t_seq[1] <- 0
  acf_xi_seq <- xi
  acf_t_seq <- t_seq
  
  # theta初始值
  theta <- rep(1,p)
  
  # 算法参数mini_batch
  B = 1
  
  # 记录theta正负号反转次数
  bouncing_time <- 0
  
  # zigzag
  for (step in 2:k){
    t_prior <- generate_t_mi_o(runif(p,0,1),xi_seq[,step-1])
    t_condition <- generate_t_mi_c(runif(p,0,1),x,y,xi_seq[step-1],theta)
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
        acf_xi_seq <- cbind(acf_xi_seq,xi_seq[,step])
        acf_t_seq[bouncing_time] <- t_seq[step]
      }
    }else{
      i <- i_c
      t <- t_condition[i]
      # 更新序列
      t_seq[step] <- t_seq[step-1]+t
      xi_seq[,step] <- xi_seq[,step-1]+t*theta
      # 更新theta
      a <- runif(1,0,1)
      if(a < mi_c(theta,y,x,xi_seq[,step],i,B)/Mi_c(x,y,theta,i,xi_seq[step-1],t) ){
        theta[i] <- theta[i]*-1
        bouncing_time <- bouncing_time+1
        acf_xi_seq <- cbind(acf_xi_seq,xi_seq[,step])
        acf_t_seq[bouncing_time] <- t_seq[step]
      }
    }
    if (bouncing_time == 4*10^4){
      # theta正负号反转4*10^4次后停止
      break
    }
  }
  plot(acf_xi_seq[1,],type='l')
  t_candidate <- vector()
  for (p_i in 1:p){
    t_candidate[p_i] <- acf_t_seq[which(abs(acf(acf_xi_seq[p_i,],lag.max=40000,plot=F)$acf)<0.02)[1]]
  }
  T_con[k_i] <- max(t_candidate)
  print('control variates slowest mixing time:')
  print(T_con[k_i])
  print('bouncing time:')
  print(bouncing_time)
}


# 实验结果
# 画图
par(mfrow=c(1,1))
plot(log(k_number,10),log(T_con/T_uncon,10), type='o', lwd=2, pch=16, cex=1.2, col='blue',
     xaxt='n', yaxt='n',xlim=c(0,4), ylim=c(-1,1), xlab='k', ylab='Efficiency gain')
axis(1,at=c(0,1,2,3,4),labels=c('10^0','10^1','10^2','10^3','10^4'))
axis(2,at=c(-1,0,1),labels=c('10^-1','10^0','10^1'))
grid()