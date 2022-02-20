### 对高维稀疏imbalance的数据采用pre-conditioning 的 importance sub-sampling
### 该文件中的内容是对论文4.3节的复现
###
### 由于原论文中采用10^4维度，10^6个观测进行实验，计算量太大，代码无法运行
### 可以考虑降低数据维度p和观测数量n后，再运行
### 该部分内容涉及论文补充材料S7部分
###
### 本文件是对论文4.3部分的实验进行复现
### pre-conditioning 的 importance sub-sampling zigzag process代码部分将会单独放在
### 4.3高维稀疏pre_conditioning.txt中,方便使用
### 若仅想查看pre-conditioning 的 importance sub-sampling算法请直接前往txt
###
###
######################以下为代码部分#######################

set.seed(114514)
# x维度10^4，样本量10^6
p = 10^4
n = 10^6

# xi真值
xi_true <- rep(-1/2,p)

# x和y
# 混合分布(1-alpha)0 + alpha标准正态 ，alpha=0.01
alpha = 0.01
# 从混合分布中抽取500个x
x <- array( rnorm(p*n,0,1) ,dim=c(p,n))
a <- array( as.integer(runif(p*n,0,1) < alpha) ,dim=c(p,n))
x <- x*a
# 根据xi真值产生逻辑回归y
y <- as.integer( runif(n,0,1) < 1/(1+exp(-xi_true%*%x)) )

# xi先验初值
xi <- rep(2,p)
# xi先验分布为正态N(0,sigma^2)，sigma=1
sigma <- 1


# 函数定义
# 抽取反弹时间t，mi_o部分
generate_t_mi_o <- function(unif,xi){
  unif <- -log(unif)
  t <- (sqrt(unif+abs(xi))-sqrt(abs(xi)))*sqrt(2*sigma^2)
  return(t)
}
# 抽取反弹时间t，mi_o部分pre-conditioning
generate_t_mi_o_p <- function(unif,xi,velocity){
  unif <- -log(unif)
  unif <- unif/velocity
  t <- (sqrt(unif+abs(xi))-sqrt(abs(xi)))*sqrt(2*sigma^2)
  return(t)
}
# 抽取反弹时间t，mi_*部分
generate_t_mi_ <- function(unif,x){
  unif <- -log(unif)
  t <- unif/rowSums(abs(x))
  return(t)
}
# 抽取反弹时间t，mi_*部分pre-conditioning
generate_t_mi_p <- function(unif,x, velocity){
  unif <- -log(unif)
  t <- unif/(rowSums(abs(x)*velocity))
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
    j <- sample(1:n,1,prob=abs(x[i,]))  # ,prob=abs(x[i,]) 和 等概 相互替换
    y_ <- y[j]
    x_ <- x[,j]
    mi_ <- x_[i]*exp(x_%*%xi)/(1+exp(x_%*%xi))-y_*x_[i]
    mi <- mi+sum(abs(x[i,]))/abs(x_[i])*mi_  # sum(abs(x[i,]))/abs(x_[i]) 和 n 替换
  }
  return(max(0,theta[i]*mi/B))
}
# mi_*上界
Mi_ <- function(x,i){
  return(sum(abs(x[i,])))  # sum(abs(x[i,])) 和 n*max(abs(x[i,])) 替换
}


# improved sub-sampling zigzag过程
# 跑k步
k = 10^5
# 初始值
xi_seq <- array(0,dim=c(p,k))
xi_seq[,1] <- xi
t_seq <- vector()
t_seq[1] <- 0
acf_xi_seq <- array(xi,dim=c(p,1))

# theta初始值
theta <- rep(1,p)

# 算法参数mini_batch
B = 10

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
      acf_xi_seq <- cbind(acf_xi_seq,xi_seq[,step])
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
      acf_xi_seq <- cbind(acf_xi_seq,xi_seq[,step])
    }
  }
}
print(max(t_seq))
t_seq_un <- t_seq
xi_seq_un <- xi_seq
acf_xi_seq_un <- acf_xi_seq

# pre-conditioning improved sub-sampling zigzag过程
# 跑k步
k = 10^5
# 初始值
xi_seq <- array(0,dim=c(p,k))
xi_seq[,1] <- xi
t_seq <- vector()
t_seq[1] <- 0
acf_xi_seq <- array(xi,dim=c(p,1))

# theta初始值
theta <- rep(1,p)

# 算法参数mini_batch
B = 10

# pre-conditioning
velocity <- rep(1/p, p)  # d=p

# 用于计算std，从而决定velocity
u_1_t <- 0
u_2_t <- 0
  
# zigzag
for (step in 2:k){
  t_prior <- generate_t_mi_o_p(runif(p,0,1),xi_seq[,step-1],velocity)
  t_condition <- generate_t_mi_p(runif(p,0,1),x,velocity)
  i_p <- which.min(t_prior)
  i_c <- which.min(t_condition)
  if(t_prior[i_p]<t_condition[i_c]){
    i <- i_p
    t <- t_prior[i]
    # 更新序列
    t_seq[step] <- t_seq[step-1]+t
    xi_seq[,step] <- xi_seq[,step-1]+t*theta*velocity
    # 更新u_12_t
    u_1_t <- u_1_t+xi_seq[,step-1]*t+1/2*t^2*theta*velocity
    u_2_t <- u_2_t+xi_seq[,step-1]^2*t+t^2*theta*velocity*xi_seq[,step-1]+1/3*t^3*(theta*velocity)^2
    # 更新velocity
    sd_t <- sqrt(u_2_t/t_seq[step]-(u_1_t/t_seq[step])^2)
    velocity <- p*sd_t/sum(sd_t)
    # 更新theta
    a <- runif(1,0,1)
    if(a < mi_o(theta,xi_seq[,step],i)/Mi_o(xi_seq[,step-1],i,t) ){
      theta[i] <- theta[i]*-1
      acf_xi_seq <- cbind(acf_xi_seq,xi_seq[,step])
    }
  }else{
    i <- i_c
    t <- t_condition[i]
    # 更新序列
    t_seq[step] <- t_seq[step-1]+t
    xi_seq[,step] <- xi_seq[,step-1]+t*theta*velocity
    # 更新u_12_t
    u_1_t <- u_1_t+xi_seq[,step-1]*t+1/2*t^2*theta*velocity
    u_2_t <- u_2_t+xi_seq[,step-1]^2*t+t^2*theta*velocity*xi_seq[,step-1]+1/3*t^3*(theta*velocity)^2
    # 更新velocity
    sd_t <- sqrt(u_2_t/t_seq[step]-(u_1_t/t_seq[step])^2)
    velocity <- p*sd_t/sum(sd_t)
    # 更新theta
    a <- runif(1,0,1)
    if(a < mi_(theta,y,x,xi_seq[,step],i,B)/Mi_(x,i) ){
      theta[i] <- theta[i]*-1
      acf_xi_seq <- cbind(acf_xi_seq,xi_seq[,step])
    }
  }
}
print(max(t_seq))

par(mfrow=c(2,2))

plot(acf_xi_seq[1,],type = 'l',main='pre-conditioning')
acf(acf_xi_seq[1,],lag.max=200)

plot(acf_xi_seq_un[1,],type = 'l',main='un pre-conditioning')
acf(acf_xi_seq_un[1,],lag.max=200)

