# ISSzigzag
R code for paper arXiv:1905.11232v2 (Efficient posterior sampling for high-dimensional imbalanced logistic regression)

ISSzigzag即importance subsampling zigzag
重要性下采样zigzag算法，用来解决zigzag算法在高维稀疏数据情况下，难以收敛到稳态分布的问题


本文件夹中为论文代码复现
使用编程语言为R

代码复现了：
1.论文所提出的算法
2.论文第4节实验内容
3.由于论文第5节实验内容所使用数据引用自其他论文，这里不予复现

共有四个文件：
算法实现.R
论文第4部分实验复现
	4_1efficiency比较.R（论文4.1节的improved与unif两种算法的efficiency比较）
	4_2控制变量.R（论文4.2节的control variates与uncontrol两种算法的efficiency比较）
	4_3高维稀疏pre_conditioning.R（论文4.3节的高维稀疏非平衡数据下采用pre_coditioning的算法）

为防止中文编码解码不匹配，
三份代码除了提供 R script 文件，也提供了txt备份。
所有代码前均有相关注释说明，请仔细阅读。


论文：Efficient posterior sampling for high-dimensional imbalanced logistic regression
