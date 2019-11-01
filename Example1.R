

rm(list = ls())

working_directory = "/Users/oscar/Documents/GitHub/MNWBS"
setwd(working_directory)
source("utils2.R")
library(ecp)
#setwd(paste(working_directory,"/Code_J",sep=""))
#source("utils_functions.R")
#source('utils.R')
library(MASS)
library(ks)
library(kdensity)
library(moments)

T_grid  = c(300,150)
#v =  c(150,300)
p_grid   = c(20,10)
#y =  matrix(0,T,p)
#z =  matrix(0,T,p)
#u = matrix(0,T,p)

NMC = 100

ind_t = 1
ind_p = 1


hausdorff_wbs =   array(0,c(NMC,  length(T_grid),  length(p_grid)))
hausdorff_energy =  array(0,c(NMC,  length(T_grid),  length(p_grid)))

hausdorff2_wbs =   array(0,c(NMC,  length(T_grid),  length(p_grid)))
hausdorff2_energy =  array(0,c(NMC,  length(T_grid),  length(p_grid)))

error_wbs =   array(0,c(NMC,  length(T_grid),length(p_grid)))
error_energy =   array(0,c(NMC,  length(T_grid),length(p_grid)))

for(ind_t in 1:length(T_grid))
{
  ####
  T =  T_grid[ind_t]
  print("T")
  print(T)
  
  for(ind_p in 1:length(p_grid)  )
  {
    p = p_grid[ind_p]
    
    print("p")
    print(p)
    
    v =  c(floor(T/3),2*floor(T/3))
    
    y =  matrix(0,T,p)
    z =  matrix(0,T,p)
    
    mu0 =   rep(0,p)
    mu1 = rep(0,p)
    #mu2 = rep(.0,p)
    mu1[1:floor(p/2)] =  1
    Sigma0 =  diag(p)
    Sigma1 = diag(p)
    #matrix(0.9,p,p)
    #  diag(Sigma1) =1    
    #Sigma1 =   diag(p)#0.5*diag(p) + 0.5*matrix(1,p,p) 
    
    
    for(iter in 1:NMC)
    {
      print("iter")
      for(t in 1:T)
      {
        if(t <v[1] ||  t > v[2])
        {
          y[t,] = mvrnorm(n = 1, mu0, Sigma0)
          z[t,] = mvrnorm(n = 1, mu0,  Sigma0)
        }
        
        if(t >=v[1] &&  t < v[2])
        {
          y[t,] = mvrnorm(n = 1, mu1,Sigma1)
          z[t,] = mvrnorm(n = 1, mu1,Sigma1)
          # y[t,] = rmvt(sigma = Sigma0, n=1,df=3)/sqrt(3)#mvrnorm(n = 1, mu1, Sigma1)
          #z[t,] = rmvt(sigma = Sigma0, n=1,df=3)/sqrt(3)#mvrnorm(n = 1, mu1, Sigma1)
        }
      }## close for generate data
      
      
      ####  
      matplot(y,type="l")
      K_max = 30
      h = 5*(K_max*log(T)/T)^{1/p}  
      #5*8^{1/p}/(T/K_max)^{1/p}
      
      
      ##
      M =   50
      alpha =  sample.int(size =M  , n = T,replace = TRUE)
      beta =   sample.int(size =M  , n = T,replace = TRUE)#alpha + floor((T- alpha)*runif(M))
      #
      for(j in 1:M)
      {
        aux =  alpha[j]
        aux2 =  beta[j]
        #
        alpha[j] = min(aux,aux2)
        beta[j] = max(aux,aux2)
      }###  close for intervals
      
      S =MNWBS_full(y,y,alpha,beta,h)
      S
      hausdorff_wbs[iter,ind_t,ind_p] =  dist_change_points(S,v)
      hausdorff2_wbs[iter,ind_t,ind_p] =  dist_change_points(v,S)
      error_wbs[iter,ind_t,ind_p] =  length(v) - length(S) 
      
      
      u = matrix(0,2*T,p)
      u[2*(1:T),  ] = y
      u[2*(1:T)-1,  ] = z
      
      temp = e.divisive(y)
      temp$estimates = setdiff(temp$estimates,c(1,T+1))
      #      temp = e.cp3o_delta(Z= u, K=10, delta=10, alpha=1, verbose=FALSE)
      #temp$estimates =  floor(temp$estimates/2)
      hausdorff_energy[iter,ind_t,ind_p] =  dist_change_points(temp$estimates,v)
      hausdorff2_energy[iter,ind_t,ind_p] =  dist_change_points(v,temp$estimates)
      error_energy[iter,ind_t,ind_p] =  length(v) - length(temp$estimates) 
      
      ####################################33
      # temp2 = kcpa(y)
      
      #############################
      print("energy")
      print(median(hausdorff_energy[iter,ind_t,ind_p]))
      print(median(hausdorff2_energy[iter,ind_t,ind_p]))
      print(mean(abs(error_energy[iter,ind_t,ind_p])))
      
      print("NP")
      print(median(hausdorff_wbs[iter,ind_t,ind_p]))
      print(median(hausdorff2_wbs[iter,ind_t,ind_p]))
      print(mean(abs(error_wbs[iter,ind_t,ind_p])))
    }### close for iter
    
    print("energy")
    print(median(hausdorff_energy[,ind_t,ind_p]))
    print(median(hausdorff2_energy[,ind_t,ind_p]))
    print(mean(abs(error_energy[,ind_t,ind_p])))
    
    print("NP")
    print(median(hausdorff_wbs[,ind_t,ind_p]))
    
    print(median(hausdorff2_wbs[,ind_t,ind_p]))
    print(mean(abs(error_wbs[,ind_t,ind_p])))
  }## close for p
}## close for T

#pdf("plot_Example_0.pdf")
#matplot(y,type="l",xlab="Time",ylab = "",main ="",cex.lab = 1.7,cex.main= 1.7)
#dev.off()
