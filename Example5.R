
rm(list = ls())

#working_directory = "/Users/oscar/Desktop/Anomaly_detection/supplementary_materials/supplementary_materials/Experiments"
working_directory = "/Users/oscar/Documents/GitHub/MNWBS"
#setwd(paste(working_directory,"/Code_J",sep=""))
source("auxiliary_functions.R")
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
    
    temp = radiological_example(1, celsium  = 1,T)
    background_train =  temp$background_train
    f0  = temp$f0
    m =  dim(f0)[2]
    
    for(iter in 1:NMC)
    {
      print("iter")
      for(t in 1:T)
      {
        if(t <v[1] ||  t > v[2])
        {
          y[t,] =   sample(x= (1:m)/m,size= p, prob = f0[1,]) 
          #mvrnorm(n = 1, mu0, Sigma0)
        }
        
        if(t >=v[1] &&  t < v[2])
        {
          #floor(p/5)
          a = 2
          y[t,1:(a)  ] = sample(x= (1:m)/m,size= a, prob = f0[2,])
          y[t,(a+1):p  ] = sample(x= (1:m)/m,size=p- a, prob = f0[1,])
          #mvrnorm(n = 1, mu1,Sigma1)
          
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
      #alpha[M] = 1
      S =MNWBS_full(y,y,alpha,beta,h)
      S
      #if(abs(lenght(S) - 2)>0)
      #{return(NULL)}
      hausdorff_wbs[iter,ind_t,ind_p] =  dist_change_points(S,v)
      hausdorff2_wbs[iter,ind_t,ind_p] =  dist_change_points(v,S)
      error_wbs[iter,ind_t,ind_p] =  length(v) - length(S) 
      
      
      u = matrix(0,2*T,p)
      u[2*(1:T),  ] = y
      u[2*(1:T)-1,  ] = z
      
      temp = e.divisive(y)
      temp$estimates = setdiff(temp$estimates,c(1,T+1))
      temp$estimates
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


pdf("plot_Example_6.pdf")
matplot(y,type="l",xlab="Time",ylab = "",main ="",cex.lab = 1.7,cex.main= 1.7)
dev.off()
