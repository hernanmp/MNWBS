
rm(list = ls())

#working_directory = "/Users/oscar/Desktop/Anomaly_detection/supplementary_materials/supplementary_materials/Experiments"
working_directory = "C:/Users/DELL/Desktop/change_point_code/MNWBS"
setwd(working_directory)

source("auxiliary_functions.R")
#working_directory = "/Users/oscar/Documents/GitHub/MNWBS"
#setwd(working_directory)
source("utils2.R")

library(ecp)
#working_directory = "/Users/oscar/Desktop/change_point_code/NWBS/"
#setwd(paste(working_directory,"/Code_J",sep=""))
#source("utils_functions.R")
#source('utils.R')
library(MASS)
library(ks)
library(kdensity)
library(moments)
#install.packages("kernseg")
library(KernSeg)
library(hdbinseg)
library(mvtnorm)

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
hausdorff_kcp =  array(0,c(NMC,  length(T_grid),  length(p_grid)))
hausdorff_sbs =  array(0,c(NMC,  length(T_grid),  length(p_grid)))
hausdorff_dcbs =  array(0,c(NMC,  length(T_grid),  length(p_grid)))

hausdorff2_wbs =   array(0,c(NMC,  length(T_grid),  length(p_grid)))
hausdorff2_energy =  array(0,c(NMC,  length(T_grid),  length(p_grid)))
hausdorff2_kcp =  array(0,c(NMC,  length(T_grid),  length(p_grid)))
hausdorff2_sbs =  array(0,c(NMC,  length(T_grid),  length(p_grid)))
hausdorff2_dcbs =  array(0,c(NMC,  length(T_grid),  length(p_grid)))

error_wbs =   array(0,c(NMC,  length(T_grid),length(p_grid)))
error_energy =   array(0,c(NMC,  length(T_grid),length(p_grid)))
error_kcp =   array(0,c(NMC,  length(T_grid),length(p_grid)))
error_sbs =   array(0,c(NMC,  length(T_grid),length(p_grid)))
error_dcbs =   array(0,c(NMC,  length(T_grid),length(p_grid)))

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
    
    v =  c(floor(T/7),2*floor(T/7),3*floor(T/7),4*floor(T/7),5*floor(T/7),6*floor(T/7))
    
    
    y =  matrix(0,T,p)
    z =  matrix(0,T,p)
    
    mu0 =   rep(0,p)
    mu1 = rep(0.2,p)
    # mu1[1:p]= 0.5 #:floor(p/5) ] = 1
    #[1] =  1
    Sigma0 =  diag(p)
    Sigma1 =  diag(p)
    
    for(iter in 1:NMC)
    {
      print("iter")
      for(t in 1:T)
      {
        if( floor(t/floor(T/7))%%2==0  || t>=7*floor(T/7))
        {
          y[t,] =  mu0 +  rmvt(n=1,sigma = Sigma0,df=3)/sqrt(3)
          #mvrnorm(n = 1, mu0, Sigma0)
          z[t,] = mu0 +  rmvt(n=1,sigma = Sigma0,df=3)/sqrt(3)
          #mvrnorm(n = 1, mu0, Sigma0)
        }
        
        
        if(floor(t/floor(T/7))%%2==1 &&  t<7*floor(T/7))
        {
          y[t,] = mu1 +  rmvt(n=1,sigma = Sigma1,df=3)/sqrt(3)
          #mvrnorm(n = 1, mu1, Sigma1)
          z[t,] =   mu1 +  rmvt(n=1,sigma = Sigma1,df=3)/sqrt(3)
          #mvrnorm(n = 1, mu1, Sigma1)
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
      #  S =  NULL
      hausdorff_wbs[iter,ind_t,ind_p] =  dist_change_points(S,v)
      hausdorff2_wbs[iter,ind_t,ind_p] =  dist_change_points(v,S)
      error_wbs[iter,ind_t,ind_p] =  length(v) - length(S) 
      
      
      u = matrix(0,2*T,p)
      u[2*(1:T),  ] = y
      u[2*(1:T)-1,  ] = z
      
      
      temp = e.divisive(y)
      temp$estimates = setdiff(temp$estimates,c(1,T+1))
      #e.cp3o_delta(Z= u, K=10, delta=10, alpha=1, verbose=FALSE)
      #temp$estimates =  floor(temp$estimates/2)
      hausdorff_energy[iter,ind_t,ind_p] =  dist_change_points(temp$estimates,v)
      hausdorff2_energy[iter,ind_t,ind_p] =  dist_change_points(v,temp$estimates)
      error_energy[iter,ind_t,ind_p] =  length(v) - length(temp$estimates) 
      
      ###########################################
      
      ####################################33
      
      
      aux = KernSeg_MultiD(y, Kmax=20, delta = 30, min.size = 2, 
                           alpha = NULL, kernel = "Gaussian", option = 0)
      
      
      best_ind = which.max(abs(diff(aux$J.est)))
      est_cp = aux$t.est[best_ind,]
      est_cp =  sort(setdiff(aux$t.est[best_ind,],c(0,1,T)))
      #est_cp
      
      hausdorff_kcp[iter,ind_t,ind_p] =  dist_change_points(est_cp,v)
      hausdorff2_kcp[iter,ind_t,ind_p] =  dist_change_points(v,est_cp)
      error_kcp[iter,ind_t,ind_p] =  length(v) - length(est_cp) 
      
      ############################################################# 
      aux = sbs.alg(t(y),do.parallel=0)
      est_sbs= sort(setdiff(aux$ecp,c(1,T) ))
      
      hausdorff_sbs[iter,ind_t,ind_p] =  dist_change_points(est_sbs,v)
      hausdorff2_sbs[iter,ind_t,ind_p] =  dist_change_points(v,est_sbs)
      error_sbs[iter,ind_t,ind_p] =  length(v) - length(est_sbs)
      ############################################################# 
      
      aux = dcbs.alg(t(y), cp.type=1, phi=-1, temporal=FALSE, do.parallel=0)$ecp
      est_dcbs = sort(setdiff(aux,c(1,T) ))
      
      hausdorff_dcbs[iter,ind_t,ind_p] =  dist_change_points(est_dcbs,v)
      hausdorff2_dcbs[iter,ind_t,ind_p] =  dist_change_points(v,est_dcbs)
      error_dcbs[iter,ind_t,ind_p] =  length(v) - length(est_dcbs)
      ############################################################# 
      # print("energy")
      # print(median(hausdorff_energy[iter,ind_t,ind_p]))
      # print(median(hausdorff2_energy[iter,ind_t,ind_p]))
      # print(mean(abs(error_energy[iter,ind_t,ind_p])))
      # 
      # print("NP")
      # print(median(hausdorff_wbs[iter,ind_t,ind_p]))
      # print(median(hausdorff2_wbs[iter,ind_t,ind_p]))
      # print(mean(abs(error_wbs[iter,ind_t,ind_p])))
    }### close for iter
    
    
    print("energy")
    print(median(hausdorff_energy[,ind_t,ind_p]))
    print(median(hausdorff2_energy[,ind_t,ind_p]))
    print(mean(abs(error_energy[,ind_t,ind_p])))
    
    print("NP")
    print(median(hausdorff_wbs[,ind_t,ind_p]))
    print(median(hausdorff2_wbs[,ind_t,ind_p]))
    print(mean(abs(error_wbs[,ind_t,ind_p])))
    
    print("kcp")
    print(median(hausdorff_kcp[,ind_t,ind_p]))
    print(median(hausdorff2_kcp[,ind_t,ind_p]))
    print(mean(abs(error_kcp[,ind_t,ind_p])))
    
    print("sbs")
    print(median(hausdorff_sbs[,ind_t,ind_p]))
    print(median(hausdorff2_sbs[,ind_t,ind_p]))
    print(mean(abs(error_sbs[,ind_t,ind_p])))
    
    print("dcbs")
    print(median(hausdorff_dcbs[,ind_t,ind_p]))
    print(median(hausdorff2_dcbs[,ind_t,ind_p]))
    print(mean(abs(error_dcbs[,ind_t,ind_p])))
  }## close for p
}## close for T
