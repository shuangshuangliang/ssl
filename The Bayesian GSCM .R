rm(list=ls())
library(MCMCpack)
library(Rsolnp)
library(gsynth)
set.seed(12345) #for reproducibility of simulation
sim.n=1000

GSCM.Mse<-rep(0,sim.n)
GSCM_r2.Mse<-rep(0,sim.n)
Bayes.Mse<-rep(0,sim.n)
Bayes1.Mse<-rep(0,sim.n)

GSCM.r<-rep(0,sim.n)
Bayes.r<-rep(0,sim.n)

Bayes.averagelength<-rep(0,sim.n)
Bayes1.averagelength<-rep(0,sim.n)
GSCM.averagelength<-rep(0,sim.n)
GSCMr2.averagelength<-rep(0,sim.n)

for (M in 1:sim.n){
  set.seed(M)
  T<-100
  T0<-60
  T11<-T0+1
  T1<-T-T0
  K<-4
  p<-4 #number of total covariates
  N<-40 #number of units
  Ntrt<-1
  Nco<-N-Ntrt
  
  theta0<-1
  
  Y<-matrix(0,N,T)
  
  alpha.i<-rnorm(N,0,sqrt(0.5))
  
  gamma0<-0
  gamma.sd<-sqrt(0.1)
  gamma.t<-numeric(T)
  gamma.t[1]<-gamma0
  for (t in 2:T){
    gamma.t[t]<-rnorm(1, gamma.t[t-1], gamma.sd)
  }
  c<-c(1,1,1,1)
  
  Z<-array(rnorm(p*N*T,1,sqrt(2)),c(N,p,T))
  True.ATT<-rep(0,T)
  
  for(i in 1:N){
    for(t in 1:T){
      D.it<-1*((i==1)&(t>T0))
      Y[i,t]<-alpha.i[i]+gamma.t[t]+sum(c*Z[i,,t])+(0.5+sqrt(t/2))*D.it+rnorm(1,0,sqrt(1))
      if(i==1){
        True.ATT[t]<-theta0*(0.5+sqrt(t/2))*D.it}
    }
  }
  
  
  
  X<-array(0,c(T,K,N))
  for (i in 1:N){
    X[,,i]<-t(Z[i,,])
  }
  
  
  Xpre<-X[1:T0,,Ntrt]
  Ypre<-Y[Ntrt,1:T0]
  
  txx<-t(X[,,2])%*%X[,,2]
  for (i in 3:N){
    txx<-txx+t(X[,,i])%*%X[,,i]
  }
  itxx<-solve(txx)
  
  
  ########################################
  ########################################
  #Bayesian R selection ##################
  ########################################
  
  xi0<-1
  v0<-1000
  C0<-1
  a0<-1
  b0<-1
  R0<-2:9
  log.f.Y<-rep(0,length(R0))
  max.BF<--Inf
  
  ##Iteration scheme
  for(i1 in 1:length(R0)){
    
    r1=R0[i1]
    if(i1==1){
      Betahat0<-matrix(rnorm(K),K,1)
      Fhat0<-matrix(0,T0,r1)
      Fstar0<-matrix(0,T0-r1,r1)
      for (i in 1:(T0-r1)){
        for (j in 1:r1){
          Fstar0[i,j]<-rnorm(1,0,1)
        }
      }
      Fhat0<-rbind(diag(r1),Fstar0)
      
      Lambdahat0<-matrix(rnorm(N*r1),N,r1)
      
      Sigmahat0<-diag(1,T0) 
      C<-Lambdahat0%*%t(Fhat0)
    }else{
      Betahat0<-Betahat
      Lambdahat0<-C[,1:r1]
      Fstar<-t(solve(t(Lambdahat0)%*%Lambdahat0+0.1*diag(r1))%*%t(Lambdahat0)%*%C[,-(1:r1)])
      Fhat0<-rbind(diag(r1),Fstar) #
      Sigmahat0<-Sigmahat
    }
    
    
    esp<-0.00001
    U<-5000
    Fhat<-Fhat0
    Lambdahat<-Lambdahat0
    Sigmahat<-Sigmahat0
    
    for(j in 1:U){
      
      #inv.Sigmahat<-diag(1/diag(Sigmahat))
      Omega<-xi0*Fhat%*%t(Fhat)+Sigmahat
      inv.Omega<-solve(Omega)
      txx<-t(X[1:T0,,1])%*%inv.Omega%*%X[1:T0,,1]
      txy<-t(X[1:T0,,1])%*%inv.Omega%*%Y[1,1:T0]
      
      for (i in 2:N){
        txx<-txx+t(X[1:T0,,i])%*%inv.Omega%*%X[1:T0,,i]
        txy<-txy+t(X[1:T0,,i])%*%inv.Omega%*%Y[i,1:T0]
      }
      sd.beta<-solve(txx+diag(K)/v0)
      Betahat<-sd.beta%*%txy #MAP of Beta
      
      for (i in 1:N){
        Lambdahat[i,]<-solve(crossprod(Fhat)+Sigmahat[1,1]*diag(r1)/xi0)%*%t(Fhat)%*%(Y[i,1:T0]-X[1:T0,,i]%*%Betahat)
      }
      
      for (t in (r1+1):T0){
        sd.ft<-solve((Sigmahat[t,t])^(-1)*t(Lambdahat)%*%Lambdahat+diag(1/C0,r1))
        Fhat[t,]<-sd.ft%*%((Sigmahat[t,t])^(-1)*t(Lambdahat)%*%(Y[,t]-t(X[t,,])%*%Betahat))
      }
      
      J0<-0
      for (t in 1:T0){
        J<-Y[,t]-t(X[t,,])%*%Betahat-Lambdahat[,]%*%Fhat[t,]
        J0<-J0+sum(J^2)
      }
      
      diag(Sigmahat)<-((J0+b0)/2)/((a0+N*T0)/2+1)
      
      if(abs(Sigmahat[1,1]-Sigmahat0[1,1])<esp&max(abs(Betahat-Betahat0))<esp){break}
      if(j==U){print("No convergence")}
      Betahat0<-Betahat
      Fhat0<-Fhat
      Lambdahat0<-Lambdahat
      Sigmahat0<-Sigmahat
    }
    
    C<-Lambdahat%*%t(Fhat) #C updated
    
    #####log(y|r)##########
    Xb0<-matrix(0,N,T0)
    for(i in 1:N){
      Xb0[i,]<-X[1:T0,,i]%*%Betahat
    }
    Y.star<-Y[1:N,1:T0]-Xb0[1:N,]-C
    #matrix.Sigma<-matrix(diag(Sigmahat),N-1,T,byrow=TRUE) #corrected
    log.f.Y[i1]<-sum(dnorm(Y.star,0,sqrt(Sigmahat[1,1]),log=TRUE))-0.5*((T0-r1)*r1+N*r1+K+1)*log(N)
    #please consider the below for the revision #
    if(max.BF<log.f.Y[i1]){
      max.BF<-log.f.Y[i1]
      ini.beta.MCMC<-Betahat
      ini.Lambda.MCMC<-Lambdahat
      ini.F.MCMC<-Fhat
      ini.Sigma.MCMC<-Sigmahat}
    #  print(paste("Current r is ",r1))
    
  }
  
  #log.f.Y
  #which.max(log.f.Y)
  #  plot(R0,log.f.Y)
  
  
  #######  Bayesian GSCM  ##########
  #############################################
  #Simulating the initial values of parameters
  tausq<-1000 #hyperparameter for delta (var)
  mu.delta<-0 #hyperparameter for delta (mean)
  r1= which.max(log.f.Y)+1 #corrected
  Bayes.r[M]<-r1
  
  
  txx<-t(X[,,2])%*%X[,,2]
  txy<-t(X[,,2])%*%Y[2,]
  for (i in 3:N){
    txx<-txx+t(X[,,i])%*%X[,,i]
    txy<-txy+t(X[,,i])%*%Y[i,]
  }
  
  hat.beta<-solve(txx)%*%txy
  
  
  Fstar <- matrix(0, nrow = T, ncol = r1)
  for (t in 1:T) {
    f_t <- matrix(0, nrow = 1, ncol = min(t, r1))
    for (i in 1:min(t, r1)) {
      f_t[i] <- rnorm(1,0,1) 
    }
    Fstar[t, 1:min(t, r1)] <- f_t
  }
  
  
  hat.Lambda<-matrix(0,N,r1)
  for (i in 1:N){
    for (j in 1:r1){
      hat.Lambda[i,j]<-rnorm(1,0,1)
    }
  }
  
  hat.Sigma<-diag(1,T)
  Psi<-diag(1,r1)
  
  v0.I.k<-diag(1/v0,K) #identity matrix of size K
  MC<-10000
  MC.beta<-matrix(0,MC,K)
  MC.Lambda<-matrix(0,MC,r1)
  MC.ATT<-matrix(0,MC,T)
  X.delta<-cbind(rep(1,T1),(T0+1):T) #M
 # T.time<-(T0+1):T
 # X.delta<-cbind(rep(1,T1),T.time-mean(T.time),(T.time-mean(T.time))^2) #M
  
  
  inv.tXX.delta<-solve(crossprod(X.delta)) #inverse of M^t M
  #beta.delta<-mvrnorm(1,mean.delta,var.delta)
  mu.delta<-rep(0,T1) #initial value of M*%*%alpha
  for (k in 1:MC){
    
    ####beta|y,Omega##########
    Omega<-xi0*Fstar%*%t(Fstar)+hat.Sigma
    inv.Omega<-solve(Omega)
    txx<-t(X[,,2])%*%(inv.Omega%*%X[,,2])
    txy<-t(X[,,2])%*%(inv.Omega%*%Y[2,])
    
    for (i in 3:N){
      txx<-txx+t(X[,,i])%*%(inv.Omega%*%X[,,i])
      txy<-txy+t(X[,,i])%*%(inv.Omega%*%Y[i,])
    }
    
    sd.beta<-solve(txx+v0.I.k) 
    mean.beta<-as.numeric(sd.beta%*%txy) 
    hat.beta<-mvrnorm(1,mean.beta,sd.beta)
    MC.beta[k,]<-hat.beta
    
    
    #####Lambda_i|F,beta,Sigma,y,Psi############
    inv.Psi<-solve(Psi)
    inv.hat.Sigma<-solve(hat.Sigma)
    sd.Lambda<-solve(t(Fstar)%*%inv.hat.Sigma%*%Fstar+inv.Psi)
    
    for (j in 2:N){
      mean.Lambda<-as.numeric(sd.Lambda%*%t(Fstar)%*%(inv.hat.Sigma%*%(Y[j,]-X[,,j]%*%hat.beta)))
      hat.Lambda[j,]<-mvrnorm(1,mean.Lambda,sd.Lambda) 
    }
    
    
    #####f_t*|Lambda,beta,Sigma,y,Psi#############
    for (i in 1:T) {
      k_t <- min(i, r1)
      I_k <- diag(k_t)
      M_t<- as.matrix(hat.Lambda[2:N, 1:k_t])
      
      sd.ft<-solve(I_k+(hat.Sigma[i,i])^(-1)*t(M_t)%*%M_t)
      mean.ft<-as.numeric(sd.ft%*%((hat.Sigma[i,i])^(-1)*t(M_t)%*%(Y[-1,i]-t(X[i,,-1])%*%hat.beta)))
      Fstar[i,1:k_t] <-mvrnorm(1,mean.ft,sd.ft)
    }
    
    ###Psi|Lambda,beta,F,y,Sigma######
    for (L in 1:r1){ 
      psi1<-(a0+N-1)/2
      psi2<-b0/2+0.5*sum(hat.Lambda[-1,L]^2)
      Psi[L,L]<-rinvgamma(1,psi1,psi2) 
    }
    
    ###sigmat|Lambda,beta,F,y,Psi###### 
    sigma1<-(a0+(N-1)*T)/2
    sigma2<-b0/2
    for (j in 1:T){
      sigma2<-sigma2+0.5*sum((Y[-1,j]-t(X[j,,-1])%*%hat.beta-hat.Lambda[-1,]%*%as.matrix(Fstar[j,]))^2)
    }
    diag(hat.Sigma)<-rinvgamma(1,sigma1,sigma2)   
    # print(paste("current iteration is",k))
    
    ########Step 2: Compute f_tl for t=1,...,T, l=1,...,r##########
    hat.F <- matrix(0, nrow = T, ncol = r1)
    for (t in 1:T){
      for (l in 1:r1){
        hat.F[t,l]<-(-1*(Fstar[l,l]<0)+1*(Fstar[l,l]>=0))*Fstar[t,l]*sqrt(Psi[l,l])
      }
    }
    
    ######Step3: i=1,and t= 1,...,T0
    #####hat\Lambda_i|F,beta,Sigma,y############
    inv.Sigma.T0<-solve(hat.Sigma[1:T0,1:T0])
    sd.Lambda1<-solve(t(hat.F[1:T0,])%*%inv.Sigma.T0%*%hat.F[1:T0,]+diag(r1)/xi0)
    mean.Lambda1<-as.numeric(sd.Lambda1%*%t(hat.F[1:T0,])%*%inv.Sigma.T0%*%(Y[1,1:T0]-X[1:T0,,1]%*%hat.beta))
    hat.Lambda[1,]<-mvrnorm(1,mean.Lambda1,sd.Lambda1)
    
    #######Step 4: i in T,and t=T0+1,...,T ####
    
    for(t in (T0+1):T){
      mean.delta<-(Y[1,t]-t(X[t,,1])%*%hat.beta-t(hat.Lambda[1,])%*%hat.F[t,]+hat.Sigma[1,1]/tausq*mu.delta[t-T0])/(1+hat.Sigma[1,1]/tausq)
      var.delta<-hat.Sigma[1,1]/(1+hat.Sigma[1,1]/tausq)
      MC.ATT[k,t]<-rnorm(1,mean.delta,sqrt(var.delta)) #delta
    } 
    # print(paste(k,"-th iteration is running"))
    a.tau<- (a0+T1)/2
    b.tau<- (b0+sum((MC.ATT[k,(T0+1):T]-mu.delta)^2))/2
    tausq<- rinvgamma(1,a.tau,b.tau)
    #mean of delta
    mean.delta<-inv.tXX.delta%*%t(X.delta)%*%MC.ATT[k,(T0+1):T]
    var.delta<- tausq*inv.tXX.delta
    beta.delta<-mvrnorm(1,mean.delta,var.delta) #alpha
    mu.delta<-X.delta%*%beta.delta #M x alpha
  } 
  
  B<-0.2*k
  Post.ATT<-MC.ATT[-(1:B),]
  
  Bayes.est.ATT<-apply(Post.ATT,2,mean)
  
  Bayes.95CI.ATT<-apply(Post.ATT,2,quantile,c(0.025,0.975))
  
########################################
  #####################################
 #####     Bayes Quadratic    ######## 
  
  MC.ATT1<-matrix(0,MC,T)
   T.time<-(T0+1):T
   X.delta1<-cbind(rep(1,T1),T.time-mean(T.time),(T.time-mean(T.time))^2) #M
  
  
  inv.tXX.delta1<-solve(crossprod(X.delta1)) #inverse of M^t M
  #beta.delta<-mvrnorm(1,mean.delta,var.delta)
  mu.delta1<-rep(0,T1) #initial value of M*%*%alpha
  for (k in 1:MC){
    
    ####beta|y,Omega##########
    Omega<-xi0*Fstar%*%t(Fstar)+hat.Sigma
    inv.Omega<-solve(Omega)
    txx<-t(X[,,2])%*%(inv.Omega%*%X[,,2])
    txy<-t(X[,,2])%*%(inv.Omega%*%Y[2,])
    
    for (i in 3:N){
      txx<-txx+t(X[,,i])%*%(inv.Omega%*%X[,,i])
      txy<-txy+t(X[,,i])%*%(inv.Omega%*%Y[i,])
    }
    
    sd.beta<-solve(txx+v0.I.k) 
    mean.beta<-as.numeric(sd.beta%*%txy) 
    hat.beta<-mvrnorm(1,mean.beta,sd.beta)
    MC.beta[k,]<-hat.beta
    
    
    #####Lambda_i|F,beta,Sigma,y,Psi############
    inv.Psi<-solve(Psi)
    inv.hat.Sigma<-solve(hat.Sigma)
    sd.Lambda<-solve(t(Fstar)%*%inv.hat.Sigma%*%Fstar+inv.Psi)
    
    for (j in 2:N){
      mean.Lambda<-as.numeric(sd.Lambda%*%t(Fstar)%*%(inv.hat.Sigma%*%(Y[j,]-X[,,j]%*%hat.beta)))
      hat.Lambda[j,]<-mvrnorm(1,mean.Lambda,sd.Lambda) 
    }
    
    
    #####f_t*|Lambda,beta,Sigma,y,Psi#############
    for (i in 1:T) {
      k_t <- min(i, r1)
      I_k <- diag(k_t)
      M_t<- as.matrix(hat.Lambda[2:N, 1:k_t])
      
      sd.ft<-solve(I_k+(hat.Sigma[i,i])^(-1)*t(M_t)%*%M_t)
      mean.ft<-as.numeric(sd.ft%*%((hat.Sigma[i,i])^(-1)*t(M_t)%*%(Y[-1,i]-t(X[i,,-1])%*%hat.beta)))
      Fstar[i,1:k_t] <-mvrnorm(1,mean.ft,sd.ft)
    }
    
    ###Psi|Lambda,beta,F,y,Sigma######
    for (L in 1:r1){ 
      psi1<-(a0+N-1)/2
      psi2<-b0/2+0.5*sum(hat.Lambda[-1,L]^2)
      Psi[L,L]<-rinvgamma(1,psi1,psi2) 
    }
    
    ###sigmat|Lambda,beta,F,y,Psi###### 
    sigma1<-(a0+(N-1)*T)/2
    sigma2<-b0/2
    for (j in 1:T){
      sigma2<-sigma2+0.5*sum((Y[-1,j]-t(X[j,,-1])%*%hat.beta-hat.Lambda[-1,]%*%as.matrix(Fstar[j,]))^2)
    }
    diag(hat.Sigma)<-rinvgamma(1,sigma1,sigma2)   
    # print(paste("current iteration is",k))
    
    ########Step 2: Compute f_tl for t=1,...,T, l=1,...,r  ##########
    hat.F <- matrix(0, nrow = T, ncol = r1)
    for (t in 1:T){
      for (l in 1:r1){
        hat.F[t,l]<-(-1*(Fstar[l,l]<0)+1*(Fstar[l,l]>=0))*Fstar[t,l]*sqrt(Psi[l,l])
      }
    }
    
    ######Step3: i=1,and t= 1,...,T0
    #####hat\Lambda_i|F,beta,Sigma,y############
    inv.Sigma.T0<-solve(hat.Sigma[1:T0,1:T0])
    sd.Lambda1<-solve(t(hat.F[1:T0,])%*%inv.Sigma.T0%*%hat.F[1:T0,]+diag(r1)/xi0)
    mean.Lambda1<-as.numeric(sd.Lambda1%*%t(hat.F[1:T0,])%*%inv.Sigma.T0%*%(Y[1,1:T0]-X[1:T0,,1]%*%hat.beta))
    hat.Lambda[1,]<-mvrnorm(1,mean.Lambda1,sd.Lambda1)
    
    #######Step 4: i in T,and t=T0+1,...,T ####
    
    for(t in (T0+1):T){
      mean.delta<-(Y[1,t]-t(X[t,,1])%*%hat.beta-t(hat.Lambda[1,])%*%hat.F[t,]+hat.Sigma[1,1]/tausq*mu.delta[t-T0])/(1+hat.Sigma[1,1]/tausq)
      var.delta<-hat.Sigma[1,1]/(1+hat.Sigma[1,1]/tausq)
      MC.ATT1[k,t]<-rnorm(1,mean.delta,sqrt(var.delta)) #delta
    } 
    # print(paste(k,"-th iteration is running"))
    a.tau<- (a0+T1)/2
    b.tau<- (b0+sum((MC.ATT1[k,(T0+1):T]-mu.delta1)^2))/2
    tausq<- rinvgamma(1,a.tau,b.tau)
    #mean of delta
    mean.delta<-inv.tXX.delta1%*%t(X.delta1)%*%MC.ATT1[k,(T0+1):T]
    var.delta<- tausq*inv.tXX.delta1
    beta.delta<-mvrnorm(1,mean.delta,var.delta) #alpha
    mu.delta1<-X.delta1%*%beta.delta #M x alpha
  } 
  
  B<-0.2*k
  Post.ATT1<-MC.ATT1[-(1:B),]
  
  Bayes.est.ATT1<-apply(Post.ATT1,2,mean)
  
  Bayes.95CI.ATT1<-apply(Post.ATT1,2,quantile,c(0.025,0.975))
  
  ########################
  ########  GSCM  ###########
  
  ###### Define the treatment assignment variable####
  D<-c(rep(0,T0),rep(1,T1),rep(0,(N*T-T)))
  ###### Define the time variable######
  time <- rep(1:T, N)
  #### Define the region ######
  unit<-rep(1:N, each = T, times = 1)
  
  # Combine all variables into a data.frame
  
  data <- data.frame(Y = as.vector(t(Y)), unit=unit, time = time, D=D,
                     Z1 = as.vector(t(Z[,1,])), Z2 = as.vector(t(Z[,2,])),
                     Z3 = as.vector(t(Z[,3,])), Z4 = as.vector(t(Z[,4,])))
  
  
  out <- gsynth(Y ~ D + Z1 + Z2 + Z3 + Z4, data = data,
                index = c("unit","time"), force = "two-way",
                CV = TRUE, r = c(2,9), se = TRUE,
                inference = "parametric", nboots = 1000,
                parallel = FALSE)
  
  GSCM.r[M]<-out$r.cv  
  
  
  out_r2 <- gsynth(Y ~ D + Z1 + Z2 + Z3 + Z4, data = data,
                   index = c("unit","time"), force = "two-way",
                   CV = FALSE, r = 2, se = TRUE,
                   inference = "parametric", nboots = 1000,
                   parallel = FALSE)
  
  
  
  #####
  t=seq(T0+1,T,1)
  TrueTrtEffect=0.5+sqrt(t/2)
  
  Bayes<-as.numeric((TrueTrtEffect-Bayes.est.ATT[T11:T])^2)
  Bayes.Mse[M]<-mean(Bayes)
  
  BayesQuadratic<-as.numeric((TrueTrtEffect-Bayes.est.ATT1[T11:T])^2)
  Bayes1.Mse[M]<-mean(BayesQuadratic)
  
  GSCM<-as.numeric((TrueTrtEffect-out$est.att[T11:T])^2)
  GSCM.Mse[M]<-mean(GSCM)
  
  GSCM_r2<-as.numeric((TrueTrtEffect-out_r2$est.att[T11:T])^2)
  GSCM_r2.Mse[M]<-mean(GSCM_r2)
  
  
  ####### average length ###########
  Bayes.linear<-as.numeric(Bayes.95CI.ATT[2,T11:T]-Bayes.95CI.ATT[1,T11:T])
  Bayes.averagelength[M]<-mean(Bayes.linear)
  
  Bayes.qua<-as.numeric(Bayes.95CI.ATT1[2,T11:T]-Bayes.95CI.ATT1[1,T11:T])
  Bayes1.averagelength[M]<-mean(Bayes.qua)
  
  GSCM1<-as.numeric(out$est.att[T11:T,4]-out$est.att[T11:T,3])
  GSCM.averagelength[M]<-mean(GSCM1)
  
  GSCMr2<-as.numeric(out_r2$est.att[T11:T,4]-out_r2$est.att[T11:T,3])
  GSCMr2.averagelength[M]<-mean(GSCMr2)
  
  print(paste("current iteration is",M))
  #  print(c("Bayes","GSCM","GSCM_r2"))
  #  print(c(Bayes.Mse[M], GSCM.Mse[M],GSCM_r2.Mse[M]))
}


plot(TrueTrtEffect,col=1,type="l",lwd=2,ylim=c(0,10))
lines(Bayes.est.ATT[T11:T],col=2,lwd=2)
lines( Bayes.95CI.ATT[1,T11:T],col=2,lwd=2,lty=2)
lines( Bayes.95CI.ATT[2,T11:T],col=2,lwd=2,lty=2)
lines(out_r2$est.att[T11:T],col=3,lwd=2)
lines(out$est.att[T11:T],col=4,lwd=2)



mean(Bayes.Mse)
mean(Bayes1.Mse)
mean(GSCM.Mse)
mean(GSCM_r2.Mse)

mean(Bayes.averagelength)
mean(Bayes1.averagelength)
mean(GSCM.averagelength)
mean(GSCMr2.averagelength)

#standard error
sd(Bayes.Mse)/sqrt(sim.n) 
sd(Bayes1.Mse)/sqrt(sim.n) 
sd(GSCM.Mse)/sqrt(sim.n) 
sd(GSCM_r2.Mse)/sqrt(sim.n) 


Bayes.r
GSCM.r

barplot(table(Bayes.r),ylim = c(0,1000),ylab = "The counts r of Bayes GSCM",xlab="r")
barplot(table(GSCM.r),ylim = c(0,1000),ylab = "The counts r of GSCM",xlab="r")

boxplot(GSCM.Mse,GSCM_r2.Mse,Bayes.Mse,Bayes1.Mse,names=c("GSCM","GSCM.r2","Bayeslinear","Bayesquadratic"))  















