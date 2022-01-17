################################ setting the progress bar
PB<-function(Methods,Rep){
	library(progress)
	pb <- progress_bar$new(
         format = paste(Methods,":completed [:bar] :percent, Execute time::elapsed",sep=""),
         total = Rep, clear = FALSE, width= 60)
	return(pb)
}
######################################

em.hsmm2hmm<-function(x, m1, m2=1, L, maxiter=100, progress=FALSE)
{

###################################################################

## USAGE
 # em.hsmm2hmm(x, m1, ...)

## ARGUMENTS
 # x=(x[1], ..., x[m]): the observed data
 # m1: the number of the state aggregates for state 0 
 # m2: the number of the state aggregates for state 1
 # L: the number of mixture components for the non-null distribution
 # maxiter: the maximum number of iterations
 # progress: the variable indicates whether to show progress

## DETAILS
 # em.hsmm2hmm calculates the MLE for a HMM which can closely approximate the HSMM.
 # the distribution of state 0 is assumed to be known as N(0, 1)
 # the distribution of state 1 is assumed to be a normal mixture with L components

## VALUES
 # dm=list(dm[[1]], dm[[2]]): the estimated p.m.f. for the dwell time distribution of state 0 and state 1
 # pc=(pc[1], ..., pc[L]): the estimated proportion of mixture components
 # f1=(mu[1], sd[1]\\...\\ mu[L], sd[L]): the estimated parameters for the non-null distribution
 # BIC: BIC value

################################################################## 

	NUM<-length(x)
# precision tolerance level

	#ptol<-1e-3
	ptol<-0

	niter<-0
	f0<-c(0, 1)

## Parameters initialization

	Omega<-matrix(c(0,1,1,0),2,2) 
      dm.new<-list(rep(0.95/m1, m1), rep(0.5/m2, m2))
      Gamma.new<-hsmm2hmm(Omega, dm.new)
	Gamma.new[which(Gamma.new!=0)]<-0.5
	
	pc.new<-rep(1/L, L)

      f1.new<-matrix(c(1:L, rep(1, L)), nrow=L, ncol=2, byrow=FALSE)



	diff<-10
	Loglikelihood.new<--10000


### The E-M Algorithm
	pb<-PB(Methods="EM_expanded_HMM",Rep=maxiter)

	while(diff>ptol && niter<maxiter)
      {
	  pb$tick()
	  niter<-niter+1
        Gamma.old<-Gamma.new
        pc.old<-pc.new
        f1.old<-f1.new
        Loglikelihood.old<-Loglikelihood.new

## updating the weights and probabilities of hidden states

        bwfw.res<-bwfw.hmm.multistates(x, Gamma.old, m1, m2=1, f0, pc.old, f1.old)

# the backward-forward variable

	  alpha<-bwfw.res$fw
        beta<-bwfw.res$bw

## Densities
	  f0x<-dnorm(x, f0[1], f0[2])      
	  f1x<-rep(0, NUM)
	  for(p in 1:L)
	  {
      	  f1x<-f1x+pc.old[p]*dnorm(x, f1.old[p, 1], f1.old[p, 2])
	  }

# gamma
      
        gamma<-matrix(0, NUM, m1+m2)
	  for (k in 1:NUM)
	  {
               gamma[k, ]<-alpha[k, ]*beta[k, ]/sum(alpha[k, ]*beta[k, ]) 
        }   
      
# eta=Pr(theta_i=p-1, theta_{i+1}=q-1 | x_1^m)
        eta<-array(0, c(NUM-1, m1+m2, m1+m2))
        b1<-rep(0, NUM-1)
        for(i in 1:(NUM-1))
        {
              for(j in 1:(m1+m2))
              {
                      eta[i, j, ]<-alpha[i, j]*Gamma.old[j, ]*(c(rep(f0x[i+1], m1), rep(f1x[i+1], m2))*beta[i+1, ])
              }
              b1[i]<-1/sum(eta[i, ,])
              eta[i, , ]<-b1[i]*eta[i, , ]
         } 


# b. transition matrix Gamma 
        for(i in 1:(m1+m2))
        {   
            for(j in 1:(m1+m2))
            {
                 q1<-sum(eta[, i, j])
                 q2<-sum(eta[, i, ])
		     if(q2==0){
				Gamma.new[i, j]<-0
		     }else{
				Gamma.new[i, j]<-q1/q2
		     }
                 
            }
         }

# c. intermediate variable
        ptheta<-rep(0, NUM)
        for(i in 1:NUM)
        {
            ptheta[i]<-sum(alpha[i, (m1+1):(m1+m2)]*beta[i, (m1+1):(m1+m2)])/sum(alpha[i, ]*beta[i, ])
        }


# c. weight variables
        omega<-matrix(0, nrow=NUM, ncol=L)
        for (c in 1:L)
        { 
          for (i in 1:NUM)
          {

                 f1x.c<-dnorm(x, mean = f1.old[c, 1], sd = f1.old[c, 2])
                 omega[i, c]<-ptheta[i]*pc.old[c]*f1x.c[i]/f1x[i]
          }
        }

# d. updating the non-null distribution 

 
        for (c in 1:L)
        {
# (i). probability weights

            q1<-sum(omega[, c])
            q2<-sum(ptheta)
            pc.new[c]<-q1/q2

# (ii). means

            q3<-sum(omega[, c]*x)
            f1.new[c, 1]<-q3/q1

# (iii). sds

            q4<-sum(omega[, c]*(x-f1.new[c, 1])*(x-f1.new[c, 1]))
            f1.new[c, 2]<-sqrt(q4/q1)

		if(progress){
			cat("The current iteration is: iter=", niter,"/", maxiter,"\n",sep="")
		}

        } 

	c0<-bwfw.res$c0
      Loglikelihood.new<--sum(log(c0))
      df1<-abs(Loglikelihood.old-Loglikelihood.new)
      diff<-df1/abs(Loglikelihood.old)

	}


	bwfw.res<-bwfw.hmm.multistates(x, Gamma.new, m1, m2=1, f0, pc.new, f1.new)

# g. return the results of the E-M algorithm
	em.var<-list(Gamma=Gamma.new, pc=pc.new, f1=f1.new, LIS=bwfw.res$LIS, ni=niter)
	return(em.var)
}
