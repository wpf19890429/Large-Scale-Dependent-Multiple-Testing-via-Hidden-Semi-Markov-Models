bwfw.hsmm2hmm<-function(x, dm, f0, pc, f1)
{

###################################################################

## USAGE
 # bwfw.hsmm2hmm(x, dm, f0, pc, f1)

## ARGUMENTS
 # x=(x[1], ..., x[m]): the observed data
 # dm=list(dm[[1]], dm[[2]])
 # dm[[1]]: the p.m.f. of the dwell time distribution of state 0 
 # dm[[2]]: the p.m.f. of the geometric distribution dm[[2]][r]=prob*(1-prob)^(r-1)
 # f0=(mu_0, sd_0): the parameters for null distribution
 # pc=(pc[1], ..., pc[L]): proportion of mixture components
 # f1=(mu[1], sd[1]\\...\\ mu[L], sd[L]): the parameters for the non-null distribution

## DETAILS
 # bwfw.hsmm2hmm calculates values for backward, forward variables, LIS variables and etc by using
 # the HMM with expanded state space to approximate the HSMM.  

## VALUES
 # alpha: rescaled backward variables
 # beta: rescaled forward variables
 # LIS_HSMM: the LIS variables
 # c0: scaling variables

################################################################## 


## Initialize

	NUM<-length(x)
      m1<-length(dm[[1]])
      m2<-length(dm[[2]])
      L<-length(pc)

## transition matrix

	Omega<-matrix(c(0, 1, 1, 0), 2, 2)
      Gamma<-hsmm2hmm(Omega, dm)

## densities

      f0x<-dnorm(x, f0[1], f0[2])      
      f1x<-rep(0, NUM)
      for(p in 1:L)
      {
      	f1x<-f1x+pc[p]*dnorm(x, f1[p, 1], f1[p, 2])
      }

# the backward-forward procedure
# a. the forward variables
# --rescaled 

	alpha<-array(0, c(NUM, m1+m2))

# scaling variable c_0

      c0<-rep(0, NUM)

# stationary distribution

	delta<-powerAx(Gamma, rep(1/nrow(Gamma), nrow(Gamma)), 100) 

# the forward variable for initial state

      alpha[1, ]<-delta*c(rep(f0x[1], m1), rep(f1x[1], m2))
      c0[1]<-1/sum(alpha[1, ])
      alpha[1, ]<-c0[1]*alpha[1, ]

# the forward variable for other states
     
      for(k in 1:(NUM-1))
      {
            alpha[k+1, ]<-t(Gamma)%*%alpha[k, ]*c(rep(f0x[k+1], m1), rep(f1x[k+1], m2))
            c0[k+1]<-1/sum(alpha[k+1, ])
            alpha[k+1, ]<-c0[k+1]*alpha[k+1, ]
      }

# b. the backward variables
# --rescaled

      b0<-rep(0, NUM)
      beta<-array(0, c(NUM, m1+m2))
      beta[NUM, ]<-1/(m1+m2)

      for (k in (NUM-1):1)
      { 
            beta[k, ]<-Gamma%*%(c(rep(f0x[k+1], m1), rep(f1x[k+1], m2))*beta[k+1, ])  
            b0[k]<-1/sum(beta[k, ])
            beta[k, ]<-b0[k]*beta[k, ]
      } 

# c. LIS_HSMM variables
# --original
# --the same formulae hold for the rescaled alpha and beta

	LIS_HSMM<-rep(0, NUM)
	for (k in 1:NUM)
	{ 
            q1<-sum(alpha[k, 1:m1]*beta[k, 1:m1])
            q2<-sum(alpha[k, (m1+1):(m1+m2)]*beta[k, (m1+1):(m1+m2)])
            LIS_HSMM[k]<-q1/(q1+q2)          
	}

# f. return the results of the bwfw proc.

	bwfw.var<-list(fw=alpha, bw=beta, LIS=LIS_HSMM, c0=c0)
      return(bwfw.var)
}


      