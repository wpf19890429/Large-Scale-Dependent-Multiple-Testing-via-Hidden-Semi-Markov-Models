rdata.shifted.Poisson.hsmm<-function(NUM, par_dwell, pro, f0, pc, f1)
{

## USAGE
 # rdata.shifted.Poisson.hsmm(NUM, lambda, prob, f0, ...)


## ARGUEMENTS
 # NUM: number of hypotheses
 # par_dwell: the parameter of the shifted Poisson distribution
 # pro: the parameter of the geometric distribution
 # Omega=(0, 1; 1, 0) is known
 # pc=(pc[1], ..., pc[L]): proportion of mixture components
 # f0: parameter of the null distribution
 # f1: parameter of the non-null mixture distribution


## DETAILS
 # rdata.shifted.Poisson.hsmm generates the underlying states and observations for HSMM 
 # where the dwell time distribution for state 0 is a shifted Poisson distribution and
 # the dwell time distribution for state 1 is a geometric distribution.  


## VALUES
 # o: continuous observed data
 # s: binary unobserved states


## Initialization

	theta<-NULL
	x<-rep(0, NUM)
      lambda<-par_dwell[[1]]

## generating the underlying states
 # Without loss of generality, assume that theta[1]=0.

	while(length(theta)<NUM)
	{
           	theta<-c(theta, rep(0, rpois(1, lambda)+1))
           	theta<-c(theta, rep(1, rgeom(1, pro)+1))
      }

## generating the observations

	nc<-length(pc)
	for (i in 1:NUM)
	{
           if (theta[i]==0)
           {
                x[i]<-rnorm(1, mean=f0[1], sd=f0[2])
           }
           else
           { 
                c<-sample(1:nc, 1, prob=pc)
                x[i]<-rnorm(1, mean=f1[c, 1], sd=f1[c, 2])
           }
	}
	data<-list(s=theta[1:NUM], o=x)
	return(data)
}
