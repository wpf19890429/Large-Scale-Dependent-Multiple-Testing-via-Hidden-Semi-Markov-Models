rdata.shifted.nbinomial.hsmm<-function(NUM, par_dwell, pro, f0, pc, f1)
{

## USAGE
 # rdata.shifted.nbinomial.hsmm(NUM, n, Pi, prob, f0, ...)


## ARGUEMENTS
 # NUM: number of hypotheses
 # par_dwell: the parameters of the shifted negative binomial distribution
 # par_dwell=list(n1, Pi)
 # pro: the parameter of the geometric distribution
 # Omega=(0, 1; 1, 0) is known
 # pc=(pc[1], ..., pc[L]): proportion of mixture components
 # f0: parameter of the null distribution
 # f1: parameter of the non-null mixture distribution


## DETAILS
 # rdata.shifted.binomial.hsmm generates the underlying states and observations for HSMM 
 # where the dwell time distribution for state 0 is a shifted negative binomial distribution and
 # the dwell time distribution for state 1 is a geometric distribution.  


## VALUES
 # o: continuous observed data
 # s: binary unobserved states


## Initialization

	theta<-NULL
	x<-rep(0, NUM)
      n1<-par_dwell[[1]]
      Pi<-par_dwell[[2]]

## generating the underlying states
 # Without loss of generality, assume that theta[1]=0.

	while(length(theta)<NUM)
	{
           	theta<-c(theta, rep(0, rnbinom(1, n1, Pi)+1))
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
                k<-sample(1:nc, 1, prob=pc)
                x[i]<-rnorm(1, mean=f1[k, 1], sd=f1[k, 2])
           }
	}
	data<-list(s=theta[1:NUM], o=x)
	return(data)
}







