rdata.unstructed.geom.hsmm<-function(NUM, par_dwell, pro, f0, pc, f1)
{

## USAGE
 # rdata.unstructed.geom.hsmm(NUM, dm1, prob, f0, ...)


## ARGUEMENTS
 # NUM: number of hypotheses
 # par_dwell: par_dwell=list(dm1)
 # dm1: the p.m.f. of the dwell time distribution of state 0 
 # pro: the parameter of the geometric distribution
 # Omega=(0, 1; 1, 0) is known
 # pc=(pc[1], ..., pc[L]): proportion of mixture components
 # f0: parameter of the null distribution
 # f1: parameter of the non-null mixture distribution


## DETAILS
 # rdata.unstructed.geom.hsmm generates the underlying states and observations for HSMM 
 # where the dwell time distribution for state 0 is a unstructed start and geometric tail distribution and
 # the dwell time distribution for state 1 is a geometric distribution.  


## VALUES
 # o: continuous observed data
 # s: binary unobserved states


## Initialization

	theta<-NULL
	x<-rep(0, NUM)
      dm1<-par_dwell[[1]]

## generating the underlying states
 # Without loss of generality, assume that theta[1]=0.
 # dm1[1:(n1-1)]: unstructed start
 # dm1[n1]:geometric tail, dm1[n1] = (1-sum(dm1[1:(n1-1)]))*Pi_1

      n1<-length(dm1)
      sum_start<-1-sum(dm1[1:(n1-1)]) 
      Pi_1<-dm1[n1]/sum_start
	while(length(theta)<NUM)
	{
            k<-sample(1:n1, 1, prob=c(dm1[1:(n1-1)], sum_start))
            if(k<n1)
            {
                  theta<-c(theta, rep(0, k))
            }
            else
            {
                  theta<-c(theta, rep(0, rgeom(1, Pi_1)+n1)) 
            }
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