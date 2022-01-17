powerAx<-function(Gamma,delta,Iter=10){
	res<-t(Gamma)%*%delta
	for(i in 1:Iter){
		res<-t(Gamma)%*%res
	}
	return(res)
}