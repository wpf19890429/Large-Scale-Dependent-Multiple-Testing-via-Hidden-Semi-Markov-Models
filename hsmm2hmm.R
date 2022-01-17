hsmm2hmm<-function(omega, dm, eps=1e-10)
{
     mv<-sapply(dm, length)
     m<-length(mv)
     G<-matrix(0, 0, sum(mv))
     for(i in 1:m)
     {
           mi<-mv[[i]]
           F<-cumsum(c(0, dm[[i]][-mi]))
           ci<-ifelse(abs(1-F)>0, dm[[i]]/(1-F), 1)
           cim<-ifelse(1-ci>0, 1-ci, 0)
           Gi<-matrix(0, mi, 0)
           for(j in 1:m)
           {
                 if(i==j)
                 {
                      if(mi==1)
                            Gi<-cbind(Gi, c(rep(0, mv[[j]]-1), cim))
                      else
                            Gi<-cbind(Gi, rbind(cbind(rep(0, mi-1), diag(cim[-mi], mi-1, mi-1)), 
                                               c(rep(0, mi-1), cim[[mi]])))
                  }
                  else
                  {
                      if(mi==1)
                            Gi<-cbind(Gi, matrix(c(omega[[i, j]]*ci, rep(0, mv[[j]]-1)), 1))
                      else
                            Gi<-cbind(Gi, cbind(omega[[i, j]]*ci, matrix(0, mv[[i]], mv[[j]]-1)))
                  }
            }
            G<-rbind(G, Gi)
       }
       return(G)
}    