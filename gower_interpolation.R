
p=5;
n=nrow(puntos_restantes)
partitions=c(0,n/p,n/p*2,n/p*3,n/p*4,n/p*5)
gower_interpolation <-function (fit1,partitions,axes){
  
  autovalores=fit1$autovalores[which(fit1$autovalores>1e-8)]
  matriz_lambda=diag(autovalores)
  m=nrow(muestra)
  matriz_y=as.matrix(fit1$points[,1:ncol(matriz_lambda)])
  matriz_G=matriz_y%*%t(matriz_y)
  g=matrix(diag(matriz_G),nrow=1,ncol=m)
  inv_lambda=solve(matriz_lambda)
  
  for (i in 1:p){
    m=nrow(puntos_restantes[(partitions[i]+1):partitions[i+1],])
    new_dist=gower.dist(puntos_restantes[(partitions[i]+1):partitions[i+1],],muestra) 
    yn=1/2*((ones(m,1)%*%g)-new_dist)%*%matriz_y%*%inv_lambda
    axes=rbind(axes,yn[,1:3])
  }
  axes
}

yn=gower_interpolation(fit1,partitions,axes)