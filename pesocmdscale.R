pesocmdscale<-function(d, k, w)
{
  n=nrow(d)
  weight.centre <- function(x, w) {
    w.c <- apply(x, 2, weighted.mean, w = w)
    x <- sweep(x, 2, w.c, "-")
    x
  }
  d=d^2
  if (missing(w))
    w <- rep(1, n)
  m <- weight.centre(d, w)
  Gw <- -1/2*t(weight.centre(t(m), w))
  Fw <- Gw * sqrt(w) %o% sqrt(w)
  e <- eigen(Fw)
  autovalores=Re(e$values)  #Tomamos la parte real como son practicamente cero las componentes complejas
  autovectores=Re(e$vectors)
  points <- sweep(autovectores, 1, sqrt(w), "/")
  points <- sweep(points, 2, sqrt(autovalores), "*") 
  #Tomas la parte real
  # Dw_pow_menos_unmedio=diag(diag(Dw)^(-1/2))
  Dw=diag(as.numeric(sweep(ones(n,1),1,sqrt(w),"/")))
  # Yw=Dw%*%autovectores%*%diag(autovalores^(1/2))
  res<-list(autovectores=autovectores,autovalores=autovalores,points=points,e=e)
  return(res)
}