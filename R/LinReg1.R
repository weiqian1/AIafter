

LinReg1 = function(X,y,Xnew,ynew=NULL,n0=5)
{
   L2R = function(X,y)
   {
      n = nrow(X); p = ncol(X); Xnew = cbind(rep(1,n),X); 
      cof = as.vector(lm.fit(Xnew,y)$coefficients); 
      cof = ifelse(1-is.na(cof)==1,cof,0); return(cof);
   }
	
   if (!is.null(ynew)) {
   X <- rbind(X,Xnew)
   y <- c(y,ynew)
   n = nrow(X); p = ncol(X);
   Frg = rep(NA,n)
   for(i in n0:(n-1))
   {
      Xtmp = X[c(1:i),]; ytmp = y[1:i]; Xnew = X[(i+1),]; 
      mtmp1 = L2R(Xtmp,ytmp);
      Frg[(i+1)] = mtmp1[1] + sum(mtmp1[-1]*Xnew);
   }
   return(Frg);
   } else {
   	X1 <- rbind(X,Xnew)
   	n1 <- nrow(X1)
   	n = nrow(X); p = ncol(X);
   Frg = rep(NA,n1)
   for(i in n0:(n-1))
   {
      Xtmp = X[c(1:i),]; ytmp = y[1:i]; Xnew = X[(i+1),]; 
      mtmp1 = L2R(Xtmp,ytmp);
      Frg[(i+1)] = mtmp1[1] + sum(mtmp1[-1]*Xnew);
   }
   Xtmp <- X[1:n,]
   ytmp <- y[1:n]
   mtmp1 = L2R(Xtmp,ytmp);
   for (i in n:(n1-1)) 
   {
      Xnew = X1[(i+1),];     
      Frg[(i+1)] = mtmp1[1] + sum(mtmp1[-1]*Xnew);
   }
   return(Frg);
   	
   	
   }
   
   
   
      
}














