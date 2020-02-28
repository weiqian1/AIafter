

CLR1 = function(X,y,Xnew,ynew=NULL,n0=5)
{
   L2CR = function(X,y)
   {
      ynew = y - X[,1]; Xnew = X[,-1]- X[,1];
      n = nrow(Xnew); p = ncol(Xnew); 
      Rinv = t(Xnew) %*% Xnew + 0.001*diag(1,p);
      C = cbind(rep(-1,p), diag(p));
      b = c(-1, rep(0.001,p));
      d = t(Xnew)%*%ynew;
      bta = solve.QP(Dmat = Rinv, dvec = d, Amat = C, bvec = b)$solution;
      bta = bta*(bta>0.001); cof = c(1-sum(bta),bta);
      return(cof);
   }
   if (!is.null(ynew)) {
   X <- rbind(X,Xnew)
   y <- c(y,ynew)
   n = nrow(X); 
   p = ncol(X);
   fcsts = rep(NA,n);

   for(i in n0:(n-1))
   {
      Xtmp = X[c(1:i),]; 
	  ytmp = y[1:i]; 
	  Xnew = X[(i+1),]; 
	  mtmp = L2CR(Xtmp,ytmp);
      fcsts[(i+1)] = sum(mtmp*Xnew);
   }
   return(fcsts); 
   } else {
   	X1 <- rbind(X,Xnew)
   	n1 <- nrow(X1)
   	n = nrow(X); p = ncol(X);
   	fcsts = rep(NA,n1);
   	
   	for(i in n0:(n-1))
   {
      Xtmp = X[c(1:i),]; 
	  ytmp = y[1:i]; 
	  Xnew = X[(i+1),]; 
	  mtmp = L2CR(Xtmp,ytmp);
      fcsts[(i+1)] = sum(mtmp*Xnew);
   }
   Xtmp <- X[1:n,]
   ytmp <- y[1:n]
   mtmp = L2CR(Xtmp,ytmp);
   for(i in n:(n1-1))
   {
	  Xnew = X1[(i+1),]; 
      fcsts[(i+1)] = sum(mtmp*Xnew);
   }
   return(fcsts); 	
   }
   
     
}

