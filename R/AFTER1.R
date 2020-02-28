
AFTER1 = function(X, y, Xnew, ynew=NULL, n0=5)
{    
	
      if (!is.null(ynew)) {
      	X <- rbind(X,Xnew)
   		y <- c(y,ynew)    	
      ##########   Part I:  initiate variables called later    ##########
       n = nrow(X);  
       p = ncol(X);

       err = X - matrix(rep(y,p),ncol=p,byrow=F); 
       W2 = W1 = err*0; 
       
       vs = ds = matrix(NA,n,p)
   
      ##########   Part II:  L1-/L2-AFTER procedure    ##########
       for(j in 1:p)
       {
          e = err[,j]; logl = matrix(0,n,2);
          
          ## get the likelihood at each point
          for(i in n0:(n-1))
          {
               ef = e[i]; 
               err_train = err[c(1:i),j];
               vs[i,j] = var(err_train);
               ds[i,j] = mean(abs(err_train));
 
               ### l2-log-likelihood
               v_tmp = vs[i,j]; l2_tmp = ef^2/(2*v_tmp) + 0.5*log(v_tmp) + 0.5*log(2*pi);

               ### l1-log-likelihood
               d_tmp = ds[i,j]; l1_tmp = abs(ef)/d_tmp + log(2*d_tmp);

               logl[(i+1),] = c(l2_tmp,l1_tmp); tmp = apply(logl,2,sum); 
               W2[(i+1),j] = tmp[1]; 
               W1[(i+1),j] = tmp[2];
          }           
      }
      
      ##########   Part III:  get the combining weights of candidate forecasts   ##########
       ### weights of L1AFTER
          W1 = W1 - apply(W1,1,min); 
          index1 = 1*(W1-20<0); 
          L1 = exp(-W1); 
          W1 = L1/apply(L1,1,sum); 
          W1 = W1*index1;
        ### weights of L2AFTER
          W2 = W2 - apply(W2,1,min); 
          index2 = 1*(W2-20<0); 
          L2 = exp(-W2); 
          W2 = L2/apply(L2,1,sum); 
          W2 = W2*index2;

     ##########   Part IV:  organize the output   ##########
         f_1 = rowSums(W1*X);
         f_2 = rowSums(W2*X);  
         output = cbind(f_1, f_2); 
      
         return(list(fcst=output, L1Wt=W1, L2Wt=W2));
         
         
       } else {
       	########## if ynew is NULL
       	
       	X1 <- rbind(X,Xnew)
   		n1 <- nrow(X1)
   		n = nrow(X); p = ncol(X);
   		err = X - matrix(rep(y,p),ncol=p,byrow=F); 
       	W2 = W1 = matrix(0,n1,p);       
       	vs = ds = matrix(NA,n,p)
   		
       for(j in 1:p)
       {
          e = err[,j]; logl = matrix(0,n+1,2);
          
          for(i in n0:n)
          {
               ef = e[i]; 
               err_train = err[c(1:i),j];
               vs[i,j] = var(err_train);
               ds[i,j] = mean(abs(err_train));
 
               v_tmp = vs[i,j]; l2_tmp = ef^2/(2*v_tmp) + 0.5*log(v_tmp) + 0.5*log(2*pi);

               d_tmp = ds[i,j]; l1_tmp = abs(ef)/d_tmp + log(2*d_tmp);

               logl[(i+1),] = c(l2_tmp,l1_tmp); tmp = apply(logl,2,sum); 
               W2[(i+1),j] = tmp[1]; 
               W1[(i+1),j] = tmp[2];
          }           
      }
      if (n1>=n+2) {
      for (i in ((n+1):(n1-1))) {
      	W2[i+1,] <- W2[n+1,]
      	W1[i+1,] <- W1[n+1,]
      }
      }
      
          W1 = W1 - apply(W1,1,min); 
          index1 = 1*(W1-20<0); 
          L1 = exp(-W1); 
          W1 = L1/apply(L1,1,sum); 
          W1 = W1*index1;
        ### weights of L2AFTER
          W2 = W2 - apply(W2,1,min); 
          index2 = 1*(W2-20<0); 
          L2 = exp(-W2); 
          W2 = L2/apply(L2,1,sum); 
          W2 = W2*index2;
      
         f_1 = rowSums(W1*X1);
         f_2 = rowSums(W2*X1);  
         output = cbind(f_1, f_2); 
      
         return(list(fcst=output, L1Wt=W1, L2Wt=W2));
   		
       	
    }          
         
}











