
BG1 = function(X,y,Xnew,ynew=NULL,n0=5,rhos=c(1,0.9)){
	### private function
	cumulative_m = function(error,rho){
		n = length(error);
		sequence = numeric(n);
		for(i in 1:n) sequence[i] = rho^(n-i);
		return(sum(error^2*sequence));
	}
	
	if (!is.null(ynew)) {
		X <- rbind(X,Xnew)
		y <- c(y,ynew)
	n = nrow(X); p = ncol(X);
	numR = length(rhos);
	err = X - matrix(rep(y,p),byrow=F,ncol=p);
	combined = matrix(NA,nrow=n,ncol=numR);
	for(i in n0:(n-1)){
		err_train = err[c(1:i),];
		for(j in 1:numR){
			rho = rhos[j];
			m_i = apply(err_train,2,function(x,rho) cumulative_m(x,rho), rho=rho);
			v_i = 1/m_i;	
			combined[(i+1),j] = sum(X[c(i+1),]*v_i/sum(v_i));					
		}
	}
	colnames(combined) = paste("BG",rhos,"ALL",sep="_");
	return(combined);
	} else {
		X1 <- rbind(X,Xnew)
		n1 <- nrow(X1)
		n = nrow(X); p = ncol(X);
		numR = length(rhos);
		err = X - matrix(rep(y,p),byrow=F,ncol=p);
		combined = matrix(NA,nrow=n1,ncol=numR);
		
		for(i in n0:n){
		err_train = err[c(1:i),];
		for(j in 1:numR){
			rho = rhos[j];
			m_i = apply(err_train,2,function(x,rho) cumulative_m(x,rho), rho=rho);
			v_i = 1/m_i;	
			combined[(i+1),j] = sum(X1[c(i+1),]*v_i/sum(v_i));					
		}
	}
	if (n1>=n+2) {
	for (i in (n+1):(n1-1)) {
		for (j in 1:numR) {
			rho = rhos[j];
			combined[(i+1),j] = sum(X1[c(i+1),]*v_i/sum(v_i));
		}
	}
	}
	colnames(combined) = paste("BG",rhos,"ALL",sep="_");
	return(combined);
		
	}
	
}
