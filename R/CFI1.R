CFI1 = function(X, y, Xnew, ynew=NULL, n0=5) {


  ### simple average
  SA = function(X,Xnew){
  	return(apply(rbind(X,Xnew),1,mean))
  	}

  ### median
  MD = function(X,Xnew){
  	return(apply(rbind(X,Xnew),1,median))
  	}

  ### simple trimmed mean 
  TM = function(X, Xnew, lb = 0.05, ub = 0.05){
	TrimMean.row = function(x){
		x_tmp = sort(x,decreasing=F);
        n = length(x);
        lb_num = floor(n*lb);
        ub_num = floor(n*ub);
		Larges = tail(x_tmp,ub_num);
		Smalls = head(x_tmp,lb_num);
		result = 1/(length(x)-lb_num-ub_num)*(sum(x) - sum(Larges) - sum(Smalls));
		return(result);
	}
	return(apply(rbind(X,Xnew),1,TrimMean.row));
  }

  
  Naive = function(X,y,Xnew,ynew=NULL,n0=5){
	BG_fcsts = BG1(X,y,Xnew,ynew,n0,c(1,0.9));
	BG_names = colnames(BG_fcsts);

	fcsts = cbind(SA(X,Xnew),MD(X,Xnew),TM(X,Xnew),BG_fcsts);
	colnames(fcsts) = c("SA","MD","TM",BG_names);	
	return(fcsts);	
  }

fwd_lr = function(X,y,size)
{
	fwd_next_lr = function(X,y,input)
	{
		L2R = function(X,y)
		{
		   X = as.matrix(X);
		   n = nrow(X); p = ncol(X); Xnew = cbind(rep(1,n),X); 
		   cof = as.vector(lm.fit(Xnew,y)$coefficients); 
		   cof = ifelse(1-is.na(cof)==1,cof,0); return(cof);
		}
	
		p = ncol(X);
		n = nrow(X);
		coefs = input[-1];
		s_inds = which(coefs!= 0);
		f_inds = which(coefs == 0);
		q = length(f_inds);
	
		coefs = list();
		length(coefs) = q;
		sse = numeric(q);
	
		for(i in 1:q)
		{
			Xnew = X[,c(s_inds,f_inds[i])];
			coefs_tmp = L2R(Xnew,y);
			coefs[[i]] = coefs_tmp;
			sse[i] = mean((y - apply(cbind(rep(1,n),Xnew),1,function(x) sum(coefs_tmp*x)))^2);	
		}
	
		f_ind = which.min(sse);
		coef_f = coefs[[f_ind]];
		output = rep(0,p+1);
		output[1] = coef_f[1];
		output[c(1+s_inds,1+f_inds[f_ind])] = coef_f[-1];
		output;
	
	}
	### end of the internal function
	
	n = nrow(X);
	p = ncol(X);
	
	size = min(p,max(size,0));
	
	mod_list = NULL;
		
	input = c(mean(y), rep(0,p));
	
	if(size == 0){
		output = input;
		break;
	}else{
		output = NULL;	
	}
		
	mod_list = rbind(mod_list,input);
	i = 0;	
	while(i<=(size-1)){
		input = fwd_next_lr(X,y,input);
		mod_list = rbind(mod_list,input);
		i = i + 1;		
	}
	
	output = mod_list;
	rownames(output) = paste("s",c(0:size),sep="");
	
	return(output);
}

L2CR = function(X,y)
{
	   ynew = y - X[,1]; Xnew = as.matrix(X[,-1]- X[,1]);
	   n = nrow(Xnew); p = ncol(Xnew); 
	   Rinv = t(Xnew) %*% Xnew + 0.001*diag(1,p);
	   C = cbind(rep(-1,p), diag(p));
	   b = c(-1, rep(0.001,p));
	   d = t(Xnew)%*%ynew;
	   bta = solve.QP(Dmat = Rinv, dvec = d, Amat = C, bvec = b)$solution;
	   bta = bta*(bta>0.001); cof = c(1-sum(bta),bta);
	   return(cof);
}


fwd_clr = function(X,y,size)
{
	### define an internal function
	fwd_next_clr = function(X,y,input)
	{
		p = ncol(X);
		coefs = input;
		s_inds = which(coefs != 0);
		f_inds = which(coefs == 0);
		q = length(f_inds);

		coefs = list();
		length(coefs) = q;
		sse = numeric(q);

		for(i in 1:q)
		{
			Xnew = X[,c(s_inds,f_inds[i])];
			coefs_tmp = L2CR(Xnew,y);
			coefs[[i]] = coefs_tmp;
			sse[i] = mean((y - apply(Xnew,1,function(x) sum(coefs_tmp*x)))^2);	
		}

		f_ind = which.min(sse);
		output = rep(0,p);
		output[c(s_inds,f_inds[f_ind])] = coefs[[f_ind]];
		output;

	}
	### end of the internal function	
	n = nrow(X);
	p = ncol(X);
	size = min(max(0,size),p);
	
	mod_list = NULL;
		
	err = X - matrix(rep(y,p),ncol=p,byrow=F);	
	sse = colMeans(err^2);
	input = rep(0,p);
	input[which.min(sse)] = 1;
	i = 1;
	
	if(size == 1){
		output = input;
		break;
	}else{
		output = NULL;	
	}
		
	mod_list = rbind(mod_list,input);
		
	while(i<=(size-1)){
		input = fwd_next_clr(X,y,input);
		mod_list = rbind(mod_list,input);
		i = i + 1;		
	}
	
	output = mod_list;
	rownames(output) = paste("s",c(1:size),sep="");
	
	return(output);
}

    X = as.matrix(X);
    n = nrow(X);
    p = ncol(X);
    
    colnames(X) = paste("X",c(1:p),sep="");
    regs.scope = paste0("y~",paste0("X", 1:p, collapse="+"));
    fcst = NULL;
    eps = 1e-10;
    
    vmax<-min(p, n0-1);

    ol_regs = fwd_lr(X[c(1:n0),],y[1:n0],vmax);
    cl_regs = fwd_clr(X[c(1:n0),],y[1:n0],vmax);

    ol_which = abs(ol_regs) > eps
    ol_formulas<-paste0("y~",as.vector(sapply(apply(ol_which[-1,-1], 1, which), function(x){paste0("X", x, collapse="+")})));  

    cl_which = cl_regs > eps
    cl_which_uniq = unique(cl_which)[-1, , drop=FALSE]       

    fcst = matrix(NA, nrow=n0, ncol=4+length(ol_formulas)+nrow(cl_which_uniq));
	
	if (!is.null(ynew)) {
	Xt <- X 
	yt <- y
	X <- rbind(X,Xnew)
	y <- c(y,ynew)
	n <- nrow(X)
	p <- ncol(X)
	
    for(i in n0:(n-1)) {
          Xtmp = X[c(1:i),];
          ytmp = y[c(1:i)];
          Xi<-matrix(as.numeric(X[i+1,]), nrow=1);
		
          ncvpath<-ncvreg(X=Xtmp, y=ytmp, penalty="lasso", lambda.min=ifelse(nrow(Xtmp) > ncol(Xtmp), 0.001, 0.01), esp = 0.01, max.iter=5000, warn = F);
          current.fcst<-predict(ncvpath, X=Xi, which=c(which.min(AIC(ncvpath)), which.min(BIC(ncvpath))));
		
          step.AIC<-stepAIC(lm(as.formula("y~1"), data=data.frame(Xtmp,y=ytmp)), scope=as.formula(regs.scope), direction="forward", trace=0, steps=min(p, i-2), k=2);
          current.fcst<-c(current.fcst, predict(step.AIC, newdata=data.frame(Xi)));
        
          step.BIC<-stepAIC(lm(as.formula("y~1"), data=data.frame(Xtmp,y=ytmp)), scope=as.formula(regs.scope), direction="forward", trace=0, steps=min(p, i-2), k=log(n))
          current.fcst<-c(current.fcst, predict(step.BIC, newdata=data.frame(Xi)));	

          
          for (j in 1:length(ol_formulas)) {
              current.fcst = c(current.fcst, predict(lm(as.formula(ol_formulas[j]), data=data.frame(Xtmp,y=ytmp)), newdata=data.frame(Xi)))
          }
          if (nrow(cl_which_uniq)!=0) { #Wei debugged 
          for (k in 1:nrow(cl_which_uniq)) {
              mtmp = L2CR(Xtmp[,cl_which_uniq[k,]],y=ytmp)
              ### CAR 11-28-17 NEED SOME CODE HERE TO SET THE COLUMNS NOT IN CL_WHICH_UNIQ TO 0
              mtmp2 = rep(0,length(Xi))
              mtmp2[cl_which_uniq[k,]] = mtmp;
              ### END CAR 11-28-17	
              current.fcst = c(current.fcst, sum(mtmp2*Xi));   ###CAR 11-28-17 replaced mtmp with mtmp2
          }
          }
          fcst = rbind(fcst,current.fcst);    
    }
    X <- Xt
    y <- yt
    
    } else {
    X1 <- rbind(X,Xnew)
	n1 <- nrow(X1)
	n = nrow(X); p = ncol(X);
	
	for(i in n0:(n-1)) {
          Xtmp = X[c(1:i),];
          ytmp = y[c(1:i)];
          Xi<-matrix(as.numeric(X[i+1,]), nrow=1);
		
          ncvpath<-ncvreg(X=Xtmp, y=ytmp, penalty="lasso", lambda.min=ifelse(nrow(Xtmp) > ncol(Xtmp), 0.001, 0.01), esp = 0.01, max.iter=5000, warn = F);
          current.fcst<-predict(ncvpath, X=Xi, which=c(which.min(AIC(ncvpath)), which.min(BIC(ncvpath))));
		
          step.AIC<-stepAIC(lm(as.formula("y~1"), data=data.frame(Xtmp,y=ytmp)), scope=as.formula(regs.scope), direction="forward", trace=0, steps=min(p, i-2), k=2);
          current.fcst<-c(current.fcst, predict(step.AIC, newdata=data.frame(Xi)));
        
          step.BIC<-stepAIC(lm(as.formula("y~1"), data=data.frame(Xtmp,y=ytmp)), scope=as.formula(regs.scope), direction="forward", trace=0, steps=min(p, i-2), k=log(n1))
          current.fcst<-c(current.fcst, predict(step.BIC, newdata=data.frame(Xi)));	

          
          for (j in 1:length(ol_formulas)) {
              current.fcst = c(current.fcst, predict(lm(as.formula(ol_formulas[j]), data=data.frame(Xtmp,y=ytmp)), newdata=data.frame(Xi)))
          }
          if (nrow(cl_which_uniq)!=0) { #Wei debugged 
          for (k in 1:nrow(cl_which_uniq)) {
              mtmp = L2CR(Xtmp[,cl_which_uniq[k,]],y=ytmp)
              
              mtmp2 = rep(0,length(Xi))
              mtmp2[cl_which_uniq[k,]] = mtmp;
              current.fcst = c(current.fcst, sum(mtmp2*Xi));   
          }
          }
          fcst = rbind(fcst,current.fcst);    
    }
	Xtmp = X[c(1:n),];
    ytmp = y[c(1:n)];
    for(i in n:(n1-1)) {
          
          Xi<-matrix(as.numeric(X1[i+1,]), nrow=1);
		
	  	  if (i==n){
          ncvpath<-ncvreg(X=Xtmp, y=ytmp, penalty="lasso", lambda.min=ifelse(nrow(Xtmp) > ncol(Xtmp), 0.001, 0.01), esp = 0.01, max.iter=5000, warn = F);
          }
          current.fcst<-predict(ncvpath, X=Xi, which=c(which.min(AIC(ncvpath)), which.min(BIC(ncvpath))));
		
	  	  if (i==n){
          step.AIC<-stepAIC(lm(as.formula("y~1"), data=data.frame(Xtmp,y=ytmp)), scope=as.formula(regs.scope), direction="forward", trace=0, steps=min(p, i-2), k=2);
          }
          current.fcst<-c(current.fcst, predict(step.AIC, newdata=data.frame(Xi)));
          if (i==n){
          step.BIC<-stepAIC(lm(as.formula("y~1"), data=data.frame(Xtmp,y=ytmp)), scope=as.formula(regs.scope), direction="forward", trace=0, steps=min(p, i-2), k=log(n))
          }
          current.fcst<-c(current.fcst, predict(step.BIC, newdata=data.frame(Xi)));	

          
          for (j in 1:length(ol_formulas)) {
              current.fcst = c(current.fcst, predict(lm(as.formula(ol_formulas[j]), data=data.frame(Xtmp,y=ytmp)), newdata=data.frame(Xi)))
          }
          if (nrow(cl_which_uniq)!=0) { #Wei debugged 
          for (k in 1:nrow(cl_which_uniq)) {
              mtmp = L2CR(Xtmp[,cl_which_uniq[k,]],y=ytmp)
              mtmp2 = rep(0,length(Xi))
              mtmp2[cl_which_uniq[k,]] = mtmp;
              current.fcst = c(current.fcst, sum(mtmp2*Xi));   
          }
          }
          fcst = rbind(fcst,current.fcst);    
    }	  	
    }
	
    naive.fcst = Naive(X,y,Xnew,ynew,n0);
    lr.fcst = LinReg1(X,y,Xnew,ynew,n0);
    clr.fcst = CLR1(X,y,Xnew,ynew,n0);
    fcst = cbind(fcst,naive.fcst,lr.fcst,clr.fcst);
	
    colnames(fcst) = paste("sf",c(1:ncol(fcst)),sep="");
    rownames(fcst) = NULL;   

    return(fcst);
}

