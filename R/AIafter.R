
AIafter <- function(X, y, Xnew, ynew=NULL, ai.test=TRUE, safe=TRUE, sig.level=0.1, cfi.nskip=35, after.nskip=5) {
		
	n = nrow(X);
    p = ncol(X);
    n1 = nrow(X)+nrow(Xnew)
    
    if (cfi.nskip>n-10) {
    	stop("cfi.nskip is too large to build stronger forecast.")
    }
    if (after.nskip>(n-cfi.nskip-5)) {
    	stop("after.nskip is too large to use AFTER effectively.")
    }
    
    fcst = NULL;
    
    # build strong forecasts
    str.fcst = CFI1(X,y,Xnew,ynew,cfi.nskip); 
    str.fcst0 <- str.fcst[1:n,]
    str.fcst1 <- str.fcst[(n+1):n1,]
    num.str.fcst<-ncol(str.fcst);

	if (ai.test==TRUE) {
		y.test <- y[-c(1:cfi.nskip)]
		X.test<-X[-c(1:cfi.nskip),]
		str.fcst.test<-str.fcst[-c(1:cfi.nskip, (n+1):n1),]
		
		s_err = matrix(rep(y.test, num.str.fcst), byrow=F, ncol=num.str.fcst) - str.fcst.test;
      	o_err = matrix(rep(y.test, p), byrow=F, ncol=p) - X.test;
      	s_best = s_err[, which.min(colMeans(s_err^2))];
      	o_best = o_err[, which.min(colMeans(o_err^2))];

      	###if s_best and o_best are too similar, skip the test
      	if (sum(abs(s_best - o_best)) < 1e-3) {
        	dm_stat<-0
        	p_dmtest<-0.5 
      	} else { 
        	dmtest<-dm.test(s_best, o_best, 
          alternative= "less");
        	dm_stat<-as.numeric(dmtest$statistic);
        	p_dmtest<-as.numeric(dmtest$p.value);
      	}

		if (p_dmtest >= sig.level)
      	GOAL="Adaptation"
      	else GOAL="Improvement"
      	
      	if(GOAL=="Adaptation"){
			fcst = AFTER1(X, y, Xnew, ynew, after.nskip)
		} else if(GOAL=="Improvement"){ 
			fcst = AFTER1(str.fcst0[-c(1:cfi.nskip),], y[-c(1:cfi.nskip)],str.fcst1,ynew,after.nskip)
			fcst$fcst = rbind(matrix(NA,nrow=cfi.nskip,ncol=ncol(fcst$fcst)), fcst$fcst);
		}
		
		if (safe==TRUE) { # compare estimation error
			both.X0<-cbind(X, str.fcst0)[-c(1:cfi.nskip),]
			both.X1 <- cbind(Xnew, str.fcst1)
      		fcst1 = AFTER1(both.X0, y[-c(1:cfi.nskip)], both.X1, ynew, after.nskip);
      		fcst1$fcst = rbind(matrix(NA,nrow=cfi.nskip,ncol=ncol(fcst1$fcst)), fcst1$fcst);
			
			# compute errors
			y.test<-y[-c(1:(cfi.nskip+after.nskip))]
			fcst1.test <- fcst1$fcst[-c(1:(cfi.nskip+after.nskip), (n+1):n1),]
			num.fcst1<-ncol(fcst1$fcst)			
			s_err1 = matrix(rep(y.test, num.fcst1), byrow=F, ncol=num.fcst1) - fcst1.test;
			fcst.test <- fcst$fcst[-c(1:(cfi.nskip+after.nskip), (n+1):n1),]
			num.fcst<-ncol(fcst$fcst)
			s_err0 = matrix(rep(y.test, num.fcst), byrow=F, ncol=num.fcst) - fcst.test;
			
			numerr1 <- min(colMeans(s_err1^2)) # MSE from combined
			numerr0 <- min(colMeans(s_err0^2)) # MSE from testing
			if (numerr1<numerr0)
				return(list(fcst=fcst1$fcst, pval_dm=p_dmtest))
			else
				return(list(fcst=fcst$fcst, pval_dm=p_dmtest))
			
		}
		return(list(fcst=fcst$fcst, pval_dm=p_dmtest))
				
	} else {
		both.X0<-cbind(X, str.fcst0)[-c(1:cfi.nskip),]
		both.X1 <- cbind(Xnew, str.fcst1)
      	fcst = AFTER1(both.X0, y[-c(1:cfi.nskip)], both.X1, ynew, after.nskip);
      	fcst$fcst = rbind(matrix(NA,nrow=cfi.nskip,ncol=ncol(fcst$fcst)), fcst$fcst);
		return(list(fcst=fcst$fcst, pval_dm=NA))	
	}
   
}



