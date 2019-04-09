

delta.estimate= function(xone,xzero, deltaone, deltazero, t, std= FALSE, conf.int = FALSE, weight.perturb = NULL) {
	delta.f = delta.estimate.RMST(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, weight = weight.perturb, delta.only=F)
	delta = delta.f$delta
	rmst.1 = delta.f$s1
	rmst.0 = delta.f$s0
	if(std | conf.int)	{
		n1 = length(xone)
		n0 = length(xzero)
		if(is.null(weight.perturb)) {weight.perturb = matrix(rexp(500*(n1+n0), rate=1), ncol = 500)} 
		delta.p.vec = apply(weight.perturb, 2, delta.estimate.RMST, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, delta.only = T)
		delta.p.vec = unlist(delta.p.vec)
		sd.delta = sd(delta.p.vec)
		mad.delta = mad(delta.p.vec)
	}
	if(conf.int){
		conf.quantile.delta = c(new.q(delta.p.vec, 0.025), new.q(delta.p.vec, 0.975))
	}
	if(!std & !conf.int) {return(list("delta" = delta, "rmst.1" = rmst.1, "rmst.0" = rmst.0))}
	if(std & !conf.int) {return(list("delta" = delta, "rmst.1" = rmst.1, "rmst.0" = rmst.0, "delta.sd" = sd.delta, "delta.mad" = mad.delta))}
	if(conf.int) {return(list("delta" = delta, "rmst.1" = rmst.1, "rmst.0" = rmst.0, "delta.sd" = sd.delta, "delta.mad" = mad.delta, "conf.int" = conf.quantile.delta))}
}

delta.estimate.RMST= function(xone,xzero, deltaone, deltazero, t, weight = NULL, delta.only=F) {
	if(is.null(weight)) {weight = rep(1,length(xone)+length(xzero))}
	censor1.t = censor.weight(xone, deltaone, t, weight = weight[1:length(xone)])
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight[(1+length(xone)):(length(xone)+length(xzero))])
	censor1.times = censor.weight(xone, deltaone, xone, weight = weight[1:length(xone)])
	censor0.times = censor.weight(xzero, deltazero, xzero, weight = weight[(1+length(xone)):(length(xone)+length(xzero))])
		M.1 = rep(0,length(xone))
		M.0 = rep(0,length(xzero))
		M.1[xone>t] = t/censor1.t
		M.1[xone<=t] = xone[xone<=t]*deltaone[xone<=t]/censor1.times[xone<=t] 	
		M.0[xzero>t] = t/censor0.t
		M.0[xzero<=t] = xzero[xzero<=t]*deltazero[xzero<=t]/censor0.times[xzero<=t] 
		delta = 1/(sum(weight[1:length(xone)])) * sum(M.1*weight[1:length(xone)]) - 1/(sum(weight[(1+length(xone)):(length(xone)+length(xzero))])) * sum(M.0*weight[(1+length(xone)):(length(xone)+length(xzero))])
		s1 = 1/(sum(weight[1:length(xone)])) * sum(M.1*weight[1:length(xone)])
		s0 = 1/(sum(weight[(1+length(xone)):(length(xone)+length(xzero))])) * sum(M.0*weight[(1+length(xone)):(length(xone)+length(xzero))])
	if(delta.only) { return(list("delta" = delta)) }
	if(!delta.only) { return(list("delta" = delta, "s1" = s1, "s0" = s0)) } 

}

R.q.event = function(xone,xzero, deltaone, deltazero, sone, szero, t, landmark, number = 40, transform = FALSE, extrapolate = TRUE, std = FALSE, conf.int = FALSE, weight.perturb = NULL, type = "np") {
	if(!(type %in% c("np","semi"))) {warning("Warning: Invalid type, default `np' for nonparametric estimator being used", call. = FALSE); type = "np"}
	warn.te = FALSE
	warn.support = FALSE
	n1 = length(xone)
	n0 = length(xzero)
	if(is.null(weight.perturb)){
		weight.perturb = matrix(rexp(500*(n1+n0), rate=1), ncol = 500)
	}
	delta = delta.estimate.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t)
	delta.estimate = delta$delta
	delta.p.vec = apply(weight.perturb, 2, delta.estimate.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, delta.only = T)
	delta.p.vec = unlist(delta.p.vec)
	sd.delta  = sd(delta.p.vec)
	mad.delta = mad(delta.p.vec)
	conf.quantile.delta = c(new.q(delta.p.vec, 0.025), new.q(delta.p.vec, 0.975))
	if(0>conf.quantile.delta[1] & 0< conf.quantile.delta[2]) {warning("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the proportion of treatment effect explained in this setting", call. = FALSE)
		warn.te = TRUE}
	if(delta.estimate < 0) {warning("Warning: it looks like you need to switch the treatment groups", call. = FALSE)}
	range.1 = range(pmin(sone,landmark))
	range.0 = range(pmin(szero,landmark))
	range.ind = (range.1[1] > range.0[1]) | (range.1[2] < range.0[2])
	if(range.ind & !extrapolate & !transform) {
		warning("Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation", call. = FALSE)
		warn.support = TRUE
	}
	
	if(type == "np"){
		delta.q.estimate = delta.q.event.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, extrapolate = extrapolate, transform = transform, number = number)$delta.q
		R.q = 1-delta.q.estimate/delta.estimate	
	}
	if(type == "semi"){
		delta.q.estimate = delta.q.event.semi.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, number = number)$delta.q
		R.q = 1-delta.q.estimate/delta.estimate	
	}
	if(std | conf.int){
		if(type == "np"){
			delta.q.p.vec.temp = t(apply(weight.perturb, 2, delta.q.event.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, extrapolate = extrapolate, transform = transform, number = number, deltaslist=FALSE, warn.extrapolate = FALSE))
			delta.q.p.vec = delta.q.p.vec.temp[,4]
			R.p = 1-delta.q.p.vec/delta.p.vec
		}
		if(type == "semi"){
			delta.q.p.vec.temp = t(apply(weight.perturb, 2, delta.q.event.semi.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, deltaslist=FALSE, number = number))
			delta.q.p.vec = delta.q.p.vec.temp[,4]
			R.p = 1-delta.q.p.vec/delta.p.vec
		}
	}
	
	if(conf.int) {	

		conf.l.quantile.delta = quantile(delta.p.vec, 0.025)
		conf.u.quantile.delta = quantile(delta.p.vec, 0.975)
		
		conf.l.quantile.delta.q = quantile(delta.q.p.vec, 0.025)
		conf.u.quantile.delta.q = quantile(delta.q.p.vec, 0.975)
		
		conf.l.quantile.R.q = quantile(R.p, 0.025)
		conf.u.quantile.R.q = quantile(R.p, 0.975)
}
	if(!std & !conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q))}
	if(std & !conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.q.sd" = sd(delta.q.p.vec),  "delta.q.mad" = mad(delta.q.p.vec), "R.q.sd" = sd(R.p), "R.q.mad" = mad(R.p)))}
	if(conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.q.sd" = sd(delta.q.p.vec),  "delta.q.mad" = mad(delta.q.p.vec), "R.q.sd" = sd(R.p), "R.q.mad" = mad(R.p), "conf.int.delta" = as.vector(c(conf.l.quantile.delta, conf.u.quantile.delta)), "conf.int.delta.q" = as.vector(c(conf.l.quantile.delta.q, conf.u.quantile.delta.q)),
"conf.int.R.q" = as.vector(c(conf.l.quantile.R.q, conf.u.quantile.R.q))))
}	
}




delta.q.event.RMST = function(xone,xzero, deltaone, deltazero, sone, szero, t, weight = NULL, landmark=landmark, deltaslist = TRUE, transform = FALSE, extrapolate = TRUE, number, warn.extrapolate = TRUE) {
	if(is.null(weight)) {weight = rep(1,length(xone)+length(xzero))}
	weight.group1 = weight[1:length(xone)]
	weight.group0 = weight[(1+length(xone)):(length(xone)+length(xzero))]
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
	censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
	censor1.t = censor.weight(xone, deltaone, t, weight = weight.group1)
	censor1.landmark = censor.weight(xone, deltaone, landmark, weight = weight.group1)
	
	#first term
	t.vec = c(landmark, landmark+c(1:(number-1))*(t-landmark)/number, t)
	#number surrogate events in control group before landmark
	number.s.o = length(szero[xzero>landmark & szero < landmark])
	phi.a.tst0 = matrix(nrow = number.s.o, ncol = number+1)
	warn.temp = FALSE
	for(i in 1:(number+1)) {
	temp.phi = pred.smooth.surv(xone.f=xone[xone>landmark & sone < landmark], deltaone.f = deltaone[xone>landmark & sone <  landmark], sone.f=log(sone[xone>landmark & sone <  landmark]), szero.one = log(szero[xzero>landmark & szero < landmark]), myt=t.vec[i], weight = weight.group1[xone>landmark & sone <  landmark], transform = transform, extrapolate = extrapolate)
	phi.a.tst0[,i] =  temp.phi$Phat.ss
	if(temp.phi$warn.flag == 1) {warn.temp = TRUE}
	}
	if(warn.temp == TRUE & warn.extrapolate == TRUE) {if(warn.extrapolate) {warning("Note: Values were extrapolated.", call. = FALSE)}}
	phi.int = vector(length = number.s.o)
	for(j in 1:number.s.o) {
		phi.int[j] = landmark + (t-landmark)/number*(phi.a.tst0[j,1]/2 + sum(phi.a.tst0[j,2:number]) + phi.a.tst0[j,number+1]/2)
	}
	first.term = sum(weight.group0[xzero>landmark & szero <landmark]*phi.int)/(sum(weight.group0)*censor0.landmark) 
	
	
	#second term
	censor1.times = censor.weight(xone, deltaone, xone, weight = weight.group1)
	M.1 = rep(0,length(xone))
	M.1[xone>t] = t/censor1.t
	M.1[xone<=t] = xone[xone<=t]*deltaone[xone<=t]/censor1.times[xone<=t] 
	denom = sum(1*(sone > landmark & xone > landmark)*weight.group1)	
	psi.a.tt0 = sum(censor1.landmark*M.1*weight.group1*(sone > landmark & xone > landmark)/(denom))

	second.term = sum(weight.group0*(szero > landmark & xzero > landmark)*psi.a.tt0)/(sum(weight.group0)*censor0.landmark) 
	
	#third.term
	censor0.times = censor.weight(xzero, deltazero, xzero, weight = weight.group0)
	M.0 = rep(0,length(xzero))
	M.0[xzero>t] = t/censor0.t
	M.0[xzero<=t] = xzero[xzero<=t]*deltazero[xzero<=t]/censor0.times[xzero<=t] 
	denom = sum(1*( xzero > landmark)*weight.group0)	
	nu.term = sum(censor0.landmark*M.0*weight.group0*(xzero > landmark)/(denom))
	third.term = (sum(weight.group0*(xzero > landmark))/(sum(weight.group0)*censor0.landmark))*nu.term

	
	delta.q = first.term+second.term - third.term
	
	if(deltaslist) {return(list("delta.q" = delta.q))}
	if(!deltaslist) {return(c(first.term, second.term, third.term, delta.q))}


}

delta.q.event.semi.RMST = function(xone,xzero, deltaone, deltazero, sone, szero, t, weight = NULL, landmark=landmark, deltaslist = TRUE, number) {
	if(is.null(weight)) {weight = rep(1,length(xone)+length(xzero))}
	weight.group1 = weight[1:length(xone)]
	weight.group0 = weight[(1+length(xone)):(length(xone)+length(xzero))]
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
	censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
	censor1.t = censor.weight(xone, deltaone, t, weight = weight.group1)
	censor1.landmark = censor.weight(xone, deltaone, landmark, weight = weight.group1)
	
	#first term
	s.predictor = sone[xone>landmark & sone <  landmark]
	x.adjust = xone[xone>landmark & sone < landmark] - landmark
	sum.model = coxph(Surv(x.adjust, deltaone[xone>landmark & sone <  landmark])~s.predictor, weights = weight.group1[xone>landmark & sone <  landmark])
	s.new = as.data.frame(szero[xzero>landmark & szero < landmark])
	names(s.new) = "s.predictor"
	predict.score = predict(sum.model, type = "risk", newdata = s.new)
	baseline.hazard = basehaz(sum.model, centered = TRUE); 
	#yes, we really do want centered = TRUE
	
	t.vec = c(landmark, landmark+c(1:(number-1))*(t-landmark)/number, t)
	#number surrogate events in control group before landmark
	number.s.o = length(szero[xzero>landmark & szero < landmark])
	phi.a.tst0 = matrix(nrow = number.s.o, ncol = number+1)
	for(i in 1:(number+1)) {
		gap.time = t.vec[i]-landmark
		baseline.t <- approx(baseline.hazard$time,baseline.hazard$hazard,gap.time, rule = 2)$y
		phi.a.tst0[,i] =  exp(-baseline.t*predict.score)
	}
	phi.int = vector(length = number.s.o)
	for(j in 1:number.s.o) {
		phi.int[j] = landmark + (t-landmark)/number*(phi.a.tst0[j,1]/2 + sum(phi.a.tst0[j,2:number]) + phi.a.tst0[j,number+1]/2)
	}
	first.term = sum(weight.group0[xzero>landmark & szero <landmark]*phi.int)/(sum(weight.group0)*censor0.landmark) 
	
	
	#second term
	censor1.times = censor.weight(xone, deltaone, xone, weight = weight.group1)
	M.1 = rep(0,length(xone))
	M.1[xone>t] = t/censor1.t
	M.1[xone<=t] = xone[xone<=t]*deltaone[xone<=t]/censor1.times[xone<=t] 
	denom = sum(1*(sone > landmark & xone > landmark)*weight.group1)	
	psi.a.tt0 = sum(censor1.landmark*M.1*weight.group1*(sone > landmark & xone > landmark)/(denom))

	second.term = sum(weight.group0*(szero > landmark & xzero > landmark)*psi.a.tt0)/(sum(weight.group0)*censor0.landmark) 
	
	#third.term
	censor0.times = censor.weight(xzero, deltazero, xzero, weight = weight.group0)
	M.0 = rep(0,length(xzero))
	M.0[xzero>t] = t/censor0.t
	M.0[xzero<=t] = xzero[xzero<=t]*deltazero[xzero<=t]/censor0.times[xzero<=t] 
	denom = sum(1*( xzero > landmark)*weight.group0)	
	nu.term = sum(censor0.landmark*M.0*weight.group0*(xzero > landmark)/(denom))
	third.term = (sum(weight.group0*(xzero > landmark))/(sum(weight.group0)*censor0.landmark))*nu.term

	
	delta.q = first.term+second.term - third.term
	
	if(deltaslist) {return(list("delta.q" = delta.q))}
	if(!deltaslist) {return(c(first.term, second.term, third.term, delta.q))}


}

censor.weight = function(data.x, data.delta, t, weight=NULL) {
	if(is.null(weight)) {weight = rep(1,length(data.x))}
	S.KM = survfit(Surv(data.x,1-data.delta)~1, weights = weight)
	S.t.KM = approx(S.KM$time,S.KM$surv,t, rule=2)$y
	return(S.t.KM)
}

pred.smooth.surv <- function(xone.f, deltaone.f, sone.f, szero.one, myt, bw = NULL, weight, transform, extrapolate = T)
  { 
  	warn.flag = 0
    if(transform){	 	
    	mean.o= mean(c(sone.f, szero.one))
  		sd.o = sd(c(szero.one, sone.f))
    	sone.f.new = pnorm((sone.f - mean.o)/sd.o)
    	szero.one.new = pnorm((szero.one - mean.o)/sd.o)
		sone.f = sone.f.new
		szero.one = szero.one.new
	} 

	if(is.null(bw))
      {
        bwini = bw.nrd(sone.f)
        n.s = length(sone.f)
        bw <- bwini/(n.s^0.11)
      }      
    kerni.ss = Kern.FUN(zz=sone.f,zi=szero.one,bw)           
    tmpind = (xone.f<=myt)&(deltaone.f==1); tj = xone.f[tmpind]; 
    kerni.1 = t(weight*t(kerni.ss))
    pihamyt0.tj.ss = helper.si(tj, "<=", xone.f, Vi=t(kerni.1)) ## n.tj x n.ss matrix ##   
    dLamhat.tj.ss = t((kerni.1[,tmpind]))/pihamyt0.tj.ss; 
	#dLamhat.tj.ss[is.na(dLamhat.tj.ss)] = 0
    ret = apply(dLamhat.tj.ss,2,sum)
    Phat.ss  =exp(-ret)
    if(sum(is.na(Phat.ss))>0 & extrapolate){
    	warn.flag = 1
    	c.mat = cbind(szero.one, Phat.ss)
    	for(o in 1:length(Phat.ss)) {
    		if(is.na(Phat.ss[o])){
    			distance = abs(szero.one - szero.one[o])
    			c.temp = cbind(c.mat, distance)
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where predication is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			Phat.ss[o] = new.est[1]   #in case there are multiple matches
    	}
    }
	}
    return(list("Phat.ss"=Phat.ss, "warn.flag" = warn.flag))
    }

VTM<-function(vc, dm){
	#takes vc and makes it the repeated row of a matrix, repeats it dm times
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }
    

Kern.FUN <- function(zz,zi,bw=NULL,kern0="gauss") ## returns an (n x nz) matrix ##
  { 
    if(is.null(bw))
      {
        bwini = bw.nrd(zz)
        n.s = length(zz)
        bw <- bwini/(n.s^0.11)
      } 
     out = (VTM(zz,length(zi))- zi)/bw
    switch(kern0,
            "epan"= 0.75*(1-out^2)*(abs(out)<=1)/bw,
            "gauss"= dnorm(out)/bw
           )
  }
  
  
  cumsum2 <- function(mydat)     #cumsum by row, col remains the same
  {
    if(is.null(dim(mydat))) return(cumsum(mydat))
    else{
      out <- matrix(cumsum(mydat), nrow=nrow(mydat))
      out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
      return(out)
    }
  }

helper.si <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
  {  
    if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
    if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
    pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
    if(is.null(Vi)){return(pos)}else{
      Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
      out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
      out[pos!=0,] <- Vi[pos,]
      if(is.null(dim(Vi))) out <- c(out)
      return(out) ## n.y x p
    }
  }
 

R.t.estimate = function(xone,xzero, deltaone, deltazero, t, landmark, std = FALSE, conf.int = FALSE, weight.perturb = NULL) {
	warn.te = FALSE
	n1 = length(xone)
	n0 = length(xzero)
	if(is.null(weight.perturb)){
		weight.perturb = matrix(rexp(500*(n1+n0), rate=1), ncol = 500)
	}
	delta = delta.estimate.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t)
	delta.estimate = delta$delta
	delta.p.vec = apply(weight.perturb, 2, delta.estimate.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, delta.only = T)
	delta.p.vec = unlist(delta.p.vec)
	sd.delta  = sd(delta.p.vec)
	mad.delta = mad(delta.p.vec)
	conf.quantile.delta = c(new.q(delta.p.vec, 0.025), new.q(delta.p.vec, 0.975))
	if(0>conf.quantile.delta[1] & 0< conf.quantile.delta[2]) {warning("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the proportion of treatment effect explained in this setting", call. = FALSE)
		warn.te = TRUE}
	if(delta.estimate < 0) {warning("Warning: it looks like you need to switch the treatment groups", call. = FALSE)}
	delta.t.estimate = delta.t.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, landmark=landmark)$delta.t
	R.t = 1-delta.t.estimate/delta.estimate	
	if(std | conf.int){
		delta.t.p.vec.temp = t(apply(weight.perturb, 2, delta.t.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, landmark=landmark))
		delta.t.p.vec = unlist(delta.t.p.vec.temp)
		R.t.p = 1- delta.t.p.vec/delta.p.vec
	}
	if(conf.int) {	

		conf.l.quantile.delta = quantile(delta.p.vec, 0.025)
		conf.u.quantile.delta = quantile(delta.p.vec, 0.975)
		
		conf.l.quantile.delta.t = quantile(delta.t.p.vec, 0.025)
		conf.u.quantile.delta.t = quantile(delta.t.p.vec, 0.975)
		
		conf.l.quantile.R.t = quantile(R.t.p, 0.025)
		conf.u.quantile.R.t = quantile(R.t.p, 0.975)
}
	if(!std & !conf.int) {return(list("delta" = delta.estimate, "delta.t" =delta.t.estimate, "R.t" = R.t))}
	if(std & !conf.int) {return(list("delta" = delta.estimate, "delta.t" =delta.t.estimate, "R.t" = R.t, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.t.sd" = sd(delta.t.p.vec),  "delta.t.mad" = mad(delta.t.p.vec), "R.t.sd" = sd(R.t.p), "R.t.mad" = mad(R.t.p)))}
	if(conf.int) {return(list("delta" = delta.estimate, "delta.t" =delta.t.estimate, "R.t" = R.t, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.t.sd" = sd(delta.t.p.vec),  "delta.t.mad" = mad(delta.t.p.vec), "R.t.sd" = sd(R.t.p), "R.t.mad" = mad(R.t.p), "conf.int.delta" = as.vector(c(conf.l.quantile.delta, conf.u.quantile.delta)), "conf.int.delta.t" = as.vector(c(conf.l.quantile.delta.t, conf.u.quantile.delta.t)),
"conf.int.R.t" = as.vector(c(conf.l.quantile.R.t, conf.u.quantile.R.t))))
}	
}

 
delta.t.RMST = function(xone,xzero, deltaone, deltazero, t, weight = NULL, landmark=landmark) {
	if(is.null(weight)) {weight = rep(1,length(xone)+length(xzero))}
	weight.group1 = weight[1:length(xone)]
	weight.group0 = weight[(1+length(xone)):(length(xone)+length(xzero))]
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
	censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
	censor1.t = censor.weight(xone, deltaone, t, weight = weight.group1)
	censor1.landmark = censor.weight(xone, deltaone, landmark, weight = weight.group1)	
	
	#a term
	censor1.times = censor.weight(xone, deltaone, xone, weight = weight.group1)
	M.1 = rep(0,length(xone))
	M.1[xone>t] = t/censor1.t
	M.1[xone<=t] = xone[xone<=t]*deltaone[xone<=t]/censor1.times[xone<=t] 
	denom = sum(1*( xone > landmark)*weight.group1)	
	nu.term.a = sum(censor1.landmark*M.1*weight.group1*(xone > landmark)/(denom))
	a.term = (sum(weight.group0*(xzero > landmark))/(sum(weight.group0)*censor0.landmark))*nu.term.a
		
	#b.term
	censor0.times = censor.weight(xzero, deltazero, xzero, weight = weight.group0)
	M.0 = rep(0,length(xzero))
	M.0[xzero>t] = t/censor0.t
	M.0[xzero<=t] = xzero[xzero<=t]*deltazero[xzero<=t]/censor0.times[xzero<=t] 
	denom = sum(1*( xzero > landmark)*weight.group0)	
	nu.term.b = sum(censor0.landmark*M.0*weight.group0*(xzero > landmark)/(denom))
	b.term = (sum(weight.group0*(xzero > landmark))/(sum(weight.group0)*censor0.landmark))*nu.term.b

	
	delta.t = a.term - b.term
	
	return(list("delta.t" = delta.t))


}


IV.event = function(xone,xzero, deltaone, deltazero, sone, szero, t, landmark, number = 40, transform = FALSE, extrapolate = TRUE, std = FALSE, conf.int = FALSE, weight.perturb = NULL, type = "np") {
	if(!(type %in% c("np","semi"))) {warning("Warning: Invalid type, default `np' for nonparametric estimator being used", call. = FALSE); type = "np"}
	warn.te = FALSE
	warn.support = FALSE
	n1 = length(xone)
	n0 = length(xzero)
	if(is.null(weight.perturb)){
		weight.perturb = matrix(rexp(500*(n1+n0), rate=1), ncol = 500)
	}
	delta = delta.estimate.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t)
	delta.estimate = delta$delta
	delta.p.vec = apply(weight.perturb, 2, delta.estimate.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, delta.only = T)
	delta.p.vec = unlist(delta.p.vec)
	sd.delta  = sd(delta.p.vec)
	mad.delta = mad(delta.p.vec)
	conf.quantile.delta = c(new.q(delta.p.vec, 0.025), new.q(delta.p.vec, 0.975))
	if(0>conf.quantile.delta[1] & 0< conf.quantile.delta[2]) {warning("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the proportion of treatment effect explained in this setting", call. = FALSE)
		warn.te = TRUE}
	if(delta.estimate < 0) {warning("Warning: it looks like you need to switch the treatment groups", call. = FALSE)}
	range.1 = range(pmin(sone,landmark))
	range.0 = range(pmin(szero,landmark))
	range.ind = (range.1[1] > range.0[1]) | (range.1[2] < range.0[2])
	if(range.ind & !extrapolate & !transform) {
		warning("Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation", call. = FALSE)
		warn.support = TRUE
	}
	
	if(type == "np"){
		delta.q.estimate = delta.q.event.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, extrapolate = extrapolate, transform = transform, number = number)$delta.q
		R.q = 1-delta.q.estimate/delta.estimate	
	}
	if(type == "semi"){
		delta.q.estimate = delta.q.event.semi.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, number = number)$delta.q
		R.q = 1-delta.q.estimate/delta.estimate	
	}
	delta.t.estimate = delta.t.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, landmark=landmark)$delta.t
	R.t = 1-delta.t.estimate/delta.estimate	
	IV = R.q-R.t
	if(std | conf.int){
		delta.t.p.vec.temp = t(apply(weight.perturb, 2, delta.t.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, landmark=landmark))
		delta.t.p.vec = unlist(delta.t.p.vec.temp)
		R.t.p = 1- delta.t.p.vec/delta.p.vec
		if(type == "np"){
			delta.q.p.vec.temp = t(apply(weight.perturb, 2, delta.q.event.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, extrapolate = extrapolate, transform = transform, number = number, deltaslist=FALSE, warn.extrapolate = FALSE))
			delta.q.p.vec = delta.q.p.vec.temp[,4]
			R.p = 1-delta.q.p.vec/delta.p.vec
		}
		if(type == "semi"){
			delta.q.p.vec.temp = t(apply(weight.perturb, 2, delta.q.event.semi.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, deltaslist=FALSE, number = number))
			delta.q.p.vec = delta.q.p.vec.temp[,4]
			R.p = 1-delta.q.p.vec/delta.p.vec
		}
		IV.p = R.p-R.t.p
	}
	
	if(conf.int) {	

		conf.l.quantile.delta = quantile(delta.p.vec, 0.025)
		conf.u.quantile.delta = quantile(delta.p.vec, 0.975)
		
		conf.l.quantile.delta.q = quantile(delta.q.p.vec, 0.025)
		conf.u.quantile.delta.q = quantile(delta.q.p.vec, 0.975)
		
		conf.l.quantile.R.q = quantile(R.p, 0.025)
		conf.u.quantile.R.q = quantile(R.p, 0.975)
		
		conf.l.quantile.delta.t = quantile(delta.t.p.vec, 0.025)
		conf.u.quantile.delta.t = quantile(delta.t.p.vec, 0.975)
		
		conf.l.quantile.R.t = quantile(R.t.p, 0.025)
		conf.u.quantile.R.t = quantile(R.t.p, 0.975)
		
		conf.l.quantile.IV = quantile(IV.p, 0.025)
		conf.u.quantile.IV = quantile(IV.p, 0.975)
}
	if(!std & !conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q, "delta.t" =delta.t.estimate, "R.t" = R.t, "IV" = IV))}
	if(std & !conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q,  "delta.t" =delta.t.estimate, "R.t" = R.t, "IV" = IV, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.q.sd" = sd(delta.q.p.vec),  "delta.q.mad" = mad(delta.q.p.vec), "R.q.sd" = sd(R.p), "R.q.mad" = mad(R.p), "delta.t.sd" = sd(delta.t.p.vec),  "delta.t.mad" = mad(delta.t.p.vec), "R.t.sd" = sd(R.t.p), "R.t.mad" = mad(R.t.p), "IV.sd" = sd(IV.p), "IV.mad" = mad(IV.p)))}
	if(conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q,  "delta.t" =delta.t.estimate, "R.t" = R.t, "IV" = IV, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.q.sd" = sd(delta.q.p.vec),  "delta.q.mad" = mad(delta.q.p.vec), "R.q.sd" = sd(R.p), "R.q.mad" = mad(R.p), "delta.t.sd" = sd(delta.t.p.vec),  "delta.t.mad" = mad(delta.t.p.vec), "R.t.sd" = sd(R.t.p), "R.t.mad" = mad(R.t.p), "IV.sd" = sd(IV.p), "IV.mad" = mad(IV.p), "conf.int.delta" = as.vector(c(conf.l.quantile.delta, conf.u.quantile.delta)), "conf.int.delta.q" = as.vector(c(conf.l.quantile.delta.q, conf.u.quantile.delta.q)),
"conf.int.R.q" = as.vector(c(conf.l.quantile.R.q, conf.u.quantile.R.q)), "conf.int.delta.t" = as.vector(c(conf.l.quantile.delta.t, conf.u.quantile.delta.t)),
"conf.int.R.t" = as.vector(c(conf.l.quantile.R.t, conf.u.quantile.R.t)), "conf.int.IV" = as.vector(c(conf.l.quantile.IV, conf.u.quantile.IV))))
}	
}





new.q = function(x, p) {
	if(sum(is.na(x))>0) {return(NA)}
	if(sum(is.na(x))==0) {return(quantile(x,p))}
}
