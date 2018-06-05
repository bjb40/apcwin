#general functions for apc project


#' @export
window = function(var,winlength,breaks){
  #this function makes windows from a continuous meaasures
  #
  #Input: var = numeric vector
  #       winlength = number identifying width of windows
  #       breaks = vector of breaks (like cut in base)
  #       can only provide 1 winlenth or breaks

  #check and throw errors
  if(missing(winlength)){
    winlength=NA
  }

  if(missing(breaks)){
    nobreak=TRUE
  } else{nobreak=FALSE}

  if(nobreak & is.na(winlength)){
    stop('
Must specify either a window length (winlength)
or vector of breaks (breaks -- works like "cut" in base R.)')
  }

  r=range(var)

if(!is.na(winlength)){
  w=winlength
  #full possibility
  vec=r[1]:r[2]
  #vector of windows
  win=1:w
  breaks = r[1]-1
  for(i in 0:(ceiling(length(vec)/w)-1)){
    p.break=vec[max(win+w*i)]
    if(is.na(p.break)){p.break=r[2]}
    breaks=append(breaks,p.break)
  }
}#end winlength

  c=cut(var,breaks=breaks)
  attr(c,'range') = r
  attr(c,'breaks') = breaks

  class(c) = append(class(c),'window')

  return(c)

}

#' @export
scopedummy=function(w,unique.vals=NULL){
  #this is a method for window object that transforms it
  #into a factor covering all possible continuous values
  #the (row) names are the values
  #i.e. 1,2,3,4,5,6 with two windows output two factors
  # 1,1,1,2,2,2
  #it can supply unique values, otherwise, it presumes continuous

  UseMethod('scopedummy',w)

}

#' @export
scopedummy.default = function(w,unique.vals=NULL) {
  cat('\nError: Must use window object.\n\n')
}

#' @export
scopedummy.window = function(w,unique.vals=NULL) {
  #w is a window object
  #unique.vals is a vector of numbers that shows unique values, if null, then it
  #presumes continuous
  r=attr(w,'range')
  if(is.null(unique.vals)){
    span=r[1]:r[2]
  } else {
    span=unique.vals
  }

  f=cut(span,breaks=attr(w,'breaks'))
  return(f)
}

#' @export
range=function(w){
  UseMethod('range',w)
}

#' @export
relevel.window = function(w,ref){
  #relevels window object ---
  #the only thing this does is to preserve the "breaks"
  #attribute so that it is still a window

  nw = relevel(w,ref)
  attr(nw,'breaks') = attr(w,'breaks')
  class(nw) = append(class(nw),'window')

  return(nw)
}

#' @export
range.window = function(w){
  #returns range object (from underlying continuous)
  return(attr(w,'range'))
}

#' @export
expand=function(w){
  UseMethod('expand',w)
}

#' @export
expand.window=function(w){
  r=range(w)
  return(r[1]:r[2])
}

##############
#gibbs sampler
###############

#' @export
lin_gibbs = function(y,x,iter=1000){
  rmse=ll=r2=s2=matrix(1,iter)
  b= matrix(0,iter,ncol(x))
  xtxi = solve(t(x)%*%x,tol=1e-22)
  #note: no intercept because i feed it one n the x matrix
  m=lm(y~x-1)
  pars=coefficients(m)
  res=residuals(m)
  n=length(y)

  #simulate sigma from inverse gamma marginal
  s2 = sqrt(1/rgamma(iter,nrow(x)-ncol(x)/2,t(res)%*%res))

  for (i in 1:iter){
    #simulate beta from mvn
    b[i,]=pars+t(rnorm(length(pars),mean=0,sd=1))%*%chol(s2[i]*xtxi)
    yhat = x %*% b[i,]
    sse = sum((y-(yhat))^2)
    sst = sum((y-mean(y))^2)
    r2[i] = 1-(sse/sst)
    rmse[i] = sqrt(sse/n)
    ll[i]=sum(dnorm(y,mean=yhat,sd=s2[i],log=TRUE))
  }

  colnames(b) = colnames(x)
  ###BIC estimate for Bayes Factor (posterior probability weight)
  #calculations based on Rafferty 1995 (SMR)
  #p. 135 explains the differnce between bic and bic_prime
  #bic is "overall model fit" where
  #bic prime "provides an assessment of whether the model is explaining enough
  #variation to justify the number of parameters it's using..."
  ###p. 135, Eq. 26

  bic_prime=n*log(1-mean(r2))+(ncol(x)-1)*log(n)
  #bic equation ...eq 23 http://www.stat.washington.edu/raftery/Research/PDF/kass1995.pdf
  bic=-2*mean(ll)+log(n)*ncol(x)
  aic=-2*mean(ll)+2*ncol(x)

  return(list(betas=b,
              y=y,
              x=x,
              sigma=s2,
              r2=r2,
              bic=bic,
              aic=aic,
              bic_prime=bic_prime,
              ll=mean(ll)
              ))

}#end linear gibbs

##############
#ML sampler --- can use to identify best fits...
#iterates faster, but not a true bayesian
###############

#' @export
lin_ml = function(y,x){

  n=length(y)
  m=lm(y~x-1)
  b=coefficients(m)
  s2=summary(m)$sigma
  r2=summary(m)$r.squared
  bic_prime=n*log(1-r2)+(ncol(x)-1)*log(n)
  ll=sum(dnorm(y,mean=predict(m),sd=s2,log=TRUE))

  #bic equation ...eq 23 http://www.stat.washington.edu/raftery/Research/PDF/kass1995.pdf
  bic=-2*ll+log(n)*ncol(x)
  aic=-2*ll+2*ncol(x)

  return(list(betas=matrix(b,1,length(b)),
              y=y,
              x=x,
              sse=sum(residuals(m)^2),
              rank = m$rank,
              sigma=s2,
              r2=r2,
              bic=bic,
              aic=aic,
              bic_prime=bic_prime,
              ll=ll
  ))
}

#@@@@@@@@@@@@@@@@@@@@@
#Misc Functions/objects
#@@@@@@@@@@@@@@@@@@@@@

#' @export
rnd = function(db,rd){
  # rounds input to preserve leading zeros
  #
  # Args:
  #   db: an object with numeric types
  #   rd: length to round (including leading zeros, default=3)
  #
  # Returns:
  #   an object of of db translated to characters with leading zeros

  if(missing(rd)){rd=3}
  rdl=paste0('%.',rd,'f')
  return(sprintf(rdl,round(db,digits=rd)))
}

#' @export
sig = function(pv){
  # returns stars based on pvalue
  #
  # Args:
  #   pv: a p-value
  #
  # Returns:
  #   a string with stars for values * <.05 **<.01 *** < .001
  s=' '
  if(length(pv)>0){
    if(pv<.001){s='***'} else if(pv<.01){s='**'} else if (pv<.05){s='*'} else if (pv<.1){s='+'}
  }
  return(s)

}




