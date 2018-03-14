#Dev R 3.3.0 "Supposedly Educational"
#general functions for apc project
#
#Bryce Bartlett

#make an S3 object:
#http://www.cyclismo.org/tutorial/R/s3Classes.html#memory-management

#print_current = function(){
#  print('\nmemoryfix branch')
#}

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
    stop("
Must specify either a window length (winlength)
or vector of breaks (breaks -- works like base cut)")
  }

  #consider warning for small cell values
  #consider error for out-of-range window or break sets

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
  #inherit the name???

  class(c) = append(class(c),'window')

  return(c)

}

scopedummy=function(w,unique.vals=NULL){
  #this is a method for window object that transforms it
  #into a factor covering all possible continuous values
  #the (row) names are the values
  #i.e. 1,2,3,4,5,6 with two windows output two factors
  # 1,1,1,2,2,2
  #it can supply unique values, otherwise, it presumes continuous

  UseMethod('scopedummy',w)

}

scopedummy.default = function(w,unique.vals=NULL) {
  cat('\nError: Must use window object.\n\n')
}

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

range=function(w){
  UseMethod('range',w)
}

relevel.window = function(w,ref){
  #relevels window object ---
  #the only thing this does is to preserve the "breaks"
  #attribute so that it is still a window

  nw = relevel(w,ref)
  attr(nw,'breaks') = attr(w,'breaks')
  class(nw) = append(class(nw),'window')

  return(nw)
}

range.window = function(w){
  #returns range object (from underlying continuous)
  return(attr(w,'range'))
}

expand=function(w){
  UseMethod('expand',w)
}

expand.window=function(w){
  r=range(w)
  return(r[1]:r[2])
}

#ols funciton for estimating an apc model
apc_lm = function(formula,data,age,per,coh,windows) {
  #runs a linear model on an apc problem (need to id window)
  #
  # Inputs
  #  formula=r formula object
  #  data = dataframe
  #  age,per,coh = named columns of dataframe for age, period, cohort
  #  need only specify 2... If all 3 are specified, model checks for
  #  linear dependence
  #  windows = list identifying window constraints across apc
  #            if empty, will pick random windows min=3,max=max(variable))
  #
  # Output
  #   list including
  #   (1) effects =  ols estimates from window model
  #   (2) smallblock = block APC estimates
  #   (3) linear = linearlized block estimates (cubic)

  #@@@@
  #input checks
  no.age=no.period=no.cohort=FALSE

  if(missing(age)){
    no.age=TRUE
    age='age'
  }
  if(missing(per)){
    no.period=TRUE
    per='period'
  }
  if(missing(coh)){
    no.cohort=TRUE
    coh='cohort'
  }

  if(sum(no.period,no.cohort,no.age)>1){
    stop('Must specify 2 of 3 required estimates: age, period, or cohort.')
  } else if(no.age){
      data[,age]=data[,per]-data[,coh]
  } else if(no.period){
      data[,per]=data[,age]+data[,coh]
  } else if(no.cohort){
      data[,coh]=data[,per]-data[,age]
  }

  #@@@@
  #window check
  if(missing(windows)){
    windows=list(age=0,period=0,cohort=0)

    #id maximum window for constraint(min is 3)
    max_w=function(var){
      r=range(var)
      return(r[2]-r[1])
      }
    windows$age=round(runif(1,3,max_w(data[,age])))
    windows$period=round(runif(1,3,max_w(data[,per])))
    windows$cohort=round(runif(1,3,max_w(data[,coh])))

  }

  #@@@@
  #build model matrix from window constraints
  wins=list(
    a=window(data[,age],winlength=windows$age),
    p=window(data[,per],winlength=windows$period),
    c=window(data[,coh],winlength=windows$cohort)
  )


  ndat=data[,!colnames(data) %in% c(age,per,coh)]
  ndat=cbind(ndat,wins$a,wins$p,wins$c)
  colnames(ndat)[-1:(3-ncol(ndat))] = c(age,per,coh)

  blockdat=lapply(wins,scopedummy);
  names(blockdat)=c(age,per,coh)

  #@@@@
  #estimate OLS model based on window constraints

  results=lm(formula,data=ndat)

  #@@@@
  #prepare small block estimates

  b=coefficients(results)
  cov=vcov(results)

  predat=lapply(blockdat,FUN=function(x)
    model.matrix(~.,data=as.data.frame(x)))

  blockeff=list(beta=list(),cov=list())
  predict=list(est=list(),cov=list())
  beta.hat=list()

  for(eff in c(age,per,coh)){
     blockeff[['beta']][[eff]] = b [grepl(paste0(eff,'|Intercept'),names(b))]
     blockeff[['cov']][[eff]] = cov[grepl(paste0(eff,'|Intercept'),rownames(cov)),
            grepl(paste0(eff,'|Intercept'),colnames(cov))]

  predict[['est']][[eff]]=predat[[eff]] %*% blockeff[['beta']][[eff]]

  #this is _assuming_ at default (left out) variables
  #should re-estimate based on means??

  predict[['cov']][[eff]] = predat[[eff]] %*% blockeff[['cov']][[eff]] %*% t(predat[[eff]])

  beta.hat[[eff]] = data.frame(cohort=expand(ndat[,eff]),
                        est=predict[['est']][[eff]],
                        se=sqrt(diag(predict[['cov']][[eff]])))
  }




  #@@@@
  #prepare linear estimates(cubic)

  #@@@@
  #return

  return(
    list(
      newdat=ndat,
      results=results,
      block.eff=beta.hat
    )
  )

}

plt = function(ols,varnames){
  #this function takes a,p,c measures
  #and plots them in a panel of 3
  #dependency -- ggplot2 (can remove for package)
  #
  #Input:
  #     ols=lm object
  #     varnames=charater vector (1-3) identifying
  #     variablenames for age period and cohort (in that order)
  #     This will work for prefixes...
  #
  #Output: list of ggplot objects

  #select coefficients from
  b = coefficients(ols)

  return(b)

}


##############
#gibbs sampler
###############

lin_gibbs = function(y,x,iter=1000){
  #iter = 1000

  rmse=ll=r2=s2=matrix(1,iter)
  b= matrix(0,iter,ncol(x))
  #ytilde=yhat=matrix(0,iter,length(y))
  xtxi = solve(t(x)%*%x,tol=1e-22)
  m=lm(y~x-1) #why no interecept???
  pars=coefficients(m)
  res=residuals(m)
  n=length(y)

  #simulate sigma from inverse gamma marginal
  s2 = 1/rgamma(iter,nrow(x)-ncol(x)/2,.5*t(res)%*%res)

  #set ppd
  ppd = matrix(0,iter,length(y))

  for (i in 1:iter){
    #print(i)
    #simulate beta from mvn
    b[i,]=pars+t(rnorm(length(pars),mean=0,sd=1))%*%chol(s2[i]*xtxi)
    yhat = x %*% b[i,]
    sse = sum((y-(yhat))^2)
    sst = sum((y-mean(y))^2)
    r2[i] = 1-(sse/sst)
    rmse[i] = sqrt(sse/n)
    ll[i]=sum(dnorm(y,mean=yhat[i,],sd=s2[i],log=TRUE))
  }

  colnames(b) = colnames(x)
  ###BIC estimate for Bayes Factor (posterior probability weight)
  #p. 135 explains the differnce between bic and bic_prime
  #bic is "overall model fit" where
  #bic prime "provides an assessment of whether the model is explaining enough
  #variation to justify the number of parameters it's using..."
  ###p. 135, Eq. 26 from Rafferty 1995 (SMR)

  bic_prime=n*log(1-mean(r2))+(ncol(x)-1)*log(n)
  #bic equation ...eq 23 http://www.stat.washington.edu/raftery/Research/PDF/kass1995.pdf
  bic=-2*mean(ll)+log(n)*ncol(x)

  #print(dim(yhat))

  #sigma is poorly named
  return(list(betas=b,
              y=y,
              x=x,
              #yhat=yhat,
              sigma=s2,
              #ytilde=yhat+rnorm(length(yhat),mean=0,sd=s2),
              r2=r2,
              #rmse=rmse,
              bic=bic,
              bic_prime=bic_prime#,
              #ll=ll
              ))

}#end linear gibbs

##############
#ML sampler --- can use to identify best fits...
#iterates faster, but not a true bayesian
###############

lin_ml = function(y,x){

  n=length(y)
  m=lm(y~x-1)
  b=coefficients(m)
  s2=summary(m)$sigma
  r2=summary(m)$r.squared
  bic_prime=n*log(r2)+(ncol(x)-1)*log(n)
  ll=sum(dnorm(y,mean=predict(m),sd=s2))

  #bic equation ...eq 23 http://www.stat.washington.edu/raftery/Research/PDF/kass1995.pdf
  bic=-2*mean(ll)+log(n)*ncol(x)

  return(list(betas=matrix(b,1,length(b)),
              y=y,
              x=x,
              #yhat=yhat,
              sigma=s2,
              #ytilde=yhat+rnorm(length(yhat),mean=0,sd=s2),
              r2=r2,
              #rmse=rmse,
              bic=bic,
              bic_prime=bic_prime#,
              #ll=ll
  ))
}

#@@@@@@@@@@@@@@@@@@@@@
#Misc Functions/objects
#@@@@@@@@@@@@@@@@@@@@@

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




