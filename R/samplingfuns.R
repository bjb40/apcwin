#funcitons for running window sampler algorithm
#
#Bryce Bartlett


window.sample=function(var,alph){
  require(MCMCpack)

  #input is continuous of a,p,c
  #alpha is a vector the length of unique entries in var that is fed to dirichelet
  #output is factor with uniform, random window constraints
  #see dirichelet, e.g. for alternative algorithms
  vals=unique(var) #center on 0; assume continuous
  len=length(vals)

  alph=unlist(alph)

  dp=rdirichlet(1,alph)
  #segment the stick
  segments=round(cumsum(dp*len))
  #identify changes in segments
  sb=segments
  for(i in seq_along(segments)){
    if(i==1){next}
    sb[i] = segments[i]-segments[i-1]
  }

  #id inclusive breaks
  breaks=vals[sb>=1]
  #print(breaks)

  #because breaks are *inclusive*, must include min-1; ensure max
  if(min(breaks)>(min(var)-1)){
    breaks=c((min(var)-1),breaks)
  }
  if(max(breaks)<max(var)){
    breaks=c(breaks,max(var))
  }

  wins=window(var,breaks=breaks)

  return(wins)
}

draw_chains = function(dat,dv='y',apc=c('a','p','c'),
                       cores=1,method='gibbs',chains=1,samples=10,draws=1000){
#   source('config~.R')

  y = dat[,dv]
  #allmods=list() #may run into size constraints/may need to limit to best mods...
  effects=xhats=ppd=list()
  tm=Sys.time()
  avtm=0

dat$a = dat[,apc[1]]
dat$p = dat[,apc[2]]
dat$c = dat[,apc[3]]

#dat = dat %>%
#  rename(a=dat[,apc[1]],
#         p=dat[,apc[2]],
#         c=dat[,apc[3]])

#set of numbers of random samples
n.samples=samples

#holder df for model summary data
win = data.frame(a=as.numeric(NA),
                 p=as.numeric(NA),
                 c=as.numeric(NA))

modsum = data.frame( r2=as.numeric(NA),
                     sigma=as.numeric(NA),
                     bic=as.numeric(NA),
                     bic_prime=as.numeric(NA))

breaks=list(a=list(),p=list(),c=list())

d = c('a','p','c')
dl = c(length(unique(dat$a)),
       length(unique(dat$p)),
       length(unique(dat$c)))
names(dl) = d

#set starting values
all.alphas = lapply(d,function(x)
  data.frame(t(rep(dl[x]/2,length(unique(dat[,x]))))))

names(all.alphas) = d #names(all.nwins) = d

#accept rate
acc=0
#count boundary conditions rejections
bound=0

#mcmc sampler (prior model probabilities are equal)
for(s in 2:(n.samples+1)){

  #print(s)
  #reset dataframe
  x=dat[,c('a','p','c')]

  all.alphas= lapply(all.alphas, function(x)
    rbind(x,x[s-1,]+rnorm(ncol(x),mean=0,sd=0.5)))

  for(d in seq_along(all.alphas)){
    rownames(all.alphas[[d]]) = 1:nrow(all.alphas[[d]])
    }

  if(any(unlist(all.alphas)<0)){

    bound=bound+1
    out.al=sum(unlist(all.alphas)<0)
    #out.wi=sum(unlist(all.nwins)<2)
    #print(c(out.al,out.wi))
    #s=s-1
    #mnum=mnum-1 #should consolidate these
    #cat('\n\nOut-of-Sample-Space Warning.\n\n')
    #acc=acc-1
    for(d in seq_along(all.alphas)){
      #all.nwins[[d]][s]=all.nwins[[d]][s-1]
      all.alphas[[d]][s,]=all.alphas[[d]][s-1,]
      #note that this samples different windows with same hyper param
    }

  }

  #skip if unideintified
  la = length(levels(x$a)) == length(unique(dat$a))
  lp = length(levels(x$p)) == length(unique(dat$p))
  lc = length(levels(x$c)) == length(unique(dat$c))
  if(all(la,lp,lc)){
    for(d in seq_along(all.alphas)){
      #all.nwins[[d]][s]=all.nwins[[d]][s-1]
      all.alphas[[d]][s,]=all.alphas[[d]][s-1,]
    }
  }

  #draw random window samples
  x$a=window.sample(x$a,all.alphas$a[s,])
  x$p=window.sample(x$p,all.alphas$p[s,])
  x$c=window.sample(x$c,all.alphas$c[s,])

  #collect model data
  nr=data.frame(a=length(levels(x$a)),
                p=length(levels(x$p)),
                c=length(levels(x$c)))
  win=rbind(win,nr)

  #collect breaks data
  breaks$a[[s]]=attr(x$a,'breaks')
  breaks$p[[s]]=attr(x$p,'breaks')
  breaks$c[[s]]=attr(x$c,'breaks')

  #add esitmate
  mnum = s

  #reassign random references to each vector
  a.lev=length(levels(x$a)); a.b = sample(1:a.lev,1)
  p.lev=length(levels(x$p)); p.b = sample(1:p.lev,1)
  c.lev=length(levels(x$c)); c.b = sample(1:c.lev,1)

  form.c = as.formula(paste("~C(a,contr.treatment(a.lev,base=a.b))+
                            C(p,contr.treatment(p.lev,base=p.b))+
                            C(c,contr.treatment(c.lev,base=c.b))"))

  #generate model matrix
  xmat = model.matrix(form.c,data=x)

  m= NULL
  if(method=='gibbs'){
    m = lin_gibbs(y=y,x=xmat)
  } else if(method=='ml'){
    m = lin_ml(y=y,x=xmat)
  }

  modsum = rbind(modsum,
                 m[c('sigma','r2','bic','bic_prime')])

  #m = allmods[[s]] = tryCatch({lin_gibbs(y=y,x=xmat)},
  #                            finally=next)


  #if(s==1){next}

  #selection criterion
  #bayes factor approximation

  bf=exp((modsum[s,'bic']-modsum[s-1,'bic'])/2)
  #print(bf)
  #print(alphas)
  R = min(1,bf,na.rm=TRUE)
  #print(R)
  if (R < runif(1)){
    acc = acc+1
  } else {

    for(d in seq_along(all.alphas)){
      all.alphas[[d]][s,]=all.alphas[[d]][s-1,]
    }


  }


}#end sampling loop

breaks = lapply(breaks, function(x)
    x[2:length(x)])


res = list(
  modsum=modsum[2:nrow(modsum),],
#  effects=effects,
#  xhats=xhats,
  breaks=breaks,
  win=win[2:nrow(modsum),],
  n.samples=n.samples,
  acc=acc,
  bound=bound,
  method=method
)

return(res)

} #end draw chains


####
#function to draw predicted values

####
#sampler for window frame models

apcsamp = function(dat,dv='y',apc=c('a','p','c'),
                   cores=1,method='gibbs',
                   chains=1,samples=100,draws=1000){
#dat is a dataframe
#y is a character vector for "y"
#apc is a character vector for age, period, and cohort

require(parallel)

times = Sys.time()

  #generate limit cases for deleted cases
  limits = list()
  for(i in seq_along(apc)){
      vec = apc[1:3 != i]

      y = dat[,dv]
      x = as.data.frame(apply(dat[,vec],2,as.factor))
      xmat = model.matrix(~., data = x)

      limits[[i]] = lin_ml(y,xmat)

  }

  names(limits) = paste0('no_',apc)


  #create holder environment for shared variables
  par=new.env()

  #assign to environemnt to pass to clusters
  #assign to environemnt to pass to clusters
  assign('args',list(
    dat=dat,
    dv=dv,
    apc=apc,
    method=method,
    samples=samples,
    draws=draws
  ),env=par)

  cores=cores

  cat('\n\nBegin drawing', cores*samples,'total samples from',cores,'parallel chains.
      This will take some time. Calculating Estimate ...')
  start=Sys.time()
  tt=do.call(draw_chains,get('args',par))
  #sz=object.size(tt)
  tot=Sys.time()-start
  cat('\nOne sample took',round(tot,1),'seconds\n')
  #    'and uses',format(sz,'Mb'),'\n')
  cat(samples*cores,'parallelized would take at least',
      round((samples*tot)/60,2),'minutes.\n\n')
  #    'and use',paste0(format(sz*cores*samples,'Gb'),'.\n\n'))

  print(Sys.time())

  cl <- makeCluster(mc <- getOption("cl.cores", cores))
  clusterExport(cl=cl,varlist='args',envir=par)
  #clusterExport(cl=cl, varlist=c("dat","apc",'method','samples','draws'))
  chains=parLapply(cl=cl,1:cores,function(...){
    args=get('args',env=par)
    require(apcwin)
    draw_chains(dat=args$dat,
                dv=args$dv,
                apc=args$apc,
                method=args$method,
                samples=args$samples,
                draws=args$draws)
  })

  #end cluster
  rm(par)
  cat('\n\nEnd drawing chains...')
  stopCluster(cl)

  modsum = extract(chains,'modsum', as.df=TRUE)

  #check against limits
  for(i in names(limits)){
    #exp((modsum[s,'bic']-modsum[s-1,'bic'])/2)
    modsum[,paste0(i,'_bf')] =
      (modsum$bic - limits[[i]]$bic) / 2
    modsum[,paste0(i,'_pfprime')] =
      (modsum$bic_prime - limits[[i]]$bic_prime) / 2
  }

  #calculate weight by bic and bic_prime
  k=min(modsum$bic)
  d=-.5*(modsum$bic-k)
  modsum$w=exp(d)/sum(exp(d))

  k = min(modsum$bic_prime)
  d=-.5*(modsum$bic_prime-k)
  modsum$w_prime=exp(d)/sum(exp(d))

  t.samples = sum(extract(chains,'n.samples'))

  breaks = list(
    a=do.call(c,lapply(chains,function(x)
      x[['breaks']][['a']])),
    p=do.call(c,lapply(chains,function(x)
      x[['breaks']][['p']])),
    c=do.call(c,lapply(chains,function(x)
      x[['breaks']][['c']]))
  )

  fullt = Sys.time() - times



  res = list(
  data = dat,
  apc = apc,
  dv = dv,
  limits = limits,
  summaries = cbind(
    extract(chains,'win',as.df=TRUE),
            modsum),
  breaks = breaks,
  chains = length(chains),
  n.samples = t.samples,
  acc = sum(extract(chains,'acc'))/t.samples,
  bound = sum(extract(chains,'bound'))/t.samples,
  method = extract(chains,'method'),
  timespent = fullt)

  class(res) = append(class(res),'apcsamp')

  return(res)

}

###helper function for chain extraction
extract = function(l,name,as.df=FALSE,span=NULL){
  #extracts and combines name object from list 'l'
  #span limits the extraction, if, for example, I need to drop first off
  if(is.null(span) & !as.df){span=1:length(l[[1]][[name]])}

  res=do.call(base::c,lapply(l,function(x) x[[name]][span]))
  if(as.df){
    res=do.call(rbind,lapply(l,function(x) x[[name]]))
  }

  return(res)
}


###########
#prepare estimated effects

###this needs work.... shoul feed it an a and a b
###will need a fully bayesian esimtator
draw_effs = function(sampobj,tol=NULL){
  #sampobj is an apcsamp object
  #tol is a number under 1 that idientifies
  #how many posterior samples
  #returns apceffects object

  #tolerance is the inverse of the number of samples
  #from sampobj

  dat = sampobj$data

  #generate holer for effects dataframe (target of function)
  apcvals = lapply(apc,function(x) unique(dat[,x]))
  names(apcvals) = apc

  effects = lapply(apcvals,function(x)
    data.frame(matrix(vector(),0,length(x))))

  fits = data.frame(matrix(vector(),0,3))
  colnames(fits) = c('r2','bic','bic_prime','sigma')

  if(is.null(tol)){tol=1/sampobj$n.samples}
  n.samples = 1/tol

  ###
  #create random index object and draw with replacement
  s.index = sample(1:nrow(sampobj$summaries),
                   n.samples,replace=TRUE,
                   prob=sampobj$summaries$w)

  ###
  #iterate through, drawing one Bayesian posterior each



  for(s in s.index){

  #sample number (easier to reference)
  snum = s.index[s]

  #reset dataframe
  x=dat[,c('a','p','c')]

  #set window breaks
  x$a = window(x$a,breaks=sampobj$breaks$a[[snum]])
  x$p = window(x$p,breaks=sampobj$breaks$p[[snum]])
  x$c = window(x$c,breaks=sampobj$breaks$c[[snum]])

  #set random reference -- query whether need to
  #save reference instead; answer is probably...
  #reassign random references to each vector
  a.lev=length(levels(x$a)); a.b = sample(1:a.lev,1)
  p.lev=length(levels(x$p)); p.b = sample(1:p.lev,1)
  c.lev=length(levels(x$c)); c.b = sample(1:c.lev,1)

  y = dat[,sampobj$dv]
  form.c = as.formula(paste("~C(a,contr.treatment(a.lev,base=a.b))+
                            C(p,contr.treatment(p.lev,base=p.b))+
                            C(c,contr.treatment(c.lev,base=c.b))"))

  #generate model matrix
  xmat = model.matrix(form.c,data=x)

  ####
  #draw bayesian posterior
  m = lin_gibbs(y=y,x=xmat,iter=1)
  fits = rbind(fits,m[c('r2','bic','bic_prime','sigma')])

  #create grand means
  grand.means = data.frame(
    a = window(mean(dat$a),breaks=attr(x$a,'breaks')),
    p = window(mean(dat$p),breaks=attr(x$p,'breaks')),
    c = window(mean(dat$c),breaks=attr(x$c,'breaks'))
  )

  grand.means$a=relevel(grand.means$a,ref=a.b)
  grand.means$p=relevel(grand.means$p,ref=p.b)
  grand.means$c=relevel(grand.means$c,ref=c.b)

  grand.means=(model.matrix(form.c,grand.means))

  #generate dummy variables to describe full range of dimension
  blockdat = list()
  blockdat$a = scopedummy(w=x$a,unique.vals=unique(dat$a))
  blockdat$p = scopedummy(w=x$p,unique.vals=unique(dat$p))
  blockdat$c = scopedummy(w=x$c,unique.vals=unique(dat$c))

  blockdat$a = relevel(blockdat$a,ref=a.b)
  blockdat$p = relevel(blockdat$p,ref=p.b)
  blockdat$c = relevel(blockdat$c,ref=c.b)

  predat=lapply(blockdat,FUN=function(x)
    model.matrix(~.,data=as.data.frame(x)))

  betas=list()
  for(eff in names(predat)){
    #fix colnames
    colnames(predat[[eff]]) = sub('x',eff,colnames(predat[[eff]]))

    #calculate means for xhat, & id effects at issue---this was replaced...
    xhat=grand.means[rep(seq(nrow(grand.means)), nrow(predat[[eff]])),]

    #replace means of effect dimensions with indicator in matrix
    calceff = grepl(paste0(eff,'.lev|Intercept'),colnames(xhat))
    xhat[,calceff] = predat[[eff]]

    effects[[eff]] = rbind(effects[[eff]],t(xhat %*% t(m$betas)))

  }


  }#end of sampling loop for "s"

  for(i in seq_along(effects)){
    colnames(effects[[i]]) = apcvals[[i]]}

  res = list(sampobj=sampobj,
             fit=fits,
             sampled=s.index,
             effects=effects)

  class(res) = 'apceffects'

  return(res)


  }



###########3
#summary function for effects class

summary.apcsamp = function(samp){

  ss=samp$summaries
  best=which(ss$bic == min(ss$bic))
  ####
  #get mean bic of APC
  cat('\n\nSamples:', samp$n.samples,
      '\tover',samp$chains,'chains.\n')
  cat('Acceptance:',samp$acc,
      '\tBoundary Faults:',samp$bound)


  cat('\n\nBest-fitting from full APC\n')
  cat('BIC:\t',ss$bic[best])
  cat('\nR-squared:\t',ss$r2[best])

  cat('\n\n\nMean Number Window Breaks:\n')
  print(apply(ss[,1:3],2,weighted.mean,w=ss$w))

  cat('Summary for limit case of null effects:\n')
  cat('BIC\n')
  print(unlist(lapply(samp$limits,function(x) x$bic)))
  cat('R-squared\n')
  cat(unlist(lapply(samp$limits,function(x) x$r2)))

  cat('\n\nBayes Factors\n')
  bf = ss[best,c(8,10,12)]
  print(bf)

}
