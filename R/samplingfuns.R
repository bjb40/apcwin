#funcitons for running window sampler algorithm
#
#Bryce Bartlett

#' @export
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

#' Draw a single chain.
#'
#' @param dat A dataframe.
#' @param dv A character variable or column number indicating the dependent variable.
#' @param apc A character vector naming the age, period, and cohort variables (defaults to 'a', 'p', and 'c').
#' @param cores Integer indicating cores to use.
#' @param method One of "gibbs" or "ml". For a fully bayesian estimator for a faster approximation.
#' @param samples Integer indicating samples; used only for "gibbs" method.
#' @param draws Integer indicating the number of models to sample in the chain.
#' @param starvals Optional list for chain starting values.
#' @return An object of apcsamp.
#' @examples
#' data(apcsamp)
#' draw_chains(apcsamp,method='ml',samples=250,cores=4,chains=4)
draw_chains = function(dat,dv='y',apc=c('a','p','c'),
                       cores=1,
                       method='gibbs',
                       chains=1,
                       samples=10,
                       draws=1000,
                       startvals=NULL){

  y = dat[,dv]
  effects=xhats=ppd=list()
  tm=Sys.time()
  avtm=0

dat$a = dat[,apc[1]]
dat$p = dat[,apc[2]]
dat$c = dat[,apc[3]]

#set of numbers of random samples
n.samples=samples

#holder df for model summary data
win = data.frame(a=as.numeric(NA),
                 p=as.numeric(NA),
                 c=as.numeric(NA))

modsum = data.frame( r2=as.numeric(NA),
                     #sse = as.numeric(NA),
                     #rank = as.numeric(NA),
                     #n=nrow(dat),
                     sigma=as.numeric(NA),
                     bic=as.numeric(NA),
                     aic=as.numeric(NA),
                     bic_prime=as.numeric(NA))

breaks=list(a=list(),p=list(),c=list())

d = c('a','p','c')
dl = c(length(unique(dat$a)),
       length(unique(dat$p)),
       length(unique(dat$c)))
names(dl) = d

#set starting values
#use ratio of a/number of windows; sample a

  all.alphas = lapply(d,function(x)
    data.frame(t(rep(dl[x]/dl[x],length(unique(dat[,x]))))))

  names(all.alphas) = d #names(all.nwins) = d

if(!is.null(startvals)){
  svals = unlist(lapply(startvals,length))
  alphvals = unlist(lapply(all.alphas,length))
  errormessage = paste(c("Problem with dimensions of supplied starting values.",
  "\n\tNumber of values supplied:",svals,
  "\n\tNumber of values needed:",alphvals),collapse=' ')


  if(!all(svals==alphvals)){
    stop(errormessage)
  } else{
    all.alphas = startvals
  }

}

#accept rate
acc=0
#count boundary conditions rejections
bound=0

#mcmc sampler (prior model probabilities are equal)
for(s in 2:(n.samples+1)){

  #print(s)
  #reset dataframe
  x=dat[,c('a','p','c')]

  #draw numerator for nwindows sd =0.5
  #divide by nwindows
  nalph = list()
  for(d in names(dl)){
    nalph[[d]] = all.alphas[[d]][s-1,]*dl[[d]]
    nalph[[d]] = nalph[[d]] +
      rnorm(ncol(nalph[[d]]),mean=0,sd=0.5)
    all.alphas[[d]] = rbind(all.alphas[[d]],
                            nalph[[d]]/dl[[d]])
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
                 m[c('sigma','r2','bic','aic','bic_prime')])

  #specialized anova approximation
  #bf=n*log(sum(residuals(m1)^2)/sum(residuals(m2)^2)) + (k1-k2)*log(n)
  #simple and old approximation

  #original--
  #bf=exp((modsum[s,'bic']-modsum[s-1,'bic'])*.5) #.5 is even prior

  #odds
  modds = exp((modsum[s,'bic']-modsum[s-1,'bic']))
  mprob = modds/(1+modds)

  #jumping kernel
  R = min(1,mprob,na.rm=TRUE)
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
  alphas=all.alphas,
  breaks=breaks,
  win=win[2:nrow(modsum),],
  n.samples=n.samples,
  acc=acc,
  bound=bound,
  method=method
)

return(res)

} #end draw chains

#helper funciton to listwise delete --
#can update to deal with missing data and controls later
#designed for use with apcsamp
checkdat = function(...){

  dat = data.frame(dat[,c(dv,apc)])
  comp = complete.cases(dat)

  if(any(comp)){
    cat('Deleting',sum(!comp), 'incomplete cases')
    dat = dat[comp,]
  }

  return(dat)

}



#' A wrapper for "draw_chains" that will use multiple search chains in the algorithm.
#'
#' @param dat A dataframe.
#' @param dv A character variable or column number indicating the dependent variable.
#' @param apc A character vector naming the age, period, and cohort variables (defaults to 'a', 'p', and 'c').
#' @param cores Integer indicating cores to use.
#' @param method One of "gibbs" or "ml". For a fully bayesian estimator for a faster approximation.
#' @param samples Integer indicating samples; used only for "gibbs" method.
#' @param draws Integer indicating the number of models to sample in the chain.
#' @param startvals Optional list for chain starting values.
#'
#' @details what does this look like.
#'
#' @return An object of apcsamp.
#'
#' @examples
#' data(apcsim)
#' apcsamp(apcsim,method='ml',samples=250,cores=4,chains=4)
#'
#' @export
apcsamp = function(dat,dv='y',apc=c('a','p','c'),
                   cores=1,method='gibbs',
                   chains=1,samples=100,draws=1000){


require(parallel)

#make sure data is in proper form, an listwise delete
  dat = checkdat()


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

  cl <- parallel::makeCluster(mc <- getOption("cl.cores", cores))
  parallel::clusterExport(cl=cl,varlist='args',envir=par)
  #clusterExport(cl=cl, varlist=c("dat","apc",'method','samples','draws'))
  chains=parallel::parLapply(cl=cl,1:cores,function(...){
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
  parallel::stopCluster(cl)

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
  #https://stats.stackexchange.com/questions/249888/use-bic-or-aic-as-approximation-for-bayesian-model-averaging
  k=min(modsum$aic)
  d=-.5*(modsum$aic-k)
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

  alphas = list(
    a = do.call(rbind,lapply(chains, function(x)
      x[['alphas']][['a']])),
    p = do.call(rbind,lapply(chains,function(x)
      x[['alphas']][['p']])),
    c = do.call(rbind,lapply(chains,function(x)
      x[['alphas']][['c']]))
  )


  fullt = Sys.time() - times

  res = list(
  data = dat,
  apc = apc,
  dv = dv,
  alphas = alphas,
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

#' @export
convergence = function(samp) {
  UseMethod('convergence',samp)
}


#' @export
convergence.default = function(samp) {
  cat('\nError: Must use apcsamp object.\n\n')
}

###helper function for convergence
#' @export
convergence.apcsamp = function(samp){
  require(coda)

  d = samp[['apc']]

  alphas = lapply(d,function(dims){
    alph = as.data.frame(samp$alpha[[dims]])
    colnames(alph) = paste0(dims,
      unique(samp$data[,dims])[order(unique(samp$data[,dims]))])
    return(alph)
  })

  alphas = do.call(cbind,alphas)
  nc = ncol(alphas)

  chainlen = nrow(alphas)/samp$chains
  chain = factor(ceiling(1:nrow(alphas)/chainlen))

  alphas =split(as.data.frame(alphas),chain)

  alphas = lapply(alphas,coda::mcmc)
  alphas = coda::as.mcmc.list(alphas)

  res = do.call(rbind,coda::gelman.diag(alphas,
                                        multivariate=FALSE,
                                        autoburnin=TRUE))

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

#' Draws effect estimates from an APC sample using Bayesian Model Averaging.
#'
#' @param sampobj An apcsamp object (\code{\link{apcsamp}})
#' @param means Logical: TRUE returns mean estimates; FALSE generates centered estimates (relatvie to grand mean).
#' @param betas Logical: TRUE returns a vector of the betas.
#' @param tol A number between 0 and 1, which indicates tolerance to rounding and the number
#' to be drawn. e.g. tol=0.001 is tolerance to the third decimal and draws 1,000 samples.
#' @return An object of apceffects.
#' @examples
#' data(apcsim)
#'
#' sampleobject = apcsamp(apcsim,method='ml',
#' samples=250,cores=4,chains=4)
#'
#' effectsobject = draw_effs(sampleobject)
#'
#' @export
draw_effs = function(sampobj,
                        means=FALSE,
                        betas=FALSE,
                        tol=NULL){

  require(dplyr)

  dat = sampobj$data

  #generate holer for effects dataframe (target of function)
  apcvals = lapply(sampobj$apc,
                   function(x) unique(dat[,x]))
  names(apcvals) = sampobj$apc

  effects = lapply(apcvals,function(x)
    data.frame(matrix(vector(),0,length(x))))

  fits = data.frame(matrix(vector(),0,4))
  colnames(fits) = c('r2','bic','bic_prime','sigma')

  if(is.null(tol)){tol=1/sampobj$n.samples}
  n.samples = floor(1/tol)

  ###
  #create random index object and draw with replacement
  s.index = sample(1:nrow(sampobj$summaries),
                   n.samples,replace=TRUE,
                   prob=sampobj$summaries$w)


  if(betas){raw.betas = list()}
  intercept=vector()

  ###
  #iterate through, drawing one Bayesian posterior each

  for(s in s.index){

    #sample number (easier to reference)
    #snum = s.index[s]

    #reset dataframe
    x=dat[,c('a','p','c')]

    breaks = sampobj$breaks
    #set window breaks
    x$a = window(x$a,breaks=breaks$a[[s]])
    x$p = window(x$p,breaks=breaks$p[[s]])
    x$c = window(x$c,breaks=breaks$c[[s]])

    #id dependent variable
    y = dat[,sampobj$dv]
    #create model matrix with "effect coding"
    xmat = model.matrix(~+a+p+c,data=x,
                        contrasts.arg=list(
                          a='contr.sum',
                          p='contr.sum',
                          c='contr.sum'
                        ))


    ####
    #draw bayesian posterior
    m = lin_gibbs(y=y,x=xmat,iter=1)
    fits = rbind(fits,m[c('r2','bic','bic_prime','sigma')])

    #caluclate marginalize effects, i.e. calcate marginal means using effects coding
    #this is just the sum of intercept and effect, ommited category is -1*sum(othercats)
    #omitted  categroy for effects coding in R is the last category
    #see Hardy 1993 and summary from APC project for examples
    marginals = lapply(c('a','p','c'), FUN = function(d){
                       e=as.data.frame(m$betas) %>% dplyr::select(starts_with(d))
                       #calculate effect ofr omitted category
                       e[,ncol(e)+1] = -1*rowSums(e)
                       #add grand mean for expected marginal
                       if(means){e = e + m$betas[,1]}
                       colnames(e) = paste0(d,levels(x[,d]))
                       return(e)

                       })#end lapply

    #save intercept value / otherwise it is in the marginal...
    intercept = c(intercept,m$betas[1])

    names(marginals) = c('a','p','c')

    #generate dummy variables to describe full range of dimension
    blockdat = list()
    blockdat$a = scopedummy(w=x$a,unique.vals=unique(dat$a))
    blockdat$p = scopedummy(w=x$p,unique.vals=unique(dat$p))
    blockdat$c = scopedummy(w=x$c,unique.vals=unique(dat$c))


    predat=lapply(c('a','p','c'), FUN=function(d)
       model.matrix(~.-1,data=as.data.frame(blockdat[[d]])) %*%
         t(as.matrix(marginals[[d]]))
       )

    names(predat) = c('a','p','c')

    #draw betas, if requested, across entire length
    if(betas){
      beta = lapply(c('a','p','c'), FUN=function(d){
        b = m$betas[grepl(d,colnames(m$betas))]
        names(b) = paste0(d,levels(x[,d])[1:length(b)])
        return(b)
        })
      beta = c(list(`(Intercept)` = m$betas[1]),beta)
      raw.betas[[length(raw.betas)+1]] = do.call(c,beta)
    }#end conditional for betas

    for(eff in names(predat)){

      effects[[eff]] =
        rbind(effects[[eff]],t(predat[[eff]]))

    }


  }#end of sampling loop for "s"

  for(i in seq_along(effects)){
    colnames(effects[[i]]) = apcvals[[i]]}

  res = list(sampobj=sampobj,
             fit=fits,
             sampled=s.index,
             effects=effects,
             intercept=intercept)

  if(betas){res[['betas']] = raw.betas}

  class(res) = 'apceffects'

  return(res)


}



###########
#summary function for effects class

#' @export
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
  cat('\n\nBayes Factor less than 6 indicates much better fit.',
      'Lowest BIC is usually best fitting.')

}


#################
#summarize apceffects object
summarize.apceffects = function(effectsobj){
  cat('coming soon...')
}

colQuant = function(df,alpha=0.05){
  #returns quantiles

  lb = alpha/2
  ub = 1-lb

  res = apply(df,2,quantile,probs=c(lb,ub))

  return(res)

}

#################
#plot apceffects object

#' @export
plot.apceffects = function(effectsobj,
                           alpha=0.05,
                           adj.se=TRUE){
  #need to add a ci
  #alpha is % level (two-tail)
  #adj.se adjusts ll and ul using mean level window breaks
  require(ggplot2)
  require(dplyr)

  eff = effectsobj[['effects']]

  #calculate mean
  preds = lapply(eff,
          function(x) colMeans(x)
    )

  #calculate intervals
  for(i in seq_along(preds)){
    preds[[i]] = rbind(preds[[i]],
          colQuant(eff[[i]],alpha=alpha))
  }

  preds = lapply(preds,t)
  preds = lapply(names(preds),
   function(nm){
    x = as.data.frame(preds[[nm]])
    colnames(x) = c('Fit','ul','ll')
    x$x = as.numeric(row.names(x))
    x$dim = nm
    return(x)
    }

  )

  preds=do.call(rbind,preds)

  #adj.se calculates an adjusted standard error where each dimension shares a variance,
  #and the pooled "n" is the weighted average of observed window breaks
  if(adj.se){
    sds = unlist(lapply(names(effectsobj$effects),function(x){
      sd(unlist(effectsobj$effects[[x]]))
    }))

    names(sds) = c('a','p','c')
    sds = as.data.frame(sds)

    #sds$n = c(length(unique(tdat$a)),length(unique(tdat$p)),length(unique(tdat$c)))
    sds$n = apply(effectsobj$sampobj$summaries[,c('a','p','c')],2,
                  weighted.mean,w=effectsobj$sampobj$summaries$w)
    sds$se = sds$sds/sqrt(sds$n)

    preds = merge(preds,as.data.frame(sds),by.x='dim',by.y='row.names')
    preds$crit = qt(1-alpha,df=preds$n)

  }

  preds = preds %>%
    arrange(dim) %>%
    mutate(dim_f = factor(dim,labels=c('Age','Cohort','Period')))

  preds$dim_f = factor(preds$dim_f,levels=c('Age','Period','Cohort'))

  #return ggplot object
  plt = ggplot(preds,
               aes(x=x,y=Fit)) +
    geom_line() +
    geom_ribbon(aes(ymax=ul,ymin=ll),alpha=0.25) +
    facet_wrap(~dim_f,scales='free_x') +
    xlab('') +
    theme_classic() +
    theme(axis.text.x = element_text(angle=45,hjust=1))


  if(adj.se == TRUE){
    plt = plt +
      geom_errorbar(aes(ymax=Fit + crit*se,ymin=Fit-crit*se))
  }


  return(plt)

}





