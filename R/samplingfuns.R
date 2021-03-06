#funcitons for running window sampler algorithm
#
#Bryce Bartlett

#' @export
window_sample=function(var,alph){
  #input is continuous of a,p,c
  #alpha is a vector the length of unique entries in var that is fed to dirichelet
  #output is factor with uniform, random window constraints
  #see dirichelet, e.g. for alternative algorithms
  vals=unique(var) #center on 0; assume continuous
  len=length(vals)

  alph=unlist(alph)

  dp=MCMCpack::rdirichlet(1,alph)
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
#' @param formula A formula object.
#' @param data A dataframe.
#' @param windowvars A character vector naming the age, period, and cohort variables (defaults to 'a', 'p', and 'c'). This can also be used to calculate breaks for other sorts of variables.
#' @param cores Integer indicating cores to use.
#' @param method One of "gibbs" or "ml". For a fully bayesian estimator for a faster approximation.
#' @param samples Integer indicating samples; used only for "gibbs" method.
#' @param draws Integer indicating the number of models to sample in the chain.
#' @param startvals Optional list for chain starting values.
#' @return An object of apcsamp.
#' @examples
#' data(apcsamp)
#' draw_chains(apcsamp,method='ml',samples=250,cores=4,chains=4)
#' @export
draw_chains = function(formula,
                       data,
                       windowvars=c('a','p','c'),
                       cores=1,
                       method='ml',
                       chains=1,
                       samples=10,
                       draws=1000,
                       startvals=NULL){

  #set up data
  frame = model.frame(formula,data)
  y = model.response(frame)

  tm=Sys.time()
  avtm=0

  #set of numbers of random samples
  n.samples=samples

  #holder df for model summary data
  win = data.frame(matrix(vector(), 0, length(windowvars),
                  dimnames=list(c(), windowvars)),
                  stringsAsFactors=FALSE)
  #generate holder data
  win[1,] = NA

  modsum = data.frame( r2=as.numeric(NA),
                       sigma=as.numeric(NA),
                       bic=as.numeric(NA),
                       aic=as.numeric(NA),
                       bic_prime=as.numeric(NA))

  breaks=lapply(as.list(win),list)

  #calculate length of variables
  if(length(windowvars)==1){
    windat = as.data.frame(data[,windowvars])
    colnames(windat) = windowvars
  } else {
    windat = data[,windowvars]
  }

  dl = apply(windat,2,function(x)
    length(unique(x)))

  #generate dataframe to hold hyperparameters
  all.alphas = lapply(names(dl),function(x)
    as.data.frame(
      matrix(rep(1,dl[x]),
             1,dl[x],
             dimnames=list(c(),paste0(x,'.',1:dl[x])))
    ))
  names(all.alphas) = names(dl)

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

  #reset model frame
  frame = model.frame(formula,data)

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

    for(d in seq_along(all.alphas)){
      all.alphas[[d]][s,]=all.alphas[[d]][s-1,]
    }

  }

  #draw random window samples across windowvars
  for(d in windowvars){
   frame[,d] <- window_sample(frame[,d],
                              all.alphas[[d]][s,])
  }

  x = model.matrix(formula,frame)

  #check whether (t(x) %*% x )' is singular
  identified = class(try(solve(t(x)%*%x),silent=TRUE))=='matrix'

  if(!identified){
    bound=bound+1

    for(d in seq_along(all.alphas)){
      all.alphas[[d]][s,]=all.alphas[[d]][s-1,]
    }
  }

  #collect model data
  if(length(windowvars)==1){
    winframe = as.data.frame(frame[,windowvars])
    colnames(winframe) = windowvars
  } else {
    winframe = frame[,windowvars]
  }

  br = lapply(winframe,function(i)
    attr(i,'breaks'))

  for(i in windowvars){
    breaks[[i]][[s]]=br[[i]]
  }

  #number of factors is the number of breaks - 1 as described
  #in paper describing model-space "G"
  win = rbind(win,as.data.frame(lapply(br,length))-1)

  #add esitmate
  mnum = s

  m= NULL
  if(method=='gibbs'){
    m = lin_gibbs(y=y,x=x,iter=draws)
    #generate means from draws
    m$sigma = mean(m$sigma);m$r2=mean(m$r2)
  } else if(method=='ml'){
    m = lin_ml(y=y,x=x)
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


#' A wrapper for "draw_chains" that will use multiple search chains in the algorithm.
#'
#' @param formula A formula for a model.
#' @param data A dataframe.
#' @param windowvars A character vector naming the age, period, and cohort variables (defaults to 'a', 'p', and 'c'). This can also be used to calculate window breaks for other sorts of continuous varaibles.
#' @param cores Integer indicating cores to use.
#' @param method One of "gibbs" or "ml". For a fully bayesian estimator for a faster approximation.
#' @param draws Integer indicating samples; used only for "gibbs" method.
#' @param samples Integer indicating the number of models to sample in the chain.
#' @param startvals Optional list for chain starting values.
#'
#' @details what does this look like.
#'
#' @return An object of apcsamp.
#'
#' @examples
#' data(apcsim)
#' apcsamp(y1~a + I(a^2) + p + c, data=apcsim,
#'   windowvars=c('p','c'),method='ml',
#'   samples=250,cores=4,chains=4)
#'
#' @export
apcsamp = function(formula,data,windowvars=c('a','p','c'),
                   cores=1,method='ml',
                   chains=1,samples=100,draws=1000,
                   ...){

  #make sure data is in proper form, an listwise delete
  obs = nrow(data)
  data = model.frame(formula,data)
  cat('Deleting',obs-nrow(data),'incomplete observations.\n')

  times = Sys.time()


  #create holder environment for shared variables
  par=new.env()

  #assign to environemnt to pass to clusters
  #assign to environemnt to pass to clusters
  assign('args',list(
    data=data,
    formula=formula,
    windowvars=windowvars,
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
    draw_chains(data=args$data,
                formula=args$formula,
                windowvars=args$windowvars,
                method=args$method,
                samples=args$samples,
                draws=args$draws)
  })

  #end cluster
  rm(par)
  cat('\n\nEnd drawing chains...')
  parallel::stopCluster(cl)

  #need to fix when windowvars is 1 dimensional
  modsum = extract(chains,'modsum', as.df=TRUE)

  #calculate weight by bic and bic_prime
  #https://stats.stackexchange.com/questions/249888/use-bic-or-aic-as-approximation-for-bayesian-model-averaging
  k=min(modsum$aic)
  d=-.5*(modsum$aic-k)
  modsum$w=exp(d)/sum(exp(d))

  k = min(modsum$bic_prime)
  d=-.5*(modsum$bic_prime-k)
  modsum$w_prime=exp(d)/sum(exp(d))

  t.samples = sum(extract(chains,'n.samples'))

  #need to fix breaks and alphas
  breaks = list()
  alphas=list()
  for(i in windowvars){
    breaks[[i]] = do.call(c,lapply(chains,function(j)
      j[['breaks']][[i]]))
    alphas[[i]] = do.call(rbind,lapply(chains,function(j)
      j[['alphas']][[i]]))
  }

  fullt = Sys.time() - times

  res = list(
  data = data,
  windowvars = windowvars,
  formula = formula,
  alphas = alphas,
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

  if(length(windowvars)==1){
    colnames(res$summaries)[1] = windowvars
  }
  class(res) = append(class(res),'apcsamp')

  return(res)

}

#' @export
convergence = function(samp) {
  UseMethod('convergence',samp)
}

#' A wrapper for "draw_effs" that will use multiple cores to speed up sampling from the posterior.
#'
#' @param sampobj A sample object.
#' @return An object of apceffects.
#' @export
postsamp = function(sampobj,tol=NULL,cores=NULL,...){

  #get total draws from posterior and divide by cores
  if(is.null(tol)){tol=1/sampobj$n.samples}
  if(is.null(cores)){cores=sampobj$chains}
  draws = ceiling(1/tol/cores) #round up

  #create holder environment for parellel chains
  par = new.env()

  #assign environment data to pass arguments
  assign('args',
         list(sampobj=sampobj,tol=1/draws,...),
         env=par)

  cl <- parallel::makeCluster(mc <- getOption("cl.cores", cores))
  parallel::clusterExport(cl=cl,varlist='args',envir=par)

  chains=parallel::parLapply(cl=cl,1:cores,function(...){
    args=get('args',env=par)
    require(apcwin)
    require(dplyr)
    draw_effs(sampobj=args$sampobj,
              tol = args$tol)
  })

  #end cluster
  parallel::stopCluster(cl)
  rm(par)

  #unzip chains and rebuild result

  res = list(
    sampobj = sampobj,
    fit = extract(chains,'fit',as.df=TRUE),
    sampled = extract(chains,'sampled'),
    intercept = extract(chains,'intercept')
  )

  for(i in sampobj$windowvars){
    res[['effects']][[i]] = do.call(rbind,
                                    lapply(chains,function(j)
                                      j[['effects']][[i]]))
  }

  class(res) = 'apceffects'

  return(res)

}


#' A wrapper that samples from model space and draws from posterior.
#'
#' @param formula a formula.
#' @return An object of apceffects.
#' @export
#swslm = function(formula,data,windowvars=c('a','p','c'),method='ml',cores=1,chains=1,samples=1000,...){
swslm = function(formula,...){

  args = list(formula=formula,...)

  #draw samples
  sampobj <- do.call(apcsamp,args)

  #draw from posterior
  args[['sampobj']] = sampobj
  effobj <- postsamp(sampobj)

  #define and return implied data breaks


  #return implied model object

#  effobj[['imp_breaks']]
#  effobj[['imp_model']]

#  return(args)
  return(effobj)


}

#' @export
convergence.default = function(samp) {
  cat('\nError: Must use apcsamp object.\n\n')
}

###helper function for convergence
#' @export
convergence.apcsamp = function(samp){
  require(coda)

  e = samp$summaries$aic

  chainlen = length(e)/samp$chains
  aic = split(e, ceiling(seq_along(e)/chainlen))
  aic = lapply(aic,coda::mcmc)
  aic = coda::as.mcmc.list(aic)

  geweke = geweke.diag(aic)
  pv = data.frame(
    zval=unlist(lapply(geweke,function(z) z$z)),
    pv=unlist(lapply(geweke,function(z)
      (1-pnorm(abs(z$z))*.5)))
)
  row.names(pv) = paste('Chain',1:nrow(pv))
  res = list(
    geweke=pv,
    rhat = gelman.diag(aic)$psrf[1],
    size = effectiveSize(aic)
  )

  return(res)

}


###helper function for chain extraction
extract = function(l,name,as.df=FALSE,span=NULL){
  #extracts and combines name object from list 'l'
  #span limits the extraction, if, for example, I need to drop first off
  if(is.null(span) & !as.df){span=1:length(l[[1]][[name]])}

  res=do.call(base::c,lapply(l,function(x) x[[name]][span]))
  if(as.df){
    res=do.call(rbind,lapply(l,function(x) as.data.frame(x[[name]])))
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
#' sampleobject = apcsamp(y1~a+p+c,data=apcsim,method='ml',
#' samples=250,cores=4,chains=4)
#'
#' effectsobject = draw_effs(sampleobject,tol=0.01)
#'
#' @export
draw_effs = function(sampobj,
                     means=FALSE,
                     betas=FALSE,
                     tol=NULL,
                     ...){

  #require(dplyr)

  dat = sampobj$data

  #generate holder for effects dataframe (target of function)
  apcvals = lapply(sampobj$windowvars,
                   function(x) unique(dat[,x]))
  names(apcvals) = sampobj$windowvars

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

    #reset dataframe
    y=model.response(dat)
    data = sampobj$data

    #set window breaks
    for(d in sampobj$windowvars){
      data[,d] = window(data[,d],breaks=sampobj$breaks[[d]][[s]])
    }

    contrasts = as.list(rep('contr.sum',length(sampobj$windowvars)))
    names(contrasts) = sampobj$windowvars

    #create model matrix with "effect coding"
    x = model.matrix(sampobj$formula,data=data,
                        contrasts.arg=contrasts)


    ####
    #draw bayesian posterior
    m = lin_gibbs(y=y,x=x,iter=1)
    fits = rbind(fits,m[c('r2','bic','bic_prime','sigma')])

    #caluclate marginalize effects, i.e. calcate marginal means using effects coding
    #this is just the sum of intercept and effect, ommited category is -1*sum(othercats)
    #omitted  categroy for effects coding in R is the last category
    #see Hardy 1993 and summary from APC project for examples
    marginals = lapply(sampobj$windowvars, FUN = function(d){
                       e=as.data.frame(m$betas) %>% dplyr::select(starts_with(d))
                       #calculate effect ofr omitted category
                       e[,ncol(e)+1] = -1*rowSums(e)
                       #add grand mean for expected marginal
                       if(means){e = e + m$betas[,1]}
                       colnames(e) = paste0(d,levels(data[,d]))
                       return(e)

                       })#end lapply

    #save intercept value / otherwise it is in the marginal...
    intercept = c(intercept,m$betas[1])

    names(marginals) = sampobj$windowvars

    #!here
    #generate dummy variables to describe full range of dimension
    blockdat = list()
    for(d in sampobj$windowvars){
      blockdat[[d]] = scopedummy(w=data[,d],
                                 unique.vals=unique(sampobj$data[,d]))
    }


    predat=lapply(sampobj$windowvars, FUN=function(d)
       model.matrix(~.-1,data=as.data.frame(blockdat[[d]])) %*%
         t(as.matrix(marginals[[d]]))
       )

    names(predat) = sampobj$windowvars

    #draw betas, if requested, across entire length
    if(betas){
      beta = lapply(sampobj$windowvars, FUN=function(d){
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

#helper function to build weights from bic or aic
calcwt = function(vector){
  #vector is of aic or bic
  k=min(vector)
  d=-.5*(vector-k)
  return(exp(d)/sum(exp(d)))
}


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


  cat('\n\nModel fit details.\n')
  cat('BIC:\t',weighted.mean(ss$bic,w=ss$w))
  cat('\nAIC:\t',weighted.mean(ss$aic,w=ss$w))
  cat('\nR-squared:\t',weighted.mean(ss$r2,w=ss$w))

  cat('\n\n\nMean Number Window Breaks:\n')
  print(apply(ss[,samp$windowvars],2,weighted.mean,w=ss$w))


  conv = convergence(samp)
  gw = paste(apply(conv$geweke,2,rnd)[,'pv'],collapse=',')

  cat('\n\nConvergence Diagnostics for sampling step:\n')
  cat('\n\tGewecke Diagnostic (Bayesian p-value for convergence--0.5 is best; tail probabilities indicate nonconvergence):\t',gw)
  cat('\n\tR-Hat (a value of 1 indicates convergence):\t',rnd(conv$rhat))

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
  #require(ggplot2)
  #require(dplyr)

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

    names(sds) = effectsobj$sampobj$windowvars
    sds = as.data.frame(sds)

    #sds$n = c(length(unique(tdat$a)),length(unique(tdat$p)),length(unique(tdat$c)))
    sds$n = apply(effectsobj$sampobj$summaries[,effectsobj$sampobj$windowvars],2,
                  weighted.mean,w=effectsobj$sampobj$summaries$w)
    sds$se = sds$sds/sqrt(sds$n)

    preds = merge(preds,as.data.frame(sds),by.x='dim',by.y='row.names')
    preds$crit = qt(1-alpha,df=preds$n)

  }

  preds = preds %>%
    arrange(dim) %>%
    mutate(dim_f = factor(dim,labels=effectsobj$sampobj$windowvars))

  preds$dim_f = factor(preds$dim_f,levels=effectsobj$sampobj$windowvars)

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





