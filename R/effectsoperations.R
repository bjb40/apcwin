####
#these functions return operations on the objects defined
#in samplingfuns

#' @export
delt = function(effectsobj){
  delt = lapply(names(effectsobj$effects),FUN=function(x){
    df = effectsobj$effects[[x]]
    df = df[,order(as.numeric(colnames(df)))]
    res=list()
    cnames=character()
    for(col in 1:(ncol(df)-1)){
      res[[col]] = df[,col+1] - df[,col]
      cnames = c(cnames,
                 paste0(colnames(df)[col+1],'-',colnames(df[col])))
    }

    res = do.call(rbind,res)
    rownames(res) = cnames
    return(res)
  })

  names(delt) = names(effectsobj$effects)
  return(delt)
}

#'Test for statistical differneces in window breaks.
#'
#' @param effectsobj An apceffects object.
#' @return Provides a list for each winow object showing whether there are statistical differences.
#' @export
tstdiff = function(effectsobj){

  delt = delt(effectsobj)
  diff = lapply(names(effectsobj$effects),FUN=function(x){
    df = length(unique(effectsobj$sampobj$data[,x]))
    diff = apply(delt[[x]],1,mean)
    se = (sd(unlist(delt[[x]]))/sqrt(df))
    tval = diff/se

    res = data.frame(
      diff = diff,
      df = df,
      se = se,
      tval = diff/se,
      pval = dt(tval,df=df),
      sig = sapply(dt(tval,df=df),sig)
    )

    return(res)
  })

  names(diff) = names(effectsobj$effects)
  return(diff)
}

#' Returns a dataframe with factor variables including window breaks at a specifie level.
#'
#' @param effectsobj: an apc effects object.
#' @param pval: a probability-value indicating a cut-off for the window breaks
#' @return A dataframe with factor varaibles appended ".win" based on the breaks.
#' @export
modeled_windows=function(effectsobj,pval=0.05){
  #pval identifies sensitivity for breakpoints

  data = effectsobj$sampobj$data

  uvals = list()
  for(d in effectsobj$sampobj$windowvars){
    uvals[[d]] = unique(data[,d])[order(unique(data[,d]))]
  }

  #draw differences from effects object
  dd = tstdiff(effectsobj)
  brk = lapply(dd,function(x){
    res = x$pval<=pval
    #last value is always true for window break
    res = c(res,TRUE)
    return(res)})

  #generate breaks from list
  for(d in effectsobj$sampobj$windowvars){
    winbrks=c(min(data[,d])-1,uvals[[d]][brk[[d]]])
    data[,paste0(d,'.win')] = window(data[,d],breaks=winbrks)
  }


  return(data)

}

