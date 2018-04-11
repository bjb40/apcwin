####
#collecting effects operations

delt = function(effectsobj){
  delt = lapply(names(effectsobj$effects),FUN=function(x){
    df = effectsobj$effects[[x]]
    df = df[,order(colnames(df))]
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
