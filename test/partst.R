
#create holder environment for shared variables
par=new.env()

#assign to environemnt to pass to clusters
assign('args',list(
dat=dat,
dv='egal',
apc=c('a','p','c'),
method='ml',
samples=3,
draws=1000
),env=par)



cores=2

cat('\n\nBegin drawing chains...')

cl <- makeCluster(mc <- getOption("cl.cores", cores))
clusterExport(cl=cl,varlist='args',envir=par)
#clusterExport(cl=cl, varlist=c("dat","apc",'method','samples','draws'))
chains=parLapply(cl=cl,1:cores,function(...){
                      require(apcwin)
                      draw_chains(dat=args$dat,
                                  dv=args$dv,
                                  apc=args$apc,
                                  method=args$method,
                                  samples=args$samples,
                                  draws=args$draws)
  })

#end cluster
cat('\n\nEnd drawing chains...')
stopCluster(cl)
