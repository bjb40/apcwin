
data(apcsim)
apcsim$c = apcsim$p-apcsim$a

#this draws 100 samples on 4 cores for 400 model samples
testsamp = apcsamp(y1~a+p+c,
                   windowvars=c('a','p','c'),
                   data=apcsim,
                   method='ml',
                   samples=100,
                   cores=4)

#this summarizes the results of the apc sampler
summary(testsamp)

#this estimates effects from the sampled models
testeff = draw_effs(testsamp,
                    tol=0.01)

#this plots the window breaks, relative to the grand mean
plot(testeff)

#this tests the equality of contiguous values
deltas = tstdiff(testeff)
print(deltas[[1]])

#this returns a new data-frame with implied winow constraints
#at the specified sensitivity (measured as a p-value of differneces)
newdata = modeled_windows(testeff,pval=0.1)

