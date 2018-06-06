data(apcsim)
apcsim$c = apcsim$p-apcsim$a

apcsim = apcsim %>%
  rename(Age=a,
         Period=p,
         Cohort=c)

#this draws 2500 samples on 4 cores for 10000 model samples
testsamp = apcsamp(y1~Age+I(Age^2)+Period+Cohort,
                   windowvars=c('Period','Cohort'),
                   data=apcsim,
                   method='ml',
                   samples=100,
                   cores=4)

#this summarizes the results of the apc sampler
summary(testsamp)

testeff = draw_effs(testsamp,
                    tol=0.01)

#plot the window breaks, relative to the grand mean
plot(testeff)

#test equality of contiguous values
deltas = tstdiff(testeff)
print(deltas[[1]])

#return a new data-frame with implied winow constraints
#at the specified sensitivity (measured as a p-value of differneces)
newdata = modeled_windows(testeff,pval=0.1)

#new model
summary(lm(y1~Age + I(Age^2) +Period.win+Cohort.win,data=newdata))
