######
#model objects tests

rm(list=ls())

###
#load test data
data(apcsim)

###
#run linear model on p and c

rmod = lm(y3~a+p,data=apcsim)

###
#print log likelihood, sigma and betas

ll = logLik(rmod)

calcll = sum(dnorm(apcsim$y3,mean=predict(rmod),
                   sd=summary(rmod)$sigma,log=TRUE))

lmobj = lin_ml(y=apcsim$y3,
               x=model.matrix(~a+p,data=apcsim))

gbobj = lin_gibbs(y=apcsim$y3,
                  x=model.matrix(~a+p,data=apcsim))

print(round(coef(rmod),2))
print(round(apply(gbobj$betas,2,mean),2))
print(round(lmobj$betas,2))

#r2 is off in lin_mod because delete intercept
print(c(summary(rmod)$r.squared,lmobj$r2,mean(gbobj$r2)))

print(c(ll,lmobj$ll,mean(gbobj$ll)))
print(c(summary(rmod)$sigma,lmobj$sigma,mean(gbobj$sigma)))

#maybe param+1 inf formula?? (because s^2? takes a df)
print(c(BIC(rmod),lmobj$bic,gbobj$bic))

#off b/c of the r2 -- al
print(c(lmobj$bic_prime,gbobj$bic_prime))
