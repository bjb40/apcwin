###
#load and clean example GSS data
library(dplyr)

load(paste0('H:/projects/apc/output/~dat.RData'))

dat = dat %>%
  mutate(happy=3-happy,
         fechld = 4-fechld,
         fefam = fefam-1,
         fepol = fepol-1,
         fepresch = fepresch-1,
         egal = fechld + fepresch + fefam + fepol,
         birthyear = year - age,
         female=sex-1,
         race=factor(race,labels=c('White','Black','Other')))

#limit to complete cases // inlcude educaton only as covariate
dat = dat %>%
  dplyr::select(egal,age,year,birthyear,female,race,educ)

t=nrow(dat)

dat = dat[complete.cases(dat),]

dat = dat %>%
  rename(a=age,
         p=year,
         c=birthyear)


#######
#compare gibs to ml

basem = lm(egal~female+a,data=dat)
mlm = lin_ml(y=dat$egal,
             x=model.matrix(~female+a,data=dat))
mgb = lin_gibbs(y=dat$egal,
                x=model.matrix(~female+a,data=dat))


#######
#testing functions

tt.gibbs=draw_chains(dat,dv='egal',samples=3)
tt.ml=draw_chains(dat,dv='egal',samples=3,method='ml')

#short test
chains=apcsamp(dat,dv='egal',method='ml',samples=3,cores=2)
chains=apcsamp(dat,dv='egal',samples=3)

#big test
chains=apcsamp(dat,dv='egal',samples=100,
               method='ml',cores=3)

###########
#testing internal ata
#load testing data

#need to make the object smaller--breaks & functions for effects

###
#load test data
#Note, you need to name your variables
#as follows a = age, p = period, and c= cohort
data(apcsim)
apcsim$c = apcsim$p-apcsim$a

#this draws 2500 samples on 4 cores for 10000 model samples
testsamp = apcwin::apcsamp(dat=apcsim,
                   dv='y1',
                   method='ml',
                   samples=50,
                   cores=2)

###
#tibble test
library(tidyverse)
tst = as_tibble(apcsim)

tst2 = apcwin::apcsamp(dat=tst,
                       dv='y1',
                       method='ml',
                       samples=10,
                       cores=2)

#this draws a posterior effect sample
#it takes a "sample" object (calculated by apcsamp)
testeff = draw_effs(testsamp,
                    tol=0.001)

#this plots the results of the dimensions
plot(testeff,alpha=0.05)

#ml.draw1 = apcsamp(apcsim,
#                  dv='y1',
#                  method='ml',
#                  samples=1000,
#                  cores=4)

#ml.draw2 = apcsamp(apcsim,
#                   dv='y2',
#                   method='ml',
#                   samples=1000,
#                   cores=4)


#ml.gss = apcsamp(dat,
#                dv='egal',
#                method='ml',
#                samples=2500,
#                cores=4)

