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

data(apcsim)

apcsim$c = apcsim$p-apcsim$a
ml.draw = apcsamp(apcsim,
                  dv='y1',
                  method='ml',
                  samples=3,
                  cores=3)

