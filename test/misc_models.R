
##
#sim data

#read gss
library(tidyverse)
library(haven)

#previously cleaned data for gender egalitarian project
dat.f = "H:/projects/gender_egal/output/private~/subdat.csv"

if(file.exists(dat.f)){
  dat = read.csv(dat.f)
  #save object for runing on the server
  save(dat,file=paste0('~dat.RData')) #~ keeps it from pushing to git
} else{
  load(paste0(outdir,'~dat.RData'))
}

dat = dat %>%
  rename(a=age,
         p=year) %>%
  mutate(c = p-a,
         y=sin(p)+(c-1950)*.01+rnorm(nrow(dat)))

noc.m = lm(y~factor(a)+factor(p),data=dat)
noa.m = lm(y~factor(p)+factor(c),data=dat)
true.m = lm(y~I(sin(p))+I(c-1950),data=dat)

###
#substitution model
#substitutions
"
https://brownmath.com/twt/sumdiff.htm#sincosAplusmnB
sin(180)
sin(120)*cos(60)+cos(120)*sin(60)

p = c+a
sin(p) = sin(c+a) =
sin(c)cos(a)+cos(c)sin(a)
"
sub.m = lm(y~I(sin(c)*cos(a)+cos(c)*sin(a))+I(c-1950),data=dat)
sub2.m = lm(y~I(sin(p))+I(p-a-1950),data=dat)
sub3.m = lm(y~I(sin(p))+p+a,data=dat)
sub4.m = lm(y~I(sin(p))+c+I(a+a),data=dat)
sub5.m = lm(y~I(sin(c)*cos(a))+I(cos(c)*sin(a))+I(c-1950),data=dat)
stargazer(list(sub4.m,sub5.m),type='text')

print(rbind(coef(true.m),coef(sub5.m)))

noc=AIC(noc.m)
noa=AIC(noa.m)
true = AIC(true.m)
sub = AIC(sub.m)

vars = dat %>%
  mutate(sinp = sin(p),
         clim = c-1950) %>%
  select(sinp,clim,y,a,c,p)

round(cor(vars,use='pairwise.complete.obs'),2)

library(ggplot2)

ages = dat %>%
  group_by(a) %>%
  summarize(p=mean(p,na.rm=TRUE),
            c=mean(c,na.rm=TRUE))

newdat = data.frame(a=18:89)
newdat = merge(newdat,ages,by='a')

pred = data.frame(predict(true.m,
               newdata=newdat %>%
                 mutate(p=mean(dat$p,na.rm=TRUE),
                        c=mean(dat$c,na.rm=TRUE)),
               se.fit=TRUE))

ggplot(cbind(pred,newdat),
       aes(x=a,
           y=fit,
           ymax=fit+1.96*se.fit,
           ymin=fit-1.96*se.fit)) +
  geom_line() +
  geom_ribbon(alpha=0.3) +
  theme_classic()

library(apcwin)
samp = apcsamp(dat = dat,
        cores=3,
        method='ml',
        samples=100)

print(summary(samp))

effs = draw_effs(samp,tol=0.01)

plot(effs)
View(tstdiff(effs)$a)

quantile(tstdiff(effs)$c$diff,probs=0.5)
mean(tstdiff(effs)$c$diff)

quantile(tstdiff(effs)$a$diff,probs=0.5)
mean(tstdiff(effs)$a$diff)

#partial derivatives
#https://cran.r-project.org/web/packages/margins/vignettes/TechnicalDetails.pdf
"
a=p-c
p=c+a
c=p-a
"
