---
author: Bryce Bartlett
title: R package for Bayesian Model Averaging of the Age-Period-Cohort identification problem.
---

# Introduction

There is a long-standing problem in demography known as the age-period-cohort identification problem. There is a fairly straightforward issue: some social characteristic, say, happiness, may be correlated with age. People generally get happier as they get older. This is the age effect. Happiness can also be associated with current events: everyone is less happy during a recession. This is a period effect. Finally, different birth cohorts (sometimes called generations) are often correlated with different characteristics. This is a cohort effect.

Unfortunately, age, period, and cohort measures are perfectly correlated. If you know any two, you also know the third. If you know how old someone is (age) and when they were born (cohort), then you know what year (period) it is ($age+cohort=period$). Because these measures are perfectly correlated, they are not identified, and inestimable. 

There are a number of proposals to address this problem, but this package deploys a novel Bayesian Model Averaging approach. (The manuscript is still in preparation, but will be linked here once complete). This works by estimating age, period, and cohort effects using piecewise constant functions. This is one of the oldest proposals. The problem is that while this breaks the statistical identification problem, solutions are sensitive to the "window" breaks.

*What does it do:* Estimate hundreds of models. Combine the using Bayesian statistics. Simulation studies show good promise that this approach makes the model less sensitive.

# Installation

This package is still in development, but you can install and use it from github using the R library devtools. Here is the code-block:

```
install.packages('devtools')
library(devtools)
install_github('bjb40/apcwin')
```

# Details


# Example



# Algorithm



