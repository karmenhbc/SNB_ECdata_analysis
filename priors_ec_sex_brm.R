# Process data
library("tidyverse") 
library("readxl") 
library("dplyr")

# Effects and meta-analysis
library("esc")
library("meta")
library("metafor")

# Distributions
library("fitdistrplus")
library("Rfast")
library("mixtools")
library("metRology")

## Data
corpus <- read_xlsx("HedgesG.xlsx",sheet=1)
## g to z-score and regression coefficients
corpus <- mutate(rowwise(corpus),z=ifelse(is.na(se),NA,es/se),
                 beta=ifelse(is.na(se),NA,convert_d2logit(d=es,se=se,totaln=sample.size)$es),
                 se.beta=ifelse(is.na(se),NA,convert_d2logit(d=es,se=se,totaln=sample.size)$se))

## priors estimate
theta.values=round(seq(-10,10,0.01),2)

# Distribution of the meta-analytic effect size (informative prior)
corpus.meta <- metagen(TE=beta,
                       seTE=se.beta,
                       data = corpus,
                       studlab = paste(reference),
                       comb.fixed = FALSE,
                       comb.random = TRUE,
                       method.tau = "REML",
                       hakn = TRUE,
                       prediction = TRUE,
                       sm = "OR")


# Prior 1: Prediction interval
# Normal distribution, with mean/sd as follows
prior1.mean <- corpus.meta$TE.random
prior1.sd <- corpus.meta$seTE.predict

# Prior 2: Chebyshev's inequality
# Normal distribution, with mean/sd as follows
beta.norm <- normal.mle(as.numeric(na.exclude(corpus$beta)))
prior2.mean <- beta.norm$param[1]
prior2.sd <- sqrt(beta.norm$param[2]*(4.47/1.96))

# Prior 3: Generalised Student's t
# Student't distribution, with mean zero and sd/df as follows
theta.t <- fitdist(data=as.numeric(na.exclude(corpus$z)),distr="t.scaled",fix.arg=list("mean"=0),start=list(df=4,sd=1),lower=c(3,0),upper=c(7,Inf))
prior3.mean <- 0
prior3.sd <- theta.t$estimate["sd"]
prior3.df <- theta.t$estimate["df"]
