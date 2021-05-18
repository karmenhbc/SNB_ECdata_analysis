# Required packages
library("tidyverse")
library("readxl")
library("rstan")
library("lme4")
library("brms")
library("tidybayes")
library("sjstats")
library("ggthemes")
library("ggmcmc")
library("RColorBrewer")
library("R2WinBUGS")
library("boa")
library("MplusAutomation")
library("psych")
library("bayesplot")
rstan_options(auto_write = TRUE)

# Data
ec.data <- read_xlsx("ECDATA.xlsx",sheet=1,na=c("","NA"))
ec.data <- mutate(ec.data,
                  ENTHESIS.ABR=recode(ENTHESIS,
                                      "Origin long head biceps brachii"="LHBB",
                                      "Origin long head triceps brachii"="LHTB",
                                      "Insertion supraspinatus"="ISS",
                                      "Insertion infraspinatus"="IIS",
                                      "Insertion teres minor"="ITM",
                                      "Insertion subscapularis"="ISK",
                                      "Origin common extensors"="COE",
                                      "Origin common flexors"="COF",
                                      "Insertion triceps brachii"="ITB",
                                      "Insertion braquialis"="IBQ",
                                      "Insertion biceps braquii"="IBB",
                                      "Insertion brachioradialis"="IBR",
                                      "Ischial tuberosity"="COIT",
                                      "Insertion gluteus medius"="IGM",
                                      "Insertion gluteus minimus"="IGN",
                                      "Insertion iliopsoas"="IIO",
                                      "Insertion obturator externus"="IOE",
                                      "Origin gastrocnemius"="OGN",
                                      "Origin Gastrocnemius"="OGN",
                                      "Insertion cuadriceps femoris (patella superior)"="IQF",
                                      "Origin patellar tendon (patella inferior)"="OPT",
                                      "Insertion patellar tendon (tibial tuberosity)"="IPT",
                                      "Insertion triceps surae"="ITS",
                                      "Insertion triceps Surae"="ITS"),
                  ENTHESIS=recode(ENTHESIS,
                                  "Origin long head biceps brachii"="Origin of the long head of the biceps brachii",
                                  "Origin long head triceps brachii"="Origin of the long head of the triceps brachii",
                                  "Insertion supraspinatus"="Insertion of the supraspinatus",
                                  "Insertion infraspinatus"="Insertion of the infraspinatus",
                                  "Insertion teres minor"="Insertion of the teres minor",
                                  "Insertion subscapularis"="Insertion of the subscapularis",
                                  "Origin common extensors"="Common origin of the extensors",
                                  "Origin common flexors"="Common origin of the flexors",
                                  "Insertion triceps brachii"="Insertion of the triceps brachii",
                                  "Insertion braquialis"="Insertion of the brachialis",
                                  "Insertion biceps braquii"="Insertion of the biceps brachii",
                                  "Insertion brachioradialis"="Insertion of the brachioradialis",
                                  "Ischial tuberosity"="Common origin of biceps femoris, semitendinosus and semimembranosus",
                                  "Insertion gluteus medius"="Insertion of the gluteus medius",
                                  "Insertion gluteus minimus"="Insertion of the gluteus minimus",
                                  "Insertion iliopsoas"="Insertion of the iliopsoas",
                                  "Insertion obturator externus"="Insertion of the obturator externus",
                                  "Origin gastrocnemius"="Origin of the medial and lateral heads of the gastrocnemius",
                                  "Origin Gastrocnemius"="Origin of the medial and lateral heads of the gastrocnemius",
                                  "Insertion cuadriceps femoris (patella superior)"="Insertion of quadriceps femoris",
                                  "Origin patellar tendon (patella inferior)"="Origin of the patellar tendon",
                                  "Insertion patellar tendon (tibial tuberosity)"="Insertion of the patellar tendon",
                                  "Insertion triceps surae"="Insertion of the triceps surae",
                                  "Insertion triceps Surae"="Insertion of the triceps surae"))


# Define predictor variables
ec.data <- mutate(ec.data,
                  SEX=factor(SEX),
                  SIDE=factor(SIDE))
# Joint order (for figures and tables)
cat.enthesis <- c("Origin of the long head of the biceps brachii","Origin of the long head of the triceps brachii","Insertion of the supraspinatus","Insertion of the infraspinatus","Insertion of the teres minor","Insertion of the subscapularis","Common origin of the extensors","Common origin of the flexors","Insertion of the triceps brachii","Insertion of the brachialis","Insertion of the biceps brachii","Insertion of the brachioradialis","Common origin of biceps femoris, semitendinosus and semimembranosus","Insertion of the gluteus medius","Insertion of the gluteus minimus","Insertion of the iliopsoas","Insertion of the obturator externus","Origin of the medial and lateral heads of the gastrocnemius","Insertion of quadriceps femoris","Origin of the patellar tendon","Insertion of the patellar tendon","Insertion of the triceps surae")
### features
bf_TC<- bf(TCZ2~SEX+me(AGEmean, AGEsd)+SIDE+(1|p|IND)+(1|q|ENTHESIS.ABR), family=bernoulli)
bf_BFZ1<- bf(BFZ1~SEX+me(AGEmean, AGEsd)+SIDE+(1|p|IND)+(1|q|ENTHESIS.ABR), family=poisson)
bf_ERZ1<- bf(ERZ1~SEX+me(AGEmean, AGEsd)+SIDE+(1|p|IND)+(1|q|ENTHESIS.ABR), family=poisson)
bf_BFZ2<- bf(BFZ2~SEX+me(AGEmean, AGEsd)+SIDE+(1|p|IND)+(1|q|ENTHESIS.ABR), family=poisson)
bf_ERZ2<- bf(ERZ2~SEX+me(AGEmean, AGEsd)+SIDE+(1|p|IND)+(1|q|ENTHESIS.ABR), family=poisson)
bf_FPOZ2<- bf(FPOZ2~SEX+me(AGEmean, AGEsd)+SIDE+(1|p|IND)+(1|q|ENTHESIS.ABR), family=poisson)
bf_MPOZ2<- bf(MPOZ2~SEX+me(AGEmean, AGEsd)+SIDE+(1|p|IND)+(1|q|ENTHESIS.ABR), family=poisson)
bf_CAZ2<- bf(CAZ2~SEX+me(AGEmean, AGEsd)+SIDE+(1|p|IND)+(1|q|ENTHESIS.ABR), family=poisson)
### set prior  ## azul oscuro - prediction interval ## beta.fit3
sexprior <- c(set_prior("normal(0.39, 0.488)", class = "b", coef = "SEXMale", resp = "BFZ1"),
              set_prior("normal(0.39, 0.488)", class = "b", coef = "SEXMale", resp = "ERZ1"),
              set_prior("normal(0.39, 0.488)", class = "b", coef = "SEXMale", resp = "TCZ2"),
              set_prior("normal(0.39, 0.488)", class = "b", coef = "SEXMale", resp = "BFZ2"),
              set_prior("normal(0.39, 0.488)", class = "b", coef = "SEXMale", resp = "ERZ2"),
              set_prior("normal(0.39, 0.488)", class = "b", coef = "SEXMale", resp = "FPOZ2"),
              set_prior("normal(0.39, 0.488)", class = "b", coef = "SEXMale", resp = "MPOZ2"),
              set_prior("normal(0.39, 0.488)", class = "b", coef = "SEXMale", resp = "CAZ2"))
### final brm model 
ec.brm.A<- brm(bf_TC + bf_BFZ1 + bf_ERZ1 + bf_BFZ2 + bf_ERZ2 + bf_FPOZ2 + bf_MPOZ2 + bf_CAZ2,
                data=ec.data, seed  = 123,
                prior = sexprior,
                chains = 4, iter=10000,
                control = list(adapt_delta = 0.98))
summary(ec.brm.A)

