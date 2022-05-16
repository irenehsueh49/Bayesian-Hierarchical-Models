---
title: "Irene Hsueh's BS 849 Homework 4"
author: "Irene Hsueh"
date: "2/21/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(rjags)
library(coda)
library(formatR)
set.seed(1234)
```

# Randomized Controlled Trial Data
```{r}
rct <- read.csv("C:/Irene Hsueh's Documents/MS Applied Biostatistics/BS 849 - Bayesian Modeling for Biomedical Research & Public Health/Class 4 - Hierarchical Models/Homework/rct_data.csv")
N <- as.matrix(rct[c("n_trt", "n_placebo")])
Y <- as.matrix(rct[c("cases_trt", "cases_placebo")])
hospital_data <- list(N=N, Y=Y, trt=c(1,0))
```


# Fixed Effects Model
```{r}
#alpha and beta are fixed effects for all hospitals
fixed_effects_model_bugs <- 
"model {
for (i in 1:22) {
for (j in 1:2) {
Y[i,j] ~ dbin(p[i,j], N[i,j])
logit(p[i,j]) <- alpha + beta*trt[j]
}
}
OR <- exp(beta)

#Prior Distribution 
alpha ~ dnorm(0, 0.001)
beta ~ dnorm(0, 0.001)
}"

fixed_effects_model <- jags.model(textConnection(fixed_effects_model_bugs), data=hospital_data, n.adapt=1000)
fixed_effects_model_gibbs <- update(fixed_effects_model, n.iter=1000)
fixed_effects_model_test <- coda.samples(fixed_effects_model, c("OR", "alpha", "beta"), n.iter=10000)
summary(fixed_effects_model_test)
```


# Hierarchical Effects and Intercept Model
```{r}
#alpha and beta are fixed effects for all hospitals
hierarchical_model_bugs <- 
"model {
for (i in 1:22) {
for (j in 1:2) {
Y[i,j] ~ dbin(p[i,j], N[i,j])
logit(p[i,j]) <- alpha[i] + beta[i]*trt[j]
}
alpha[i] ~ dnorm(alpha0, tau_alpha)
beta[i] ~ dnorm(beta0, tau_beta)
OR[i] <- exp(beta[i])
}

rank_OR <- rank(OR[])

beta_predict ~ dnorm(beta0, tau_beta)
OR_predict <- exp(beta_predict)

hospital1_beta_predict ~ dnorm(beta[1], tau_beta)
hospital1_OR_predict <- exp(hospital1_beta_predict)

#Prior Distribution 
alpha0 ~ dnorm(0, 0.001)
beta0 ~ dnorm(0, 0.001)
tau_alpha ~ dgamma(1, 1)
tau_beta ~ dgamma(1, 1)
}"

hierarchical_model <- jags.model(textConnection(hierarchical_model_bugs), data=hospital_data, n.adapt=1000)
hierarchical_model_gibbs <- update(hierarchical_model, n.iter=1000)
hierarchical_model_test <- coda.samples(hierarchical_model, c("OR", "OR_predict", "hospital1_OR_predict"), n.iter=10000)
summary(hierarchical_model_test)
```



# Independent Model
```{r}
#alpha and beta are fixed effects for all hospitals
independent_model_bugs <- 
"model {
for (i in 1:22) {
  for (j in 1:2) {
  Y[i,j] ~ dbin(p[i,j], N[i,j])
  logit(p[i,j]) <- alpha[i] + beta[i]*trt[j]
  }
OR[i] <- exp(beta[i])

#Prior Distribution 
alpha[i] ~ dnorm(0, 0.001)
beta[i] ~ dnorm(0, 0.001)
}

rank_OR <- rank(OR[])
}"

independent_model <- jags.model(textConnection(independent_model_bugs), data=hospital_data, n.adapt=1000)
independent_model_gibbs <- update(independent_model, n.iter=1000)
independent_model_test <- coda.samples(independent_model, c("OR"), n.iter=10000)
summary(independent_model_test)
```





