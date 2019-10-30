
# BinGSD

<!-- badges: start -->
<!-- badges: end -->

The goal of BinGSD is to calculate boundaries and conditional power for single-arm group sequential test with binary endpoint. Two kinds of tests are available in this package: asymptotic test based on normal distribution; exact test based on binomial distribution. BinGSD also provides functions to compute boundary crossing probabilities given a specific design.


## Installation

The released version of BinGSD can be downloaded from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("BinGSD")
```

## Architecture 

There are six functions in BinGSD to satisfy the needs of user to design and analyse the single-arm group-sequential trial with binary endpoints:

* asymdesign: compute sample size and boundaries for each analysis, under asymptotic test.
* exactdesign: compute sample size and boundaries for each analysis, under exact test.
* asymprob: given the design, compute boundary crossing probabilities under asymptotic test.
* exactprob: given the design, compute boundary crossing probabilities under exact test.
* asymcp: given the design and interim statistic, compute conditional power under asymptotic test.
* exactcp: given the design and interim statistic, compute conditional power under exact test.

## Example

To design a single-arm group-sequential trial via BinGSD, user should first define some terms:

``` r
library(BinGSD) ##load package
I=c(0.2,0.4,0.6,0.8,1) ##information fractions
beta=0.2   ##desired overall type II error rate
betaspend=c(0.1,0.2,0.3,0.2,0.2) ##proportions of type II error spent at each analysis
alpha=0.05  ##desired overall type I error
p_0=0.3   #response rate under null hypothesis
p_1=0.5   #response rate under alternative hypothesis
K=5     #number of planned analyses
tol=1e-6 #tolerance level

```
Call function `asymdesign` to obtain sample sizes and boundaries based on asymptotic test, under the settings above:

``` r
asymdesign(I,beta,betaspend,alpha,p_0,p_1,K,tol) 
```
The output is an object of class `asymdesign`, including the last and only upper bound, lower boundaries, sample sizes, boundary crossing probabilities, actual overall type I error, power of the test, etc.

For usage of other functions in BinGSD, please refer to the manual or start with the examples at the last section of vignette.
