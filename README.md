This repository provides the R routines from the article [Flexible robust mixture regression modeling](https://doi.org/10.57805/revstat.v20i1.365). It is a joint work with Carlos Antonio Abanto-Valle, published in _REVSTAT - Statistical Journal_.

The paper provides a flexible methodology for the class of finite mixture of regressions with a scale mixture of skew-normal errors (SMSN-FMRM), relaxing the constraints previously imposed during the estimation process[^1]. A Bayesian inference procedure is developed based on the data augmentation principle and Markov chain Monte Carlo (MCMC) algorithms. A simulation study is implemented to understand the possible effects caused by the restrictions, and an example with a well-known dataset illustrates the performance of the proposed methods.

The repo includes:

- **MCMC_MixReg_SkewNormal.R** : MCMC routine for the mixture of regressions based on the skew-normal distribution;
- **MCMC_MixReg_SkewT.R** : MCMC routine for the mixture of regressions based on the skew-t distribution;
- **MCMC_MixReg_SkewSlash.R** : MCMC routine for the mixture of regressions based on the skew-slash distribution;
- **SkewSlash_density.R** : Density function for the skew-slash distribution;
- Directory **./Full_Conditionals** contains the full conditional distributions for the above models.

[^1]: Zeller, C.B., Cabral, C.R.B. & Lachos, V.H. [Robust mixture regression modeling based on scale mixtures of skew-normal distributions](https://doi.org/10.1007/s11749-015-0460-4). TEST 25, 375–396 (2016).
