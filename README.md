In this repository, one can find the R routines used in the article [Flexible robust mixture regression modeling](https://doi.org/10.57805/revstat.v20i1.365). It is a joint work with Carlos Antonio Abanto-Valle, published in _REVSTAT - Statistical Journal_.

The paper provides a flexible methodology for the class of finite mixture of regressions with scale mixture of skew-normal errors (SMSN-FMRM), relaxing the constraints previously imposed during the estimation process^[[1]](https://doi.org/10.1007/s11749-015-0460-4). A Bayesian inference procedure is developed based on the data augmentation principle and Markov chain Monte Carlo (MCMC) algorithms. A simulation study is implemented to understand the possible effects caused by the restrictions, and an example with a well-known dataset illustrates the performance of the proposed methods.

This repo includes:

- **MCMC_MixReg_SkewNormal.R** : MCMC routine for the mixture of regressions based on the skew-normal distribution;
- **MCMC_MixReg_SkewT.R** : MCMC routine for the mixture of regressions based on the skew-t distribution;
- **MCMC_MixReg_SkewSlash.R** : MCMC routine for the mixture of regressions based on the skew-slash distribution;
- Directory **./Full_Conditionals** contains the full conditional distributions for the above models.
