CKMFit: circadian mRNA kinetic model fitting 
====================
This R package is fitting circadian (temproal) profiles of pre-mRNA and mRNA with a kineticc model,
dissecting the contributions of rhythmic transcription and mRNA degradation 
and finally inferring mRNA kineitc parameters, e.g. half-life, rhythmic amplitude and phase of rhythmic degradation.

## Installation
#### Prerequisites
* R 3.4.1 (currently tested by JW), >= R 3.0.0 should work as well (but to check) 

```
DESeq2, limma, emdbook, deSolve, fdrtool, circular, preprocessCore, 
gtools, biomaRt, numDeriv, Matrix, graphics, stats, utils
```
#### Cloning the git repository
```
cd dir_to_place_repository
git clone https://github.com/lengfei5/mRNA_degradation_kinetics_fitting
```

#### Installing the R package
```
install.packages("devtools")
library(devtools)
install_github("lengfei5/CKMFit")
```

## Getting Started
```{r}
rm(list=ls())

## import example and create an object containing data table, geneNames, geneLengths, sizeFactors, dispersion estiamtion and variance estimation 
dataDir = "data/"
load(file = paste0(dataDir, "fitting_degradation_all_data_example_readCount_rpkm.Rdata"))

zt = seq(0,94,by = 2)
ZT.int = grep('.count.premRNA', colnames(T))
ZT.ex = grep('.count.mRNA', colnames(T))
length.int = which(colnames(T) == "length.premRNA")
length.ex = which(colnames(T) == "length.mRNA")

## creat a MDfitDataSet object (a S3 class)
TEST.readCount.NB = FALSE

if(TEST.readCount.NB){
  mds = MDfitDataSet(P = T[, ZT.int], M = T[, ZT.ex], length.P = T[, length.int], length.M = T[, length.ex], zt=zt,
                     mode = "NB", fitType.dispersion = "local")
}else{
  mds = MDfitDataSet(P = T[, ZT.int], M = T[, ZT.ex], zt=zt, mode = "logNormal", fitType.var = "pool")
  
}

## Specify required parameter for the main function and test 
outliers.removal = TRUE
debug = TRUE
identifiablity.analysis.gamma = TRUE
gg = "Per3"
gene.index = which(T$gene==gg)

ptm <- proc.time()
res.fit = make.fits.with.all.models.for.one.gene.remove.outliers(mds, gene.index = gene.index, debug = debug,
                                                                            outliers.removal = outliers.removal,
                                                                            identifiablity.analysis.gamma = identifiablity.analysis.gamma);
cat("------------- time required ---------------\n")
proc.time() - ptm
```

## Improvements
- [x] Since we are also planning to add the option for "Gaussian noise", we should think how to design the fucntions in such way
  that they can be easily to be adapted to do it. 
  And also keep in mind that JW has done it, at least partially, 
  in the inital effort (in the folder`origin/`).   
- [x] The empirical Bayes for the variance in gaussian noise have been implemented in limma pacakge for microarray; however, it can not be  
    directly borrowed, because the limma estimate first the gene-wide variance by fitting a GLM, which is applicable in our case. 
    we need to understand how it works in some detailed steps and to ajust it for our case.
- [x] In the limma package, two main papers were done for the variance estiamtion with EB shrinkage  
    Smyth (2004) and Phipson et al. (2016), the latter addes a robust option to deal with outliers;  
    And "limma-trend" option was also added, which is relevant to our case and   
    integrated the idea from the paper Sartor et al. BMC (2006) which proposed the intensity-based EB method for variance estimation
    
- [x] Gaussian mode is implemented in the parameter optimizaiotn function
- [x] Implement Gaussian mode for outlier detection
- [x] Implement Gaussian mdoe for identifiability analysis; Some inspiration could come from the bbmle package in which a profile-likelihood was   
        implemented. 
        

## TO-Do list
- [ ] Now the code is designed just for fitting one gene. 
  Ideally the code can easily fit all genes in the data in parallel.
  Thus the parallization should be taken into consideration now. 

- [ ] Since the Gaussian need to calculate error function in log scale, probably need to impute the data if there are zeros
  
- [ ] Revise the hessian funciton for SE calculation, because either optim function ofr hessina function can yield NA for some parameters

- [ ] Not sure we should change S3 class to S4 (more strict in the definition and less error-prone in usage)

- [ ] Not sure we should do something similar to limma or DESeq2, wrapping data and funciton in one object; and extracting function will show the resutls
  because I think this will somehow a easy solution for all data and function dependencies. 

- [ ] Parameter cleaning in the last step is not clear how to integrate from the origin code

- [ ] Headers in all scripts should probably removed or modified

- [ ] Connect the general parameter boundaries (modifiable by used) and gene-specific boundaries (refine the boundaries by the gene data) 
