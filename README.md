## Exploring Resampling with Neighborhood Bias on Imbalanced Regression Problems - EPIA 2017

This repository has all the code used in the experiments carried out in the paper *"Exploring Resampling with Neighborhood Bias on Imbalanced Regression Problems"* [1].


This repository is organized as follows:

* **R_Code** folder - contains all the code for reproducing the experiments presented in the paper;
* **Data** folder - contains the 18 regression data sets used in the experiments carried out;
* **Figures** folder - contains all the extra figures obtained from the experimental evaluation carried out on 18 real world data sets;
* **Tables** folder - contains all the tables obtained from the experimental evaluation carried out on 18 real world data sets;


### Requirements

The experimental design was implemented in R language. Both code and data are in a format suitable for R environment.

In order to replicate these experiments you will need a working installation
  of R. Check [https://www.r-project.org/] if you need to download and install it.

In your R installation you also need to install the following additional R packages:

  - e1071
  - randomForest
  - earth
  - performanceEstimation
  - UBL
  - uba

  All the above packages with exception of uba, can be installed from CRAN Repository directly as any "normal" R package. Essentially you need to issue the following command within R:

```r
install.packages(c("e1071"", "randomForest", "earth", "performanceEstimation", "UBL"))
```

Additionally, you will need to install uba package from a tar.gz file that you can download from [http://www.dcc.fc.up.pt/~rpribeiro/uba/]. 

For installing this package issue the following command within R:
```r
install.packages("uba_0.7.7.tar.gz",repos=NULL,dependencies=T)
```

To replicate the figures in this repository you will also need to install the package:

  - ggplot2

As with any R package, we only need to issue the following command:

```r
install.packages("ggplot2")
```

Check the other README files in each folder to see more detailed instructions on how to run the experiments.

*****

### References
[1] Paula Branco, Luís Torgo and Rita P. Ribeiro (2017). *"Exploring Resampling with Neighborhood Bias on Imbalanced Regression Problems"* 18th Portuguese Conference on Artificial Intelligence, EPIA 2017 (to appear).
