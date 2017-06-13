# Regression Data Sets Used #

### Load the Data Sets in R ###
For loading all the 18 regression data sets in R you need to issue the following command:

```r
load("DataSets18.RData")
```

This will load into R an object named `DSs` which contains a list with 18 objects of class `dataset` from package DMwR.


### Data Sets Description ###
The main characteristics of the 18 regression data sets included in this folder are as follows:


| Data Set   | N    | tpred | p.nom | p.num | nRare | % Rare |
|------------|------|-------|-------|-------|-------|--------|
| servo      | 167  | 4     | 2     | 2     | 34    | 0.204  |
| a6         | 198  | 11    | 3     | 8     | 33    | 0.167  |
| Abalone    | 4177 | 8     | 1     | 7     | 679   | 0.163  |
| machineCpu | 209  | 6     | 0     | 6     | 34    | 0.163  |
| a3         | 198  | 11    | 3     | 8     | 32    | 0.162  |
| a4         | 198  | 11    | 3     | 8     | 31    | 0.157  |
| a1         | 198  | 11    | 3     | 8     | 28    | 0.141  |
| a7         | 198  | 11    | 3     | 8     | 27    | 0.136  |
| boston     | 506  | 13    | 0     | 13    | 65    | 0.128  |
| a2         | 198  | 11    | 3     | 8     | 22    | 0.111  |
| fuelCons   | 1764 | 38    | 12    | 26    | 164   | 0.093  |
| availPwr   | 1802 | 16    | 7     | 9     | 157   | 0.087  |
| cpuSm      | 8192 | 13    | 0     | 13    | 713   | 0.087  |
| maxTorq    |1802  | 33    | 13    | 20    | 129   | 0.072  |
| bank8FM    | 4499 | 9     | 0     | 9     | 288   | 0.064  |
| ConcrStr   | 1030 | 8     | 0     | 8     | 55    | 0.053  |
| Accel      | 1732 | 15    | 3     | 12    | 89    | 0.051  |
| airfoild   | 1503 | 5     | 0     | 5     | 62    | 0.041  |

where,

N represents the total number of cases;

tpred represents the number of predictors;

p.nom representd the number of nominal predictors;

p.num represents the number of numeric predictors;

nRare represents the number of cases with target variable relevance above 0.8; and

% Rare=nRare/N.

