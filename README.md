Y. Bai, A. Carriero, T. E. Clark, and M. Marcellino. (2022). Macroeconomic forecasting in a multi-country context. Journal of Applied Econometrics, *forthcoming*.

This file details the data and part of the computer programs used to produce the results in the published paper.

Regarding replication and data, please note the following:

(1)	 We are not able to distribute the short-term interest rate data under the terms of our license with GFD. The interest rate data in “data.xlsx” are obtained by adding a random number to the original series, generated by function RAND in Excel. However, we have provided the GFD tickers (3rd column) in the first line of each sheet. Researchers with access to GFD can download the series with these tickers.

(2)	 Compared to the original file we have used. the main “MCVAR.m” file is restructured to facilitate replication. However, estimation and forecasting computations require significant CPU time. Most of the paper’s results are obtained over Bocconi BIDSA server (12/88 cores used) with parallel computing toolbox. 

1.	Data folder files:
data.xlsx: Quarterly time series data of real GDP growth, inflation, and a short-term interest rate (different from the original series, as explained before). Except GDP growth, the first line of each sheet includes tickers of the data sources.

datapreparation.m: The main file used to construct the dataset. dataQ.mat is the final transformed dataset we use for estimation and forecasting computations.

Other m files: These files are used for data transformation and cleaning, which are needed to run datapreparation.m.


2.	Function folder files
Programs are written in Matlab, and models are estimated with Bayesian methods (Gibbs samplers).

There are three subfolders:
MCVAR: M files used to run Gibbs samplers for multi-country VAR-SV models with various prior specifications. We have also included some files to produce results in the robustness section (Prior grouping of coefficients, Stochastic volatility, and Alternative hierarchical shrinkage in country-specific VAR-SV).

Prediction: M files used to run predictive simulations and compute CRPS for density forecast evaluation.

Util: Auxiliary files used for estimation and forecasting computations.

3.	Main file: 
MCVAR.m: The main file to reproduce forecasting results. There are five options available: Minnesota prior, Horseshoe prior, Normal-Gamma prior, Normal-Gamma-Gamma prior and SSSS prior. You can change the relevant lines to reproduce the results in the robustness check sections. Please note that since for SSSS prior we cannot use the triangular algorithm, estimation is not feasible if lag length is larger than 1. We set L=1.

Contact information:  yu.bai [AT] unibocconi.it
