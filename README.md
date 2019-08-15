# copula-dyn-pred
R code for "A Gaussian copula approach for dynamic prediction of survival with a longitudinal biomarker" 

This repository contains R code for the simulation and real data analysis. The version of R used to perform the analyses was 3.4.3.

The R script "SimulationCode.R" contains the code for generating data under Simulation Scenario 1a and fitting the joint, landmark, and copula models (including the code for maximum likelihood estimation under the Gaussian and Student t copula). It requires loading "PredictionFunctions.R", which contains the functions to obtain the dynamic predictions for the joint, landmark, and copula models.

For the real data analysis, the R data set heart.valve was used from the joineR package. 
