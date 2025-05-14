# Bayesian Machine Learning approach for modelling gene expression time series

## Overview
Here, we compared two time-series differential expression tools – GPcounts, which utilises a negative-binomial Gaussian Process framework, and maSigPro, which applies polynomial regression – to analyse noisy RNA-seq data from mice under stress conditions. GPcounts is implemented in Python using TensorFlow and GPflow. maSigPro is implemented using R. 
Further details for each package can be found in their respective repositories:
https://github.com/ManchesterBioinference/GPcounts

https://github.com/mjnueda/maSigPro

## Repository Structure
- **GPcounts Report.pdf**: The report covers the two regression methods used for gene series analysis.
- **GPcounts Poster.pdf**: Academic poster for the report.
- **Example_Normalised_Counts.csv**: Count data for 100 representative genes.
- **Example_Sampletable.csv**: Sample table for representative genes.
- **maSigPro.R**: R script implementing maSigPro regression analysis using a polynomial regression model. 
- **GPcounts.py**: Python script implementing GPcounts regression analysis.
- **Comparisons.ipynb**: Jupyter Notebook identifying the model- and sex-specific genes identified using each model. 
- **GO_Enrichment_Heatmpa.R**: Example R script for generating a heatmap of the top 30 GO enriched terms. 
