import numpy as np
import pandas as pd
import gpflow
from IPython.display import display
import tensorflow as tf
from GPcounts.RNA_seq_GP import rna_seq_gp
import scipy
import scipy.stats as stats
import statsmodels
from statsmodels.stats.multitest import multipletests
import os 

file_path = r"Example_Normalised_Counts.csv"
df = pd.read_csv(file_path)
# Filter out genes with 15 or fewer observations (nonzero counts)
filtered_df = df.loc[(df.iloc[:, 1:] > 0).sum(axis=1) > 15]

filtered_df.set_index(filtered_df.columns[0], inplace = True)
Y = filtered_df 
gene_names = Y.index.tolist()

#Y = Y.drop(["Name"], axis=1) #If gene name is in it 
#X = pd.read_csv(r"Example_Sampletable.csv",index_col=[0], header=[0])
X = pd.read_csv(r"/mnt/iusers01/fatpou01/bmh01/msc-bioinf-2024-2025/d20849jl/All sampletable.csv",index_col=[0], header=[0])
X = X[["Time"]]

sparse = True
gp_counts = rna_seq_gp(X,Y.loc[gene_names],M=10, sparse=sparse) 
lik_name = "Negative_binomial" 
results = gp_counts.Two_samples_test(lik_name)

results_df = pd.DataFrame(results)

llr_values = results.loc[:,"log_likelihood_ratio"]

degrees_freedom = 1
p_values = [1 - stats.chi2.cdf(abs(llr), degrees_freedom) for llr in llr_values]

results_df_chisq = 1
results_df["p_value"] = 1 - stats.chi2.cdf(abs(results_df["log_likelihood_ratio"]), results_df_chisq)

p_value = results_df["p_value"].values
rejected, q_values, _, _ = multipletests(p_values, method='fdr_bh')
results_df["q_value"] = q_values
results_df["significant"] = rejected

results_df.to_csv("Example_Results.csv", index=True)