#import sys
#sys.path.append('/home/mccuem/norto253/mypython/lib/python2.7/site-packages')
import numpy as np
import rpy2.robjects as ro
#from scipy import stats
import pandas as pd


####load GRM, covariates, and pheno
G0 = np.load("/home/mccuem/norto253/WP_2Million/Elaine_v3/R_scripts/G0_mat.npy")
G1 = np.load("/home/mccuem/norto253/WP_2Million/Elaine_v3/FastLMM_input_and_output/G1_mat.npy")
test_snps = np.load("/home/mccuem/norto253/WP_2Million/Elaine_v3/FastLMM_input_and_output/test_snps_mat.npy")
covar = np.load("/home/mccuem/norto253/WP_2Million/Elaine_v3/R_scripts/covar_mat.npy")
y = np.load("/home/mccuem/norto253/WP_2Million/Elaine_v3/R_scripts/y_mat.npy")
from pysnptools.standardizer import DiagKtoN
#G0_standardized_val = DiagKtoN(G0.shape[0]).standardize(G0)
G0_standardized_val = G0
G1_standardized_val = DiagKtoN(G1.shape[0]).standardize(G1)
#G1_standardized_val = G1
from fastlmm.inference.lmm_cov import LMM as fastLMM
lmm = fastLMM(X=covar, Y=y, G=G0_standardized_val)

W = G1_standardized_val.copy()
UGup,UUGup = lmm.rotate(W)
i_up = np.zeros((W.shape[1]), dtype=np.bool)
i_G1 = np.ones((W.shape[1]), dtype=np.bool)
result = lmm.findH2_2K(nGridH2=10, minH2=0.0, maxH2=0.99999, i_up=i_up, i_G1=i_G1, UW=UGup, UUW=UUGup)
if result['h2'] > -1:
    h2 = result['h2']
else:
    h2 = result['h2'][0]
if result['h2_1'] > -1:
    h2_1 = result['h2_1']
else:
    h2_1 = result['h2_1'][0]

res = lmm.nLLeval_2K(h2=h2, h2_1=h2_1, dof=None, scale=1.0, penalty=0.0, snps=test_snps, UW=UGup, UUW=UUGup, i_up=i_up, i_G1=i_G1)

beta = res['beta']

items = [
         ('SnpWeight', beta[:,0]),
         ('SnpWeightSE', np.sqrt(res['variance_beta'][:,0])),
         ('SnpFractVarExpl', np.sqrt(res['fraction_variance_explained_beta'][:,0]))
         ]
frame = pd.DataFrame.from_items(items)
frame.to_csv("gwas_result.txt", sep="\t", index=False)


