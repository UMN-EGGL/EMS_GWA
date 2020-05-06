#import sys
#sys.path.append('/home/mccuem/norto253/mypython/lib/python2.7/site-packages')

import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from scipy import stats
import pandas as pd


def _mix(G, G0_standardized_val, G1_standardized_val, mixing):
    #logging.info("concat G1, mixing {0}".format(mixing))
    G[:,0:G0_standardized_val.shape[1]] = G0_standardized_val
    G[:,0:G0_standardized_val.shape[1]] *= (np.sqrt(1.0-mixing))
    G[:,G0_standardized_val.shape[1]:] = G1_standardized_val
    G[:,G0_standardized_val.shape[1]:] *= np.sqrt(mixing)

def _find_mixing(G, covar, G0_standardized_val, G1_standardized_val, h2, y):
    import fastlmm.util.mingrid as mingrid
    assert h2 is None, "if mixing is None, expect h2 to also be None"
    resmin=[None]
    def f(mixing,G0_standardized_val=G0_standardized_val,G1_standardized_val=G1_standardized_val,covar=covar,y=y,**kwargs):
        _mix(G, G0_standardized_val,G1_standardized_val,mixing)
        lmm = fastLMM(X=covar, Y=y, G=G, K=None, inplace=True)
        result = lmm.findH2()
        if (resmin[0] is None) or (result['nLL']<resmin[0]['nLL']):
            resmin[0]=result
        return result['nLL']
    mixing,nLL = mingrid.minimize1D(f=f, nGrid=10, minval=0.0, maxval=1.0,verbose=False)
    h2 = resmin[0]['h2']
    return mixing, h2

RcppCNPy = importr('RcppCNPy')
lrgpr = importr('lrgpr')
corpcor = importr('corpcor')

ro.r('source("/home/mccuem/norto253/lrgpr/R/lrgpr.R")')
ro.r('source("/home/mccuem/norto253/lrgpr/R/genericFunctions.R")')



####load GRM, covariates, and pheno
G0 = np.array(ro.r('G0<-npyLoad("/home/mccuem/norto253/WP_2Million/Elaine_v3/R_scripts/G0_mat.npy")'))
covar = np.array(ro.r('covar<-npyLoad("/home/mccuem/norto253/WP_2Million/Elaine_v3/R_scripts/covar_mat.npy")'))
y = np.array(ro.r('y<-npyLoad("/home/mccuem/norto253/WP_2Million/Elaine_v3/R_scripts/y_mat.npy")'))


#norm_factor = 1./np.sqrt((G0**2).sum() / float(G0.shape[0]))
#G0_standardized_val = norm_factor * G0
from pysnptools.standardizer import DiagKtoN
#G0_standardized_val = DiagKtoN(G0.shape[0]).standardize(G0)
G0_standardized_val = G0

from fastlmm.inference.lmm_cov import LMM as fastLMM
lmm = fastLMM(X=covar, Y=y, G=G0_standardized_val)
result = lmm.findH2()
#dir(lmm) #lists attributes

if result['h2'] > -1:
    h2 = result['h2']
else:
    h2 = result['h2'][0]

residual_var = 1-h2
delta = residual_var/h2
m1_df = sum(lmm.S/(lmm.S + delta))
m1 = result['nLL'][0]*-1

#nr, nc = G0_standardized_val.shape
#G0_standardized_val_vec = ro.FloatVector(G0_standardized_val.transpose().reshape((G0_standardized_val.size)))
#G0r = ro.r.matrix(G0_standardized_val_vec, nrow=nr, ncol=nc)
#ro.globalenv['G0'] = G0r
#ro.r('dcmp<-fast.svd(G0)')
#ro.r('fitA<-lrgpr(y ~ covar[,1:2], decomp=dcmp, diagnostic=TRUE)')
#m1_df = np.array(ro.r('fitA_df<-sum(fitA$hii)'))
#m1 = np.array(ro.r('fitA$logLik'))


keep_list = ro.r('keep_list<-list()')
ro.r('load("/home/mccuem/norto253/WP_2Million/Elaine_v3/FastLMM_input_and_output/snp_list.Rdat")')
#ro.r('load("/Users/Nichol/Google Drive/GWASems/FastLMM_input_and_output/snp_list2.Rdat")')
ro.r('load("/home/mccuem/norto253/WP_2Million/Elaine_v3/FastLMM_input_and_output/X_data.Rdat")')

res = np.array(ro.r('res<-matrix(NA, length(snp_list), 13)'))


f_handle = open('result.txt', 'w')

items = [
         (0),
         (0),
         (m1),
         (m1_df),
         (0),
         (0),
         (delta),
         (h2),
         (0),
         (residual_var),
         (1),
        ]
with open('result.txt','a') as f_handle:
    np.savetxt(f_handle, np.array(items), newline=" ",   fmt='%.3f')
    f_handle.write("\n")


#for i in range(1, 5):

for i in range(1, res.shape[0]+1):
    ro.globalenv['i'] = i
    #keep_list2 = ro.r('keep_list2<-c(snp_list[i], keep_list)')
    #keep_list2 = ro.r('keep_list2<-c(snp_list2[snp_list2[,1]==snp_list2[snp_list2[,3]==snp_list[i],][1],3], keep_list)')
    keep_list2 = ro.r('keep_list2<-c(snp_list[i], colnames(X_data)[which(colnames(X_data)==snp_list[i])+1],colnames(X_data)[which(colnames(X_data)==snp_list[i])-1], keep_list)')
    G1 = np.array(ro.r('XX<-as.matrix(X_data[,colnames(X_data)%in%keep_list2])/(sqrt(length(keep_list2)))'))
    
    #norm_factor = 1./np.sqrt((G1**2).sum() / float(G1.shape[0]))
    #G1_standardized_val = norm_factor * G1
    from pysnptools.standardizer import DiagKtoN
    G1_standardized_val = DiagKtoN(G1.shape[0]).standardize(G1)
    #G1_standardized_val = G1
    
    lmmB = lmm
    W = G1_standardized_val.copy()
    UGup,UUGup = lmmB.rotate(W)
    i_up = np.zeros((W.shape[1]), dtype=np.bool)
    i_G1 = np.ones((W.shape[1]), dtype=np.bool)
    result = lmmB.findH2_2K(nGridH2=10, minH2=0.0, maxH2=0.99999, i_up=i_up, i_G1=i_G1, UW=UGup, UUW=UUGup)
    m2 = result['nLL'][0]*-1
    if result['h2'] > -1:
        h2 = result['h2']
    else:
        h2 = result['h2'][0]
    if result['h2_1'] > -1:
        h2_1 = result['h2_1']
    else:
        h2_1 = result['h2_1'][0]
    h2_all = h2 + h2_1
    residual_var = 1-(h2_all)
    delta = (1 - h2_all)/h2_all
    mixing = h2_1/h2_all
    G = np.empty((G0.shape[0],G0.shape[1]+G1.shape[1]))
    G[:,0:G0_standardized_val.shape[1]] = G0_standardized_val
    G[:,0:G0_standardized_val.shape[1]] *= (np.sqrt(1.0-mixing))
    G[:,G0_standardized_val.shape[1]:] = G1_standardized_val
    G[:,G0_standardized_val.shape[1]:] *= np.sqrt(mixing)

    dcmp_u, dcmp_s, dcmp_v = np.linalg.svd(G, full_matrices=False)
    s2 = dcmp_s**2
    m2_df = sum(s2/(s2 + delta))
    m2 = result['nLL'][0]*-1

    #lmmC = fastLMM(X=covar, Y=y, G=G)
    #result = lmmC.findH2()
    #s2 = lmmC.S
    #m2_df = sum(s2/(s2 + delta))
    
    #nr, nc = G.shape
    #Gvec = ro.FloatVector(G.transpose().reshape((G.size)))
    #Gr = ro.r.matrix(Gvec, nrow=nr, ncol=nc)
    #ro.globalenv['G'] = Gr
    #ro.r('dcmp<-fast.svd(G)')
    #ro.r('fitA<-lrgpr(y ~ covar[,1:2], decomp=dcmp, diagnostic=TRUE)')
    #m2_df = np.array(ro.r('fitA_df<-sum(fitA$hii)'))
    #m2 = np.array(ro.r('fitA$logLik'))
    
    df_diff =  m2_df-m1_df
    if(df_diff<0):
        df_diff = 0.001
    pval = 1 - stats.chi2.cdf(-2*(m1 - m2), df_diff)
    #pval = 1 - stats.chi2.cdf(-2*(m1 - m2), 1)
    if(m2<m1):
        pval = 1
    if(m2>m1 and pval<0.05):
        m1 = m2
        m1_df = m2_df
        #keep_list = np.array(ro.r('keep_list<-c(keep_list, snp_list[i])'))
        #keep_list = np.array(ro.r('keep_list<-c(keep_list, snp_list2[snp_list2[,1]==snp_list2[snp_list2[,3]==snp_list[i],][1],3])'))
        keep_list = np.array(ro.r('keep_list<-c(keep_list, snp_list[i], colnames(X_data)[which(colnames(X_data)==snp_list[i])+1],colnames(X_data)[which(colnames(X_data)==snp_list[i])-1])'))
    items = [
             (i),
             (len(keep_list)),
             (m1),
             (m1_df),
             (m2),
             (m2_df),
             (delta),
             (h2),
             (h2_1),
             (residual_var),
             (pval),
             ]
    with open('result.txt','a') as f_handle:
        np.savetxt(f_handle, np.array(items), newline=" ",   fmt='%.3f')
        f_handle.write("\n")
