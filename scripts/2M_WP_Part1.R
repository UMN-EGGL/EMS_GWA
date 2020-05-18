DIR_working<-"~/WP_2Million/Elaine_v3/R_scripts"
DIR_data<-"~/WP_2Million/Elaine_v3/Data"
DIR_fastLMM<-"~/WP_2Million/Elaine_v3/FastLMM input and output"

setwd(DIR_working)

library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
# Here
library(inline)
library(rhdf5)
library(corpcor)
library(RcppCNPy)
library(plyr)
library(synbreed)
library(GenABEL)
library(gap)
library(Matrix)
library(aod)
library(gdata)
library(itertools)
#library(asreml)
library(nadiv)
#library(gap)
library(sfsmisc)
library(lmtest)
library(CGEN)
library(HardyWeinberg)
#library(itertools)
library(plotrix)
library(doParallel)
library(abind)  
registerDoParallel(cores=3)

library(extrafont)
font_import()
fonts()
library(extrafontdb)
#fonttable()
#font_install("fontcm")
fonts()
loadfonts(device = "pdf", quiet=TRUE)
library(shape)

library(lrgpr)
source('~/lrgpr/R/lrgpr.R')
source('~/lrgpr/R/genericFunctions.R')
source('~/lrgpr/R/plink.R')
source('~/lrgpr/R/plots.R')
library(formula.tools)
as.character.formula<-formula.tools:::as.character.formula


source("~/WP_2Million/Elaine_v3/R_scripts/helper functions.R")

center_mean <- function(x) {
         ones = rep(1, nrow(x))
         x_mean = ones %*% t(colMeans(x))
         x - x_mean
    }

relmat_centered<-function(data){
    a<-center_mean(data)
    b<-tcrossprod(a)
    c<-as.matrix(b/ncol(data))
    rownames(c)<-rownames(data)
    colnames(c)<-rownames(data)
    return(c)
}

###read pheno files
p_welsh <- read.table("~/WP_2Million/Elaine_v3/Data/WP_2015/WP_2M_Imp_Combo_Pheno.txt", header=T)
rownames(p_welsh)<- p_welsh$VCF_ID

p_welsh$Breed<-as.factor('Welsh')
p_welsh[1:10,1:10]

#Reduce the dataframe to only account for specific sections of horses and then create the text file to keep only those horses for plink###
Section_Index <- which((p_welsh$Section=='A')|(p_welsh$Section=='B')|(p_welsh$Section=='C')|(p_welsh$Section=='D'))
p_welsh <- p_welsh[Section_Index,]
VCF_ID_Table <-data.frame(p_welsh$VCF_ID,p_welsh$VCF_ID)
write.table(VCF_ID_Table,file="Section_A_B_C_D_VCFID.txt",quote=F,col.names=F,row.name=F)

p<-p_welsh
rownames(p)
rm(p_welsh)


p$Height_s<-scale(p$Height)
p$NH_s<-scale(p$NH)
p$GH_s<-scale(p$GH)
p$GLU_s<-scale(p$GLU)
p$INS[p$INS<1.53& !is.na(p$INS)] <- 1.5 
p$logINS <- log10(p$INS)
p$logINS_s<-scale(p$logINS)
p$GLU_OGT_s<-scale(p$GLU_OGT)
p$INS_OGT[p$INS_OGT<1.53& !is.na(p$INS_OGT)] <- 1.5 
p$logINS_OGT <- log10(p$INS_OGT)
p$logINS_OGT_s<-scale(p$logINS_OGT)
p$logTG <- log10(p$TG)
p$logTG_s<-scale(p$logTG)
p$sqrtNEFA <- sqrt(p$NEFA)
p$sqrtNEFA_s<-scale(p$sqrtNEFA)
p$logACTH2 <- log10(p$ACTH_2)
p$logACTH2_s<-scale(p$logACTH2)
p$sqrtLeptin<-sqrt(p$Leptin)
p$sqrtLeptin_s<-scale(p$sqrtLeptin)
p$sqrtAdiponectin<-sqrt(p$Adiponectin)
p$sqrtAdiponectin_s<-scale(p$sqrtAdiponectin)
p$Female <- ifelse(p$Female=='y', c(1), c(0))
p$Female_s <- scale(p$Female)
p$AGE_s<-scale(p$AGE)
p$AGE_s[is.na(p$AGE_s)]<-mean(p$AGE_s, na.rm=T)
p$NH_s<-as.numeric(p$NH_s)
p$GH_s<-as.numeric(p$GH_s)
p$GLU_s<-as.numeric(p$GLU_s)
p$GLU_OGT_s<-as.numeric(p$GLU_OGT_s)
p$logINS_s<-as.numeric(p$logINS_s)
p$logINS_OGT_s<-as.numeric(p$logINS_OGT_s)
p$logTG_s<-as.numeric(p$logTG_s)
p$sqrtNEFA_s<-as.numeric(p$sqrtNEFA_s)
p$logACTH2_s<-as.numeric(p$logACTH2_s)
p$sqrtLeptin_s<-as.numeric(p$sqrtLeptin_s)
p$sqrtAdiponectin_s<-as.numeric(p$sqrtAdiponectin_s)
p$AGE_s<-as.numeric(p$AGE_s)
p$Month<-relevel(p$Month, ref='DEC')
p$LAM_recode<-0
p$LAM_recode[p$LAM=='y']<-1
p$obese<-0
p$obese[p$BCS>=7 & !is.na(p$BCS)]<-1
p$Month<-reorder(p$Month, new.order=c('JAN','MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'DEC'))
p$Number_ID<-as.factor(p$Number_ID)
p$Owner<-as.factor(p$Owner2)
unique(p$Owner)


p$cluster<-p$Owner
p$cluster2<-as.factor(paste(as.character(p$Owner), as.character(p$Month)))

#load family data
fam_file_Welsh <- "~/WP_2Million/Elaine_v3/Data/WP_2015/WP_Imputed_2Million_Combined2.tfam"
FAM_Welsh<- read.fam(fam_file_Welsh)
head(FAM_Welsh)
p<-p[match(FAM_Welsh[,2], rownames(p)),]


###load Welsh genotype data
tped_file_Welsh<- "~/WP_2Million/Elaine_v3/Data/WP_2015/WP_Imputed_2Million_Combined2.tped"
convertToBinary(tped_file_Welsh, "./test_Welsh.binary", "TPED", onlyCheckFormat = FALSE,simpleAnnotation = FALSE)

M<- attach.big.matrix("./test_Welsh.binary_descr", readonly=TRUE)
M[1:10,1:10]

alleles <- readLines("./test_Welsh.binary_alleles")
alleles <- alleles[-1]
xx<-ldply(strsplit(alleles, " "))
snp_info<- cbind(xx[[2]], xx[[1]],xx[[4]],xx[[5]],xx[[6]])
rm(xx, alleles)
colnames(snp_info)<-c('SNP', 'Chr', 'Pos', 'REF_allele', 'ALT_allele')
rownames(snp_info)<-seq(1:nrow(snp_info))
head(snp_info)

#clean for non-bi-allelic snps, fixed alleles, snps with missing genotypes
NA_geno_sum<-colSums(is.na(M[,1:ncol(M)]))>0
snp_info<-snp_info[snp_info[,'REF_allele']%in%c('A', 'G', 'C', 'T') & snp_info[,'ALT_allele']%in%c('A', 'G', 'C', 'T')  & !snp_info[,'SNP']%in%colnames(M)[NA_geno_sum],]


ALT_allele_freq<-getAlleleFreq(M[,match(as.character(snp_info[,'SNP']), colnames(M))])
snp_info<-cbind(snp_info, ALT_allele_freq)
MAF<-as.numeric(snp_info[,'ALT_allele_freq'])
MAF[MAF>0.5]<-1-(MAF[MAF>0.5])
snp_info<-cbind(snp_info, MAF)
rm(MAF)
hist(as.numeric(snp_info[,'MAF']))
table(as.numeric(snp_info[,'MAF'])==0)



z <- filebacked.big.matrix(nrow(M), nrow(snp_info), 
                           backingfile="M.bin",
                          descriptorfile="M.desc",
                          dimnames=list(rownames(p), colnames(M)[match(as.character(snp_info[,'SNP']), colnames(M))]))                         
z[1:nrow(M),1:nrow(snp_info)]<- M[,match(as.character(snp_info[,'SNP']), colnames(M))]  
rm(z)  
gc()
M <- attach.big.matrix("M.desc")




##calculate hardy weinberg
HWE<-array(1, nrow(snp_info))
names(HWE)<-snp_info[,'SNP']


mylm<-function (i){
	xx<-M[,grep(paste('chr_', i, '_', sep=""), colnames(M))]	
	xx[,colMeans(xx)/2>0.5]<-(xx[,colMeans(xx)/2>0.5]*-1)+2

	hwe<-apply(xx, 2, function (j) {
		df<-matrix(0,1,3)
		colnames(df)<-c(0,1,2)
		x<-table(j)
		df[,match(names(x), colnames(df))]<-x
		HWChisq(df)$pval
		})
	return(hwe)
}

registerDoParallel(cores=3)
jobs<- foreach(i=1:31, .combine=abind) %dopar% mylm(i)
HWE[match(names(jobs), names(HWE))]<-jobs
rm(jobs)
snp_info<-cbind(snp_info, HWE)
rm(HWE)

save(snp_info, file="snp_info_Welsh2M.Rdat")
