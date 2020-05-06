install.packages('tidyverse')

# install some packages from CRAN
install.packages(
    c('Rcpp','RcppArmadillo','RcppEigen','inline',
      'corpcor','RcppCNPy','plyr',
      'gap','Matrix','aod','gdata','itertools','nadiv',
      'sfsmisc','lmtest','HardyWeinberg','plotrix',
      'doParallel','abind','extrafont','extrafontdb','shape'
    )
)

# install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}

BiocManager::install(
    pkgs=c('rhdf5','LDheatmap','snpStats','CGEN'),
    ask=FALSE
)

# These are needed for synbreed
install.packages(c('BGLR', 'doBy', 'igraph', 'qtl'))

install.packages('regress_1.3-15.tar.gz',repo=NULL)
install.packages('synbreed_0.12-9.tar.gz',repo=NULL)

# GenAbel
install.packages("GenABEL.data_1.0.0.tar.gz",repo=NULL)
install.packages("GenABEL_1.8-0.tar.gz",repo=NULL)


# Now install lrgpr
install.packages(
    c("Rcpp", "RcppGSL", "RcppProgress", "MASS", "formula.tools",
      "BH", "doParallel", "bigmemory", "bigmemory.sri", "aod")
)
install.packages("lrgpr_0.1.12.tar.gz",repo=NULL)
