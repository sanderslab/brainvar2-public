install.packages("BiocManager", repos='http://cran.us.r-project.org')

# BiocManager::install("GenomicFeatures", dependencies=TRUE, ask = FALSE)
BiocManager::install("txdbmaker", ask = FALSE)
BiocManager::install("tximport")
install.packages("data.table", repos='http://cran.us.r-project.org')