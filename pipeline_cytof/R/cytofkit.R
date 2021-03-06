################################################################################################################
# Wrapper for cytofkit
# Purpose: clustering of cytof data based on customisable marker selection
# Output: Rdata object containing dataframe of dimensionality reduced data and cluster annotations
# 	  png of tSNE with phenograph and clusterX clustering
################################################################################################################

### packages
stopifnot(
	require(Rtsne),
	require(optparse)
	)

### command line args
option_list <- list(
	    make_option(c("--infile", "-i"), default=NULL, help="FCS file(s), comma seperated"),
   	    make_option(c("--out_tsne"), default=NULL, help="RData file - tSNE data"),	
	    make_option(c("--markers", "-m"), default=NULL, help="Comma seperated list of markers (in caps)"),
	    make_option(c("--perplexity"), default=30, help="perplexity value for tSNE"),
	    make_option(c("--iterations"), default=1000, help="tSNE iterations"),
	    )

opt <- parse_args(OptionParser(option_list=option_list))

# import data from FCS file(s)
mat <- readRDS(opt$infile)
mat <- as.matrix(mat)

#tSNE
# ma_tsne <- cytof_dimReduction(data=mat, method = "tsne")

mat_tsne_out <- Rtsne(mat,
		      initial_dims=ncol(mat),
		      dims=2,
		      pca=TRUE,
		      check_duplicates=FALSE,
		      perplexity=as.numeric(opt$perplexity),
		      max_iter=as.numeric(opt$iterations))

mat_tsne <- mat
mat_tsne$tsne_1 <- mat_tsne_out$Y[,1]
mat_tsne$tsne_2 <- mat_tsne_out$Y[,2]
mat_tsne$perplexity <- opt$perplexity
mat_tsne$iterations <- opt$iterations

# save Rdata file
saveRDS(mat_tsne, opt$out_tsne)
saveRDS(mat, opt$out_norm)
