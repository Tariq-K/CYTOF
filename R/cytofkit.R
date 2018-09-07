################################################################################################################
# Wrapper for cytofkit
# Purpose: clustering of cytof data based on customisable marker selection
# Output: Rdata object containing dataframe of dimensionality reduced data and cluster annotations
# 	  png of tSNE with phenograph and clusterX clustering
################################################################################################################

### packages
stopifnot(
	require(cytofkit),
	require(gridExtra),
	require(grid),
	require(ggplot2),
	require(lattice),
	require(stringr),
	require(Rtsne),
	require(optparse)
	)

### command line args
option_list <- list(
	    make_option(c("--infiles", "-i"), default=NULL, help="FCS file(s), comma seperated"),
	    make_option(c("--out_norm"), default=NULL, help="RData file - normalised data"),
   	    make_option(c("--out_tsne"), default=NULL, help="RData file - tSNE data"),	
	    make_option(c("--markers", "-m"), default=NULL, help="Comma seperated list of markers (in caps)"),
	    make_option(c("--perplexity"), default=30, help="perplexity value for tSNE"),
	    make_option(c("--iterations"), default=1000, help="tSNE iterations"),
	    make_option(c("--dataTransform"), default="cytofAsinh", help="Data transformation method: cytofAsinh, arcsinh, or none"),
	    make_option(c("--no_events"), default=90000, help="No. of events to sample from each sample. Default=90,000, for all (1 donor only)set as 'all' ")
	    )

opt <- parse_args(OptionParser(option_list=option_list))

# import data from FCS file(s)
files <- as.list(str_split(opt$infiles, ",", simplify=TRUE))

if (length(files) >1){
   # merge x events (specified by opt$no_events) randomly sampled without replacement from multiple donors
   mat <- cytof_exprsMerge(files, comp=FALSE, transformMethod=opt$dataTransform, mergeMethod="ceil", fixedNum=opt$no_events)
   mat <- as.data.frame(mat)
} else {
  # limit events from 1 donor by opt$no_events
  mat <- cytof_exprsExtract(fcsFile=opt$infile, comp=FALSE, transformMethod=opt$dataTransform)
  mat <- as.data.frame(mat)
  if (!opt$no_events == "all"){
    subsample <- sample(1:nrow(mat), opt$no_events)
    mat <- mat[subsample, ]
    }
  }

# now remove markers used for gating & any redundant markers
print(head(mat))

names <- str_split(names(mat), "_", simplify=T)[, 2] # first rename cols
names <- str_replace(names, ">", "")
colnames(mat) <- str_to_upper(names)
clust_markers <- strsplit(opt$markers, ",")[[1]] # get markers for clustering as vector
print(str(clust_markers))
mat <- mat[clust_markers] # subset data

# remove any duplicate rows before running tSNE
dups <- duplicated(mat)
mat <- mat[!dups, ]

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
