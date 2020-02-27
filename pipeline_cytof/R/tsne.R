################################################################################################################
# Purpose: tSNE dimensionality reduction of normalised CYTOF data
# Output: normalised data with tSNE dimensions & parameters
################################################################################################################

### packages
stopifnot(
	require(Rtsne),
	require(optparse)
	)

### command line args
option_list <- list(
	    make_option(c("--infile", "-i"), default=NULL, help="FCS file(s), comma seperated"),
   	    make_option(c("--outfile"), default=NULL, help="RData file - tSNE data"),	
	    make_option(c("--markers", "-m"), default=NULL, help="Comma seperated list of markers (in caps)"),
	    make_option(c("--perplexity"), default=30, help="perplexity value for tSNE"),
	    make_option(c("--iterations"), default=1000, help="tSNE iterations")
	    )

opt <- parse_args(OptionParser(option_list=option_list))

# import data from FCS file(s)
print("Importing data, and subsetting on specified markers for tSNE")
print(paste("Markers:", opt$markers))

df <- read.table(opt$infile, sep="\t")

clust_markers <- strsplit(opt$markers, ",")[[1]] # get markers for clustering as vector
mat <- df[unique(clust_markers)] # subset data on markers for tSNE

mat <- as.matrix(mat)

#tSNE
print("Running tSNE")
mat_tsne_out <- Rtsne(mat,
		      initial_dims=ncol(mat),
		      dims=2,
		      pca=TRUE,
		      check_duplicates=FALSE,
		      perplexity=as.numeric(opt$perplexity),
		      max_iter=as.numeric(opt$iterations))

mat_tsne <- as.data.frame(mat_tsne_out$Y[,1])
colnames(mat_tsne) <- c("tsne_1")
mat_tsne$tsne_2 <- mat_tsne_out$Y[,2]
mat_tsne$perplexity <- as.factor(opt$perplexity)
mat_tsne$iterations <- as.factor(opt$iterations)
mat_tsne$id <- rownames(df)

# save table
print("Done!")
write.table(mat_tsne, file=opt$outfile, sep="\t", row.names=FALSE)
