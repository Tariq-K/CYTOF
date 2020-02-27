# packages
stopifnot(
	require(cytofkit),
	require(optparse),
	require(Rphenograph)
	)
	
# options
option_list <- list(
	       make_option(c("--infile", "-i"), help="normalised CYTOF data"),
	       make_option(c("--outfile", "-o"), help="RData object containing clustering info"),
	       make_option(c("--k"), default=30, help="k parameter for Rphenograph"),
   	       make_option(c("--markers", "-m"), default=NULL, help="Comma seperated list of markers (in caps)"),
  	       make_option(c("--noClusters"), default=NULL, help="No clusters parameter for FlowSOM")
	       )

opts <- parse_args(OptionParser(option_list=option_list))

# load data
data <- read.table(opts$infile, sep="\t")

clust_markers <- strsplit(opts$markers, ",")[[1]] # get markers for clustering as vector
data <- data[unique(clust_markers)] # subset data on markers for tSNE

# clustering
clust <- Rphenograph(data, k=as.integer(opts$k))
clusters <- as.data.frame(factor(membership(clust[[2]])))

# annotate df
colnames(clusters) <- c(paste0("phenograph_k", opts$k))
rownames(clusters) <- rownames(data) # add index 

# save df
write.table(clusters, file=opts$outfile, sep="\t")
