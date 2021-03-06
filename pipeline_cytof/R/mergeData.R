# packages
stopifnot(
	require(optparse),
	require(stringr)
	)
	
# options
option_list <- list(
	       make_option(c("--infile", "-i"), help="tSNE reduced CYTOF RData file"),
	       make_option(c("--outfile", "-o"), help="RData object containing normalised data, tsne dimensions, and clusters"),
               make_option(c("--ClusterX"), default="none", help="ClusterX clusters"),
               make_option(c("--Rphenograph"), default="none", help="Rphenograph clusters"),
               make_option(c("--FlowSOM"), default="none", help="FlowSOM clusters")
	       )

opts <- parse_args(OptionParser(option_list=option_list))

# load data
tsne <- readRDS(opts$infile)

# clusters
# set cluster value to "none" if it was not performed. Then filter on clust methods != "none" before merging
if (opts$Rphenograph != "none") {
   rphenograph <- readRDS(opts$Rphenograph)
   } else {
   rphenograph <- NULL
   }

if (opts$ClusterX != "none") {
   clusterx <- readRDS(opts$ClusterX)
   } else {
   clusterx <- NULL
   }

if (opts$FlowSOM != "none") {
   flowsom <- readRDS(opts$FlowSOM)
   } else {
   flowsom <- NULL
   }

clusters <- list()
cols = c("rphenograph", "clusterx", "flowsom")
i = 0
for (x in list(rphenograph, clusterx, flowsom)) {
    i = i + 1
    if (!is.null(x)) {
       col <- cols[[i]]
       clusters[[col]] <- x
       }
    } # clusters is list of named clustering results


# merge data in clusters list with tSNE data
n = 0 
for (c in clusters){
    n = n + 1
    if (n == 1) {
       name = names(clusters)[n]
       clust <- cbind(tsne, col=c)
       colnames(clust)[dim(clust)[2]] <- paste(name)
       } else {
       	 name = names(clusters)[n]
	 df <- cbind(clust, col=c)
	 clust <- df
	 colnames(clust)[dim(clust)[2]] <- paste(name)
	 }
}

get_donor_id <- function(x){
    d_id <- str_split(x, "\\.", simplify=TRUE)[1]
    return(d_id)	 
}

clust$donor_id <- sapply(rownames(clust), get_donor_id) # donor_id col can be used for identification of event source downstream

print(head(clust))

saveRDS(clust, opts$outfile)
