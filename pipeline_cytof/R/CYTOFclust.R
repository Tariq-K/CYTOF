# packages
stopifnot(
	require(cytofkit),
	require(optparse)
	)
	
# options
option_list <- list(
	       make_option(c("--infile", "-i"), help="tSNE reduced or high dimensionality CYTOF RData file"),
	       make_option(c("--outfile", "-o"), help="RData object containing clustering info"),
               make_option(c("--clusterMethod", "-m"), help="Rphenograph, FlowSOM, or ClusterX"),
	       make_option(c("--noClusters"), default=NULL, help="No clusters parameter for FlowSOM")
	       )

opts <- parse_args(OptionParser(option_list=option_list))

# load data
data <- readRDS(opts$infile)
print(head(data))

# clustering
if (opts$clusterMethod == "FlowSOM"){
   clust <- cytof_cluster(xdata=data, method=opts$clusterMethod, FlowSOM_k=opts$clust_no)
   }
if (opts$clusterMethod == "Rphenograph"){
   clust <- cytof_cluster(xdata=data, method=opts$clusterMethod) # Phenograph has an internal parameter k, 
   }	    			      				 # by default this is set to 30 in cytofkit
if (opts$clusterMethod == "ClusterX"){				 
   clust <- cytof_cluster(ydata=data, method=opts$clusterMethod)
   }

saveRDS(clust, opts$outfile)
