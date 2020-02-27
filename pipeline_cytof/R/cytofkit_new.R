################################################################################################################
# Wrapper for cytofkit
# Purpose: subsetting and normalisation of cytof data based on customisable marker selection
# Output: normalised FCS data
################################################################################################################

### packages
stopifnot(
	require(cytofkit),
	require(stringr),
	require(optparse)
	)

### command line args
option_list <- list(
	    make_option(c("--infiles", "-i"), default=NULL, help="FCS file(s), comma seperated"),
	    make_option(c("--outfile"), default=NULL, help="RData file - normalised data"),
	    make_option(c("--dataTransform"), default="cytofAsinh", help="Data transformation method: cytofAsinh, arcsinh, or none"),
   	    make_option(c("--markers", "-m"), default=NULL, help="Comma seperated list of markers (in caps)"),		
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

# save table
write.table(mat, file=opt$outfile, sep="\t")
