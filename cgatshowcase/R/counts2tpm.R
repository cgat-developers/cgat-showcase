library(optparse)

option_list <- list(
					make_option(c("--counts"), default="must_specify",
					help="Specify a counts table to convert to tpm's"),
					make_option(c("--genome"), default="must_specify",
					help="Specify a genome name in ensembl convention"),
					make_option(c("--meanfraglength"), default="must_specify",
					help="Specify a mean fragment length of your libaries"),
					make_option(c("--effectivelength"), default="must_specify",
					help="The effective lengths from kallisto"))


opt <- parse_args(OptionParser(option_list=option_list))

print("Running with the following options:")
print(opt)

effectivelength <- read.csv(opt$effectivelength)
rownames(effectivelength) <- effectivelength$X
effectivelength$X <- NULL
effectivelength <- rowMeans(effectivelength)


################################
# tpm function
################################
# function from github gist https://gist.github.com/slowkow/c6ab0348747f86e2748b

counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

counts <- read.csv(opt$counts)
rownames(counts) <- counts$X
counts$X <- NULL

# calculating tpm values assuming a 
tpm <- counts_to_tpm(counts, effectivelength, rep(as.numeric(opt$meanfraglength), length(colnames(counts))))

# output the tpm as a table
write.table(tpm, "DEresults.dir/tpm.tsv", sep="\t")
