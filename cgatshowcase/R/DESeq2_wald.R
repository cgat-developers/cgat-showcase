# This script will process the ouput of kallisto and then perofrm Deseq2 using the wald test.
# Returns a counts table, MA and dispersion plot. As an additional object it also outputs a
# RDS data object for importing into an R environment for further exploration.

setwd(".")

library(tximport)
library(DESeq2)
library(optparse)
library(biomaRt)


option_list <- list(
					make_option(c("--design"), default="must_specify",
					help="To run DESEq2 you need to specify the location of a design.tsv file according to the pipeline documentation"),
					make_option(c("--contrast"), default="must_specify",
					help="must specify a contrast in the pipeline.yml file"),
					make_option(c("--refgroup"), default="must_specify",
					help="must specify a reference group to compare against in the pipeline.yml file"),
					make_option(c("--fdr"), default=0.05,
					help="set an optional fdr, will default to 0.05"),
					make_option(c("--biomart"), default="must_specify",
					help="must specify a biomart reference dataset"))

opt <- parse_args(OptionParser(option_list=option_list))

print("Running with the following options:")
print(opt)

# Reads in the design file
design = read.table(opt$design, header=TRUE, fill=TRUE)


dir <- "kallisto.dir"
sample_track <- design$track
files <- file.path(dir, sample_track, "abundance.tsv")

# Convert the transcripts to gene ids
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = opt$biomart,
                         host="www.ensembl.org")

t2g <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id","ensembl_gene_id"), mart = mart)

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = t2g)

# Set up the Deseq2 object using tximport package
rownames(design) <- design$track
dds <- DESeqDataSetFromTximport(txi.kallisto, design, ~ group)

# run the deseq2 model using the wald test
dds = suppressMessages(
  DESeq(dds, test="Wald", fitType="parametric"))

dds$group <- relevel(dds$group, ref = opt$refgroup)


# return the results table and output as a dataframe
res = suppressMessages(results(dds))
res = as.data.frame(res)

# write the counts table to file
colnames(txi.kallisto$counts) <- design$track 
write.csv(txi.kallisto$counts, file="DEresults.dir/counts.csv")

# write the results dataframe
write.csv(res, file="DEresults.dir/res.csv")

# write the effective length
colnames(txi.kallisto$length) <- design$track 
write.csv(txi.kallisto$length, "DEresults.dir/length.csv")

# save the RDS for processing outside of the pipeline
saveRDS(dds, "DEresults.dir/dds.rds")

# Generate plotting of the Deseq2 objects
dir.create("plots.dir/", showWarnings = FALSE)
png(paste0(c("plots.dir/", "MA.png"), collapse="_"))
plotMA(dds, alpha=opt$fdr)
dev.off()

png("plots.dir/dispersion.png")
plotDispEsts(dds)
dev.off()
