setwd(".")

library(tximport)
library(DESeq2)
library(optparse)
library(biomaRt)

# Need a specific way of making sure the pipeline can run with different TxDb datasets
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


option_list <- list(
					make_option(c("--design"), default="must_specify",
					help="To run DESEq2 you need to specify the location of a design.tsv file according to the pipeline documentation"),
					make_option(c("--contrast"), default="must_specify",
					help="must specify a contrast in the pipeline.yml file"),
					make_option(c("--fdr"), default=0.05,
					help="set an optional fdr, will default to 0.05"))

opt <- parse_args(OptionParser(option_list=option_list))

print("Running with the following options:")
print(opt)


design = read.table(opt$design, header=TRUE, fill=TRUE)


dir <- "kallisto.dir"

files <- file.path(dir, samples$track, "abundance.h5")


mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         #dataset = "hsapiens_gene_ensembl",
                         dataset = "hsapiens_gene_ensembl",
                         host="www.ensembl.org")

t2g <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id","ensembl_gene_id"), mart = mart)

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = t2g)

rownames(design) <- design$track
dds <- DESeqDataSetFromTximport(txi.kallisto, design, ~ group)

dds = suppressMessages(
  DESeq(dds, test="Wald", fitType="parametric"))

dds$group <- relevel(dds$group, ref = opt$contrast)



res = suppressMessages(results(dds))
res = as.data.frame(res)

dir.create("plots.dir/")
png(paste0(c("plots.dir", "MA.png"), collapse="_"))
plotMA(res, alpha=opt$fdr)
dev.off()

png("plots.dir/dispersion.png")
plotDispEsts(dds)
dev.off()
