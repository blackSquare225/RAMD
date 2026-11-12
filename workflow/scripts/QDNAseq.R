args = commandArgs(trailingOnly = TRUE)

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (!requireNamespace("QDNAseq", quietly = TRUE)) {
    system("sudo apt-get update")
    system("sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev build-essential")
    install.packages("curl")
    BiocManager::install("QDNAseq")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", dependencies = TRUE)
}

if (!requireNamespace("QDNAseq.hg38", quietly = TRUE)) {
  devtools::install_github("asntech/QDNAseq.hg38@main")
}

library(QDNAseq)
library(QDNAseq.hg38)

binSize <- as.numeric(args[1])  # Bin size (e.g., 100)
bamFile <- args[2]              # Path to BAM file
outdir <- args[3] 

# Print the arguments to confirm
cat("Bin size:", binSize, "\n")
cat("BAM file:", bamFile, "\n")
cat("Output directory:", outdir, "\n")


# Create the output directory if it doesn't exist
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
  cat("Output directory created:", outdir, "\n")
}

#provide QDNAseq bin annotations of size 1, 5, 10, 15, 30, 50, 100, 500 and 1000 kbp
#for the human genome build hg38
#see https://github.com/asntech/QDNAseq.hg38 for more details

bins <- getBinAnnotations(binSize=binSize, genome="hg38")

#load files.bam
readCounts <- binReadCounts(bins, bamfiles=bamFile, ext = "bam")
name <-gsub(".minimap2.wgs.srt.bam","", basename(bamFile))

#Plot a raw copy number profile (read counts across the genome)
pdf(file.path(outdir,paste0(name, ".QDNAseq.",binSize,"kbp.raw.count.pdf")))
plot(readCounts, logTransform=F, ylim=c(-50, 200))
dev.off()

#highlight bins that will be removed with default filtering
##N.B. see https://rdrr.io/bioc/QDNAseq/man/highlightFilters.html for more details about filtering (i.e. mappability, bases, ect..)
pdf(file.path(outdir,paste0(name, ".QDNAseq.",binSize,"kbp.raw.count.filter.pdf")))
plot(readCounts, logTransform=F, ylim=c(-50, 200))
highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE,
                 mappability = NA, bases = NA)
dev.off()

#Apply filters chosen above and plot median read counts as a function of GC content and mappability
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

pdf(file.path(outdir,paste0(name, ".QDNAseq.",binSize,"kbp.isobar.pdf")))
isobarPlot(readCountsFiltered)
dev.off()

#Chromosomes: character vector specifying which chromosomes to filter out (default=c("X", "Y", "MT"))

#Estimate the correction for GC content and mappability
readCountsFiltered <- estimateCorrection(readCountsFiltered)

#plot for the relationship between the observed standard deviation in the data
#and its read depth.
#The theoretical expectation is a linear relationship, which is shown in the plot
#with a black line. Samples with low-quality DNA will be noisier than expected
#and appear further above the line than good-quality samples.

pdf(file.path(outdir,paste0(name, ".QDNAseq.",binSize,"kbp.noiseplot.pdf")))
noisePlot(readCountsFiltered)
dev.off()

#apply the correction for GC content and mappability
copyNumbers <- correctBins(readCountsFiltered)

#median normalization
copyNumbersNormalized <- normalizeBins(copyNumbers)

#Smoothing outliers
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

#Plotting
pdf(file.path(outdir,paste0(name, ".QDNAseq.",binSize,"kbp.cn.smooth.pdf")))
plot(copyNumbersSmooth)
dev.off()


#DOWNSTREM ANALYSIS
#apply the segmentation algorithm (Circular binary segmetation) and normalize the data
# see https://www.rdocumentation.org/packages/QDNAseq/versions/1.8.0/topics/segmentBins for deteails regarding transforming data function
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")

copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

#Copy number profile after segmenting
pdf(file.path(outdir,paste0(name, ".QDNAseq.",binSize,"kbp.cn.segmented.pdf")))
plot(copyNumbersSegmented)
dev.off()

#call aberrations, and plot the final results.
copyNumbersCalled <- callBins(copyNumbersSegmented)

pdf(file.path(outdir,paste0(name, ".QDNAseq.",binSize,"kbp.cn.called.pdf")))
plot(copyNumbersCalled)
dev.off()

pdf(file.path(outdir,paste0(name, ".QDNAseq.",binSize,"kbp.cn.frequency.pdf")))
frequencyPlot(copyNumbersCalled)
dev.off()


#Called data can be exported as VCF file or SEG for further downstream analysis.
exportBins(copyNumbersCalled, file=file.path(outdir,paste0(name, ".QDNAseq.",binSize,"kbp.vcf")), format = "vcf")
exportBins(copyNumbersCalled, file=file.path(outdir,paste0(name, ".QDNAseq.",binSize,"kbp.seg")), format="seg")
exportBins(copyNumbersCalled, file=file.path(outdir,paste0(name, ".QDNAseq.",binSize,"kbp.bed")), format="bed")
