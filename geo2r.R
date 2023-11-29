# Geo2R from website

# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Data plots for selected GEO samples

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE239282", "file=GSE239282_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# pre-filter low count genes
# keep genes with at least 2 counts > 10
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]

# log transform raw counts
# instead of raw counts can display vst(as.matrix(tbl)) i.e. variance stabilized counts
dat <- log10(tbl + 1)

# box-and-whisker plot
dev.new(width=3+ncol(tbl)/6, height=5)
par(mar=c(7,4,2,1))
boxplot(dat, boxwex=0.7, notch=T, main="GSE239282", ylab="lg(cnt + 1)", outline=F, las=2)
dev.off()

# UMAP plot (dimensionality reduction)
library(umap)
dat <- dat[!duplicated(dat), ] # first remove duplicates
ump <- umap(t(dat), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="GSE239282 UMAP plot, nbrs =15", xlab="", ylab="", pch=20, cex=1.5)
library(car)
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

