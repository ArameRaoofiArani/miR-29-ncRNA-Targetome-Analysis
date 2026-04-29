
```r
library(DESeq2)

# Read data
counts <- read.delim("GSE270877_Raw_Counts.txt", header=TRUE, sep="\t")

# Candidate lncRNAs
lnc_ids <- c(
  "ENSG00000130600", # H19
  "ENSG00000233016", # SNHG7
  "ENSG00000227036"  # LINC00511
)

lnc_labels <- c("H19","SNHG7","LINC00511")

# Expression matrix
expr_counts <- counts[,2:19]
rownames(expr_counts) <- counts$gene_id
expr_counts <- round(as.matrix(expr_counts))

# Metadata
samples <- colnames(expr_counts)

celltype <- ifelse(grepl("^T", samples), "TIGK", "HaCaT")
time <- ifelse(grepl("0", samples),"0h",
        ifelse(grepl("6", samples),"6h","24h"))

meta <- data.frame(
  row.names=samples,
  celltype=factor(celltype),
  time=factor(time)
)

# DESeq2
dds <- DESeqDataSetFromMatrix(expr_counts, meta, ~celltype + time)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
expr <- assay(vsd)

# Split
expr_hacat <- expr[, meta$celltype=="HaCaT"]
expr_tigk  <- expr[, meta$celltype=="TIGK"]

# Save
saveRDS(expr,"expr_all.rds")
saveRDS(expr_hacat,"expr_hacat.rds")
saveRDS(expr_tigk,"expr_tigk.rds")
saveRDS(counts,"counts.rds")
saveRDS(lnc_ids,"lnc_ids.rds")
saveRDS(lnc_labels,"lnc_labels.rds")
