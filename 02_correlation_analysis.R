# Load data
expr <- readRDS("expr_all.rds")
expr_hacat <- readRDS("expr_hacat.rds")
expr_tigk <- readRDS("expr_tigk.rds")
counts <- readRDS("counts.rds")
lnc_ids <- readRDS("lnc_ids.rds")
lnc_labels <- readRDS("lnc_labels.rds")

# Target genes (final set used in paper)
targets <- c(
"ADAM12","LOX","COL1A1","FN1","COL5A1","FBN1","MMP2","SPARC","LAMA2",
"WNT4","PDGFRB","VEGFA","SMAD3","TGFB1","TGFB2",
"LAMC2","DSC2","DSG2","ITGB1","ITGA6",
"MCL1","BCL2","CASP8",
"TNFAIP3","REL","AHR",
"HDAC9","HDAC4","DNMT3A","TET2","SMARCC1"
)

# Match genes
targets_present <- targets[targets %in% counts$gene_name]
target_ids <- counts$gene_id[match(targets_present, counts$gene_name)]

# Correlation function
make_cor <- function(expr_mat){
  mat <- matrix(NA, nrow=length(target_ids), ncol=length(lnc_ids))
  for(i in 1:length(target_ids)){
    for(j in 1:length(lnc_ids)){
      mat[i,j] <- cor(
        expr_mat[lnc_ids[j],],
        expr_mat[target_ids[i],],
        use="pairwise.complete.obs",
        method="pearson"
      )
    }
  }
  rownames(mat) <- targets_present
  colnames(mat) <- lnc_labels
  mat
}

cor_all <- make_cor(expr)
cor_hacat <- make_cor(expr_hacat)
cor_tigk <- make_cor(expr_tigk)

saveRDS(cor_all,"cor_all.rds")
saveRDS(cor_hacat,"cor_hacat.rds")
saveRDS(cor_tigk,"cor_tigk.rds")
