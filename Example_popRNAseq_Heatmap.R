#Heatmap for popRNAseq Data
library('ggplot2')
library('RColorBrewer')
library('dplyr')
pal <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

# rnaseq is the matrix with normalized counts to use to generate the Heatmap
rnaseq.data <- read.csv('Normalized_counts_freshy_isolated.csv')
rnaseq <- as.matrix(rnaseq.data)

# popRNAseq.DGE contains all genes with DGE with an adj_pValue less than 0.05
popRNAseq.DGE <- read.csv('CD66P_vs_CD66N_adj_pValue_less_than_0.05.csv')
genes.dge.popRNAseq <- as.character(popRNAseq.DGE$gene)
gene.keep <- row.names(rnaseq) %in% genes.dge.popRNAseq

# Filter the rnaseq matrix to keep only genes with sigificant Adj_pValue in DGE of popRNAseq analysis
mat.use <- rnaseq[gene.keep,]

# markers.spb.vs.mpb contains DGE of MPB vs SPB with an adj_pValue less than 0.05
markers.spb.vs.mpb <- read.csv('markers.spb.vs.mpb_adj.pVal_less_than_005.csv')
genes.dge.spb.vs.mpb <- as.character(markers.spb.vs.mpb$gene)

# Create a Matrix of normalized counts with genes from DGE of popRNAseq that are common with genes of DGE of MPB and SPB cells
gene.sc <- row.names(mat.use) %in% genes.dge.spb.vs.mpb
mat.use2 <- mat.use[gene.sc,]

# popRNAseq.DGE.filtered contains all genes with DGE with an adj_pValue less than 0.05 and that are common with the DGE of MPB and SPB cells
popRNAseq.DGE.filtered <- popRNAseq.DGE[row.names(mat.use2),]
top <- popRNAseq.DGE.filtered %>% top_n(20,log2fc)
bottom <- popRNAseq.DGE.filtered %>% top_n(-20,log2fc)
minmax <- bind_rows(top,bottom) 
minmax.ordered <- minmax[with(minmax, order(-log2fc)), ] 
to.plot <- as.character(minmax.ordered$gene)
mat.heatmap <- mat.use2[to.plot,]
heatmap.2(mat.heatmap,trace="none",main="Top 40 most variable genes up and down",col=pal)





