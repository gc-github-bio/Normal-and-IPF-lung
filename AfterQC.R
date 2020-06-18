library('Seurat')
library('ggplot2')
library('dplyr')

# Import datasets after qc steps
# 'Normal' datasets 1-6, IPF datasets 7-13

dataDir = c('/media/gc/GC_2T/qc/dataset1/',
            '/media/gc/GC_2T/qc/dataset2/',
            '/media/gc/GC_2T/qc/dataset3/',
            '/media/gc/GC_2T/qc/dataset4/',
            '/media/gc/GC_2T/qc/dataset5/',
            '/media/gc/GC_2T/qc/dataset6/',
            '/media/gc/GC_2T/qc/dataset7/',
            '/media/gc/GC_2T/qc/dataset8/',
            '/media/gc/GC_2T/qc/dataset9/',
            '/media/gc/GC_2T/qc/dataset10/',
            '/media/gc/GC_2T/qc/dataset11/',
            '/media/gc/GC_2T/qc/dataset12/',
            '/media/gc/GC_2T/qc/dataset13/')


seu1 <- readRDS(dataDir[1])
seu2 <- readRDS(dataDir[2])
seu3 <- readRDS(dataDir[3])
seu4 <- readRDS(dataDir[4])
seu5 <- readRDS(dataDir[5])
seu6 <- readRDS(dataDir[6])
seu7 <- readRDS(dataDir[7])
seu8 <- readRDS(dataDir[8])
seu9 <- readRDS(dataDir[9])
seu10 <- readRDS(dataDir[10])
seu11 <- readRDS(dataDir[11])
seu12 <- readRDS(dataDir[12])
seu13 <- readRDS(dataDir[13])

#set setwd() to desired project: "normal", "IPF" , or "all" 

#Integration of QC processed normal datasets
#merge all normal datasets and set up an object list
seu <- merge(x= seu1, y = c(seu2, seu3, seu4, seu5, seu6),project = "Normal")
seu.list <- SplitObject(seu, split.by = "id")

#Integration of QC processed IPF datasets
#merge all IPF datasets and set up an object list (Note: seu obj name need to be changed if processed at the same time with Normal)
seu <- merge(x= seu7, y = c(seu8, seu9, seu10, seu11, seu12, seu13),project = "IPF")
seu.list <- SplitObject(seu, split.by = "id")

#For all datasets, normal and IPF together (Note: seu obj name need to be changed if processed at the same time with previous)
seu <- merge(x= seu1, y = c(seu2, seu3, seu4, seu5, seu6, seu7, seu8, seu9, seu10, seu11, seu12, seu13),project = "All")
seu.list <- SplitObject(seu, split.by = "id")

#Run SCTranform on each object separately. 
library('future')
plan("multiprocess", workers = 10)

for (i in 1:length(seu.list)) {
  seu.list[[i]] <- SCTransform(seu.list[[i]], vars.to.regress = c("percent.mt","percent.rp"),variable.features.n = 5000,verbose = FALSE)
}

seu.features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = 5000)
seu.list <- PrepSCTIntegration(object.list = seu.list, anchor.features = seu.features, verbose = FALSE)
seu.anchors <- FindIntegrationAnchors(object.list = seu.list, normalization.method = "SCT", anchor.features = seu.features, verbose = FALSE)
seu.integrated <- IntegrateData(anchorset = seu.anchors, normalization.method = "SCT", verbose = FALSE)
DefaultAssay(seu.integrated) <- 'integrated'
seu.integrated <- RunPCA(seu.integrated, verbose = FALSE)

#Verify appropriate number of PC to use
ElbowPlot(seu.integrated, ndims = 50)
seu.integrated <- RunUMAP(seu.integrated, dims = 1:30, verbose = FALSE)
seu.integrated <- FindNeighbors(seu.integrated, dims = 1:30, force.recalc = TRUE, verbose = FALSE)
seu.integrated <- FindClusters(seu.integrated, resolution = 3, algorithm = 4, verbose = FALSE)

#DEG for all clusters
Idents(seu.integrated) <- 'seurat_clusters'
all.markers <- FindAllMarkers(seu.integrated, test.use = "MAST")
write.csv(all.markers,'all.markers.csv')

#Assign clusters using published signatures
Jaffe.Signature <- readRDS('rds/JaffeKleinSignature.rds')

# Create a score for previously published cell types: Basal, Secretory, Ciliated, Ionocyte
seu.integrated <- AddModuleScore(
  object = seu.integrated,
  features = Jaffe.Signature$'.',
  ctrl = 50,
  name = '.')

#Verify which cluster is identified by the signature
DefaultAssay(seu.integrated) <- 'RNA'
seu.integrated <- NormalizeData(seu.integrated, verbose = FALSE)
p1 <- FeaturePlot(object = seu.integrated, features = c("."), min.cutoff = "q10", max.cutoff = "q90")
p2 <- DimPlot(object = seu.integrated, label = TRUE)
CombinePlots(plots = list(p1,p2),ncol = 1)

#Assign clusters to specific cell types
seu.integrated <- RenameIdents(seu.integrated,
                               '1' = 'j.bas', '2' = '2', '3' = '3', '4' = 'j.sec', '5' = 'j.sec',
                               '6' = 'j.cil','7' = '7','8' = 'j.bas', '9' = '9','10' = '10',
                               '11' = '11','12' = 'j.cil','13' = 'j.bas','14' = 'j.cil', '15' = 'j.cil',
                               '16' = 'j.bas','17' = 'j.cil','18' = '18','19' = '19', '20' = 'j.bas',
                               '21' = 'j.bas','22' = 'j.cil','23' = 'j.sec','24' = '24', '25' = '25',
                               '26' = '26','27' = '27','28' = '28','29' = 'j.sec', '30' = 'j.cil',
                               '31' = '31','32' = '32','33' = '33','34' = 'j.bas', '35' = '35',
                               '36' = '36','37' = 'j.ion')

seu.integrated$j.first.iter <- Idents(seu.integrated)

#Find markers for 1st iteration 
j.sec.markers <- FindMarkers(seu.integrated, ident.1 = "j.sec", ident.2 = NULL, test.use = "MAST")
j.bas.markers <- FindMarkers(seu.integrated, ident.1 = "j.bas", ident.2 = NULL, test.use = "MAST")
j.cil.markers <- FindMarkers(seu.integrated, ident.1 = "j.cil", ident.2 = NULL, test.use = "MAST")
j.ion.markers <- FindMarkers(seu.integrated, ident.1 = "j.ion", ident.2 = NULL, test.use = "MAST")

#save tables
write.csv(j.sec.markers,'j.sec.markers.csv')
write.csv(j.bas.markers,'j.bas.markers.csv')
write.csv(j.cil.markers,'j.cil.markers.csv')
write.csv(j.ion.markers,'j.ion.markers.csv')

# Create a score for each of the cell types, Basal, Secretory, Ciliated, and Ionocyte after first iteration
top <- j.sec.markers %>% top_n(25,avg_logFC)
sign <- list(row.names(top))

seu.integrated <- AddModuleScore(
  object = seu.integrated,
  features = sign,
  ctrl = 50,
  name = 'sign_')
#Repeat for ech cell type

#Verify cell assignment
p1 <- FeaturePlot(object = seu.integrated, features = c("sign_"), min.cutoff = "q10",max.cutoff = "q90", pt.size = 1)
p2 <- DimPlot(object = seu.integrated, label = TRUE)
CombinePlots(plots = list(p1,p2),ncol = 1)

#Repeat iteration as needed

#Assign clusters to specific cell types after final iteration
Idents(seu.integrated) <- 'seurat_clusters'
seu.integrated <- RenameIdents(seu.integrated,
                               '1' = 'Basal', '2' = 'Ciliated', '3' = 'Basal', '4' = 'Secretory', '5' = 'Secretory',
                               '6' = 'Ciliated','7' = 'Basal','8' = 'Basal', '9' = 'Basal','10' = 'Secretory',
                               '11' = 'Secretory','12' = 'Ciliated','13' = 'Basal','14' = 'Ciliated', '15' = 'Ciliated',
                               '16' = 'Basal','17' = 'Ciliated','18' = 'Secretory','19' = 'Ciliated', '20' = 'Basal',
                               '21' = 'Basal','22' = 'Ciliated','23' = 'Secretory','24' = 'Secretory', '25' = 'Basal',
                               '26' = 'Basal','27' = 'Basal','28' = 'Secretory','29' = 'Secretory', '30' = 'Ciliated',
                               '31' = 'Secretory','32' = 'Basal','33' = 'Ciliated','34' = 'Basal', '35' = 'Secretory',
                               '36' = 'Secretory','37' = 'Ionocyte')

seu.integrated$major <- Idents(seu.integrated)

#Repeat iteration within major cell types to define cell type subclustering

#Assign cell types subclustering
Idents(seu.integrated) <- 'seurat_clusters'
seu.integrated <- RenameIdents(seu.integrated,
                               '1' = 'MPB', '2' = 'Ciliated2', '3' = 'MPB', '4' = 'Goblet', '5' = 'Goblet',
                               '6' = 'Ciliated2','7' = 'SPB','8' = 'MPB', '9' = 'AB','10' = 'Club',
                               '11' = 'Goblet','12' = 'Ciliated1','13' = 'PB','14' = 'Ciliated1', '15' = 'Ciliated1',
                               '16' = 'SPB','17' = 'Ciliated1','18' = 'Serous2','19' = 'Ciliated1', '20' = 'MPB',
                               '21' = 'MPB','22' = 'Ciliated1','23' = 'Club','24' = 'Club', '25' = 'MPB',
                               '26' = 'MPB','27' = 'AB','28' = 'Mucous','29' = 'Club', '30' = 'Ciliated1',
                               '31' = 'Serous1','32' = 'AB','33' = 'Ciliated1','34' = 'PB', '35' = 'Club',
                               '36' = 'Club','37' = 'Ionocyte')

seu.integrated$minor <- Idents(seu.integrated)

#Fig 1A
DimPlot(object = seu.integrated, label = TRUE)

#Create separe objects for basal, secretory, and ciliated cell types

seu.integrated.bas <- subset(x = seu.integrated, idents = c("Basal"), invert = FALSE)
seu.integrated.sec <- subset(x = seu.integrated, idents = c("Secretory"), invert = FALSE)
seu.integrated.cil <- subset(x = seu.integrated, idents = c("Ciliated"), invert = FALSE)

#DEG

Idents(seu.integrated) <- 'major'

markers.Basal <- FindMarkers(seu.integrated, ident.1 = "Basal", ident.2 = NULL, test.use = 'MAST')
markers.Secretory <- FindMarkers(seu.integrated, ident.1 = "Secretory", ident.2 = NULL, test.use = 'MAST')
markers.Ciliated <- FindMarkers(seu.integrated, ident.1 = "Ciliated", ident.2 = NULL, test.use = 'MAST')
markers.Ionocyte <- FindMarkers(seu.integrated, ident.1 = "Ionocyte", ident.2 = NULL, test.use = 'MAST')

write.csv(markers.Basal,'markers.Basal.csv')
write.csv(markers.Secretory,'markers.Secretory.csv')
write.csv(markers.Ciliated,'markers.Ciliated.csv')
write.csv(markers.Ionocyte,'markers.Ionocyte.csv')

Idents(seu.integrated) <- 'minor'

markers.MPB <- FindMarkers(seu.integrated, ident.1 = "MPB", ident.2 = NULL, test.use = 'MAST')
markers.PB <- FindMarkers(seu.integrated, ident.1 = "PB", ident.2 = NULL, test.use = 'MAST')
markers.SPB <- FindMarkers(seu.integrated, ident.1 = "SPB", ident.2 = NULL, test.use = 'MAST')
markers.AB <- FindMarkers(seu.integrated, ident.1 = "AB", ident.2 = NULL, test.use = 'MAST')

write.csv(markers.MPB,'markers.MPB.csv')
write.csv(markers.PB,'markers.PB.csv')
write.csv(markers.SPB,'markers.SPB.csv')
write.csv(markers.AB,'markers.AB.csv')

markers.Club <- FindMarkers(seu.integrated, ident.1 = "Club", ident.2 = NULL, test.use = 'MAST')
markers.Goblet <- FindMarkers(seu.integrated, ident.1 = "Goblet", ident.2 = NULL, test.use = 'MAST')
markers.Mucous <- FindMarkers(seu.integrated, ident.1 = "Mucous", ident.2 = NULL, test.use = 'MAST')
markers.Serous1 <- FindMarkers(seu.integrated, ident.1 = "Serous1", ident.2 = NULL, test.use = 'MAST')
markers.Serous2 <- FindMarkers(seu.integrated, ident.1 = "Serous2", ident.2 = NULL, test.use = 'MAST')

write.csv(markers.Club,'markers.Club.csv')
write.csv(markers.Goblet,'markers.Goblet.csv')
write.csv(markers.Mucous,'markers.Mucous.csv')
write.csv(markers.Serous1,'markers.Serous1.csv')
write.csv(markers.Serous2,'markers.Serous2.csv')

markers.Ciliated1 <- FindMarkers(seu.integrated, ident.1 = "Ciliated1", ident.2 = NULL, test.use = 'MAST')
markers.Ciliated2 <- FindMarkers(seu.integrated, ident.1 = "Ciliated2", ident.2 = NULL, test.use = 'MAST')

write.csv(markers.Ciliated1,'markers.Ciliated1.csv')
write.csv(markers.Ciliated2,'markers.Ciliated2.csv')


#within

Idents(seu.integrated) <- 'minor'
minor.markers <- FindAllMarkers(seu.integrated, test.use = "MAST")
write.csv(minor.markers,'minor.markers.csv')

Idents(seu.integrated.bas) <- 'minor'
minor.bas.markers <- FindAllMarkers(seu.integrated.bas, test.use = "MAST")
write.csv(minor.bas.markers,'minor.bas.markers.csv')

Idents(seu.integrated.sec) <- 'minor'
minor.sec.markers <- FindAllMarkers(seu.integrated.sec, test.use = "MAST")
write.csv(minor.sec.markers,'minor.sec.markers.csv')

Idents(seu.integrated.cil) <- 'minor'
minor.cil.markers <- FindAllMarkers(seu.integrated.cil, test.use = "MAST")
write.csv(minor.cil.markers,'minor.cil.markers.csv')


#Fig 1B
gene.plot <- c("KRT15","KRT17","S100A2","KRT5","MMP10","KRT14","KRT6A","AREG","PHLDA1","SFN","DST","ADIRF","BCAM","FHL2","MYC","ATP1B3","IFITM3",
               "SCGB1A1","MSMB","TFF3","MUC5B","BPIFB1","BPIFA1","PI3","SCGB3A1","LYZ","PRB3","PRB4","PRR4","LTF","ZG16B","AZGP1","AC020656.1","PIP","PRH2",
               "CAPS","TMEM190","TPPP3","C9orf24","C20orf85","RSPH1","OMG","FAM183A","C1orf194","C9orf116","TUBA1A","MORN2","ODF3B","PIFO","AL357093.2","SNTN","CETN2","DNAAF1","C5orf49","CCDC78",
               "RARRES2","ASCL3","ID3","TMEM61","SCNN1B","AKR1B1","BPIFA2","ATP6V1G3","CLCNKB","HES6","HEPACAM2","FAM24B","FOXI1","PCP4")

seu.integrated$id.type <- paste(Idents(seu.integrated), seu.integrated$id, sep = "_")
Idents(seu.integrated) <- "id.type"

#adjust levels if necessary

DoHeatmap(seu.integrated, features = gene.plot ,assay = "RNA", slot ="scale.data",
          disp.min = -2.5,draw.lines = FALSE,raster = FALSE,
          disp.max = NULL) + scale_fill_gradientn(colors = colorRampPalette(c("#ffffff", 
                                                                              "#bfbfbf", 
                                                                              "#000000"))(100)) + NoLegend()

#Fig1C-G
fig1.sign <- readRDS('/fig1.sign.RDS')

#add score for each signature (hsa04310, KW-0346, GO:0051301, GO:0007050, HSA-3700989)

seu.integrated <- AddModuleScore(
  object = seu.integrated,
  features = fig1.sign$.,
  ctrl = 50,
  name = 'sign_')

#plot each signature
FeaturePlot(object = seu.integrated, features = c("sign_"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")

#Fig 1H, J, L
gene.plot.bas <- c("MMP10","CTGF","PHLDA1","CYR61","TNFRSF12A","DST","DKK1","STK17A","TSLP","KRT17","BCAM","C1orf56","POSTN","G0S2","MMP13","RGCC","CAV1","ARID5B","FST",
                   "KRT15","TUBA1B","HIST1H4C","STMN1","PTTG1","UBE2C","PCLAF","TUBB","TK1","CENPW","CDC20","TOP2A","BIRC5","MKI67","SERPINB3","CLDN4","KRT16","PLAC8",
                   "PLAUR","SLPI","GDF15","S100P","SERPINB1","RAB11FIP1","FAM3D","GSN","CXCL17","LMO7","WFDC2","MDK","APOBEC3A","FOS","JUN","ZFP36","JUNB","ID1","IER2",
                   "ATF3","DUSP1","EGR1","GADD45B","HSPA1A","HSPA1B","NFKBIA","FOSB","HES1","DDIT4","SOCS3")

Idents(seu.integrated.bas) <- 'minor'
seu.integrated.bas$id.basal <- paste(Idents(seu.integrated.bas), seu.integrated.bas$id, sep = "_")
Idents(seu.integrated.bas) <- "id.basal"

#adjust levels if necessary

DoHeatmap(seu.integrated.bas, features = gene.plot.bas, assay = "RNA", slot ="scale.data",
          disp.min = -2.5,draw.lines = FALSE,raster = FALSE,
          disp.max = NULL) + scale_fill_gradientn(colors = colorRampPalette(c("#ffffff", 
                                                                              "#bfbfbf", 
                                                                              "#000000"))(100)) + NoLegend()


gene.plot.sec <- c("KRT17","SFN","S100A2","AQP3","SERPINB1","TACSTD2","ALDH3A1","KRT6A","SLC25A5","HMGA1","PLAC8","SERPINB4","GLUL","GPX2","S100A10","EMP1","SCGB1A1","VMO1",
                   "BPIFA1","MUC5AC","TFF3","MUC5B","BPIFB2","GLYATL2","AZGP1","TCN1","LYZ","PIP","LTF","AC020656.1","DMBT1","C6orf58","HP","CHI3L1","PRB4","PRB3","SCGB3A2","S100A1")

Idents(seu.integrated.sec) <- 'minor'
seu.integrated.sec$id.sec <- paste(Idents(seu.integrated.sec), seu.integrated.sec$id, sep = "_")
Idents(seu.integrated.sec) <- "id.sec"

#adjust levels if necessary

DoHeatmap(seu.integrated.sec, features = gene.plot.sec, assay = "RNA", slot ="scale.data",
          disp.min = -2.5,draw.lines = FALSE,raster = FALSE,
          disp.max = NULL) + scale_fill_gradientn(colors = colorRampPalette(c("#ffffff", 
                                                                              "#bfbfbf", 
                                                                              "#000000"))(100)) + NoLegend()

gene.plot.cil <- c("CD59","SARAF","PTGES3","TMEM59","CD24","C12orf75","B2M","PFN2","CALM2","SUMO2","H3F3A","FTH1","CABCOCO1","DSTN","FTL","ALDH1A1","SRP9","CETN2","OMG",
                   "LRRIQ1","MALAT1","WDR60","PTRH1","CFAP157","CCDC189","DRC3","CCDC17","AC009502.1","TTLL10","DNAH12","TOGARAM2","JUN","HES1","ANKRD18A","AHI1","VWA3A",
                   "LRRC74B","HIST1H4C","ADPRHL2","HEXIM1","CSKMT","AC023157.3","IER2","PLCG2","SOX4","AL355075.4","SERPINF1","S100A2")

Idents(seu.integrated.cil) <- 'minor'
seu.integrated.cil$id.cil <- paste(Idents(seu.integrated.cil), seu.integrated.cil$id, sep = "_")
Idents(seu.integrated.cil) <- "id.cil"

#adjust levels if necessary

DoHeatmap(seu.integrated.cil, features = gene.plot.cil, assay = "RNA", slot ="scale.data", 
          disp.min = -2.5,draw.lines = FALSE,raster = FALSE,
          disp.max = NULL) + scale_fill_gradientn(colors = colorRampPalette(c("#ffffff", 
                                                                              "#bfbfbf", 
                                                                              "#000000"))(100)) + NoLegend()

#FigI, K, M

fig1.vln.sign <- readRDS('/fig1.vln.sign.RDS')

#add score for each signature (MPB, PB, SPB, AB, Club, Goblet, Mucous, Serous1, Serous2, CIliatd1, Ciliated2)

seu.integrated.bas <- AddModuleScore(
  object = seu.integrated.bas,
  features = fig1.vln.sign$.,
  ctrl = 50,
  name = 'sign_')

seu.integrated.sec <- AddModuleScore(
  object = seu.integrated.sec,
  features = fig1.vln.sign$.,
  ctrl = 50,
  name = 'sign_')

seu.integrated.cil <- AddModuleScore(
  object = seu.integrated.cil,
  features = fig1.vln.sign$.,
  ctrl = 50,
  name = 'sign_')

#Violin plot for each signature 

VlnPlot(object=seu.integrated.bas,features=c("sign_"), pt.size = 0, combine =TRUE)+
  geom_boxplot(width=0.1,fill="gray", alpha=0.6)

VlnPlot(object=seu.integrated.sec,features=c("sign_"), pt.size = 0, combine =TRUE)+
  geom_boxplot(width=0.1,fill="gray", alpha=0.6)

VlnPlot(object=seu.integrated.cil,features=c("sign_"), pt.size = 0, combine =TRUE)+
  geom_boxplot(width=0.1,fill="gray", alpha=0.6)


#Fig2A, G-J

fig2.umap.sign <- readRDS('/fig2.umap.sign.RDS') # o-glycan biosynthesis and CD66 transcripts

seu.integrated <- AddModuleScore(
  object = seu.integrated,
  features = fig2.umap.sign$.,
  ctrl = 50,
  name = 'sign_')

FeaturePlot(object = seu.integrated, features = c("sign_oGlyc", "sign_cd66", "CEACAM1", "CEACAM5", "CACAM6"), min.cutoff = "q10", max.cutoff = "q90")

#Fig2B
#Diffusion Map

library('destiny')
library('slingshot')
library('SingleCellExperiment')
library('RColorBrewer')
library('rgl')

seu <- readRDS('/normal.datasets.rds')
seu.destiny <- subset(x = seu, idents = c("PB","AB","MPB","SPB","club_goblet"), invert = FALSE)
colourCount = length(unique(Idents(seu.destiny)))
pal <- getPalette(colourCount)
palette(pal)

# prepare for Destiny and Slingshot
object.sce <- as.SingleCellExperiment(seu.destiny)
object.matrix <- as.matrix(logcounts(object.sce))
use.genes <- VariableFeatures(seu.destiny)
object.matrix <- object.matrix[use.genes,]
object.dataframe <- as.data.frame(t(object.matrix))
object.dataframe$subtype <- Idents(seu.destiny)
ct<-as.ExpressionSet(object.dataframe)

# Prepare Diffusion Map
dif <- DiffusionMap(ct, k = 2000)
plot(dif, col_by ='subtype')

#plot 3D
rd3 <- cbind(coord1 = dif$coord1, coord2 = dif$coord2, coord3 = dif$coord3)
reducedDims(object.sce) <- SimpleList(DiffMap3D = rd3)
traj <- slingshot(reducedDims(object.sce)$DiffMap3D, clusterLabels = Idents(seu.destiny))
plot3d(reducedDims(object.sce)$DiffMap3D, col = pal[Idents(seu.destiny)], aspect = 'iso', alpha = 0.3)
plot3d(traj, lwd = 3, add = TRUE)
rgl.postscript("destiny_Normal_plot.pdf",fmt= "pdf") 
rgl.close()


#Fig2F

gene.plot.CD <- c("BCAM","TNFRSF12A","ITGB1","PRNP","IFITM1","NGFR","ATP1B3","F3","HMMR","BSG","PLAUR","ALCAM","CEACAM5","CDH1","CEACAM6","F11R","MUC1","SDC1","ITGA2","CD55","CD24")

Idents(seu.integrated.bas) <- 'minor'

DoHeatmap(seu.integrated.bas, features = gene.plot.CD, assay = "RNA", slot ="scale.data",
          disp.min = -2.5,draw.lines = FALSE,raster = FALSE,
          disp.max = NULL) 


#FIG3C, (D same way)

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

# popRNAseq.DGE contains all genes with DGE with an adj_pValue less than 0.05

top <- popRNAseq.DGE %>% top_n(20,log2fc)
bottom <- popRNAseq.DGE %>% top_n(-20,log2fc)
minmax <- bind_rows(top,bottom) 
minmax.ordered <- minmax[with(minmax, order(-log2fc)), ] 
to.plot <- as.character(minmax.ordered$gene)
mat.heatmap <- mat.use[to.plot,]
heatmap.2(mat.heatmap,trace="none",main="Top 40 most variable genes up and down",col=pal)


#FIGS4A, (B same way)

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


#Fig4A
# UMAP of IPF integrated datasets

Idents(seu.integrated) <- 'minor'
DimPlot(object = seu.integrated, label = TRUE)

#Fig4B
Idents(seu.integrated) <- 'major'
seu.integrated$id.all <- paste(Idents(seu.integrated), seu.integrated$id, sep = "_")
Idents(seu.integrated) <- "id.all"

gene.plot.all <- (c("S100A2","KRT17","KRT15","KRT5","KRT6A","SFN","AQP3","S100A9","ALDH3A1","KRT19",
                "BPIFA1","MSMB","LTF","SCGB1A1","SCGB3A1","BPIFB1","SLPI","PIGR","CXCL17","CEACAM6",
                "TPPP3","CAPS","TMEM190","C9orf24","C20orf85",
                "C1orf194","FAM183A","RSPH1","CETN2","PIFO",
                "ASCL3","CLCNKB","ATP6V1G3","HEPACAM2","FOXI1","SCG2","STAP1","CEL","IGF1","CLNK"))

cols = c("#99d8c9","#66c2a4","#41ae76","#238b45","#006d2c","#00441b","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026",
         "#99d8c9","#66c2a4","#41ae76","#238b45","#006d2c","#00441b","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026",
         "#99d8c9","#66c2a4","#41ae76","#238b45","#006d2c","#00441b","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026",
         "#99d8c9","#66c2a4","#41ae76","#238b45","#006d2c","#00441b","#fed976","#feb24c","#fc4e2a","#bd0026")

#adjust levels if necessary

DoHeatmap(seu2, features = gene.plot.all, assay = "RNA", slot ="scale.data", group.colors =cols,
          disp.min = -2.5,draw.lines = FALSE,raster = FALSE,
          disp.max = NULL) + scale_fill_gradientn(colors = colorRampPalette(c("#ffffff", 
                                                                              "#bfbfbf", 
                                                                              "#000000"))(100)) + NoLegend()

#Fig4C-F

fig4.UMAP.sign <- readRDS('/fig4.UMAP.sign.RDS')

#add score for CD66 transcripts

seu.integrated <- AddModuleScore(
  object = seu.integrated,
  features = fig4.UMAP.sign$.,
  ctrl = 50,
  name = '.')

FeaturePlot(object = seu.integrated, features = c("CEACAM1", "CEACAM5", "CEACAM6", "sign_cd66"), min.cutoff = "q10", max.cutoff = "q90")

#Fig4G
#Use Normal + IPF integrated object

Idents(seu.integrated) <- 'minor'

all.spb <- subset(x = seu.integrated, idents = c("SPB"), invert = FALSE)

markers.to.plot <- c("KRT5","CEACAM6", "CEACAM5", "CEACAM1", "SERPINB3", "S100A8", "S100A9", "RACK1", "LCN2", "CSTA", "S100P", "SERPINB4",
                     "MMP7", "MMP1", "SFTPB", "IGFBP7", "FOS","JUN", "SCGB3A1", "JUNB")

# type: _CO and _IPF
DotPlot(all.spb, features = rev(markers.to.plot), cols = c("#3db85b","#b30269"), dot.min = 0.1,dot.scale = 15,  split.by = 'type') + RotatedAxis()

#Fig4H 

seu <- readRDS('/ipf.datasets.rds')
seu.destiny <- subset(x = seu, idents = c("AB","MPB","SPB","SCGB3A2+","SCGB3A1+/3A2+", "Club", "Goblet"), invert = FALSE)
colourCount = length(unique(Idents(seu.destiny)))
pal <- getPalette(colourCount)
palette(pal)

# prepare for Destiny and Slingshot
object.sce <- as.SingleCellExperiment(seu.destiny)
object.matrix <- as.matrix(logcounts(object.sce))
use.genes <- VariableFeatures(seu.destiny)
object.matrix <- object.matrix[use.genes,]
object.dataframe <- as.data.frame(t(object.matrix))
object.dataframe$subtype <- Idents(seu.destiny)
ct<-as.ExpressionSet(object.dataframe)

# Prepare Diffusion Map
dif <- DiffusionMap(ct, k = 2000)
plot(dif, col_by ='subtype')

#plot 3D
rd3 <- cbind(coord1 = dif$coord1, coord2 = dif$coord2, coord3 = dif$coord3)
reducedDims(object.sce) <- SimpleList(DiffMap3D = rd3)
traj <- slingshot(reducedDims(object.sce)$DiffMap3D, clusterLabels = Idents(seu.destiny))
plot3d(reducedDims(object.sce)$DiffMap3D, col = pal[Idents(seu.destiny)], aspect = 'iso', alpha = 0.3)
plot3d(traj, lwd = 3, add = TRUE)
rgl.postscript("destiny_IPF_plot.pdf",fmt= "pdf") 
rgl.close()


#Fig5A, B
#Analysis NOTCH signaling
library(iTALK)

#Import L-R Database for Notch signaling
pathway.genes <- as.data.frame(read.table("database.txt", header = TRUE))
dataitalk <- as.data.frame(t(as.matrix(read.csv('basal.cells.csv', header = TRUE))))

Index.ligands <- database$Ligand.ApprovedSymbol %in% pathway.genes$HALLMARK_PATHWAY
Index.Receptors <- database$Receptor.ApprovedSymbol %in% pathway.genes$HALLMARK_PATHWAY
database$pathwayLigand <- Index.ligands
database$pathwayReceptor <- Index.Receptors
database.pathway <- database %>% filter(pathwayLigand == TRUE | pathwayReceptor == TRUE)

#top highly expressed genes
highly_exprs_genes<-rawParse(dataitalk,top_genes=20,stats='mean')
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_col<-structure(c('#f0ee95','#9df79c','#b9d3ed','#f5e105'),names=unique(dataitalk$cell_type))
par(mfrow=c(1,2))
res<-NULL

for(comm_type in comm_list){
  
  res_cat<-FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm_type, database = database.pathway)
  res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=F),]
  NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=0.5,edge.max.width=6)
  LRPlot(res_cat[1:50,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:50],link.arr.width=res_cat$cell_to_mean_exprs[1:50])
  title(comm_type)
  res<-rbind(res,res_cat)
}

res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:50,]
par(mfrow=c(1,2))
NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=0.5,edge.max.width=6)
LRPlot(res[1:50,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res$cell_from_mean_exprs[1:50],link.arr.width=res$cell_to_mean_exprs[1:50])







