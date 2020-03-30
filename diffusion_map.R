library('Seurat')
library('destiny')
library('slingshot')
library('SingleCellExperiment')
library('ggplot2')
library('RColorBrewer')
library('rgl')

seu <- readRDS('/normal.datasets.rds')
seu.destiny <- subset(x = seu, idents = c("PB","AB","MPB","SPB","club_goblet"), invert = FALSE)
colourCount = length(unique(Idents(seu.destiny)))
pal <- getPalette(colourCount)
palette(pal)

# prepare for Destiny and Slingshot
object.sce <- as.SingleCellExperiment(seu.destiny)
cellLabels <- Idents(seu.destiny)
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
rgl.postscript("destiny_plot.pdf",fmt= "pdf") 
rgl.close()












