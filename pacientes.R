# Creado por Omar Vargas Hernández
# Octubre 2024

setwd("~/Proyecto_INMEGEN")

# cargar paquete pacman
library( "pacman")

# cargar otros paquetes
p_load( "Seurat", "harmony", "scIntegrationMetrics", "openxlsx", "ggplot2")

# Leemos los datos de cada muestra

# Data de los pacientes
cell_pt17 <- Read10X( data.dir = "data/UMN-017-FMAC_ileon_cellrangercounts/outs/filtered_feature_bc_matrix/" )
cell_pt19 <- Read10X( data.dir = "data/UMN-019-ICC_ileon_cellrangercounts/outs/filtered_feature_bc_matrix/" )

# Creamos objetos seurat para cada uno de los pacientes
seurat_pt17 <- CreateSeuratObject( counts = cell_pt17,
                                   project = "pt17" )
seurat_pt19 <- CreateSeuratObject( counts = cell_pt19,
                                   project = "pt19" )

# Limpiamos objetos
rm( list =  ls( pattern = "cell_") )      # se borran los objetos
gc( )                                     # se limpia la memoria restante

# Filtrar datos seurat pt17

# Revisamos el porcentaje de conteos que corresponden a genes mitocondriales
seurat_pt17[[ "percent.mt" ]] <- PercentageFeatureSet( object = seurat_pt17,
                                                         pattern = "^MT-" )

# Procedemos a filtrar las celulas, de acuerdo a las metricas de control de calidad comunes

seurat_pt17_T <- subset( x = seurat_pt17,
                         subset = nFeature_RNA > 200 &
                           nFeature_RNA < 2500 &
                           percent.mt < 5)
seurat_pt17_T@project.name <- "pt17T"
levels(seurat_pt17_T$orig.ident) <- c("pt17T")
Idents(seurat_pt17_T) <- seurat_pt17_T$orig.ident

#para el otro paciente

# Revisamos el porcentaje de conteos que corresponden a genes mitocondriales
seurat_pt19[[ "percent.mt" ]] <- PercentageFeatureSet( object = seurat_pt19,
                                                         pattern = "^MT-" )

seurat_pt19_T <- subset( x = seurat_pt19,
                         subset = nFeature_RNA > 200 &
                           nFeature_RNA < 2500 &
                           percent.mt < 5)

seurat_pt19_T@project.name <- "pt19T"
levels(seurat_pt19_T$orig.ident) <- c("pt19T")
Idents(seurat_pt19_T) <- seurat_pt19_T$orig.ident

rm(seurat_pt17)
rm(seurat_pt19)

##Normalizamos

seurat_pt17_T <- NormalizeData(seurat_pt17_T, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_pt19_T <- NormalizeData(seurat_pt19_T, normalization.method = "LogNormalize", scale.factor = 10000)

#Detección de genes variables
seurat_pt17_T <- FindVariableFeatures(seurat_pt17_T, selection.method = "vst", nfeatures = 2000)
seurat_pt19_T <- FindVariableFeatures(seurat_pt19_T, selection.method = "vst", nfeatures = 2000)

# Combinar matrices de conteo de ambos pacientes
seurat_pt17y19_merged <- merge(seurat_pt17_T, y = seurat_pt19_T, 
                            add.cell.ids = c("S1", "S2"), 
                            project = "pt17y19_merged")

# Identificar genes variables
seurat_pt17y19_merged <- FindVariableFeatures(seurat_pt17y19_merged)

# Preprocesamiento antes de Harmony
seurat_pt17y19_merged <- ScaleData(seurat_pt17y19_merged)
seurat_pt17y19_merged <- RunPCA(seurat_pt17y19_merged)
seurat_pt17y19_merged <- FindNeighbors(seurat_pt17y19_merged)
seurat_pt17y19_merged <- RunUMAP(seurat_pt17y19_merged, dims = 1:30)
seurat_pt17y19_merged <- FindClusters(seurat_pt17y19_merged, resolution = 0.5)
DimPlot(seurat_pt17y19_merged, group.by = "orig.ident")

# Métricas ANTES de la integración 

integration_metrics0 <- getIntegrationMetrics(
  object = seurat_pt17y19_merged, 
  metrics = NULL,
  meta.label = "seurat_clusters",
  meta.batch = "orig.ident",
  method.reduction = "umap",
  ndim = 30 
)

print(integration_metrics0)

# Convertir la lista de métricas a un data frame
metrics_df0 <- as.data.frame(integration_metrics0)


# Integración con Harmony 
seurat_pt17y19_harmony <- RunHarmony(seurat_pt17y19_merged, 
                                  group.by.vars = "orig.ident", 
                                  plot_convergence = TRUE)

seurat_pt17y19_harmony <- RunUMAP(seurat_pt17y19_harmony, reduction = "harmony", dims = 1:30)
DimPlot(seurat_pt17y19_harmony, group.by = "orig.ident")

# Clustering DESPUÉS de Harmony
seurat_pt17y19_harmony <- FindNeighbors(seurat_pt17y19_harmony, reduction = "harmony")
seurat_pt17y19_harmony <- FindClusters(seurat_pt17y19_harmony, resolution = 0.5)

# Métricas DESPUÉS de la integración

integration_metrics1 <- getIntegrationMetrics(
  object = seurat_pt17y19_harmony, 
  metrics = NULL,
  meta.label = "seurat_clusters",
  meta.batch = "orig.ident",
  method.reduction = "umap",
  ndim = 30 
)

# Imprime los resultados
print(integration_metrics1)

# Convertir la lista de métricas a un data frame
metrics_df1 <- as.data.frame(integration_metrics1)

# Guardar los 2 data frames en un archivo XLSX
write.xlsx(rbind(metrics_df0,metrics_df1), file = "pt17ypt19-metricas.xlsx")
