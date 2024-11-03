# Creado por Omar Vargas Hernandez
# Octubre 2024
# Texto sin acentos

##########OBJETIVO##########
#Crear 2 subconjuntos de una sola matriz de conteos, para evaluar las metricas de integracion

setwd("~/Proyecto_INMEGEN")

# cargar paquete pacman
library( "pacman")

# cargar otros paquetes
p_load( "Seurat", "harmony", "scIntegrationMetrics", "openxlsx", "ggplot2")

# Leemos los datos del paciente

cell_pt17 <- Read10X( data.dir = "data/UMN-017-FMAC_ileon_cellrangercounts/outs/filtered_feature_bc_matrix/")

# Creamos el objeto seurat para el paciente
seurat_pt17 <- CreateSeuratObject( counts = cell_pt17, project = "pt17" )

# Limpiamos
rm( list =  ls( pattern = "cell_") )      # se borran los primeros objetos
gc( )                                     # se limpia la memoria restante

# Revisamos el porcentaje de conteos que corresponden a genes mitocondriales
seurat_pt17[[ "percent.mt" ]] <- PercentageFeatureSet( object = seurat_pt17, pattern = "^MT-" )

#Hacemos los subconjuntos
seurat_pt17$subset <- sample(c(1, 2), size = ncol(seurat_pt17), replace = TRUE)

# Procedemos a filtrar las celulas, de acuerdo a las metricas de control de calidad comunes

seurat_pt17T  <- subset( x = seurat_pt17, subset = nFeature_RNA > 200 &
                                                    nFeature_RNA < 2500 &
                                                    percent.mt < 5 )
seurat_pt17s1 <- subset( x = seurat_pt17, subset = nFeature_RNA > 200 &
                                                    nFeature_RNA < 2500 &
                                                    percent.mt < 5 &
                                                    subset == 1 )
seurat_pt17s2 <- subset( x = seurat_pt17, subset = nFeature_RNA > 200 &
                                                    nFeature_RNA < 2500 &
                                                    percent.mt < 5 &
                                                    subset == 2 )
#Arreglamos rotulos
seurat_pt17T@project.name  <- "pt17T"
seurat_pt17s1@project.name <- "pt17s1"
seurat_pt17s2@project.name <- "pt17s2"

levels(seurat_pt17T$orig.ident)  <- c("pt17T")
levels(seurat_pt17s1$orig.ident) <- c("pt17s1")
levels(seurat_pt17s2$orig.ident) <- c("pt17s2")

Idents(seurat_pt17T)  <- seurat_pt17T$orig.ident
Idents(seurat_pt17s1) <- seurat_pt17s1$orig.ident
Idents(seurat_pt17s2) <- seurat_pt17s2$orig.ident

# guardamos los objetos
#saveRDS( object = seurat_pt17,   file = "seurat_pt17.rds")

##Normalizamos

seurat_pt17s1 <- NormalizeData(seurat_pt17s1, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_pt17s2 <- NormalizeData(seurat_pt17s2, normalization.method = "LogNormalize", scale.factor = 10000)

#Detección de genes variables
seurat_pt17s1 <- FindVariableFeatures(seurat_pt17s1, selection.method = "vst", nfeatures = 2000)
seurat_pt17s2 <- FindVariableFeatures(seurat_pt17s2, selection.method = "vst", nfeatures = 2000)

# Combinar subconjuntos
seurat_pt17_merged <- merge(seurat_pt17s1, y = seurat_pt17s2, 
                               add.cell.ids = c("s1", "s2"), 
                               project = "pt17_merged")

# Identificar genes variables
seurat_pt17_merged <- FindVariableFeatures(seurat_pt17_merged)

# Preprocesamiento antes de Harmony
seurat_pt17_merged <- ScaleData(seurat_pt17_merged)
seurat_pt17_merged <- RunPCA(seurat_pt17_merged)
seurat_pt17_merged <- FindNeighbors(seurat_pt17_merged)
seurat_pt17_merged <- RunUMAP(seurat_pt17_merged, dims = 1:30)
seurat_pt17_merged <- FindClusters(seurat_pt17_merged, resolution = 0.5)
DimPlot(seurat_pt17_merged, group.by = "orig.ident") #graficos por separado
#split.by = "orig.ident"

# Métricas ANTES de la integración 

integration_metricsA <- getIntegrationMetrics(
  object = seurat_pt17_merged, 
  metrics = NULL,
  meta.label = "seurat_clusters",
  meta.batch = "orig.ident",
  method.reduction = "umap",
  ndim = 30 
)

print(integration_metricsA)

# Convertir la lista de métricas a un data frame
metrics_dfA <- as.data.frame(integration_metricsA)

# Integración con Harmony 
seurat_pt17_harmony <- RunHarmony(seurat_pt17_merged, 
                                     group.by.vars = "orig.ident", 
                                     plot_convergence = TRUE)

seurat_pt17_harmony <- RunUMAP(seurat_pt17_harmony, reduction = "harmony", dims = 1:30)
DimPlot(seurat_pt17_harmony, group.by = "orig.ident")

# Clustering DESPUÉS de Harmony
seurat_pt17_harmony <- FindNeighbors(seurat_pt17_harmony, reduction = "harmony")
seurat_pt17_harmony <- FindClusters(seurat_pt17_harmony, resolution = 0.5)

#Métricas de integración

# Métricas de integración
integration_metricsB <- getIntegrationMetrics(
  object = seurat_pt17_harmony, 
  metrics = NULL,
  meta.label = "seurat_clusters",
  meta.batch = "orig.ident",
  method.reduction = "umap",
  ndim = 30 
)

# Imprime los resultados
print(integration_metricsB)

# Convertir la lista de métricas a un data frame
metrics_dfB <- as.data.frame(integration_metricsB)

# Guardar los 2 data frames en un archivo XLSX

write.xlsx(rbind(metrics_dfA,metrics_dfB), file = "pt17-5.xlsx")
