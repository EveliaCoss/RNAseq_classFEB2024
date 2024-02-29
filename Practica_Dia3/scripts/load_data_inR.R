######
# Script : Importar datos de cuentas en R
# Author: Sofia Salazar, Diego Ramirez y Evelia Coss
# Date: 27/02/2024
# Description: El siguiente script nos permite importar los datos provenientes del alineamiento de STAR a R,
# para el posterior analisis de Expresion diferencial con DESEq2.
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: metadata.csv, cuentas de STAR (Terminacion ReadsPerGene.out.tab)
#   - Output: Matriz de cuentas (CSV y RData)
#######

# qlogin 
# module load r/4.0.2
# R

# --- Load data -----
# Cargar archivos
#indir <- "/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/STAR_output"
indir <- "/mnt/Guanina/bioinfo24/data/STAR_output/"
outdir <- "/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/results/"

# Opcion A - moverme a la carpeta y buscar
setwd(indir)
files <- dir(pattern = "ReadsPerGene.out.tab")

# Opcion B -  sin movernos de carpeta
files <- dir(indir, pattern = "ReadsPerGene.out.tab")

# crear matriz de cuentas
counts <- c() # esta sera la matriz
for(i in seq_along(files)){
  x <- read.table(file = files[i], sep = "\t", header = F, as.is = T)
  # as.is para no convertir tipo de datos
  counts <- cbind(counts, x[,2])
}

# Cargar Metadatos
metadata <- read.csv("/mnt/Guanina/bioinfo24/data/metadata.csv", header = F)
# Renombrar columnas en la metadata
colnames(metadata) <- c("sample_id", "type")
# Convertir a formato dataframe
counts <- as.data.frame(counts)
rownames(counts) <- x[,1] # Renombrar las filas con el nombre de los genes
colnames(counts) <- sub("_ReadsPerGene.out.tab", "", files)

# Eliminar las 4 primeras filas
# counts <- counts[5:129239, ] # Filtramos los rows con informacion general sobre el mapeo
counts <- counts[-c(1:4),]

# Almacenar metadata y matriz de cuentas
save(metadata, counts, file = paste0(outdir, "counts/raw_counts.RData"))
write.csv(counts, file = paste0(outdir,"counts/raw_counts.csv"))

# Guardar informacion de ejecucion
sessionInfo()

# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRblas.so
# LAPACK: /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
# 
# loaded via a namespace (and not attached):
#   [1] compiler_4.0.2 tools_4.0.2
