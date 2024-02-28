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
indir <- "/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/STAR_output"
outdir <- "/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/results/"
setwd(indir)
files <- dir(pattern = "ReadsPerGene.out.tab")
counts <- c() # esta sera la matriz
for(i in seq_along(files)){
  x <- read.table(file = files[i], sep = "\t", header = F, as.is = T)
  # as.is para no convertir tipo de datos
  counts <- cbind(counts, x[,2])
}

# Cargar Metadatos
metadata <- read.csv("/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/metadata.csv", header = F)
# Renombrar columnas con el ID de los transcriptomas
colnames(metadata) <- c("sample_id", "type")
# Convertir a formato dataframe
counts <- as.data.frame(counts)
rownames(counts) <- x[,1] # Renombrar las filas con el nombre de los genes
colnames(counts) <- sub("_ReadsPerGene.out.tab", "", files)

# almacenar metadata y matriz de cuentas
save(metadata, counts, file = paste0(outdir, 'counts/STAR_counts.RData'))

# Guardar informacion de ejecucion
sessionInfo()