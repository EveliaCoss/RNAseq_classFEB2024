######
# Script : Analisis de expresion diferencial
# Author: Sofia Salazar, Diego Ramirez y Evelia Coss
# Date: 27/02/2024
# Description: El siguiente script nos permite realiza el Analisis de expresion Diferencial
# a partir de los datos provenientes del alineamiento de STAR a R,
# Primero correr el script "load_data_inR.R"
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: Cargar la variable raw_counts.RData que contiene la matriz de cuentas y la metadata
#   - Output: DEG
#######

# qlogin 
# module load r/4.0.2
# R

# --- Load packages ----------
library(DESeq2)

# --- Load data -----
# Cargar archivos
outdir <- "/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/results/"
figdir <- '/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/results/figures/'

#Cargar variable "counts", proveniente del script "load_data_inR.R"
load("/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/results/counts/raw_counts.RData") 
samples <- metadata$sample_id # Extraer los nombres de los Transcriptomas
metadata$type <- as.factor(metadata$type) # convertir a factor

# --- DEG ----
counts <- counts[which(rowSums(counts) > 10),] #Seleccionamos genes con mas de 10 cuentas

# Convertir al formato dds
dds <- DESeqDataSetFromMatrix(countData =  counts, 
            colData = metadata, design = ~type) #Se hace un DESeqDataSet para realizar un analisis

dim(dds) # checar las dimensiones
#[1] 33470     8

##  -- Asignar la referencia y generar contrastes -----
# Las comparaciones se realizan por pares
#Si no se indica de manera explicita que se va a comparara, lo va a tomar de manera alfabetica, 
# en este caso se indica que control es la referencia, 
dds$type <- relevel(dds$type, ref = "CONTROL") 

## --- Obtener archivo dds ----

dds <- DESeq(dds)

# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

# Obtener la lista de coeficientes o contrastes
resultsNames(dds)

# [1] "Intercept"                 "type_PLS_15min_vs_CONTROL"
# [3] "type_PLS_30min_vs_CONTROL" "type_PLS_4h_vs_CONTROL"

# Guardar la salida del diseno
save(metadata, dds, file = paste0(outdir, 'dds_Times_vs_control.RData'))

## --- Normalizacion de los datos ---------
# Opcion 1. log2(n + 1)
ntd <- normTransform(dds)

# Opcion 2. regularized logarithm or rlog
# Normalizacion de las cuentas por logaritmo y podrias hacer el analisis usando este objeto en lugar del dds
ddslog <- rlog(dds, blind = F) 

# Opcion 3. vsd
# Estima la tendencia de dispersion de los datos y calcula la varianza, hace una normalizacion de las 
# cuentas con respecto al tamaño de la libreria
vsdata <- vst(dds, blind = F) 

## --- Deteccion de batch effect ----

# Almacenar la grafica
png(file = paste0(figdir, "PCA_rlog.png"))
plt <- plotPCA(ddslog, intgroup = "type")
print(plt)
dev.off()

# Almacenar la grafica
png(file = paste0(figdir, "PCA_vsd.png"))
plt <- plotPCA(vsdata, intgroup = "type")
print(plt)
dev.off()

# Guardar la salida del diseno (vsdata)
save(metadata, vsdata, file = paste0(outdir, 'vst_Times_vs_control.RData'))

# En la grafica de las primeras dos componentes principales son notorias las diferencias 
# entre tipos de muestras con respecto a las componente principales que capturan su varianza, 
# cada componente principal representa una combinacion lineal de las variables (en este caso genes) 
# que explican la mayor cantidad de varianza en nuestros datos (las cuentas).


## ---- Obtener informacion del contraste 1 ----
# results(dds, contrast=c("condition","treated","untreated"))
res_15t <- results(dds, name = "type_PLS_15min_vs_CONTROL")
res_15t

summary(res_15t)

# out of 33470 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1657, 5%
# LFC < 0 (down)     : 811, 2.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 18169, 54%
# (mean count < 12)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Guardar los resultados
write.csv(res_15t, file=paste0(outdir, 'DE_15min_vs_control.csv'))

## ---- Obtener informacion del contraste 2 ----
res_30t <- results(dds, name = "type_PLS_30min_vs_CONTROL")
res_30t

summary(res_30t)

# out of 33470 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1309, 3.9%
# LFC < 0 (down)     : 1300, 3.9%
# outliers [1]       : 0, 0%
# low counts [2]     : 19467, 58%
# (mean count < 15)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Guardar los resultados
write.csv(res_30t, file=paste0(outdir, 'DE_30min_vs_control.csv'))

## ---- Obtener informacion del contraste 3 ----
res_4t <- results(dds, name = "type_PLS_4h_vs_CONTROL")
res_4t

summary(res_4t)

# out of 33470 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 4289, 13%
# LFC < 0 (down)     : 3409, 10%
# outliers [1]       : 0, 0%
# low counts [2]     : 6489, 19%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Guardar los resultados
write.csv(res_4t, file=paste0(outdir, 'DE_4h_vs_control.csv'))
