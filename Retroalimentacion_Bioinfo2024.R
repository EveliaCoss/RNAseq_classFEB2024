# 
# Author: Michael Love
# Modified: Evelia Coss
# R version 4.4.0

# --- Paquetes ----
# Tuve que reinstalar Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Instalar paquete "pasilla"
BiocManager::install("pasilla")
# https://bioconductor.org/packages/release/data/experiment/html/pasilla.html

# Cargar paquete
library("pasilla")
library("DESeq2")

# the pasilla data constructed from the count matrix method above. 
# This data set is from an experiment on Drosophila melanogaster cell cultures and investigated the effect of RNAi knock-down of 
# the splicing factor pasilla (Brooks et al. 2011). The detailed transcript of the production of the pasilla data is provided 
# in the vignette of the data package pasilla.

# --- Cargar datos ----
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)

# --- Acomodar informacion ----
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# Visualizar primeross datos 
head(cts,2)

# metadata (toda la informacion de las variables)
coldata

#            condition        type
# treated1     treated single-read
# treated2     treated  paired-end
# treated3     treated  paired-end
# untreated1 untreated single-read
# untreated2 untreated single-read
# untreated3 untreated  paired-end
# untreated4 untreated  paired-end

# Ordenar filas y columnas con el mismo orden
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))

# Salvar nueva variable ordenada
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# --- Crear objeto DESEq2 ----
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds

# class: DESeqDataSet 
# dim: 14599 7 
# metadata(1): version
# assays(1): counts
# rownames(14599): FBgn0000003 FBgn0000008 ... FBgn0261574 FBgn0261575
# rowData names(0):
#   colnames(7): treated1 treated2 ... untreated3 untreated4
# colData names(2): condition type

# --- Filtro (opcional) ----
# Podemos hacer un filtro de limpieza, 
# Eliminar genes que no tienen por lo menos 10 reads en 3 muestras
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds

# class: DESeqDataSet 
# dim: 8148 7 
# metadata(1): version
# assays(1): counts
# rownames(8148): FBgn0000008 FBgn0000017 ... FBgn0261573 FBgn0261574
# rowData names(0):
#   colnames(7): treated1 treated2 ... untreated3 untreated4
# colData names(2): condition type

# --- Filtro (opcional) ----

# cuales son los grupos?
levels(dds$condition)

# condition treated vs untreated
# Para asignar que condiciones es la referencia se hace esto:
dds$condition <- relevel(dds$condition, ref = "untreated")

# NOTA: si tienen replicas tecnicas, existe la funcion "collapseReplicates"

# --- Expresion diferencial ------

dds <- DESeq(dds)

# Obtener los resultados de una comparacion
# Opcion A (aqui solo estamos comparando una variable)
res <- results(dds)
summary(res)

# out of 8148 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 533, 6.5%
# LFC < 0 (down)     : 536, 6.6%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 5)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

res

# log2 fold change (MLE): condition treated vs untreated 
# Wald test p-value: condition treated vs untreated 
# DataFrame with 8148 rows and 6 columns
# baseMean log2FoldChange     lfcSE       stat    pvalue      padj
# <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
#   FBgn0000008   95.28865     0.00399148  0.225010  0.0177391 0.9858470  0.996699
# FBgn0000017 4359.09632    -0.23842494  0.127094 -1.8759764 0.0606585  0.289604
# FBgn0000018  419.06811    -0.10185506  0.146568 -0.6949338 0.4870968  0.822681
# FBgn0000024    6.41105     0.21429657  0.691557  0.3098756 0.7566555  0.939146
# FBgn0000032  990.79225    -0.08896298  0.146253 -0.6082822 0.5430003  0.848881
# ...                ...            ...       ...        ...       ...       ...
# FBgn0261564   1160.028     -0.0857255  0.108354 -0.7911643 0.4288481  0.789246
# FBgn0261565    620.388     -0.2943294  0.140496 -2.0949303 0.0361772  0.206423
# FBgn0261570   3212.969      0.2971841  0.126742  2.3447877 0.0190379  0.133380
# FBgn0261573   2243.936      0.0146611  0.111365  0.1316493 0.8952617  0.977565
# FBgn0261574   4863.807      0.0179729  0.194137  0.0925784 0.9262385  0.986726


# --- Obtener los resultados --------

# Ver todas las comparaciones
resultsNames(dds)

# [1] "Intercept"                      "condition_treated_vs_untreated"

# Seleccionar una de las comparaciones
# Opcion B
res <- results(dds, name="condition_treated_vs_untreated")
summary(res)

# out of 8148 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 533, 6.5%
# LFC < 0 (down)     : 536, 6.6%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 5)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# NOTA: Da lo mismo que la Opcion A, porque no hay mas variables

# Opcion C
res <- results(dds, contrast=c("condition","treated","untreated"))
summary(res)

# out of 8148 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 533, 6.5%
# LFC < 0 (down)     : 536, 6.6%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 5)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# NOTA: Da lo mismo que la Opcion A y B, porque no hay mas variables

# Si quieres seleccionar por pvalue
res05 <- results(dds, alpha=0.05)
summary(res05)

# out of 8148 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 416, 5.1%
# LFC < 0 (down)     : 437, 5.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 5)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# ---  Ejemplo con 2 variables  ----
# Tutorial completo: https://github.com/tavareshugo/tutorial_DESeq2_contrasts
# Seccion que vamos a ver: https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md

# simulate data
dds <- makeExampleDESeqDataSet(n = 1000, m = 12, betaSD = 2)
dds$colour <- factor(rep(c("pink", "white"), each = 6))
dds$colour <- relevel(dds$colour, "white")
dds$condition <- factor(rep(c("sun", "shade"), 6))
dds <- dds[, order(dds$colour, dds$condition)]
colnames(dds) <- paste0("sample", 1:ncol(dds))

colData(dds)

## DataFrame with 12 rows and 2 columns
##          condition   colour
##           <factor> <factor>
## sample1      shade    white
## sample2      shade    white
## sample3      shade    white
## sample4      sun      white
## sample5      sun      white
## ...            ...      ...
## sample8      shade     pink
## sample9      shade     pink
## sample10     sun       pink
## sample11     sun       pink
## sample12     sun       pink

# Modelo
design(dds) <- ~ colour + condition + colour:condition

# Reasignar referencias
dds$colour <- relevel(dds$colour, ref = "white")
dds$condition <- relevel(dds$condition, ref = "shade")


dds <- DESeq(dds) # Crear el objeto de DESEQ
resultsNames(dds) # Observar contrastes

# [1] "Intercept"               "colour_pink_vs_white"    "condition_sun_vs_shade"  "colourpink.conditionsun"

# NOTA: colour:condition el efecto del sol o la sombra atraves de las condiciones de colores.

# Extraer la informacion del primer contraste
res1 <- results(dds, name="colour_pink_vs_white")
summary(res1)

# out of 998 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 212, 21%
# LFC < 0 (down)     : 199, 20%
# outliers [1]       : 1, 0.1%
# low counts [2]     : 20, 2%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

res1

# Extraer la informacion del segundo contraste
res2 <- results(dds, name="condition_sun_vs_shade")
summary(res2)

# out of 998 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 1, 0.1%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

res2


# ---  Ajustar el modelo  ----

# Podemos modificar nuestro modelo

# get the model matrix
mod_mat <- model.matrix(design(dds), colData(dds))
mod_mat

# Define coefficient vectors for each condition
pink_shade <- colMeans(mod_mat[dds$colour == "pink" & dds$condition == "shade", ])
pink_sun <- colMeans(mod_mat[dds$colour == "pink" & dds$condition == "sun", ])
white_shade <- colMeans(mod_mat[dds$colour == "white" & dds$condition == "shade", ])
white_sun <- colMeans(mod_mat[dds$colour == "white" & dds$condition == "sun", ])

# We are now ready to define any contrast of interest from these vectors (for completeness we show the 
# equivalent syntax using the coefficient's names from DESeq).

## Pink vs White (in the shade) ----
res1 <- results(dds, contrast = pink_shade - white_shade)

# or equivalently
res1 <- results(dds, contrast = list("colour_pink_vs_white"))

summary(res1)

# out of 998 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 209, 21%
# LFC < 0 (down)     : 201, 20%
# outliers [1]       : 2, 0.2%
# low counts [2]     : 58, 5.8%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

## Pink vs White (in the sun) ----

res2 <- results(dds, contrast = pink_sun - white_sun)

# or equivalently
res2 <- results(dds, contrast = list(c("colour_pink_vs_white",
                                       "colourpink.conditionsun")))
summary(res2)

# out of 998 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 227, 23%
# LFC < 0 (down)     : 200, 20%
# outliers [1]       : 2, 0.2%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

## Sun vs Shade (for whites) ----

res3 <- results(dds, contrast = white_sun - white_shade)

# or equivalently
res3 <- results(dds, contrast = list(c("condition_sun_vs_shade")))

summary(res3)

# out of 998 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2, 0.2%
# LFC < 0 (down)     : 1, 0.1%
# outliers [1]       : 2, 0.2%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

#  Sun vs Shade (for pinks) ----

res4 <- results(dds, contrast = pink_sun - pink_shade)
# or equivalently
res4 <- results(dds, contrast = list(c("condition_sun_vs_shade", 
                                       "colourpink.conditionsun")))

summary(res4)

# out of 998 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 2, 0.2%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

## Interaction between colour and condition (i.e. do pinks and whites respond differently to the sun?): ----

res5 <- results(dds, 
                contrast = (pink_sun - pink_shade) - (white_sun - white_shade))

# or equivalently
res5 <- results(dds, contrast = list("colourpink.conditionsun"))
summary(res5)

# out of 998 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.1%
# outliers [1]       : 2, 0.2%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# In conclusion, although we can define these contrasts using DESeq coefficient names, 
# it is somewhat more explicit (and perhaps intuitive?) what it is we're comparing using matrix-based contrasts. 

