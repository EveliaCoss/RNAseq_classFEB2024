#######
# Script : Analisis de terminos GO
# Author: Sofia Salazar, Diego Ramirez y Evelia Coss
# Date: 01/03/2024
# Description: El siguiente script nos permite realiza la Determinacion funcional de los genes diferencialmente expresados
# a partir de los datos provenientes del alineamiento de STAR a R,
# Primero correr el script "load_data_inR.R"
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: metadata.csv, cuentas de STAR (Terminacion ReadsPerGene.out.tab)
#   - Output: Matriz de cuentas (CSV y RData)
#######

# qlogin 
# module load r/4.0.2
# R

# --- Load packages ----------
library(gprofiler2)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(dplyr)

# --- Load data -----
# Cargar archivos
indir <- "/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/results/"
outdir <- "/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/results/"
figdir <- '/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/results/figures/'

# ---- Analisis de terminos Go ----
# Seleccionar bases de datos
sources_db <- c("GO:BP", "KEGG", "REAC", "TF", "MIRNA", "CORUM", "HP", "HPA", "WP")

# Seleccionar solo archivos CSV
files <- dir(indir, pattern = "^DE_(.+)\\.csv$") 

# ---- Ejemplo de UN SOLO ARCHIVO --------
# Extraer el nombre del primer archivo
plot_name <- gsub("^DE_(.+)\\.csv$", "\\1",  files[1]) #name

# Cargar archivo
df <- read.csv(file = paste0(indir, files[1]), row.names = 'X')
head(df)

# baseMean log2FoldChange     lfcSE       stat      pvalue
# Q5GH67                 2.946080     -1.4509823 1.7966778 -0.8075918 0.419325606
# ENSMUST00000192692.1   7.679819     -1.3231310 1.0757085 -1.2300089 0.218693775
# ENSMUST00000193244.1   1.548829      1.9331254 2.5641372  0.7539087 0.450904040
# A0A140LHJ6            29.519295      0.8629139 0.5619564  1.5355531 0.124648017
# P56716                11.956961     -0.8251543 0.7594324 -1.0865408 0.277239787
# Q61473               232.679365     -0.5870415 0.2042779 -2.8737397 0.004056432
# padj
# Q5GH67                       NA
# ENSMUST00000192692.1         NA
# ENSMUST00000193244.1         NA
# A0A140LHJ6           0.39585706
# P56716               0.59994644
# Q61473               0.03436736


# Agregar informacion sobre la expresion
abslogFC <- 2 # Corte de 2 log2FoldChange
df <- df %>% 
  dplyr::mutate(Expression = case_when(log2FoldChange >= abslogFC & padj < 0.05 ~ "Up-regulated",
                                       log2FoldChange <= -(abslogFC) & padj < 0.05 ~ "Down-regulated",
                                       TRUE ~ "Unchanged")) 

# Obtener los nombres de los genes
# > UP 
up_genes <- df %>% filter(Expression == 'Up-regulated') %>% 
  arrange(padj, desc(abs(log2FoldChange)))
# Extraer solo el nombre de los genes
up_genes <- rownames(up_genes) 

head(up_genes)
# [1] "Q60765" "P17879" "Q9Z0F5" "Q3TX21" "P01101" "Q544D6"

# > Down
down_genes <- df %>% filter(Expression == 'Down-regulated') %>% 
  arrange(padj, desc(abs(log2FoldChange)))
# Extraer solo el nombre de los genes
down_genes <- rownames(down_genes) 

head(down_genes)
# [1] "ENSMUST00000127786.3" "ENSMUST00000197196.1" "ENSMUST00000180445.2"
# [4] "ENSMUST00000128339.7" "Q8BLQ0"               "ENSMUST00000200199.1"

# 
multi_gp <- gost(list("Upregulated" = up_genes, 
                      "Downregulated" = down_genes), 
                 correction_method = "fdr", user_threshold = 0.05,
                 multi_query = F, ordered_query = T, 
                 sources = sources_db, 
                 evcodes = TRUE,  # intersection = intersection - a comma separated list of genes 
                 # from the query that are annotated to the corresponding term
                 organism = 'mmusculus') # para humano es hsapiens


## ---- colors ---
# paleta de colores
Category_colors <- data.frame(
  category = c("GO:BP", "GO:CC", "GO:MF", "KEGG",
               'REAC', 'TF', 'MIRNA', 'HPA', 'CORUM', 'HP', 'WP'), 
  label = c('Biological Process', 'Cellular Component', 'Molecular Function',  "KEGG",
            'REAC', 'TF', 'MIRNA', 'HPA', 'CORUM', 'HP', 'WP'),
  colors =  c('#FF9900', '#109618','#DC3912', '#DD4477',
              '#3366CC','#5574A6', '#22AA99', '#6633CC', '#66AA00', '#990099', '#0099C6'))

## ----manhattan plot--------
gostp1 <- gostplot(multi_gp, interactive = FALSE)

# Guardar grafica
ggsave(paste0(figdir, "ManhattanGO_", plot_name, ".png"),
       plot = gostp1, dpi = 300)

## ----Dataframe de todos los datos --------
# Convertir a dataframe
gost_query <- as.data.frame(multi_gp$result)

# Extarer informacion en modo matriz de todos los resultados
bar_data <- data.frame("term" = as.factor(gost_query$term_name), "condition" = gost_query$query, 
                       "count" = gost_query$term_size, "p.adjust" = gost_query$p_value, 
                       'category' = as.factor(gost_query$source), "go_id" = as.factor(gost_query$term_id),
                       'geneNames' = gost_query$intersection
)


## ---- DOWN genes ----
bar_data_down <- subset(bar_data, condition == 'Downregulated')

# Ordenar datos y seleccion por pvalue
bar_data_down <-head(bar_data_down[order(bar_data_down$p.adjust),],40) # order by pvalue
bar_data_down_ordered <- bar_data_down[order(bar_data_down$p.adjust),] # order by pvalue
bar_data_down_ordered<- bar_data_down_ordered[order(bar_data_down_ordered$category),] # order by category
bar_data_down_ordered$p.val <- round(-log10(bar_data_down_ordered$p.adjust), 2)
bar_data_down_ordered$num <- seq(1:nrow(bar_data_down_ordered)) # num category for plot

# Guardar dataset
save(bar_data_down_ordered, file = paste0(outdir, "DOWN_GO_", plot_name, ".RData"))

# agregar colores para la grafica
bar_data_down_ordered_mod <- left_join(bar_data_down_ordered, Category_colors, by= "category")

### ---- DOWN genes (barplot) ----
# Generar la grafica
g.down <- ggplot(bar_data_down_ordered_mod, aes(p.val, reorder(term, -num), fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = p.val),
    color = "black",
    hjust = 0,
    size = 2.2,
    position = position_dodge(0)
  ) +
  labs(x = "-log10(p-value)" , y = NULL) +
  scale_fill_manual(name='Category', 
                    labels = unique(bar_data_down_ordered_mod$label), 
                    values = unique(bar_data_down_ordered_mod$colors)) +
  theme(
    legend.position = "right",
    axis.title.y = element_blank(),
    strip.text.x = element_text(size = 11, face = "bold"),
    strip.background = element_blank()
  )+ theme_classic()


# Guardar la figura
ggsave(paste0(figdir,"barplotDOWN_GO_", plot_name, ".png"),
       plot = g.down + theme_classic(), dpi = 600, width = 10, height = 5)



## ---- UP genes ----
bar_data_up <- subset(bar_data, condition == 'Upregulated')

# Ordenar datos y seleccion por pvalue
bar_data_up <-head(bar_data_up[order(bar_data_up$p.adjust),],40) # order by pvalue
bar_data_up_ordered <- bar_data_up[order(bar_data_up$p.adjust),] # order by pvalue
bar_data_up_ordered<- bar_data_up_ordered[order(bar_data_up_ordered$category),] # order by category
bar_data_up_ordered$p.val <- round(-log10(bar_data_up_ordered$p.adjust), 2)
bar_data_up_ordered$num <- seq(1:nrow(bar_data_up_ordered)) # num category for plot

# Guardar dataset
save(bar_data_up_ordered, file = paste0(outdir, "UP_GO_", plot_name, ".RData"))

# agregar colores para la grafica
bar_data_up_ordered_mod <- left_join(bar_data_up_ordered, Category_colors, by= "category")

### ---- UP genes (barplot) ----
# Generar la grafica
g.up <- ggplot(bar_data_up_ordered_mod, aes(p.val, reorder(term, -num), fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = p.val),
    color = "black",
    hjust = 0,
    size = 2.2,
    position = position_dodge(0)
  ) +
  labs(x = "-log10(p-value)" , y = NULL) +
  scale_fill_manual(name='Category', 
                    labels = unique(bar_data_up_ordered_mod$label), 
                    values = unique(bar_data_up_ordered_mod$colors)) +
  theme(
    legend.position = "right",
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x = element_text(size = 11, face = "bold"),
    strip.background = element_blank() 
  ) + theme_classic()

# Guardar la figura
ggsave(paste0(figdir, "barplotUP_GO_", plot_name, ".png"),
       plot = g.up + theme_classic(), dpi = 600, width = 10, height = 5)



