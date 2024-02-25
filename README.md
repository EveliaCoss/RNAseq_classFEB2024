# Análisis de datos de RNA-Seq

Instructora: Dra. Evelia Coss, Posdoc de la Dra. Alejandra Medina, LIIGH-UNAM

Clases para los alumnos de Ciencias Genomicas de 4to semestre de la ENES, Juriquilla (27 Feb - 1 Marzo 2024). Formando parte de la clase de Bioinformática y Estadística 2. 

## Descripción

El módulo consta de sesiones teóricas y prácticas impartidas de forma presencial, que cubrirán aspectos básicos del tópico como:

- Calidad y limpieza de archivos fastq
- Alineamiento y ensamblaje con el genoma de referencia usando STAR
- Generación del archivo de cuentas crudas
- Importar datos en R
- Normalización y corrección por batch
- Expresión diferencial con DESEq2 y edgeR
- Análisis funcional de los genes detectados
- Visualización grafica de los resultados

Se darán presentaciones detalladas del uso de programas clave, todos de código fuente abierto, usando datos tomados de las bases de datos. También se presentará el uso de algunos scripts de Bash y R muy sencillos, con el objetivo de aprender los aspectos básicos de estos lenguajes para el análisis de datos transcriptómico.

## Contenido 📌

- Dia 1. Aspectos generales de RNA-Seq / Control de calidad de los datos
- Dia 2. Diversos pipeline para Alineamiento, ensamblaje y conteo de reads
- Dia 3. Importar datos en R, Normalización y Corrección por batch / DEG con DESeq2
- Dia 4. GSEA - Análisis funcional

### Dia 1. Aspectos generales de RNA-Seq / Control de calidad de los datos

- Fecha: martes 27 de Febrero 2024
- Presentación:
    - [Aspectos generales de RNA-Seq](https://eveliacoss.github.io/RNAseq_classFEB2024/Presentaciones/Dia1_AspectosGenerales.html#1)
    - [Control de calidad de los datos](https://eveliacoss.github.io/RNAseq_classFEB2024/Presentaciones/Dia1_AspectosGenerales.html#43)
    - [Mis primeros pasos en Bash](https://eveliacoss.github.io/RNAseq_classFEB2024/Presentaciones/Dia1_AspectosGenerales.html#62)
- Datos: `/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/`
- Lecturas y cursos recomendados:

### Dia 2. Diversos pipeline para Alineamiento, ensamblaje y conteo de reads

- Fecha: miércoles 28 de Febrero 2024
- Presentación:
- Lecturas y cursos recomendados:

### Dia 3. Importar datos en R, Normalización y Corrección por batch / DEG con DESeq2

- Fecha: jueves 29 de Febrero 2024
- Presentación:
- Lecturas y cursos recomendados:

### Dia 4. GSEA - Análisis funcional

- Fecha: viernes 1 de marzo 2024
- Presentación:
- Lecturas y cursos recomendados:

## Requisitos

- Contar con una terminal en tu sistema operativo
  - Si cuentas con Windows tener una terminal como [MobaXTerm](https://mobaxterm.mobatek.net) o descargar y acceder a la terminal de [Visual Studio Code](https://code.visualstudio.com/)
  - Si cuentas con una Mac o Linux, ya tienes una terminal incluida.
- Tener conocimientos básicos del uso de R y bash.
- Tener instalado R version 4.3.2 y RStudio
- Paquetes de R con Bioconductor
  - Bioconductor
  - DESEq2
  - tximport
  - topGO
  - biomaRT

```
# Instalar Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")

# Paquetes / librerias
paquetes = c("DESeq2", "tximport", "topGO", "biomaRT")
BiocManager::install(paquetes)
```

- Paquetes de R (CRAN)
  - RColorBrewer (opcional)
  - tidyverse
  - dplyr
  - stringr
  - ggrepel
  - ggplot2

```
install.packages("tidyverse")
#install.packages("RColorBrewer")
install.packages("dplyr")
#install.packages("stringr")
install.packages("ggrepel")
install.packages("ggplot2")
```

## Clases previas

- Introduccion a Rmarkdown
  
En la clase previa les enseñe a crear reportes en Rmarkdown, si necesitan revisar la clase les dejo el link: https://github.com/EveliaCoss/RmarkdownGraphs_notes

- Crear llaves y alias
 
También aprendimos a crear llaves (ssh-keygen) y alias para acceder a los servidores de una manera segura y rápida: https://github.com/EveliaCoss/keygen

## Cursos para practicar 

- [VieRnes de Bioinformática 2023](https://github.com/EveliaCoss/ViernesBioinfo2023)
- [VieRnes de Bioinformática 2024](https://github.com/EveliaCoss/ViernesBioinfo2024)

## Referencias 📚
- [GOterms en S. cereviacae](https://www.yeastgenome.org/goSlimMapper)
- [Go Term finder](https://go.princeton.edu/cgi-bin/GOTermFinder?)
- [REVIGO - pagina principal](http://revigo.irb.hr/FAQ)
- [REVIGO - Reducir terminos GO, ejemplos](https://www.bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html)
- [ggprofiler2](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html)
- [Pathway enrichment analysis and visualization of omics data](https://cytoscape.org/cytoscape-tutorials/protocols/enrichmentmap-pipeline/#/)
- [Biomedical knowledge mining using GOSemSim and Clusterprofiler](https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html)
- [Pathview - Pagina principal](https://pathview.r-forge.r-project.org/)
- [Pathview - Manual](https://pathview.r-forge.r-project.org/pathview.pdf)
- [KEGG - Pathway ID](https://www.genome.jp/kegg/pathway.html)
