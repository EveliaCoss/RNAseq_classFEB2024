# sustituye los paths con los tuyos
for file in /mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/data/*.fastq.gz; do fastqc $file -o /mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/quality1; done
cd /mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/quality1
multiqc .
