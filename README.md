# 2025_PRC2_Development -- Scripts used for data analysis in "Polycomb repressive complex 2 regulates sexual development in Neurospora crassa"

1. Scripts used for downloading RNA-seq files from the SRA can be found at the following repository : https://github.com/UGALewisLab/downloadSRA. This will create a directory (./FastqFiles) with raw fastQ files. The accession numbers used for this study can be found in Table S1. 
   
3. Script1-MappingRNAseq: This will iterate through the FastqFiles directory created in the prior step and map RNA-seq files using STAR. This will output sorted .bam files, bigWig files, and Featurecounts files. The featurecounts output files should be moved to a local directory for analysis in R. Script1a is used to launch the mapping script (Script 1b).
   
4. Script2-CompilingRNAseqCounts: This R script should be run from a local machine and will require the following file organization:
      ./MatrixPlotting/CountFiles/BatchX -- directory containing batch files
      ./Datasets/RNAseq_Sample_File.txt -- file containing all metadata for RNAseq experiments
       For this analysis, different samples were organized into batches, so each Batch folder contains a countfiles for different categories. These count files were read in, processed to a single matrix, and then averaged.
  
5. Script3a-IdentifyingH3K27me3_Genes: bash Script to identify peaks and consensus peaks in 4 published H3K27me2_me3 and H3K36me3 datasets. SRA files were downloaded the same way as RNAseq files, and were mapped via this repository: https://github.com/UGALewisLab/FungiCutAndRun.
   
7. Script3b-Finalizing_H3K27me3_Genes: R script for identifying genes that are statistically signficant outliers for gene expression in the H3K27me3-H3K36me3 consensus dataset generated above. This was used to create Table S3 - Silent H3K27me3 Genes

8. Script4-DifferentialExpression: Script for identifying genes that show statistically significant gene expression changes (Table S2)

9. Script5-Figures: Scripts used for all main text figures

10. Script6-Supplement: Scripts used for all supplementary figures.
