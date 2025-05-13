#!/bin/bash
#SBATCH --job-name=ad_MACS
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ad45368@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=1:00:00
#SBATCH --output=../OutputSE/logs/MACS3.%j.out
#SBATCH --error=../OutputSE/logs/MACS3.%j.err

#### ChIP-seq output from https://github.com/UGALewisLab/FungiCutAndRun should be organized in a folder ./Output/SortedBamFiles

OUTDIR=/PATH/TO/OUTPUT/DIRECTORY
bam="$OUTDIR/SortedBamFiles"
PLOTDIR="$OUTDIR/Metaplots"
BWDIR="$OUTDIR/BigWigs"

mkdir $PLOTDIR

###This is a non-IP ChIP sample (input)
PeakInput="${bam}/SRR7690277.bam"

module load MACS3/3.0.0b1-foss-2022a-Python-3.10.4
#
macs3 callpeak -t "${bam}/SRR7690281.bam" -c ${PeakInput} -f BAM -n "SRR7690281"  --broad -g 41037538 --broad-cutoff 0.1 --outdir "${OUTDIR}/Peaks" --min-length 800 --max-gap 500 --nomodel
macs3 callpeak -t "${bam}/SRR7690285.bam" -c ${PeakInput} -f BAM -n "SRR7690285"  --broad -g 41037538 --broad-cutoff 0.1 --outdir "${OUTDIR}/Peaks" --min-length 800 --max-gap 500 --nomodel
macs3 callpeak -t "${bam}/SRR12229309.bam" -c ${PeakInput} -f BAMPE -n "SRR12229309"  --broad -g 41037538 --broad-cutoff 0.1 --outdir "${OUTDIR}/Peaks" --min-length 800 --max-gap 500 --nomodel
macs3 callpeak -t "${bam}/SRR12229314.bam" -c ${PeakInput} -f BAM -n "SRR12229314"  --broad -g 41037538 --broad-cutoff 0.1 --outdir "${OUTDIR}/Peaks" --min-length 800 --max-gap 500 --nomodel
#

##### Intersect peaksets with genes ####
module load BEDTools/2.30.0-GCC-12.2.0
GENES=/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed

bedtools intersect  -wa -f 0.5 -e -header -a ${GENES} -b ${OUTDIR}/Peaks/SRR7690281_peaks.broadPeak > ${OUTDIR}/Peaks/SRR7690281_peaks_annotated_50.bed
bedtools intersect  -wa -f 0.5 -e -header -a ${GENES} -b ${OUTDIR}/Peaks/SRR7690285_peaks.broadPeak > ${OUTDIR}/Peaks/SRR7690285_peaks_annotated_50.bed
bedtools intersect  -wa -f 0.5 -e -header -a ${GENES} -b ${OUTDIR}/Peaks/SRR12229309_peaks.broadPeak > ${OUTDIR}/Peaks/SRR12229309_peaks_annotated_50.bed
bedtools intersect  -wa -f 0.5 -e -header -a ${GENES} -b ${OUTDIR}/Peaks/SRR12229314_peaks.broadPeak > ${OUTDIR}/Peaks/SRR12229314_peaks_annotated_50.bed


### next step is to make merged peaksets for each modification. h3k36me2 and h3k36me3 will be combined to a 2/3 gene set
# SampleID	Tissue	Factor	Strain	Replicate
# SRR7690285	Conidia	H3K36me3 	WT 	1
# SRR12229309	Mycelia	H3K36me3 	WT	2

# SRR7690281	Conidia	H3K27me2/3 	WT 	1
# SRR12229314	Mycelia	H3K27me2/3 	WT	2

### Find consensus peaksets for H3K36me3 and H3K27me3
bedtools intersect -wa -f 0.5 -e -header -a ${GENES} -b ${OUTDIR}/Peaks/SRR12229309_peaks.broadPeak ${OUTDIR}/Peaks/SRR7690285_peaks.broadPeak >  ${OUTDIR}/Peaks/H3K36me3_peakset.bed
bedtools intersect -wa -f 0.5 -e -header -a ${GENES} -b  ${OUTDIR}/Peaks/SRR7690281_peaks.broadPeak ${OUTDIR}/Peaks/SRR12229314_peaks.broadPeak >  ${OUTDIR}/Peaks/H3K27me2_me3_peakset.bed

### HEATMAP PLOTTING - checking expression and enrichment at GENES and comparing with RNA-seq (SRR5320488)
# module load deepTools/3.5.2-foss-2022a
#
computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R ${OUTDIR}/Peaks/H3K36me3_peakset.bed \
    -S ${BWDIR}/SRR7690285.bin_25.smooth_75Bulk.bw  ${BWDIR}/SRR12229309.bin_25.smooth_75Bulk.bw ${BWDIR}/SRR5320488.bw \
    -o H3K36peakset.gz --samplesLabel "WT_H3K36me3_1"  "WT_H3K36me3_2" "RNA-seq" \
    --outFileSortedRegions H3K36regions.bed --skipZeros
plotHeatmap --matrixFile H3K36peakset.gz -o H3K36peakset.png --outFileSortedRegions H3K36me3_regions.bed --outFileNameMatrix H3K36me3_regions.gz --linesAtTickMarks --colorMap PuOr_r --sortUsingSamples 5  --kmeans 4 -
#
### Plotting for H3K27me2_3
#
computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R ${OUTDIR}/Peaks/H3K27me2_me3_peakset.bed \
    -S ${BWDIR}/SRR7690281.bin_25.smooth_75Bulk.bw ${BWDIR}/SRR12229314.bin_25.smooth_75Bulk.bw  ${BWDIR}/SRR5320488.bw \
    -o ${PLOTDIR}/H3K27peakset.gz --samplesLabel "WT_H3K27me2_3_1"  "WT_H3K27me2_3_2"  "RNA-seq" \
    --outFileSortedRegions ${PLOTDIR}/H3K27regions.bed --skipZeros

    plotHeatmap --matrixFile  ${PLOTDIR}/H3K27peakset.gz -o  ${PLOTDIR}/H3K27peakset.png --outFileSortedRegions  ${PLOTDIR}/H3K27peakset_REGIONS.bed --outFileNameMatrix  ${PLOTDIR}/H3K27peakset_MATRIX.gz --linesAtTickMarks --colorMap PuOr_r --clusterUsingSamples 1 2 --sortUsingSamples 4  --kmeans 4
#

#ran back through BEDTools to find consensus peakset
bedtools intersect -wa -u -header -a ${OUTDIR}/Peaks/H3K36me3_peakset.bed -b ${OUTDIR}/Peaks/H3K27me2_me3_peakset.bed  >  ${OUTDIR}/Peaks/H3K27_H3K36_shared_regions.bed

#### identifying H3K36 vs H3K27me3 overlaps --- doing this as reciprocals ##
bedtools intersect -wa -header -a ${GENES} -b ${OUTDIR}/Peaks/H3K27_H3K36_shared_regions.bed  >  ${OUTDIR}/Peaks/H3K27me3_H3K36me3_comarked_genes.bed


###Performed final cleanup in Excel - extracted NCU numbers from each gene, removed duplicates, and removaled all non-protein coding genes.also removed mitochondrial genes.
#Final numbers:
# H3K27me3/H3K36 shared = 659 genes
# H3K27me3 alone = 25 GENES
# ash1 dep H3K36me3 only= 2164
