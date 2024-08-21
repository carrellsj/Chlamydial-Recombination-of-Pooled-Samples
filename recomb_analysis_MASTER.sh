#!/bin/bash
# version 2017.10.12
# version 2017.10.24 - added code to compute depth of coverage of recombination reads
# version 2017.12.13 - added code to genotype indel variants; some general code cleanup
# version 2018.01.13 - added mapping quality filter
# version 2019.11.22 - added Host read filter; added parental plasmid filter
# version 2019.11.23 - MOSAIK for hybrid alignment
# version 2019.11.26 - remove host .sam file 

FASTQ_BASE_FILENAME=DxL2-1_paired_trimmed
PARENTAL_STRAIN_1=Drif
PARENTAL_STRAIN_2=L2ofl
PARENTAL_PLASMID=Parent_Plasmid
HYBRID_GENOME_NAME=DxL2_Hybrid
CROSS_STRAIN_SNP_GENOTYPE_FILE=DxL2_SNPs.txt
CROSS_STRAIN_INDEL_GENOTYPE_FILE=DxL2_indels.txt
CROSS_STRAIN_VARIANT_FEATURE_FILE=DxL2_variants.bed
HOST_GENOME_FILE=C57BL_mus_musculus
MYCOPLASMA_GENOME_FILE=Mycoplasma

MIN_MAPPING_QUALITY_SCORE=9
MIN_NUM_ALLELES_PER_PARENT_STRAIN=2

NUM_CORES=32                      ## for Steve's mac, change to:  8          for Biomed, change to: 32
UNIX_SORT_BIN=sort                ## for Steve's mac, change to:  gsort      for Biomed, change to: sort
LEGACY_SAMTOOLS_SORT=false        ## for Steve's mac, change to:  true       for Biomed, change to: false

## exit immediately if any command produces an error condition
set -e

## get a list of all read IDs for which either read in the pair contains a 3' homopolymer repeat (AAAAA or TTTTT)
egrep -B1 '^[^\@][ACGT]*(TTTTT|AAAAA)$' ${FASTQ_BASE_FILENAME}_R1.fq | egrep '^\@' | cut -f1 -d\ > ${FASTQ_BASE_FILENAME}_read_ids_homopolymers.txt
egrep -B1 '^[^\@][ACGT]*(TTTTT|AAAAA)$' ${FASTQ_BASE_FILENAME}_R2.fq | egrep '^\@' | cut -f1 -d\ >> ${FASTQ_BASE_FILENAME}_read_ids_homopolymers.txt
${UNIX_SORT_BIN} ${FASTQ_BASE_FILENAME}_read_ids_homopolymers.txt | uniq > ${FASTQ_BASE_FILENAME}_read_ids_homopolymers_unique.txt
grep  '^\@' ${FASTQ_BASE_FILENAME}_R1.fq | cut -f1 -d\  > ${FASTQ_BASE_FILENAME}_all_read_ids.txt
${UNIX_SORT_BIN} --parallel=${NUM_CORES} -S 50% ${FASTQ_BASE_FILENAME}_all_read_ids.txt > ${FASTQ_BASE_FILENAME}_all_read_ids_sorted.txt
comm -2 -3 ${FASTQ_BASE_FILENAME}_all_read_ids_sorted.txt ${FASTQ_BASE_FILENAME}_read_ids_homopolymers_unique.txt | cut -f2 -d\@ > ${FASTQ_BASE_FILENAME}_read_ids_no_homopolymers.txt

## make fastq files of only the read pairs without 3' homopolymers
seqtk subseq ${FASTQ_BASE_FILENAME}_R1.fq ${FASTQ_BASE_FILENAME}_read_ids_no_homopolymers.txt > \
                           ${FASTQ_BASE_FILENAME}_R1_no_homopolymers.fq
seqtk subseq ${FASTQ_BASE_FILENAME}_R2.fq ${FASTQ_BASE_FILENAME}_read_ids_no_homopolymers.txt > \
                           ${FASTQ_BASE_FILENAME}_R2_no_homopolymers.fq
echo 'Plasmid'

## align the reads to Parent_Plasmid
bowtie2 -p ${NUM_CORES} -x ${PARENTAL_PLASMID} \
                           -1 ${FASTQ_BASE_FILENAME}_R1_no_homopolymers.fq -2 ${FASTQ_BASE_FILENAME}_R2_no_homopolymers.fq > \
                           ${PARENTAL_PLASMID}.sam

## get the read IDs of the read pairs for which there was not a match to Parent_Plasmid
## Note:  the SAM flag field value 0x0008 means "mate is unmapped"
samtools view -f 0x0008 ${PARENTAL_PLASMID}.sam | cut -f1 | uniq > ${FASTQ_BASE_FILENAME}_read_ids_PLASMID_mismatches.txt

## make fastq files of read pairs not matched to Host strain
seqtk subseq ${FASTQ_BASE_FILENAME}_R1_no_homopolymers.fq ${FASTQ_BASE_FILENAME}_read_ids_PLASMID_mismatches.txt > \
                           ${FASTQ_BASE_FILENAME}_R1_PLASMID_mismatches.fq
seqtk subseq ${FASTQ_BASE_FILENAME}_R2_no_homopolymers.fq ${FASTQ_BASE_FILENAME}_read_ids_PLASMID_mismatches.txt > \
                           ${FASTQ_BASE_FILENAME}_R2_PLASMID_mismatches.fq

echo 'Mycoplasma'

## align the reads to Mycoplasma contamination
bowtie2 -p ${NUM_CORES} -x ${MYCOPLASMA_GENOME_FILE} \
                           -1 ${FASTQ_BASE_FILENAME}_R1_PLASMID_mismatches.fq -2 ${FASTQ_BASE_FILENAME}_R2_PLASMID_mismatches.fq > \
                           ${MYCOPLASMA_GENOME_FILE}.sam

## get the read IDs of the read pairs for which there was not a match to Parent_Plasmid
## Note:  the SAM flag field value 0x0008 means "mate is unmapped"
samtools view -f 0x0008 ${MYCOPLASMA_GENOME_FILE}.sam | cut -f1 | uniq > ${FASTQ_BASE_FILENAME}_read_ids_MYCOPLASMA_mismatches.txt

## make fastq files of read pairs not matched to Host strain
seqtk subseq ${FASTQ_BASE_FILENAME}_R1_PLASMID_mismatches.fq ${FASTQ_BASE_FILENAME}_read_ids_MYCOPLASMA_mismatches.txt > \
                           ${FASTQ_BASE_FILENAME}_R1_MYCOPLASMA_mismatches.fq
seqtk subseq ${FASTQ_BASE_FILENAME}_R2_PLASMID_mismatches.fq ${FASTQ_BASE_FILENAME}_read_ids_MYCOPLASMA_mismatches.txt > \
                           ${FASTQ_BASE_FILENAME}_R2_MYCOPLASMA_mismatches.fq

#remove previous read file 
rm -f ${FASTQ_BASE_FILENAME}_R1_PLASMID_mismatches.fq
rm -f ${FASTQ_BASE_FILENAME}_R2_PLASMID_mismatches.fq
rm -f ${MYCOPLASMA_GENOME_FILE}.sam

echo ‘Host’

## align the reads to HOST_GENOME_FILE 
bowtie2 -p ${NUM_CORES} -x ${HOST_GENOME_FILE} \
                           -1 ${FASTQ_BASE_FILENAME}_R1_MYCOPLASMA_mismatches.fq -2 ${FASTQ_BASE_FILENAME}_R2_MYCOPLASMA_mismatches.fq > \
                           ${HOST_GENOME_FILE}.sam

## get the read IDs of the read pairs for which there was not a match to Host strain
## Note:  the SAM flag field value 0x0008 means "mate is unmapped"
samtools view -f 0x0008 ${HOST_GENOME_FILE}.sam | cut -f1 | uniq > ${FASTQ_BASE_FILENAME}_read_ids_HOST_mismatches.txt

## make fastq files of read pairs not matched to Host strain
seqtk subseq ${FASTQ_BASE_FILENAME}_R1_MYCOPLASMA_mismatches.fq ${FASTQ_BASE_FILENAME}_read_ids_HOST_mismatches.txt > \
                           ${FASTQ_BASE_FILENAME}_R1_HOST_mismatches.fq
seqtk subseq ${FASTQ_BASE_FILENAME}_R2_MYCOPLASMA_mismatches.fq ${FASTQ_BASE_FILENAME}_read_ids_HOST_mismatches.txt > \
                           ${FASTQ_BASE_FILENAME}_R2_HOST_mismatches.fq

#remove previous read files and host sam file
rm -f ${FASTQ_BASE_FILENAME}_R1_MYCOPLASMA_mismatches.fq
rm -f ${FASTQ_BASE_FILENAME}_R2_MYCOPLASMA_mismatches.fq
rm -f ${HOST_GENOME_FILE}.sam

echo ‘Mapping reads to parent strain 1 - serovars D’

## align the reads to PARENTAL_STRAIN_1 requiring perfect match
bowtie2 -p ${NUM_CORES} --score-min 'C,0,-1' -x ${PARENTAL_STRAIN_1} \
                           -1 ${FASTQ_BASE_FILENAME}_R1_HOST_mismatches.fq -2 ${FASTQ_BASE_FILENAME}_R2_HOST_mismatches.fq > \
                           ${PARENTAL_STRAIN_1}_perfect.sam
##Sam to Bam D
#samtools view -@ ${NUM_CORES} ${PARENTAL_STRAIN_1}_perfect.sam > ${PARENTAL_STRAIN_1}_perfect.bam

echo ‘Mapping reads to parent strain 2 - serovars L’

## align the reads to PARENTAL_STRAIN_2 requiring perfect match
bowtie2 -p ${NUM_CORES} --score-min 'C,0,-1' -x ${PARENTAL_STRAIN_2} \
                           -1 ${FASTQ_BASE_FILENAME}_R1_HOST_mismatches.fq -2 ${FASTQ_BASE_FILENAME}_R2_HOST_mismatches.fq > \
                           ${PARENTAL_STRAIN_2}_perfect.sam

##Sam to Bam L2
#samtools view -@ ${NUM_CORES} ${PARENTAL_STRAIN_2}_perfect.sam > ${PARENTAL_STRAIN_2}_perfect.bam

## get the read IDs of the read pairs for fasttrack
## Note:  the SAM flag field value 0x0008 means "mate is unmapped"
samtools view -f 0x0049 ${PARENTAL_STRAIN_1}_perfect.sam | cut -f1 | uniq > ${FASTQ_BASE_FILENAME}_read_ids_onesidemismatches_option1.txt
samtools view -f 0x0085 ${PARENTAL_STRAIN_2}_perfect.sam | cut -f1 | uniq >> ${FASTQ_BASE_FILENAME}_read_ids_onesidemismatches_option1.txt
${UNIX_SORT_BIN} --parallel=${NUM_CORES} -S 50% ${FASTQ_BASE_FILENAME}_read_ids_onesidemismatches_option1.txt | \
                           uniq -c | awk '{ if ($1 > 1) {print $2;} }' > ${FASTQ_BASE_FILENAME}_read_ids_mismatches_fasttrack_option1.txt

## get the read IDs of the read pairs for fasttrack
## Note:  the SAM flag field value 0x0008 means "mate is unmapped"
samtools view -f 0x0085 ${PARENTAL_STRAIN_1}_perfect.sam | cut -f1 | uniq > ${FASTQ_BASE_FILENAME}_read_ids_onesidemismatches_option2.txt
samtools view -f 0x0049 ${PARENTAL_STRAIN_2}_perfect.sam | cut -f1 | uniq >> ${FASTQ_BASE_FILENAME}_read_ids_onesidemismatches_option2.txt
${UNIX_SORT_BIN} --parallel=${NUM_CORES} -S 50% ${FASTQ_BASE_FILENAME}_read_ids_onesidemismatches_option2.txt | \
                           uniq -c | awk '{ if ($1 > 1) {print $2;} }' > ${FASTQ_BASE_FILENAME}_read_ids_mismatches_fasttrack_option2.txt
#echo 'Check Fasttrack'

## get the read IDs of the read pairs for which there was not a perfect match to one or both parental strains
## Note:  the SAM flag field value 0x0008 means "mate is unmapped"
samtools view -f 0x000c ${PARENTAL_STRAIN_1}_perfect.sam | cut -f1 | uniq > ${FASTQ_BASE_FILENAME}_read_ids_mismatches.txt
samtools view -f 0x000c ${PARENTAL_STRAIN_2}_perfect.sam | cut -f1 | uniq >> ${FASTQ_BASE_FILENAME}_read_ids_mismatches.txt
${UNIX_SORT_BIN} --parallel=${NUM_CORES} -S 50% ${FASTQ_BASE_FILENAME}_read_ids_mismatches.txt | \
                           uniq -c | awk '{ if ($1 > 1) {print $2;} }' > ${FASTQ_BASE_FILENAME}_read_ids_mismatches_unique.txt

cat ${FASTQ_BASE_FILENAME}_read_ids_mismatches_fasttrack_option1.txt > ${FASTQ_BASE_FILENAME}_read_ids_mismatches_fasttrack_withdups.txt
cat ${FASTQ_BASE_FILENAME}_read_ids_mismatches_fasttrack_option2.txt >> ${FASTQ_BASE_FILENAME}_read_ids_mismatches_fasttrack_withdups.txt
cat ${FASTQ_BASE_FILENAME}_read_ids_mismatches_unique.txt >> ${FASTQ_BASE_FILENAME}_read_ids_mismatches_fasttrack_withdups.txt

${UNIX_SORT_BIN} --parallel=${NUM_CORES} -S 50% ${FASTQ_BASE_FILENAME}_read_ids_mismatches_fasttrack_withdups.txt | uniq > \
                                                ${FASTQ_BASE_FILENAME}_read_ids_mismatches_fasttrack.txt

##------------- Make FastQ file of reads
## make fastq files of only the read pairs with imperfect matches
seqtk subseq ${FASTQ_BASE_FILENAME}_R1_HOST_mismatches.fq ${FASTQ_BASE_FILENAME}_read_ids_mismatches_fasttrack.txt > \
                           ${FASTQ_BASE_FILENAME}_R1_mismatches.fq
seqtk subseq ${FASTQ_BASE_FILENAME}_R2_HOST_mismatches.fq ${FASTQ_BASE_FILENAME}_read_ids_mismatches_fasttrack.txt > \
                           ${FASTQ_BASE_FILENAME}_R2_mismatches.fq

##Bowtie2

## run bowtie2 to align the imperfect-matching reads against the hybrid genome assembly
bowtie2 -p ${NUM_CORES} -L 15 --score-min L,-1,-1 -N 1 --n-ceil L,0,1 -k 1 --rfg 50,50 -x ${HYBRID_GENOME_NAME} \
                           -1 ${FASTQ_BASE_FILENAME}_R1_mismatches.fq -2 ${FASTQ_BASE_FILENAME}_R2_mismatches.fq > \
                           ${FASTQ_BASE_FILENAME}.sam

## make a BAM file from the SAM file of aligned reads
samtools view -bS ${FASTQ_BASE_FILENAME}.sam > ${FASTQ_BASE_FILENAME}.bam

# ## prep the read data into the format required by MOSAIK
# ${MOSAIK_BIN}/MosaikBuild -q ${FASTQ_BASE_FILENAME}_R1_mismatches.fq \
#                           -q2 ${FASTQ_BASE_FILENAME}_R2_mismatches.fq \
#                           -out ${FASTQ_BASE_FILENAME}.mkb -st illumina

# ## this is the alignment step
# ${MOSAIK_BIN}/MosaikAligner -annpe ${MOSAIK_ANN}/2.1.26.pe.100.0065.ann \
#                             -annse ${MOSAIK_ANN}/2.1.26.se.100.005.ann \
#                             -j ${HYBRID_GENOME_NAME}.15 \
#                             -gop ${MOSAIK_GAP_OPEN_PENALTY} \
#                             -in ${FASTQ_BASE_FILENAME}.mkb \
#                             -ia ${HYBRID_GENOME_NAME}.dat \
#                             -out ${FASTQ_BASE_FILENAME}.mka

#samtools sort -n s001_paired_trimmed.mka.bam > s001_paired_trimmed.mka.sorted.bam

#pairToBed -type either -abam -bedpe \
#          -a ${FASTQ_BASE_FILENAME}.mka.sorted.bam \
#          -b DxL2_SNPs.bed > s001_paired_trimmed.mka.sorted_hitsnps.bed


## sort reads by read ID, in preparation for running pairToBed
if [ "$LEGACY_SAMTOOLS_SORT" = false ] ; then
    samtools sort -n ${FASTQ_BASE_FILENAME}.bam > ${FASTQ_BASE_FILENAME}.sorted.bam  ## this is the new syntax
else
    samtools sort -n ${FASTQ_BASE_FILENAME}.bam ${FASTQ_BASE_FILENAME}.sorted     ## this is the OLD syntax :LEGACY:
fi

## find read pairs for which both reads overlap at least one SNP
pairToBed -type either \
          -abam ${FASTQ_BASE_FILENAME}.sorted.bam \
          -b ${CROSS_STRAIN_VARIANT_FEATURE_FILE} > ${FASTQ_BASE_FILENAME}.sorted_hit_variants.bam

# find candidate recombination events and save to a tab-delimited text file
/local/cluster/R-4.0.3/bin/Rscript genotype_reads_with_no_indels_coords_2Snps_2021.R ${HYBRID_GENOME_NAME}.fasta \
                                               ${CROSS_STRAIN_SNP_GENOTYPE_FILE} \
                                               ${FASTQ_BASE_FILENAME}.sorted_hit_variants.bam \
                                               ${MIN_MAPPING_QUALITY_SCORE} \
                                               ${MIN_NUM_ALLELES_PER_PARENT_STRAIN}

# convert that tab-delimited text file to a SAM file
samtools view -H ${FASTQ_BASE_FILENAME}.sorted_hit_variants.bam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sam

# make a SAM file from the "_candidate_recombinations.txt" file produced by the R genotyping script:
cat ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.txt | \
            awk 'NR >= 2 { print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $7 "\t" $8 "\t=\t" $10 "\t" $11 "\t" $12 "\t" $13 }' >> \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sam

# convert that SAM file to a BAM file
samtools view -bS ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.bam

# convert that BAM file to a BED file
bamToBed -i ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.bam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.bed

# make a coverage file
if [ "$LEGACY_SAMTOOLS_SORT" = false ] ; then
    samtools sort ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.bam > \
                  ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sorted.bam  ## this is the new syntax
else
    samtools sort ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.bam \
                  ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sorted ## this is the OLD syntax :LEGACY:
fi

samtools depth ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sorted.bam > \
               ${FASTQ_BASE_FILENAME}_recombinations_depth_coverage.txt

# convert that tab-delimited text file to a SAM file
samtools view -H ${FASTQ_BASE_FILENAME}.sorted_hit_variants.bam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sam

# make a SAM file from the "_candidate_recombinations.txt" file produced by the R genotyping script:
cat ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.txt | \
            awk 'NR >= 2 { print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $7 "\t" $8 "\t=\t" $10 "\t" $11 "\t" $12 "\t" $13 }' >> \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sam

# convert that SAM file to a BAM file
samtools view -bS ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.bam

# convert that BAM file to a BED file
bamToBed -i ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.bam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.bed

# make a coverage file
if [ "$LEGACY_SAMTOOLS_SORT" = false ] ; then
    samtools sort ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.bam > \
                  ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sorted.bam  ## this is the new syntax
else
    samtools sort ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.bam \
                  ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sorted ## this is the OLD syntax :LEGACY:
fi

samtools depth ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sorted.bam > \
               ${FASTQ_BASE_FILENAME}_NOT_recombinations_depth_coverage.txt


