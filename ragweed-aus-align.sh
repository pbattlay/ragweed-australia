#####
# RAGWEED 2024
# ALIGNMENT AND SNP CALLING
#####

# Identify adapter sequences
cd ~/ha22_scratch/ragweed-aus/

ls fastq/*.fq.gz | awk -F "/" '{print $2}' | awk -F "_" '{print $1}' | sort -u > samples.txt

# adapters.sh
#!/bin/bash
#SBATCH --job-name=adapters
#SBATCH --account=om62
#SBATCH --time=00:30:00
#SBATCH --output=adapters-%a.out
#SBATCH --error=adapters-%a.err
#SBATCH --array=1-101

N=$SLURM_ARRAY_TASK_ID
SAMP=$(cat ~/ha22_scratch/ragweed-aus/samples.txt | head -n $N | tail -n 1)

module load adapterremoval
AdapterRemoval --file1 ~/ha22_scratch/ragweed-aus/fastq/${SAMP}_L2_1.fq.gz --file2 ~/ha22_scratch/ragweed-aus/fastq/${SAMP}_L2_2.fq.gz --identify-adapters > ~/ha22_scratch/ragweed-aus/fastq/${SAMP}.adapters

###


# Align FASTQ files to reference
# a template makefile
# -*- mode: Yaml; -*-
# Timestamp: 2017-02-20T16:32:12.844855
#
# Default options.
# Can also be specific for a set of samples, libraries, and lanes,
# by including the "Options" hierarchy at the same level as those
# samples, libraries, or lanes below. This does not include
# "Features", which may only be specific globally.
Options:
  # Sequencing platform, see SAM/BAM reference for valid values
  Platform: Illumina
  # Quality offset for Phred scores, either 33 (Sanger/Illumina 1.8+)
  # or 64 (Illumina 1.3+ / 1.5+). For Bowtie2 it is also possible to
  # specify 'Solexa', to handle reads on the Solexa scale. This is
  # used during adapter-trimming and sequence alignment
  QualityOffset: 33
  # Split a lane into multiple entries, one for each (pair of) file(s)
  # found using the search-string specified for a given lane. Each
  # lane is named by adding a number to the end of the given barcode.
  SplitLanesByFilenames: yes
  # Compression format for FASTQ reads; 'gz' for GZip, 'bz2' for BZip2
  CompressionFormat: gz

  # Settings for trimming of reads, see AdapterRemoval man-page
  AdapterRemoval:
     # Adapter sequences, set and uncomment to override defaults
     --adapter1: ADAPTER1
     --adapter2: ADAPTER2
     # Some BAM pipeline defaults differ from AR defaults;
     # To override, change these value(s):
     --mm: 3
     --minlength: 25
     # Extra features enabled by default; change 'yes' to 'no' to disable
     --collapse: yes
     --trimns: yes
     --trimqualities: yes
     --mate-separator: "."

  # Settings for aligners supported by the pipeline
  Aligners:
    # Choice of aligner software to use, either "BWA" or "Bowtie2"
    Program: BWA

    # Settings for mappings performed using BWA
    BWA:
      # One of "backtrack", "bwasw", or "mem"; see the BWA documentation
      # for a description of each algorithm (defaults to 'backtrack')
      Algorithm: mem
      # Filter aligned reads with a mapping quality (Phred) below this value
      MinQuality: 0
      # Filter reads that did not map to the reference sequence
      FilterUnmappedReads: no
      # May be disabled ("no") for aDNA alignments, as post-mortem damage
      # localizes to the seed region, which BWA expects to have few
      # errors (sets "-l"). See http:/pmid.us/22574660
      UseSeed: yes
      # Additional command-line options may be specified for the "aln"
      # call(s), as described below for Bowtie2 below.

    # Settings for mappings performed using Bowtie2
    Bowtie2:
      # Filter aligned reads with a mapping quality (Phred) below this value
      MinQuality: 0
      # Filter reads that did not map to the reference sequence
      FilterUnmappedReads: yes
      # Examples of how to add additional command-line options
#      --trim5: 5
#      --trim3: 5
      # Note that the colon is required, even if no value is specified
      --very-sensitive:
      # Example of how to specify multiple values for an option
#      --rg:
#        - CN:SequencingCenterNameHere
#        - DS:DescriptionOfReadGroup

  # Mark / filter PCR duplicates. If set to 'filter', PCR duplicates are
  # removed from the output files; if set to 'mark', PCR duplicates are
  # flagged with bit 0x400, and not removed from the output files; if set to
  # 'no', the reads are assumed to not have been amplified. Collapsed reads
  # are filtered using the command 'paleomix rmdup_duplicates', while "normal"
  # reads are filtered using Picard MarkDuplicates.
  PCRDuplicates: mark

  # Command-line options for mapDamage; note that the long-form
  # options are expected; --length, not -l, etc. Uncomment the
  # "mapDamage" line adding command-line options below.
  mapDamage:
    # By default, the pipeline will downsample the input to 100k hits
    # when running mapDamage; remove to use all hits
    --downsample: 100000

  # Set to 'yes' exclude a type of trimmed reads from alignment / analysis;
  # possible read-types reflect the output of AdapterRemoval
  ExcludeReads:
    # Exclude single-end reads (yes / no)?
    Single: no
    # Exclude non-collapsed paired-end reads (yes / no)?
    Paired: no
    # Exclude paired-end reads for which the mate was discarded (yes / no)?
    Singleton: no
    # Exclude overlapping paired-ended reads collapsed into a single sequence
    # by AdapterRemoval (yes / no)?
    Collapsed: no
    # Like 'Collapsed', but only for collapsed reads truncated due to the
    # presence of ambiguous or low quality bases at read termini (yes / no).
    CollapsedTruncated: no

  # Optional steps to perform during processing.
  Features:
    # Generate BAM without realignment around indels (yes / no)
    RawBAM: no
    # Generate indel-realigned BAM using the GATK Indel realigner (yes / no)
    RealignedBAM: yes
    # To disable mapDamage, write 'no'; to generate basic mapDamage plots,
    # write 'plot'; to build post-mortem damage models, write 'model',
    # and to produce rescaled BAMs, write 'rescale'. The 'model' option
    # includes the 'plot' output, and the 'rescale' option includes both
    # 'plot' and 'model' results. All analyses are carried out per library.
    mapDamage: plot
    # Generate coverage information for the raw BAM (wo/ indel realignment).
    # If one or more 'RegionsOfInterest' have been specified for a prefix,
    # additional coverage files are generated for each alignment (yes / no)
    Coverage: yes
    # Generate histogram of number of sites with a given read-depth, from 0
    # to 200. If one or more 'RegionsOfInterest' have been specified for a
    # prefix, additional histograms are generated for each alignment (yes / no)
    Depths: no
    # Generate summary table for each target (yes / no)
    Summary: yes
    # Generate histogram of PCR duplicates, for use with PreSeq (yes / no)
    DuplicateHist: no


# Map of prefixes by name, each having a Path key, which specifies the
# location of the BWA/Bowtie2 index, and optional label, and an option
# set of regions for which additional statistics are produced.
Prefixes:
  # Replace 'NAME_OF_PREFIX' with name of the prefix; this name
  # is used in summary statistics and as part of output filenames.
#  Alternaria_alternata:
#    # Replace 'PATH_TO_PREFIX' with the path to .fasta file containing the
#    # references against which reads are to be mapped. Using the same name
#    # as filename is strongly recommended (e.g. /path/to/Human_g1k_v37.fasta
#    # should be named 'Human_g1k_v37').
#    Path: "/home/ahomet/v/vanessab/data/bigdata/Reference_genomes/Alternaria_alternata/LMXP01.1.fsa_nt.fasta"
  ragweed2021:
    # Replace 'PATH_TO_PREFIX' with the path to .fasta file containing the
    # references against which reads are to be mapped. Using the same name
    # as filename is strongly recommended (e.g. /path/to/Human_g1k_v37.fasta
    # should be named 'Human_g1k_v37').
    Path: /home/pbattlay/om62_scratch/ragweed2022/ref/ragweed-dipasm-hap1.fasta
#    RegionsOfInterest: 
#      Path:

        
  ragweed_chloroplast_MG019037:
    Path: /home/pbattlay/om62_scratch/ragweed2022/ref/aarte-chlor.fasta    

    # (Optional) Uncomment and replace 'PATH_TO_BEDFILE' with the path to a
    # .bed file listing extra regions for which coverage / depth statistics
    # should be calculated; if no names are specified for the BED records,
    # results are named after the chromosome / contig. Change 'NAME' to the
    # name to be used in summary statistics and output filenames.
#    RegionsOfInterest:



# Mapping targets are specified using the following structure. Uncomment and
# replace 'NAME_OF_TARGET' with the desired prefix for filenames.
#NAME_OF_TARGET:
   #  Uncomment and replace 'NAME_OF_SAMPLE' with the name of this sample.
#  NAME_OF_SAMPLE:
     #  Uncomment and replace 'NAME_OF_LIBRARY' with the name of this sample.
#    NAME_OF_LIBRARY:
       # Uncomment and replace 'NAME_OF_LANE' with the name of this lane,
       # and replace 'PATH_WITH_WILDCARDS' with the path to the FASTQ files
       # to be trimmed and mapped for this lane (may include wildcards).
#      NAME_OF_LANE: PATH_WITH_WILDCARDS

SAMPLENAME:
  SAMPLENAME:
    SAMPLENAME:
      Lane2_NovaSeq: /home/pbattlay/ha22_scratch/ragweed-aus/fastq/SAMPLENAME_L2_{Pair}.fq.gz

###



# make a makefile for each sample using sample name and adapter sequences
cd ~/ha22_scratch/ragweed-aus/

cat samples.txt | while read SAMP
do
A1=$(cat fastq/$SAMP.adapters | grep 'adapter1:' | sed 's/  --adapter1:  //')
A2=$(cat fastq/$SAMP.adapters | grep 'adapter2:' | sed 's/  --adapter2:  //')
cat makefile-template.yaml | sed "s/ADAPTER1/$A1/" | sed "s/ADAPTER2/$A2/" | sed "s/SAMPLENAME/$SAMP/" > bam/$SAMP.yaml
done



# an array job to run paleomix for each sample
# paleomix-aus-array.sh
#!/bin/bash
#SBATCH --job-name=paleomix-aus
#SBATCH --account=om62
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=16G
#SBATCH --time=12:00:00
#SBATCH --output=paleomix-aus-%a.out
#SBATCH --error=paleomix-aus-%a.err
#SBATCH --array=1-101

module load paleomix/1.2.13.4-python2
module load adapterremoval/2.3.1
module load mapdamage/2.0.9-u1

source /usr/local2/bioinformatics/bioansible_env.sh
module load bwa/v0.7.15
module load samtools/1.8

cd ~/ha22_scratch/ragweed-aus/bam/

N=$SLURM_ARRAY_TASK_ID
SAMP=$(cat ~/ha22_scratch/ragweed-aus/samples.txt | head -n $N | tail -n 1)

paleomix bam_pipeline run /home/pbattlay/ha22_scratch/ragweed-aus/bam/$SAMP.yaml \
	--jar-root=/home/pbattlay/om62_scratch/jar/ \
	--max-threads 12 \
	--bwa-max-threads 12 \
	--adapterremoval-max-threads 12 \
	--jre-option -Xmx16g \
	--temp-root /home/pbattlay/ha22_scratch/ragweed-aus/bam/temp/

###

# extract mean coverage for each alignment
cd ~/ha22_scratch/ragweed-aus/bam/

for i in *.summary
do
cat $i | grep "hits_coverage(ragweed2021)" | awk '{print $1, $5}' | head -n 1
done








####
# VARIANT CALLING
####

# in R
options(scipen=999)

# read in a fasta index and make a list of 1Mb windows across the genome
dat = read.table("~/scratch/ragweed/ref/ragweed-dipasm-hap1.fasta.fai")[, 1:2]

# only keep contigs > 100kb in length
dat = dat[dat$V2 > 100000, ]

# a vector of all contigs
ctgs = dat$V1

for (c in 1:length(ctgs)){
clen = dat[dat$V1 == ctgs[c], 2]
dat.c = data.frame(ctgs[c], seq(1, clen, 1000000), seq(1, clen, 1000000) + 999999)
dat.c[nrow(dat.c), 3] = clen
if (c == 1){ dat.out = dat.c } else { dat.out = rbind(dat.out, dat.c) }
}

write.table(dat.out, file = "~/scratch/ragweed/windows.list", col.names = F, row.names = F, quote = F)
# 1124 windows
###

# original 311 samples (North America; Europe) + 95 Australian samples (exclude AU18-22A; very low coverage)
samples.txt
# 406 samples

# make bam list from sample files
while read SAMP; do ls bam/$SAMP.ragweed2021.realigned.bam; done < samples.txt > bam/bam-UG.list


# UG-array.sh
#!/bin/bash
#SBATCH --job-name=UG-array
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=24:00:00
#SBATCH --output=UG-array-%a.out
#SBATCH --error=UG-array-%a.err
#SBATCH --array=1-1124%200

cd ~/scratch/ragweed/

N=$SLURM_ARRAY_TASK_ID
cat windows.list | head -n $N | tail -n 1 | awk '{print $1 ":" $2 "-" $3}' > window-$N.intervals

module load nixpkgs/16.09
module load gatk/3.8

JAVA_TOOL_OPTIONS="-Xmx8g"

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
	-T UnifiedGenotyper \
	-nt 8 \
	-R ref/ragweed-dipasm-hap1.fasta \
	-I bam/bam-UG.list \
	-L window-$N.intervals \
	-glm BOTH \
	-o vcf/window-$N.vcf

###


#####
# VARIANT FILTERING
#####

# extract statistics for each 1Mb chunk
# vcf-stats.sh
#!/bin/bash
#SBATCH --job-name=vcf-stats
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=00:20:00
#SBATCH --output=vcf-stats-%a.out
#SBATCH --error=vcf-stats-%a.err
#SBATCH --array=1-1124

cd ~/scratch/ragweed/

N=$SLURM_ARRAY_TASK_ID
vcf=$(echo "vcf/window-"$N".vcf")

module load vcftools/0.1.16

# SNPs
vcftools --vcf $vcf --remove-indels --recode --recode-INFO-all --out ${vcf/.vcf/.SNPs}
# QD
cat ${vcf/.vcf/.SNPs}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^QD=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.SNPs}.QD
# DP
cat ${vcf/.vcf/.SNPs}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^DP=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.SNPs}.DP
# MQ
cat ${vcf/.vcf/.SNPs}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^MQ=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.SNPs}.MQ
# FS
cat ${vcf/.vcf/.SNPs}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^FS=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.SNPs}.FS
# SOR
cat ${vcf/.vcf/.SNPs}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^SOR=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.SNPs}.SOR
# MQRankSum
cat ${vcf/.vcf/.SNPs}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^MQRankSum=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.SNPs}.MQRankSum
# ReadPosRankSum
cat ${vcf/.vcf/.SNPs}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^ReadPosRankSum=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.SNPs}.ReadPosRankSum

# indels
vcftools --vcf $vcf --keep-only-indels --recode --recode-INFO-all --out ${vcf/.vcf/.indels}
# QD
cat ${vcf/.vcf/.indels}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^QD=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.indels}.QD
# DP
cat ${vcf/.vcf/.indels}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^DP=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.indels}.DP
# FS
cat ${vcf/.vcf/.indels}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^FS=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.indels}.FS
# SOR
cat ${vcf/.vcf/.indels}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^SOR=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.indels}.SOR
# ReadPosRankSum
cat ${vcf/.vcf/.indels}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^ReadPosRankSum=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.indels}.ReadPosRankSum
# InbreedingCoeff
cat ${vcf/.vcf/.indels}.recode.vcf | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^InbreedingCoeff=" | awk -F "=" '{print $2}' > ${vcf/.vcf/.indels}.InbreedingCoeff

###

# concatenate statistic lists
cat vcf/*.SNPs.QD > vcf/stats/SNPs.QD.txt
cat vcf/*.SNPs.DP > vcf/stats/SNPs.DP.txt
cat vcf/*.SNPs.MQ > vcf/stats/SNPs.MQ.txt
cat vcf/*.SNPs.FS > vcf/stats/SNPs.FS.txt
cat vcf/*.SNPs.SOR > vcf/stats/SNPs.SOR.txt
cat vcf/*.SNPs.MQRankSum > vcf/stats/SNPs.MQRankSum.txt
cat vcf/*.SNPs.ReadPosRankSum > vcf/stats/SNPs.ReadPosRankSum.txt

cat vcf/*.indels.QD > vcf/stats/indels.QD.txt
cat vcf/*.indels.DP > vcf/stats/indels.DP.txt
cat vcf/*.indels.FS > vcf/stats/indels.FS.txt
cat vcf/*.indels.SOR > vcf/stats/indels.SOR.txt
cat vcf/*.indels.ReadPosRankSum > vcf/stats/indels.ReadPosRankSum.txt
cat vcf/*.indels.InbreedingCoeff > vcf/stats/indels.InbreedingCoeff.txt

# clean up
rm vcf/*.log
rm vcf/*.QD
rm vcf/*.DP
rm vcf/*.MQ
rm vcf/*.FS
rm vcf/*.SOR
rm vcf/*.MQRankSum
rm vcf/*.ReadPosRankSum
rm vcf/*.InbreedingCoeff

# statistic distribution plots in R with GATK-recommended thresholds marked
library(ggplot2)
library(gridExtra)

setwd("~/scratch/ragweed/vcf/stats/")

# SNPs
snps.QD = as.data.frame(sample(as.numeric(scan("SNPs.QD.txt", what = "vector")), 100000))
colnames(snps.QD) = "V1"
snps.MQ = as.data.frame(sample(as.numeric(scan("SNPs.MQ.txt", what = "vector")), 100000))
colnames(snps.MQ) = "V1"
snps.FS = as.data.frame(sample(as.numeric(scan("SNPs.FS.txt", what = "vector")), 100000))
colnames(snps.FS) = "V1"
snps.SOR = as.data.frame(sample(as.numeric(scan("SNPs.SOR.txt", what = "vector")), 100000))
colnames(snps.SOR) = "V1"
snps.MQRankSum = as.data.frame(sample(as.numeric(scan("SNPs.MQRankSum.txt", what = "vector")), 100000))
colnames(snps.MQRankSum) = "V1"
snps.ReadPosRankSum = as.data.frame(sample(as.numeric(scan("SNPs.ReadPosRankSum.txt", what = "vector")), 100000))
colnames(snps.ReadPosRankSum) = "V1"
snps.DP = as.data.frame(sample(as.numeric(scan("SNPs.DP.txt", what = "vector")), 100000))
colnames(snps.DP) = "V1"

snps.QD.plot = ggplot(snps.QD, aes(x = log10(V1))) + 
	geom_density() + 
	geom_vline(aes(xintercept = log10(2)), color = "red") +
	labs(x = "log10(SNP QD)") +
	theme_classic()

snps.MQ.plot = ggplot(snps.MQ, aes(x = V1)) + 
	geom_density() + 
	geom_vline(aes(xintercept = 40), color = "red") +
	labs(x = "SNP MQ") +
	theme_classic()

snps.FS.plot = ggplot(snps.FS, aes(x = log10(V1))) + 
	geom_density() + 
	geom_vline(aes(xintercept = log10(60)), color = "red") +
	labs(x = "log10(SNP FS)") +
	theme_classic()

snps.SOR.plot = ggplot(snps.SOR, aes(x = log10(V1))) + 
	geom_density() + 
	geom_vline(aes(xintercept = log10(3)), color = "red") +
	labs(x = "log10(SNP SOR)") +
	theme_classic()

snps.MQRankSum.plot = ggplot(snps.MQRankSum, aes(x = V1)) + 
	geom_density() + 
	geom_vline(aes(xintercept = -12.5), color = "red") +
	labs(x = "SNP MQRankSum") +
	theme_classic()

snps.ReadPosRankSum.plot = ggplot(snps.ReadPosRankSum, aes(x = V1)) + 
	geom_density() + 
	geom_vline(aes(xintercept = -8), color = "red") +
	labs(x = "SNP ReadPosRankSum") +
	theme_classic()

snps.DP.plot = ggplot(snps.DP, aes(x = V1)) + 
	geom_density() +
	geom_vline(aes(xintercept = mean(V1) + 1.5 * sd(V1)), color = "red") +
	geom_vline(aes(xintercept = mean(V1) - 1 * sd(V1)), color = "red") +
	labs(x = "SNP DP") +
	theme_classic()

pdf(file = "SNPs-filter-406.pdf",
  width = 14,
  height = 10.5,
  useDingbats=FALSE)
grid.arrange(snps.QD.plot, snps.DP.plot, snps.MQ.plot, snps.FS.plot, snps.SOR.plot, snps.MQRankSum.plot, snps.ReadPosRankSum.plot, nrow = 3)
dev.off()

# indels
indels.QD = as.data.frame(sample(as.numeric(scan("indels.QD.txt", what = "vector")), 100000))
colnames(indels.QD) = "V1"
indels.FS = as.data.frame(sample(as.numeric(scan("indels.FS.txt", what = "vector")), 100000))
colnames(indels.FS) = "V1"
indels.SOR = as.data.frame(sample(as.numeric(scan("indels.SOR.txt", what = "vector")), 100000))
colnames(indels.SOR) = "V1"
indels.ReadPosRankSum = as.data.frame(sample(as.numeric(scan("indels.ReadPosRankSum.txt", what = "vector")), 100000))
colnames(indels.ReadPosRankSum) = "V1"
indels.InbreedingCoeff = as.data.frame(sample(as.numeric(scan("indels.InbreedingCoeff.txt", what = "vector")), 100000))
colnames(indels.InbreedingCoeff) = "V1"
indels.DP = as.data.frame(sample(as.numeric(scan("indels.DP.txt", what = "vector")), 100000))
colnames(indels.DP) = "V1"

indels.QD.plot = ggplot(indels.QD, aes(x = log10(V1))) + 
	geom_density() + 
	geom_vline(aes(xintercept = log10(2)), color = "red") +
	labs(x = "log10(indel QD)") +
	theme_classic()

indels.FS.plot = ggplot(indels.FS, aes(x = log10(V1))) + 
	geom_density() + 
	geom_vline(aes(xintercept = log10(200)), color = "red") +
	labs(x = "log10(indel FS)") +
	theme_classic()

indels.SOR.plot = ggplot(indels.SOR, aes(x = log10(V1))) + 
	geom_density() + 
	geom_vline(aes(xintercept = log10(10)), color = "red") +
	labs(x = "log10(indel SOR)") +
	theme_classic()

indels.ReadPosRankSum.plot = ggplot(indels.ReadPosRankSum, aes(x = V1)) + 
	geom_density() + 
	geom_vline(aes(xintercept = -20), color = "red") +
	labs(x = "indel ReadPosRankSum") +
	theme_classic()

indels.InbreedingCoeff.plot = ggplot(indels.InbreedingCoeff, aes(x = V1)) + 
	geom_density() + 
	geom_vline(aes(xintercept = -0.8), color = "red") +
	labs(x = "indel InbreedingCoeff") +
	theme_classic()

indels.DP.plot = ggplot(indels.DP, aes(x = V1)) + 
	geom_density() +
	geom_vline(aes(xintercept = mean(V1) + 1.5 * sd(V1)), color = "red") +
	geom_vline(aes(xintercept = mean(V1) - 1 * sd(V1)), color = "red") +
	labs(x = "indel DP") +
	theme_classic()

pdf(file = "indels-filter-406.pdf",
  width = 14,
  height = 7,
  useDingbats=FALSE)
grid.arrange(indels.QD.plot, indels.DP.plot, indels.FS.plot, indels.SOR.plot, indels.ReadPosRankSum.plot, indels.InbreedingCoeff.plot, nrow = 2)
dev.off()

###

# array job to tag variants for filtering
# vcf-tag.sh
#!/bin/bash
#SBATCH --job-name=vcf-tag
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=00:20:00
#SBATCH --output=vcf-tag-%a.out
#SBATCH --error=vcf-tag-%a.err
#SBATCH --array=1-1124

cd ~/scratch/ragweed/

N=$SLURM_ARRAY_TASK_ID

module load nixpkgs/16.09
module load gatk/3.8

java \
-jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ref/ragweed-dipasm-hap1.fasta \
-V vcf/window-$N.SNPs.recode.vcf \
-filter "QD < 2.0" --filterName "QD2" \
-filter "QUAL < 30.0" --filterName "QUAL30" \
-filter "SOR > 3.0" --filterName "SOR3" \
-filter "FS > 60.0" --filterName "FS60" \
-filter "MQ < 40.0" --filterName "MQ40" \
-filter "MQRankSum < -12.5" --filterName "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum-8" \
-filter "DP > 5776.69" --filterName "DP-5776.69" \
-filter "DP < 551.93" --filterName "DP-551.93" \
-o vcf/window-$N.SNPs.HardFilterTag.vcf

java \
-jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ref/ragweed-dipasm-hap1.fasta \
-V vcf/window-$N.indels.recode.vcf \
-filter "QD < 2.0" --filterName "QD2" \
-filter "QUAL < 30.0" --filterName "QUAL30" \
-filter "FS > 200.0" --filterName "FS200" \
-filter "ReadPosRankSum < -20.0" --filterName "ReadPosRankSum-20" \
-filter "SOR > 10.0" --filterName "SOR-10" \
-filter "InbreedingCoeff < -0.8" --filterName "InbreedingCoeff--0.8" \
-filter "DP > 4470.71" --filterName "DP-4470.71" \
-filter "DP < 743.20" --filterName "DP-743.20" \
-o vcf/window-$N.indels.HardFilterTag.vcf

###

# array job to remove filtered variants plus additional filters
# remove variants called in less than 80% of samples
# remove variants with a genotype depth < 3

# vcf-filter.sh
#!/bin/bash
#SBATCH --job-name=vcf-filter
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=00:05:00
#SBATCH --output=vcf-filter-%a.out
#SBATCH --error=vcf-filter-%a.err
#SBATCH --array=1-1124

cd ~/scratch/ragweed/

N=$SLURM_ARRAY_TASK_ID

module load StdEnv/2020
module load vcftools

vcftools --vcf vcf/window-$N.SNPs.HardFilterTag.vcf \
--max-missing 0.8 \
--minDP 3 \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out vcf/window-$N.SNPs.HardFilter

vcftools --vcf vcf/window-$N.indels.HardFilterTag.vcf \
--max-missing 0.8 \
--minDP 3 \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out vcf/window-$N.indels.HardFilter

###

# merge filtered VCF chunks into single a VCF
ls vcf/*.HardFilter.recode.vcf > vcf-HardFilter.list

# vcf-merge.sh
#!/bin/bash
#SBATCH --job-name=vcf-merge
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=04:00:00
#SBATCH --output=vcf-merge.out
#SBATCH --error=vcf-merge.err

cd ~/scratch/ragweed/

module load StdEnv/2020
module load picard

JAVA_TOOL_OPTIONS="-Xmx8g"

java -Xmx8g -jar $EBROOTPICARD/picard.jar MergeVcfs \
	I=vcf-HardFilter.list \
	O=aa406.HardFilter.vcf

###

# list missingness for each sample
vcftools --vcf aa406.HardFilter.vcf \
--out aa406.HardFilter \
--missing-indv

# exclude samples with > 60% missing
cat aa406.HardFilter.imiss | awk '$5 > 0.6 {print $1}' | tail -n +2 > ifilter.txt

vcftools --vcf aa406.HardFilter.vcf \
--remove ifilter.txt \
--mac 1 \
--recode \
--recode-INFO-all \
--out aa398.HardFilter

###

# fix mislabeled samples 
# change MO-18 to MP-18
# change MOB-1 to MOB-8
cat aa398.HardFilter.recode.vcf | grep "^##" > aa398.unimputed.vcf
cat aa398.HardFilter.recode.vcf | grep "^#CHROM" | head -n 1 | sed 's/MO-18/MP-18/' | sed 's/MOB-1/MOB-8/' >> aa398.unimputed.vcf
cat aa398.HardFilter.recode.vcf | grep -v "^#" >> aa398.unimputed.vcf

# get list of scaffols with variants called
cat aa398.unimputed.vcf | grep -v "^#" | awk '{print $1}' | sort -u > vcontiglist.txt

# an array job to impute each scaffold separately with beagle
# impute.sh
#!/bin/bash
#SBATCH --job-name=impute
#SBATCH --account=def-rieseber
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --output=impute-%a.out
#SBATCH --error=impute-%a.err
#SBATCH --array=1-24

cd ~/scratch/ragweed/

module load beagle/5.4

N=$SLURM_ARRAY_TASK_ID
SCAF=$(cat vcontiglist.txt | head -n $N | tail -n 1)

java -Xmx32g -jar ${EBROOTBEAGLE}/beagle.22Jul22.46e.jar \
out=$SCAF.imputed \
gt=aa398.unimputed.vcf \
nthreads=1 \
window=20.0 \
chrom=$SCAF

###

# index imputed contigs
for i in h1s*.vcf.gz; do echo $i; bcftools index $i; done

# merge imputed VCFs
bcftools concat h1s*.vcf.gz > aa398.imputed.vcf

# MAF < 0.05
# SNPs only
vcftools --vcf aa398.imputed.vcf \
--maf 0.05 \
--remove-indels \
--recode \
--recode-INFO-all \
--out aa398.imputed.maf.0.05


# remove EU2-14-1 and refilter
echo EU2-14-1 >> ifilter2.txt

vcftools --vcf aa398.imputed.maf.0.05.recode.vcf \
--remove ifilter2.txt \
--mac 1 \
--maf 0.05 \
--recode \
--recode-INFO-all \
--out aa397.imputed.maf.0.05




