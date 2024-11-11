# BAYPASS

# get WorldClim variables for each population
# locally
library(raster)
library(dplyr)
library(tidyr)

setwd("~/Dropbox/HODGINSLAB/RAGWEED4")

meta = read.table("samples-meta-443.txt", header = F)
colnames(meta) = c("sample", "pop", "subrange", "range", "lat", "lon")

wcbio = getData("worldclim", var = "bio", res = 2.5)
meta = cbind(meta, extract(wcbio, meta[, c("lon", "lat")]))
meta = meta[-which(meta$pop == "-"), ]
env = meta %>% group_by(pop) %>% summarise(
	bio1 = mean(bio1),
	bio2 = mean(bio2),
	bio3 = mean(bio3),
	bio4 = mean(bio4),
	bio5 = mean(bio5),
	bio6 = mean(bio6),
	bio7 = mean(bio7),
	bio8 = mean(bio8),
	bio9 = mean(bio9),
	bio10 = mean(bio10),
	bio11 = mean(bio11),
	bio12 = mean(bio12),
	bio13 = mean(bio13),
	bio14 = mean(bio14),
	bio15 = mean(bio15),
	bio16 = mean(bio16),
	bio17 = mean(bio17),
	bio18 = mean(bio18),
	bio19 = mean(bio19))

write.table(env, file = "pops-env-443.txt", col.names = T, row.names = F, quote = F)

# select variables
pdf(file = "pairs.panels-nam.pdf",
	width = 7,
	height = 7)
pairs.panels(env[which(env$pop %in% meta[which(meta$range == "nam"), "pop"]), 2:20])
dev.off()

pdf(file = "pairs.panels-eur.pdf",
	width = 7,
	height = 7)
pairs.panels(env[which(env$pop %in% meta[which(meta$range == "eur"), "pop"]), 2:20])
dev.off()

pdf(file = "pairs.panels-aus.pdf",
	width = 7,
	height = 7)
pairs.panels(env[which(env$pop %in% meta[which(meta$range == "aus"), "pop"]), 2:20])
dev.off()

# keep BIO1, BIO2, BIO12, BIO15
# r < 0.7 in all three ranges



# format metadata
# filter for populations N >= 2
# subset by individual ranges and native-invasive contrasts

# in R
#library(tidyverse)
library(dplyr)

setwd("~/scratch/ragweed/")

meta = read.table("samples-meta-443.txt")
colnames(meta) = c("sample", "pop", "subrange", "range", "lat", "lon")

pops = meta %>% group_by(pop) %>% summarize(n = n())

# remove pop = "-" (no population)
pops = pops[pops$pop != "-", ]

# keep only samples which are in populations
meta = meta[which(meta$pop %in% pops$pop), ]

# subset and write out metadata for each baypass run
nam = subset(meta, range == "nam")
write.table(nam, file = "baypass-nam-meta.txt", col.names = F, row.names = F, quote = F)

eur = subset(meta, range == "eur")
write.table(eur, file = "baypass-eur-meta.txt", col.names = F, row.names = F, quote = F)

aus = subset(meta, range == "aus")
write.table(aus, file = "baypass-aus-meta.txt", col.names = F, row.names = F, quote = F)

#nam_eur = subset(meta, range == "nam" | range == "eur")
#write.table(nam_eur, file = "baypass-nam_eur-meta.txt", col.names = F, row.names = F, quote = F)

#nam_aus = subset(meta, range == "nam" | range == "aus")
#write.table(nam_aus, file = "baypass-nam_aus-meta.txt", col.names = F, row.names = F, quote = F)

nam_eur_2 = subset(meta, subrange != "nam.S" & range != "aus")
write.table(nam_eur_2, file = "baypass-nam_eur_2-meta.txt", col.names = F, row.names = F, quote = F)

#nam_aus_2 = subset(meta, subrange == "nam.S" | range == "aus" & pop != "AU01")
#write.table(nam_aus_2, file = "baypass-nam_aus_2-meta.txt", col.names = F, row.names = F, quote = F)

### ISSUE WITH CLUSTER SPLITTING POPULATIONS!
# FIXED HERE
nam_aus_3 = subset(meta, subrange == "nam.S" | subrange == "nam.ME" | range == "aus" & pop != "AU01")

# keep only samples which are in populations with n > 2 AFTER SUBSETTING
nam_aus_3_pops = nam_aus_3 %>% group_by(pop) %>% summarize(n = n())
nam_aus_3_pops = nam_aus_3_pops[which(nam_aus_3_pops$n > 1), ]
nam_aus_3 = nam_aus_3[which(nam_aus_3$pop %in% nam_aus_3_pops$pop), ]

write.table(nam_aus_3, file = "baypass-nam_aus_3-meta.txt", col.names = F, row.names = F, quote = F)

# make a vector of which populations are in North America (-1) and Australia (1)
nam_aus_3_contrast = rep(NA, nrow(nam_aus_3_pops))
nam_aus_3_contrast[which(nam_aus_3_pops$pop %in% nam_aus_3[nam_aus_3$range == "nam", ]$pop)] = -1
nam_aus_3_contrast[which(nam_aus_3_pops$pop %in% nam_aus_3[nam_aus_3$range == "aus", ]$pop)] = 1

# write out
cat(nam_aus_3_contrast, file = "baypass2/nam_aus_3.contrast")


###

# make bamlists
cd ~/scratch/ragweed/

#for range in nam eur aus nam_eur nam_aus
#for range in nam_eur_2 nam_aus_2
#for range in nam_aus_3
for range in nam eur aus
do
cat baypass-$range-meta.txt | awk '{print $1}' \
| while read samp; do ls bam/$samp.ragweed2021.realigned.bam; done > bamlist-baypass-$range.txt
done

# make population lists
cd ~/scratch/ragweed/

#for range in nam eur aus nam_eur nam_aus
#for range in nam_eur_2 nam_aus_2
#for range in nam_aus_3
for range in nam eur aus
do
cat baypass-$range-meta.txt | awk '{print $2}' | sort -u > poplist-baypass-$range.txt
done

wc -l bamlist-baypass-*.txt
#   95 bamlist-baypass-aus.txt
#  160 bamlist-baypass-eur.txt
#  109 bamlist-baypass-nam_aus_2.txt
#  157 bamlist-baypass-nam_aus_3.txt
#  270 bamlist-baypass-nam_aus.txt
#  314 bamlist-baypass-nam_eur_2.txt
#  335 bamlist-baypass-nam_eur.txt
#  175 bamlist-baypass-nam.txt

wc -l poplist-baypass-*.txt
#  11 poplist-baypass-aus.txt
#  39 poplist-baypass-eur.txt
#  16 poplist-baypass-nam_aus_2.txt
#  37 poplist-baypass-nam_aus_3.txt
#  65 poplist-baypass-nam_aus.txt
#  87 poplist-baypass-nam_eur_2.txt
#  93 poplist-baypass-nam_eur.txt
#  54 poplist-baypass-nam.txt

# an R script to filter beagle files for samples in individual populations
# baypass-pop-separate.R
library(data.table)

setwd("~/scratch/ragweed/")

i = commandArgs(trailingOnly = TRUE)[1]
range = commandArgs(trailingOnly = TRUE)[2]

# read in sample metadata
meta = read.table(paste0("~/scratch/ragweed/baypass-", range, "-meta.txt"))
colnames(meta) = c("sample", "pop", "subrange", "range", "lat", "lon")

# read in pops list
pops = scan(paste0("poplist-baypass-", range, ".txt"), what = 'list')

# read in beagle file
beagle = fread(paste0("baypass/", range, "-", i, ".beagle.gz"), header = T)

for (j in 1:length(pops)){

keep.cols = sort(c(1:3,
	which(meta$pop %in% pops[j]) * 3 + 1,
	which(meta$pop %in% pops[j]) * 3 + 2,
	which(meta$pop %in% pops[j]) * 3 + 3))

pop.beagle = beagle[, ..keep.cols]

# write out population beagle file
fwrite(pop.beagle, file = paste0("baypass/", range, "-", pops[j], "-", i, ".beagle.gz"), col.names = F, row.names = F, quote = F, sep = "\t", compress = "gzip")

}

###

# an R script to merge populations and make BayPass input files
# baypass-pop-merge.R
library(data.table)

setwd("~/scratch/ragweed/")

i = commandArgs(trailingOnly = TRUE)[1]
range = commandArgs(trailingOnly = TRUE)[2]

# read in pops list
pops = scan(paste0("poplist-baypass-", range, ".txt"), what = 'list')

# read in sample metadata
meta = read.table(paste0("~/scratch/ragweed/baypass-", range, "-meta.txt"))
colnames(meta) = c("sample", "pop", "subrange", "range", "lat", "lon")

# loop over populations
for (j in 1:length(pops)) {

# get sample size for population
pop.n = sum(meta$pop %in% pops[j])

# read in allele frequencies
mafs = fread(paste0("baypass/", range, "-", pops[j], "-", i, ".mafs.gz"), header = T)

mafs = as.data.frame(cbind(paste0(mafs$chromo, ":", mafs$position), mafs$PPmaf))

# get allele counts based on population sample size
mafs$V2 = round(as.numeric(mafs$V2) * pop.n * 2)

# get alternate allele count
mafs$V3 = pop.n * 2 - mafs$V2

# name columns
colnames(mafs)[2:3] = c(paste0(pops[j], "-G1"), paste0(pops[j], "-G2"))

# merge with other populations
if (j == 1) {
mafs.merged = mafs } else { mafs.merged = cbind(mafs.merged, mafs[, 2:3]) }

write.table(mafs.merged[, 1], file = paste0("baypass2/", range, "-", i, ".snps"), col.names = F, row.names = F, quote = F)
write.table(mafs.merged[, -1], file = paste0("baypass2/", range, "-", i, ".gt"), col.names = F, row.names = F, quote = F)

}

###

# RUN THIS SCRIPT FOR EACH OF nam, eur, aus, nam_eur, nam_aus
# MANUALLY SET minInd flag to 0.75 * N

# a-bp-aus.sh
#!/bin/bash
#SBATCH --job-name=a-bp-aus
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=03:00:00
#SBATCH --output=a-bp-aus-%a.out
#SBATCH --error=a-bp-aus-%a.err
#SBATCH --array=1-113

cd ~/scratch/ragweed/

N=$SLURM_ARRAY_TASK_ID
WIND=$(cat windows-chrs-only-10Mb.list | head -n $N | tail -n 1 | awk '{print $1 ":" $2 "-" $3}')

# run angsd to get genotype probabilities for each sample
module load nixpkgs/16.09
module load intel/2018.3
module load angsd/0.929

angsd \
-nThreads 8 \
-bam bamlist-baypass-aus.txt \
-out baypass/aus-$N \
-r $WIND -GL 2 -doMajorMinor 1 -doCounts 1 -doGLF 2 -SNP_pval 1e-6 -doMaf 2 -doGeno -1 -doPost 1 -minMapQ 30 -minQ 20 -trim 5 -minMaf 0.05 -geno_minDepth 2 -setMinDepthInd 2 -uniqueOnly 1 \
-minInd 71

rm baypass/aus-$N.mafs.gz
rm baypass/aus-$N.arg

# make beagle files for each population
module load StdEnv/2020
module load r/4.3.1

Rscript baypass-pop-separate.R $N aus

# get allele frequencies for each population
module load nixpkgs/16.09
module load intel/2018.3
module load angsd/0.929

cat poplist-baypass-aus.txt | while read pop
do
angsd \
-nThreads 8 \
-beagle baypass/aus-$pop-$N.beagle.gz \
-out baypass/aus-$pop-$N \
-fai ref/ragweed-dipasm-hap1.fasta.fai \
-doMaf 4 \
-minInd 2

rm baypass/aus-$pop-$N.arg
rm baypass/aus-$pop-$N.beagle.gz

done

# merge populations and make BayPass input files
module load StdEnv/2020
module load r/4.3.1

Rscript baypass-pop-merge.R $N aus

cat poplist-baypass-aus.txt | while read pop
do
rm baypass/aus-$pop-$N.mafs.gz
done

###

# Omega matrix

# an R script to downsample BayPass inputs for Omega matrix
# omega-filter.R
library(data.table)

setwd("~/scratch/ragweed/")

i = commandArgs(trailingOnly = TRUE)[1]
range = commandArgs(trailingOnly = TRUE)[2]

# read in genotype file
gt = fread(paste0("baypass2/", range, "-", i, ".gt"), header = F)
snps = scan(paste0("baypass2/", range, "-", i, ".snps"), what = 'list')

# read in omega matrix SNP list
filt.in = scan(paste0("baypass2/", range, "-omega-snps.txt"), what = 'list')

# filter genotype file
gt.filt = gt[which(snps %in% filt.in), ]

# write out filtered genotype file
write.table(gt.filt, file = paste0("baypass2/", range, "-", i, "-omega.gt"), col.names = F, row.names = F, quote = F)

###

module load r/4.3.1

#for range in nam eur aus nam_eur nam_aus
#for range in nam_eur_2 nam_aus_2
#for range in nam_aus_3
for range in nam eur aus
do

# make bedfile of sites
for N in {1..113}
do cat baypass2/$range-$N.snps | awk -F ":" '{print $1 "\t" $2 "\t" $2}'
done > baypass2/$range-snps.bed

# genes.bed
# a list of genes and haploblock regions

# select 10000 random sites outside of genes and haploblocks
module load StdEnv/2020
module load bedtools/2.30.0

bedtools intersect -a baypass2/$range-snps.bed -b genes.bed -v | awk -F "\t" '{print $1 ":" $2}' | shuf -n 10000 > baypass2/$range-omega-snps.txt

for N in {1..113}
do Rscript omega-filter.R $N $range
done

for N in {1..113}
do cat baypass2/$range-$N-omega.gt
done > baypass2/$range-omega.gt

done

#for range in nam eur aus nam_eur nam_aus
#for range in nam_eur_2 nam_aus_2
for range in nam_aus_3
do

module load StdEnv/2020 
module load gcc/9.3.0
module load baypass/2.2

g_baypass -gfile baypass2/$range-omega.gt -outprefix baypass2/$range-omega

done

# run BayPass
# RUN THIS SCRIPT FOR EACH OF nam, eur, aus

# bp-nam.sh
#!/bin/bash
#SBATCH --job-name=bp-nam
#SBATCH --account=def-rieseber
#SBATCH --mem-per-cpu=10G
#SBATCH --time=48:00:00
#SBATCH --output=bp-nam-%a.out
#SBATCH --error=bp-nam-%a.err
#SBATCH --array=1-113

cd ~/scratch/ragweed/

N=$SLURM_ARRAY_TASK_ID

module load StdEnv/2020 
module load gcc/9.3.0
module load baypass/2.2

g_baypass -gfile baypass2/nam-$N.gt -omegafile baypass2/nam-omega_mat_omega.out -outprefix baypass2/nam-$N

###

# an R script to run EAAs for each chunk and prepare WZA input files
# eaa-wza.R
library(data.table)
library(tidyverse)

setwd("~/scratch/ragweed/")

i = commandArgs(trailingOnly = TRUE)[1]
range = commandArgs(trailingOnly = TRUE)[2]

# read in pops list
pops = scan(paste0("poplist-baypass-", range, ".txt"), what = 'list')

# read in sample metadata
meta = read.table(paste0("~/scratch/ragweed/baypass-", range, "-meta.txt"))
colnames(meta) = c("sample", "pop", "subrange", "range", "lat", "lon")

# read in environmental variables for each population
env = read.table("pops-env-443.txt", header = T)

# get only populations in this range
env = env[which(env$pop %in% pops), -1]

# summarize metadata by population
meta.sum = as.data.frame(meta %>% group_by(pop) %>% summarize(n = n(), lat = mean(lat), lon = mean(lon)))

# read in BayPass input
gt = as.data.frame(fread(paste0("baypass2/", range, "-", i, ".gt"), header = F))

# loop over populations
for (j in 1:length(pops)) {

# select the first of two genotypes for a population and calculate allele frequency
afs = gt[, 2 * j - 1] / (2 * meta.sum[meta.sum$pop == pops[j], "n"])

# merge with other populations
if (j == 1) {
afs.merged = afs 
} else {
afs.merged = cbind(afs.merged, afs)
}

}

# WZA data frame
# read in SNP names
snps = as.data.frame(fread(paste0("baypass2/", range, "-", i, ".snps"), header = F))

# divide into 10kb windows
snps = snps %>% separate(V1, sep = ":", into = c("CHR", "LOC"), remove = F)
colnames(snps)[1] = "SNP"
snps$WIND = paste0(snps$CHR, ":", floor(as.numeric(snps$LOC) / 10000) * 10000 + 5000)

# minor allele frequency for WZA
mafs = rowMeans(afs.merged)
mafs[which(mafs > 0.5)] = 1 - mafs[which(mafs > 0.5)]

# read in XTX
xtx = as.data.frame(fread(paste0("baypass2/", range, "-", i, "_summary_pi_xtx.out"), header = T))$XtXst

# write out XTX data for WZA
write.table(as.data.frame(cbind(snps$SNP, snps$WIND, xtx, mafs)), file = paste0("baypass3/", range, "-xtx-", i, ".wzain"), col.names = F, row.names = F, quote = F)

# run EAA for each WorldClim variable and write out WZA input data
for (v in 1:19){
cors = apply(afs.merged, 1, function(x) cor.test(as.numeric(x), as.numeric(env[, v]), method = "kendall", exact = FALSE)$estimate)
write.table(as.data.frame(cbind(snps$SNP, snps$WIND, abs(cors), mafs)), file = paste0("baypass3/", range, "-BIO", v, "-", i, ".wzain"), col.names = F, row.names = F, quote = F)
}

###

# run EAAs and prepare data for WZA
# RUN THIS SCRIPT FOR EACH OF nam, eur, aus

# eaa-aus.sh
#!/bin/bash
#SBATCH --job-name=eaa-aus
#SBATCH --account=def-rieseber
#SBATCH --mem-per-cpu=2G
#SBATCH --time=3:00:00
#SBATCH --output=eaa-aus-%a.out
#SBATCH --error=eaa-aus-%a.err
#SBATCH --array=1-113

cd ~/scratch/ragweed/

N=$SLURM_ARRAY_TASK_ID

module load StdEnv/2020 
module load r/4.3.1

Rscript eaa-wza.R $N aus

###

# merge WZA input parts
for range in nam eur aus
do

# XtX
echo "SNP,WIND,XTX,MAF" > baypass3/${range}-xtx.wzain.csv
for N in {1..113}
do
cat baypass3/${range}-xtx-${N}.wzain | sed 's/ /,/g' >> baypass3/${range}-xtx.wzain.csv
done

# EAA
for v in {1..19}
do
echo "SNP,WIND,TAU,MAF" > baypass3/${range}-BIO${v}.wzain.csv
for N in {1..113}
do
cat baypass3/${range}-BIO${v}-${N}.wzain | sed 's/ /,/g' >> baypass3/${range}-BIO${v}.wzain.csv
done
done

# run the WZA
module load scipy-stack/2023a

python3 general_WZA_script.py \
--correlations baypass3/${range}-xtx.wzain.csv \
--summary_stat XTX \
--large_i_small_p \
--window WIND \
--MAF MAF \
--output baypass3/${range}-xtx.wzaout.csv \
--sep ","

for v in {1..19}
do
python3 general_WZA_script.py \
--correlations baypass3/${range}-BIO${v}.wzain.csv \
--summary_stat TAU \
--large_i_small_p \
--window WIND \
--MAF MAF \
--output baypass3/${range}-BIO${v}.wzaout.csv \
--sep ","
done

done



# CONTRASTS
# run BayPass contrast
# RUN THIS SCRIPT FOR EACH OF nam_eur, nam_aus, nam_eur_2, nam_aus_2, nam_aus_3

# bp-nam_aus_3.sh
#!/bin/bash
#SBATCH --job-name=bp-nam_aus_3
#SBATCH --account=def-rieseber
#SBATCH --mem-per-cpu=10G
#SBATCH --time=92:00:00
#SBATCH --output=bp-nam_aus_3-%a.out
#SBATCH --error=bp-nam_aus_3-%a.err
#SBATCH --array=1-113

cd ~/scratch/ragweed/

N=$SLURM_ARRAY_TASK_ID

module load StdEnv/2020 
module load gcc/9.3.0
module load baypass/2.2

g_baypass -gfile baypass2/nam_aus_3-$N.gt -omegafile baypass2/nam_aus_3-omega_mat_omega.out -contrastfile baypass2/nam_aus_3.contrast -outprefix baypass2/nam_aus_3-$N

###

# an R script to prepare WZA input files for contrast
# con-wza.R
library(data.table)
library(tidyverse)

setwd("~/scratch/ragweed/")

i = commandArgs(trailingOnly = TRUE)[1]
range = commandArgs(trailingOnly = TRUE)[2]

# read in pops list
pops = scan(paste0("poplist-baypass-", range, ".txt"), what = 'list')

# read in sample metadata
meta = read.table(paste0("~/scratch/ragweed/baypass-", range, "-meta.txt"))
colnames(meta) = c("sample", "pop", "subrange", "range", "lat", "lon")

# summarize metadata by population
meta.sum = as.data.frame(meta %>% group_by(pop) %>% summarize(n = n()))

# read in BayPass input
gt = as.data.frame(fread(paste0("baypass2/", range, "-", i, ".gt"), header = F))

# loop over populations
for (j in 1:length(pops)) {

# select the first of two genotypes for a population and calculate allele frequency
afs = gt[, 2 * j - 1] / (2 * meta.sum[meta.sum$pop == pops[j], "n"])

# merge with other populations
if (j == 1) {
afs.merged = afs 
} else {
afs.merged = cbind(afs.merged, afs)
}

}

# WZA data frame
# read in SNP names
snps = as.data.frame(fread(paste0("baypass2/", range, "-", i, ".snps"), header = F))

# divide into 10kb windows
snps = snps %>% separate(V1, sep = ":", into = c("CHR", "LOC"), remove = F)
colnames(snps)[1] = "SNP"
snps$WIND = paste0(snps$CHR, ":", floor(as.numeric(snps$LOC) / 10000) * 10000 + 5000)

# minor allele frequency for WZA
mafs = rowMeans(afs.merged)
mafs[which(mafs > 0.5)] = 1 - mafs[which(mafs > 0.5)]

# read in contrast
con = as.data.frame(fread(paste0("baypass2/", range, "-", i, "_summary_contrast.out"), header = T))$C2_std

# write out contrast data for WZA
write.table(as.data.frame(cbind(snps$SNP, snps$WIND, con, mafs)), file = paste0("baypass3/", range, "-c-", i, ".wzain"), col.names = F, row.names = F, quote = F)

###

# make WZA input parts and merge
module load StdEnv/2020 
module load r/4.3.1

for range in nam_eur nam_aus nam_eur_2 nam_aus_2 nam_aus_3
do
echo "SNP,WIND,C,MAF" > baypass3/${range}-c.wzain.csv
for N in {1..113}
do
Rscript con-wza.R $N $range
cat baypass3/${range}-c-${N}.wzain | sed 's/ /,/g' >> baypass3/${range}-c.wzain.csv
done


# run the WZA
module load scipy-stack/2023a

python3 general_WZA_script.py \
--correlations baypass3/${range}-c.wzain.csv \
--summary_stat C \
--large_i_small_p \
--window WIND \
--MAF MAF \
--output baypass3/${range}-c.wzaout.csv \
--sep ","

done


