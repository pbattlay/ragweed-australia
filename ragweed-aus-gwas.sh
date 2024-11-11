#####
# RAGWEED 2024
# GWAS

# PHENOTYPES
# check phenotypes for normality; outliers
# locally in R
setwd("~/Dropbox/HODGINSLAB/RAGWEED4/")

lotte = read.csv("~/Dropbox/HODGINSLAB/RAGWEED4/lottedat.csv", header = T)
sampmatch = read.csv("samplematch-227.csv", header = T)

phen = merge(sampmatch, lotte, by.x = "lotteSampleID", by.y = "SampleID")

# write out sample names
write.table(phen$paulSampleID, file = "NEWNEW-samps-gwas-226.txt", col.names = F, row.names = F, quote = F)

# phenotypes to keep
# T_die_day					Number of days after second transplant no racemes produced pollen
# T_Racemes					Number of racemes
# T_Raceme_length_longest	Length of longest raceme
# T_Longest_leaf			Length of longest leaf
# T_seed_20					Seed weight
# T_shoot_biom				Shoot biomass
# T_tot_biom				Total biomass
# T_fl_end					Flowering end
# T_fl_start				Flowering onset
# T_fit_all					Total reproductive biomass
# T_sexrat_weight			Floral sex allocation
# T_rootshoot_rat			Root/shoot ratio
# T_maxHeight_cm			Maximum height

keep = c("T_die_day", "T_Racemes", "T_Raceme_length_longest", "T_Longest_leaf", "T_seed_20", "T_shoot_biom", "T_tot_biom", "T_fl_end", "T_fl_start", "T_fit_all", "T_sexrat_weight", "T_rootshoot_rat", "T_maxHeight_cm")
phen = phen[, c(which(colnames(phen) %in% keep))]

# remove outliers
phen$T_Longest_leaf[which(phen$T_Longest_leaf <= 20)] = NA

# transform
phen$T_sexrat_weight = log10(phen$T_sexrat_weight)
#phen$T_rootshoot_rat = log10(phen$T_rootshoot_rat)
phen$T_Racemes = phen$T_Racemes ^ 0.5

# set NAs to -999
phen[is.na(phen)] = -999

# write out phenotype files
write.table(phen, file = "NEWNEW-phenos-gwas-226.txt", col.names = F, row.names = F, quote = F)
write.table(colnames(phen), file = "NEWNEW-names-gwas-226.txt", col.names = F, row.names = F, quote = F)

# GWAS
cd ~/scratch/ragweed/

# make bam list
cat NEWNEW-samps-gwas-226.txt | while read samp; do ls bam/$samp.ragweed2021.realigned.bam; done > NEWNEW-bamlist-gwas-226.txt

# a-gwas1.sh
#!/bin/bash
#SBATCH --job-name=a-gwas1
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=03:00:00
#SBATCH --output=a-gwas1-%a.out
#SBATCH --error=a-gwas1-%a.err
#SBATCH --array=1-113

cd ~/scratch/ragweed/

N=$SLURM_ARRAY_TASK_ID
WIND=$(cat windows-chrs-only-10Mb.list | head -n $N | tail -n 1 | awk '{print $1 ":" $2 "-" $3}')

module load nixpkgs/16.09
module load intel/2018.3
module load angsd/0.929
module load plink/1.9b_5.2-x86_64

angsd \
-nThreads 8 \
-bam NEWNEW-bamlist-gwas-226.txt \
-r $WIND \
-out gwas1/gwas-$N \
-GL 2 \
-doMajorMinor 4 \
-ref ref/ragweed-dipasm-hap1.fasta \
-doCounts 1 \
-doGLF 2 \
-SNP_pval 1e-6 \
-doMaf 2 \
-doGeno -1 \
-doPost 1 \
-minMapQ 30 \
-minQ 20 \
-trim 5 \
-minMaf 0.05 \
-minInd 170 \
-geno_minDepth 2 \
-setMinDepthInd 2 \
-uniqueOnly 1 \
-doPlink 2

# LD filter with plink
plink \
--tfile gwas1/gwas-$N \
--allow-extra-chr \
--indep-pairwise 50 5 0.5 \
--out gwas1/gwas-$N

###

# clean up
for N in {1..113}
do
rm gwas1/gwas-$N.arg
rm gwas1/gwas-$N.log
rm gwas1/gwas-$N.prune.out
rm gwas1/gwas-$N.nosex
done

# get a list of sites that are not in genes
# make bedfile of sites
for N in {1..113}
do zcat gwas1/gwas-$N.mafs.gz | awk '{print $1 "\t" $2 "\t" $2}' | tail -n +2
done > gwas1.bed

# already done for PCA
# make bedfile of genes
#cat ref/genes_filtered.gff \
#| awk '$3 == "gene" {print $1 "\t" $4 "\t" $5}' > genes.bed

# add haploblocks
#cat haploblocks-24-8-23.txt >> genes.bed

module load StdEnv/2020
module load bedtools/2.30.0

bedtools intersect -a gwas1.bed -b genes.bed -v | awk -F "\t" '{print $1 "_" $2}' > gwas1.nogenes.txt

###

# an R script to prune for LD and remove sites in genes
# prune-remove-genes-gwas1.R
library(data.table)

setwd("~/scratch/ragweed/gwas1/")

N = commandArgs(trailingOnly = TRUE)[1]

# read in beagle file
beagle = fread(paste0("gwas-", N, ".beagle.gz"), header = T)

# read in LD prune list
prune.in = scan(paste0("gwas-", N, ".prune.in"), what = 'list')

# prune beagle file for LD
beagle.pruned = beagle[which(beagle$marker %in% prune.in), ]

# read in gene filter list
nogenes = scan("~/scratch/ragweed/gwas1.nogenes.txt", what = 'list')

# apply gene filter
beagle.pruned.nogenes = beagle.pruned[which(beagle.pruned$marker %in% nogenes), ]

# write out pruned and filtered beagle file
fwrite(beagle.pruned.nogenes, file = paste0("gwas1-", N, ".pruned.nogenes.beagle.gz"),
	col.names = F, row.names = F, quote = F, sep = "\t", compress = "gzip", nThread = 36)

###

# prune for LD and remove sites in genes and haploblocks
# prune-remove-genes-gwas1.sh
#!/bin/bash
#SBATCH --job-name=prune-remove-genes-gwas1
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --time=03:00:00
#SBATCH --output=prune-remove-genes-gwas1.out
#SBATCH --error=prune-remove-genes-gwas1.err

module load StdEnv/2020
module load r/4.3.1

cd ~/scratch/ragweed/

for N in {1..113}; do Rscript prune-remove-genes-gwas1.R $N; done

###

# concatenate filtered beagle files and downsample to 100k sites
# NOTE: some files will be empty due to haploblocks overlapping the window
for N in {1..113}; do zcat gwas1/gwas1-$N.pruned.nogenes.beagle.gz; done \
| shuf -n 100000 | gzip > gwas1.pruned.nogenes.100k.beagle.gz

# run PCAngsd
pcangsd -t 8 -b gwas1.pruned.nogenes.100k.beagle.gz -o gwas1.pruned.nogenes.100k

# get first two PCs of covariance matrix
# in R
pca = as.data.frame(eigen(read.table("gwas1.pruned.nogenes.100k.cov"))$vectors[, c(1, 2)])
write.table(pca, file = "gwas-226-pc12.cov", col.names = F, row.names = F, quote = F)

# GWAS
wc -l NEWNEW-names-gwas-226.txt
# 13

# an R script to prepare a WZA input file from GWAS output
# gwas1-wzain.R
library(data.table)

setwd("~/scratch/ragweed/gwas1/")

i = commandArgs(trailingOnly = TRUE)[1]
N = commandArgs(trailingOnly = TRUE)[2]

# read in GWAS results
gwas = fread(paste0("NEWNEW-", i, "-", N, ".lrt0.gz"), header = T)
gwas = gwas[, c("Chromosome", "Position", "Frequency", "P")]

# make SNP column
gwas$SNP = paste0(gwas$Chromosome, ":", gwas$Position)

# make WIND column (10kbp windows)
gwas$WIND = paste0(gwas$Chromosome, ":", floor(as.numeric(gwas$Position) / 10000) * 10000 + 5000)

# convert frequency to minor allele frequency
gwas$Frequency[which(gwas$Frequency > 0.5)] = 1 - gwas$Frequency[which(gwas$Frequency > 0.5)]

# write out
write.table(gwas[, c("SNP", "WIND", "P", "Frequency")], file = paste0("NEWNEW-", i, "-", N, ".wzain.csv"), sep = ",", col.names = F, row.names = F, quote = F)

###

# run GWAS in chunks for each phenotype
# concatenate results
# prepare input for and run the WZA
# gwas1.sh
#!/bin/bash
#SBATCH --job-name=gwas1
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=02:00:00
#SBATCH --output=gwas1-%a.out
#SBATCH --error=gwas1-%a.err
#SBATCH --array=1-13

cd ~/scratch/ragweed/

i=$SLURM_ARRAY_TASK_ID

module load StdEnv/2020
module load angsd/0.936
module load r/4.3.1
module load scipy-stack/2023a

cat NEWNEW-phenos-gwas-226.txt | awk -v i=$i '{print $i}' > gwas1/NEWNEW-$i.pheno
pheno=$(cat NEWNEW-names-gwas-226.txt | head -n $i | tail -n 1)

# run GWAS in 10Mb chunks
for N in {1..113}
do
angsd \
-doMaf 4 \
-beagle gwas1/gwas-$N.beagle.gz \
-fai ref/ragweed-dipasm-hap1.fasta.fai \
-yQuant gwas1/NEWNEW-$i.pheno \
-cov gwas-226-pc12.cov \
-Pvalue 1 \
-doAsso 4 \
-out gwas1/NEWNEW-$i-$N \
-P 8

# remove arg files
rm gwas1/NEWNEW-$i-$N.arg

# prepare WZA input in R
Rscript gwas1-wzain.R $i $N

done

# concatenate results
zcat gwas1/NEWNEW-$i-1.lrt0.gz | head -n 1 > gwas1/NEWNEW-$pheno.txt
for N in {1..113}; do zcat gwas1/NEWNEW-$i-$N.lrt0.gz | tail -n +2; done >> gwas1/NEWNEW-$pheno.txt

echo "SNP,WIND,P,MAF" > gwas1/NEWNEW-$pheno.wzain.csv
for N in {1..113}; do cat gwas1/NEWNEW-$i-$N.wzain.csv >> gwas1/NEWNEW-$pheno.wzain.csv; done

# run the WZA
python3 general_WZA_script.py \
--correlations gwas1/NEWNEW-$pheno.wzain.csv \
--summary_stat P \
--window WIND \
--MAF MAF \
--output gwas1/NEWNEW-$pheno.wzaout.csv \
--sep ","

###


# WZA HAPLOBLOCK ENRICHMENT AND MANHATTAN PLOTS
library(tidyverse)

setwd("~/Dropbox/HODGINSLAB/RAGWEED4/")

# read in chromosome data for manhattan plots
contigs = read.table("RESULTS/chromo-length.txt", header = F)
contigs$V3 = c(1:nrow(contigs))
colnames(contigs) = c("CONTIG", "LENGTH", "CHR")
contigs$CHR = as.numeric(contigs$CHR)

# get breaks for manhattan plots
manbreaks = c()
s = 0
for (i in 1:nrow(contigs)) {
	nbp = contigs[contigs$CHR == i,]$LENGTH
	manbreaks = c(manbreaks, round(nbp/2) + s)
	s = s + nbp
}
manbreaksnames = 1:nrow(contigs)

# read in haploblock data; name haploblocks and get cumulative counts
hbs = as.data.frame(read.table("~/Dropbox/HODGINSLAB/RAGWEED4/haploblocks/haploblocks-24-8-23.txt", header = F))
hbs$HB = paste0(hbs$V1, "_", hbs$V2, "_", hbs$V3)
hbc = merge(contigs[, c("CONTIG", "CHR")], hbs, by.x = "CONTIG", by.y = "V1")
s = 0
for (i in 1:length(unique(contigs$CHR))) {
	nbp = contigs[contigs$CHR == i,]$LENGTH
	hbc[hbc$CHR == i, c("V2", "V3")] = hbc[hbc$CHR == i, c("V2", "V3")] + s
	s = s + nbp
}

# read in list of phenotypes
phenos = scan("NEWNEW-names-gwas-226.txt", what = 'list')

# set significance threshold for WZA window p-values
c = 0.001

for (i in 1:length(phenos)){
wza = read.csv(paste0("RESULTS/NEWNEW-", phenos[i], ".wzaout.csv"))

# filter for minimum five SNPs per window
# this is the same for each phenotype, so no need for merging by window
wza = subset(wza, SNPs >= 5)

# get empirical p-value from Z-score
wza$P = rank(-wza$Z) / nrow(wza)

wza = wza[, c("gene", "P")] %>% separate(gene, sep = ":", into = c("CONTIG", "LOC"), remove = F)
colnames(wza)[1] = "WIND"
wza$LOC = as.numeric(wza$LOC)
#wza$PHENO = phenos[i]

# test for enrichment of outlier windows (empirical p < c) in each haploblock
# keep a vector of bonferroni significant haploblocks
enriched = c()

for (h in 1:nrow(hbs)){
hbwinds = subset(wza, CONTIG == hbs[h, 1] & LOC > hbs[h, 2] & LOC < hbs[h, 3])

overlap = nrow(subset(hbwinds, P < c))
r1.out.wind = nrow(hbwinds)
r2.out.wind = nrow(subset(wza, P < c))
wind.tot = nrow(wza)

res = phyper(overlap, r1.out.wind, wind.tot - r1.out.wind, r2.out.wind, lower.tail = FALSE, log.p = FALSE)

# only record bonferroni-significant associations
if (res < 0.05 / 37) {enriched = c(enriched, hbs$HB[h])}

# add results to a vector
if (h == 1) {res.h = res} else {res.h = c(res.h, res)}

}

# add results vectors to a dataframe
if (i == 1) {res.df = res.h} else {res.df = cbind(res.df, res.h)}

# plot out WZA manhattan
# prepare XtX data frame for Manhattan plot
wza = merge(contigs, wza, by = "CONTIG")

# add a cumulative BP count
wza$LOCc = NA
s = 0
for (j in 1:length(unique(wza$CHR))) {
	nbp = max(wza[wza$CHR == j,]$LENGTH)
	wza[wza$CHR == j, "LOCc"] = wza[wza$CHR == j,"LOC"] + s
	s = s + nbp
}

# get a dataframe for enriched haploblocks
hbc.e = hbc[which(hbc$HB %in% enriched), ]

# should be the same ymax for all plots (empirical p-value)
y.max = 5

# GWAS plot
gwas.man = ggplot() +
	geom_point(data = wza,
		aes(x = LOCc, y = -log10(P), color = as.factor(CHR)), size = 0.5, shape = 19) +
	geom_point(data = subset(wza, P < c),
		aes(x = LOCc, y = -log10(P)), color = "red", size = 0.5, shape = 19) +
#	geom_text(aes(x = 675550157, y = -log10(2.8744e-10), label = " ELF3"), size = 7, vjust = "middle", hjust = "left") +
	ggtitle(phenos[i]) +
	scale_color_manual(values = rep(c("#bdbdbd", "#636363"), 500)) +
	scale_x_continuous(name = "chromosome", breaks = manbreaks, labels = manbreaksnames) +
	labs(y = expression(-log[10]*italic(p))) +
	geom_segment(aes(x = hbc[,3], y = y.max, xend = hbc[,4], yend = y.max), colour = "#a6cee3", size = 4) +
	geom_segment(aes(x = hbc.e[,3], y = y.max, xend = hbc.e[,4], yend = y.max), colour = "red", size = 4) +
	geom_hline(yintercept = -log10(c), color = "black", linetype = "dashed") +
	ylim(0, y.max) +
	theme_classic() +
	theme(legend.position = "none",
		plot.title = element_text(size = 28),
		axis.text = element_text(size = 24),
		axis.title = element_text(size = 28)
	)

ggsave(gwas.man,
	file = paste0("junk/", phenos[i], "-wza.man.png"),
	device = "png",
	width = 1920,
	height = 240,
	units = "px",
	dpi = 72)

}

colnames(res.df) = phenos
rownames(res.df) = hbs$HB

write.table(res.df, file = "RESULTS/hb-pheno-wza-enrich.csv", sep = ",", col.names = T, row.names = T, quote = F)

######