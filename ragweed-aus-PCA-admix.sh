#####
# RAGWEED 2024
# PCA AND ADMIXTURE ANALYSES
#####

cd ~/scratch/ragweed/

# install PCAngsd v1.2 locally
module load python
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd/
pip install --user -r requirements.txt
python setup.py build_ext --inplace
pip3 install -e .

# samples-meta-444.txt
# a table with sample name, population, range, lat and lon

# remove EU2-14-1 (EU sample that clusters with AUS)
cat samples-meta-444.txt | awk -F "\t" '$1 != "EU2-14-1" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' \
> samples-meta-443.txt

# make bam list
cat samples-meta-443.txt | awk -F "\t" '{print $1}' \
| while read samp; do ls bam/$samp.ragweed2021.realigned.bam; done > bamlist-443.txt

# make window list
# in R
options(scipen=999)

# read in a fasta index and make a list of 1Mb windows across the genome
dat = read.table("~/scratch/ragweed/ref/ragweed-dipasm-hap1.fasta.fai")[, 1:2]

# only keep chromosomes
dat = dat[dat$V2 > 10000000, ]

# a vector of all contigs
ctgs = dat$V1

for (c in 1:length(ctgs)){
clen = dat[dat$V1 == ctgs[c], 2]
dat.c = data.frame(ctgs[c], seq(1, clen, 10000000), seq(1, clen, 10000000) + 9999999)
dat.c[nrow(dat.c), 3] = clen
if (c == 1){ dat.out = dat.c } else { dat.out = rbind(dat.out, dat.c) }
}

write.table(dat.out, file = "~/scratch/ragweed/windows-chrs-only-10Mb.list", col.names = F, row.names = F, quote = F)
# 113 windows
###

# a-PCA-admix.sh
#!/bin/bash
#SBATCH --job-name=a-PCA-admix
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=03:00:00
#SBATCH --output=a-PCA-admix-%a.out
#SBATCH --error=a-PCA-admix-%a.err
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
-bam bamlist-443.txt \
-r $WIND \
-out angsd/PCA-admix-$N \
-GL 2 \
-doMajorMinor 1 \
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
-minInd 333 \
-geno_minDepth 2 \
-setMinDepthInd 2 \
-uniqueOnly 1 \
-doPlink 2

# LD filter with plink
plink \
--tfile angsd/PCA-admix-$N \
--allow-extra-chr \
--indep-pairwise 50 5 0.5 \
--out angsd/PCA-admix-$N

###

# clean up
for N in {1..113}
do
rm angsd/PCA-admix-$N.arg
rm angsd/PCA-admix-$N.log
rm angsd/PCA-admix-$N.prune.out
rm angsd/PCA-admix-$N.nosex
done

# get a list of sites that are not in genes
# make bedfile of sites
for N in {1..113}
do zcat angsd/PCA-admix-$N.mafs.gz | awk '{print $1 "\t" $2 "\t" $2}' | tail -n +2
done > PCA-admix.bed

# make bedfile of genes
cat ref/genes_filtered.gff \
| awk '$3 == "gene" {print $1 "\t" $4 "\t" $5}' > genes.bed

# add haploblocks
cat haploblocks-24-8-23.txt >> genes.bed

module load StdEnv/2020
module load bedtools/2.30.0

bedtools intersect -a PCA-admix.bed -b genes.bed -v | awk -F "\t" '{print $1 "_" $2}' > PCA-admix.nogenes.txt

###

# an R script to prune for LD and remove sites in genes and haploblocks
# prune-remove-genes.R
library(data.table)

setwd("~/scratch/ragweed/angsd/")

N = commandArgs(trailingOnly = TRUE)[1]

# read in beagle file
beagle = fread(paste0("PCA-admix-", N, ".beagle.gz"), header = T)

# read in LD prune list
prune.in = scan(paste0("PCA-admix-", N, ".prune.in"), what = 'list')

# prune beagle file for LD
beagle.pruned = beagle[which(beagle$marker %in% prune.in), ]

# read in gene filter list
nogenes = scan("~/scratch/ragweed/PCA-admix.nogenes.txt", what = 'list')

# apply gene filter
beagle.pruned.nogenes = beagle.pruned[which(beagle.pruned$marker %in% nogenes), ]

# write out pruned and filtered beagle file
fwrite(beagle.pruned.nogenes, file = paste0("PCA-admix-", N, ".pruned.nogenes.beagle.gz"),
	col.names = F, row.names = F, quote = F, sep = "\t", compress = "gzip", nThread = 36)

###

# prune for LD and remove sites in genes

# prune-remove-genes-PCA-admix.sh
#!/bin/bash
#SBATCH --job-name=prune-remove-genes-PCA-admix
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --time=03:00:00
#SBATCH --output=prune-remove-genes-PCA-admix.out
#SBATCH --error=prune-remove-genes-PCA-admix.err

module load StdEnv/2020
module load r/4.3.1

cd ~/scratch/ragweed/

for N in {1..113}; do Rscript prune-remove-genes.R $N; done

###

# concatenate filtered beagle files and downsample to 100k sites
for N in {1..113}; do zcat angsd/PCA-admix-$N.pruned.nogenes.beagle.gz; done \
| shuf -n 100000 | gzip > PCA-admix.pruned.nogenes.100k.beagle.gz

# run PCAngsd
pcangsd -t 8 -b PCA-admix.pruned.nogenes.100k.beagle.gz -o PCA-admix.pruned.nogenes.100k

###

# plot in R
module load r/4.3.1

library(ggplot2)
library(cowplot)

setwd("~/scratch/ragweed/")

meta = read.table("samples-meta-443.txt", header = F)
colnames(meta) = c("sample", "pop", "subrange", "range", "lat", "lon")

pca = eigen(read.table("PCA-admix.pruned.nogenes.100k.cov"))

pcs = as.data.frame(pca$vectors[, c(1, 2)])

pcs = cbind(meta[, 1:4], pcs)

# percentage of variance explained by each PC
var.ex = round((pca$values / sum(pca$values)) * 100, 2)

pca.subrange = ggplot() +
	geom_point(data = pcs,
		aes(x = V1, y = -V2, fill = subrange), size = 3, shape = 21) +
	ggtitle("") +
	labs(x = paste0("PC1 (", var.ex[1], "%)"),
		y = paste0("PC2 (", var.ex[2], "%)")) +
	scale_fill_manual(name = "cluster",
		values = c("#e31a1c", "#b2df8a", "#33a02c", "#fdbf6f", "#ff7f00", "#a6cee3", "#1f78b4", "#cab2d6", "#6a3d9a"),
		labels = c("Australia", "Europe", "Europe (east)", "Europe (mid east)", "Europe (west)", "North America (east)", "North America (midwest)", "North America (south)", "North America (west)")) +
	theme_classic() +
	theme(legend.position = c(0.85, 0.2),
		plot.title = element_text(size = 28),
		axis.text = element_text(size = 18),
		axis.title = element_text(size = 28),
		legend.title = element_text(size = 28),
		legend.text = element_text(size = 18))

pca.auspops = ggplot() +
	geom_point(data = pcs,
		aes(x = V1, y = -V2), color = "#d9d9d9", size = 3, shape = 19) +
	geom_point(data = subset(pcs, range == "aus"),
		aes(x = V1, y = -V2, fill = pop), size = 3, shape = 21) +
	ggtitle("") +
	labs(x = paste0("PC1 (", var.ex[1], "%)"),
		y = paste0("PC2 (", var.ex[2], "%)")) +
	scale_fill_manual(name = "Australian\npopulation",
		values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99")) +
	theme_classic() +
	theme(legend.position = c(0.85, 0.2),
		plot.title = element_text(size = 28),
		axis.text = element_text(size = 18),
		axis.title = element_text(size = 28),
		legend.title = element_text(size = 28),
		legend.text = element_text(size = 18))

png(file = "~/scratch/ragweed/figs/pca.png",
	width = 1920,
	height = 960)
plot_grid(pca.subrange, pca.auspops, ncol = 2, align = "hv")
dev.off()

###

# NGS admix
# ngsadmix-PCA-admix.sh
#!/bin/bash
#SBATCH --job-name=ngsadmix-PCA-admix
#SBATCH --account=def-rieseber
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=72:00:00
#SBATCH --output=ngsadmix-PCA-admix-%a.out
#SBATCH --error=ngsadmix-PCA-admix-%a.err
#SBATCH --array=1-10

cd ~/scratch/ragweed/

module load StdEnv/2020
module load angsd/0.939

N=$SLURM_ARRAY_TASK_ID

for K in {2..15}
do NGSadmix -likes PCA-admix.pruned.nogenes.100k.beagle.gz -K $K -o PCA-admix.pruned.nogenes.100k-$K-rep$N -P 8
done

###

# get likelihood estimates for each value of K
for N in {1..10}
do
for k in {2..15}
do echo $(cat PCA-admix.pruned.nogenes.100k-$k-rep$N.log | grep "best like" | awk -F "=" '{print $2}' | awk '{print $1}')
done > admix-likes-rep$N.txt
done


# plots in R
module load r/4.3.1

library(ggplot2)
library(forcats)
library(pals)
library(ggthemes)
library(cowplot)

setwd("~/scratch/ragweed/")

# admixture plots

# first determine which rep to plot
kreps = as.data.frame(cbind(
	scan("admix-likes-rep1.txt", what = 'list'),
	scan("admix-likes-rep2.txt", what = 'list'),
	scan("admix-likes-rep3.txt", what = 'list'),
	scan("admix-likes-rep4.txt", what = 'list'),
	scan("admix-likes-rep5.txt", what = 'list'),
	scan("admix-likes-rep6.txt", what = 'list'),
	scan("admix-likes-rep7.txt", what = 'list'),
	scan("admix-likes-rep8.txt", what = 'list'),
	scan("admix-likes-rep9.txt", what = 'list'),
	scan("admix-likes-rep10.txt", what = 'list')
	))

# max likelihood
kmaxlik = c()
for (i in 1:nrow(kreps)){
kmaxlik = c(kmaxlik, max(as.numeric(kreps[i, ])))
}

# rep with max likelihood
kmaxlikrep = c()
for (i in 1:nrow(kreps)){
kmaxlikrep = c(kmaxlikrep, which(as.numeric(kreps[i, ]) %in% kmaxlik[i])[1])
}

meta = read.table("samples-meta-443.txt", header = F)
colnames(meta) = c("sample", "pop", "subrange", "range", "lat", "lon")

for (k in 2:15){

q = read.table(paste0("PCA-admix.pruned.nogenes.100k-", k, "-rep", kmaxlikrep[k-1], ".qopt"))

ks = as.data.frame(matrix(ncol = 6, nrow = 0))

for (i in 1:k){
ktemp = cbind(meta[, 1:4], i, q[, i])
colnames(ktemp)[5:6] = c("grp", "prop")
if (i == 1){ ks = ktemp } else {ks = rbind(ks, ktemp)}
}

k.plot = ggplot(ks[order(ks$subrange), ], aes(factor(sample), prop, fill = factor(grp))) +
	geom_col(color = "gray", size = 0.1) +
	facet_grid(~fct_inorder(subrange), switch = "x", scales = "free", space = "free") +
	theme_minimal() + labs(y = paste0("K=", k)) +
	scale_y_continuous(expand = c(0, 0)) +
	scale_x_discrete(expand = expand_scale(add = 1)) +
	theme(
		panel.spacing.x = unit(0.1, "lines"),
		axis.ticks = element_blank(),
		axis.text.x = element_blank(),
		strip.text.x = element_blank(),
		axis.text.y = element_blank(),
		panel.grid = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 28)
		) +
  	scale_fill_manual(values = glasbey(), guide = 'none')

assign(paste0("k", k, ".plot"), k.plot)

}

png(file = "~/scratch/ragweed/figs/admix.k2-15.png",
	width = 1920,
	height = 120 * 14)
plot_grid(k2.plot, k3.plot, k4.plot, k5.plot, k6.plot, k7.plot, k8.plot, k9.plot, k10.plot, k11.plot, k12.plot, k13.plot, k14.plot, k15.plot, 
		ncol = 1, align = "hv")
dev.off()







#####


# likelihood plots
likes = read.table("PCA-admix.likes.txt", header = F)
colnames(likes) = c("K", "likelihood")

likes$delta = 0
for (i in 2:nrow(likes)){
likes$delta[i] = likes$likelihood[i-1] - likes$likelihood[i]
}

likes.plot = ggplot() +
	geom_point(data = likes,
		aes(x = K, y = likelihood), color = "black", size = 5, shape = 19) +
	ggtitle("") +
	labs(x = "K", 
		y = "likelihood estimate") +
	theme_classic() +
	theme(plot.title = element_text(size = 28),
		axis.text = element_text(size = 14),
		axis.title = element_text(size = 28))

delta.plot = ggplot() +
	geom_point(data = likes,
		aes(x = K, y = delta), color = "black", size = 5, shape = 19) +
	ggtitle("") +
	labs(x = "K", 
		y = "delta") +
	theme_classic() +
	theme(plot.title = element_text(size = 28),
		axis.text = element_text(size = 14),
		axis.title = element_text(size = 28))

png(file = "~/scratch/ragweed/figs/klikes.png",
	width = 1920,
	height = 960)
plot_grid(likes.plot, delta.plot, ncol = 2, align = "hv")
dev.off()























