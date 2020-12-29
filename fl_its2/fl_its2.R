#Nicola's Florida/Irma ITS2 analysis
#Almost entirely based on DADA2 ITS Pipeline 1.8 Walkthrough:
#https://benjjneb.github.io/dada2/ITS_workflow.html
#with edits by Carly D. Kenkel and modifications for my data by Nicola Kriefall
#6/1/20

#~########################~#
##### PRE-PROCESSING #######
#~########################~#

#fastq files should have R1 & R2 designations for PE reads
#Also - some pre-trimming. Retain only PE reads that match amplicon primer. Remove reads containing Illumina sequencing adapters

##in Terminal home directory:
##following instructions of installing BBtools from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/
##1. download BBMap package, sftp to installation directory
##2. untar: 
#tar -xvzf BBMap_(version).tar.gz
##3. test package:
#cd bbmap
#~/bin/bbmap/stats.sh in=~/bin/bbmap/resources/phix174_ill.ref.fa.gz

## my adaptors, which I saved as "adaptors.fasta"
# >forward
# AATGATACGGCGACCAC
# >forwardrc
# GTGGTCGCCGTATCATT
# >reverse
# CAAGCAGAAGACGGCATAC
# >reverserc
# GTATGCCGTCTTCTGCTTG

##Note: Illumina should have cut these out already, normal if you don't get any

##primers for ITS:
# >forward
# GTGAATTGCAGAACTCCGTG
# >reverse
# CCTCCGCTTACTTATATGCTT

##Still in terminal - making a sample list based on the second phrase before the underscore in the .fastq name
##for my 2015-2017 time points:
#ls *R1_001.fastq | cut -d '_' -f 1 > samples.list
##for my 2018 time point: 
#ls *R1_001.fastq | cut -d '_' -f 2 > samples.list
##(change -f to whatever phrase your sample name appears in)

##cuts off the extra words in the .fastq files
##for my 2015-2017 time points:
#for file in $(cat samples.list); do  mv ${file}_*R1*.fastq ${file}_R1.fastq; mv ${file}_*R2*.fastq ${file}_R2.fastq; done 
##for my 2018 time point: 
#for file in $(cat samples.list); do  mv FL18_${file}_*R1*.fastq ${file}_R1.fastq; mv FL18_${file}_*R2*.fastq ${file}_R2.fastq; done 

##gets rid of reads that still have the adaptor sequence, shouldn't be there, I didn't have any
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1.fastq in2=${file}_R2.fastq ref=adaptors.fasta out1=${file}_R1_NoIll.fastq out2=${file}_R2_NoIll.fastq; done &>bbduk_NoIll.log

##getting rid of first 4 bases (degenerate primers created them)
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1_NoIll.fastq in2=${file}_R2_NoIll.fastq ftl=4 out1=${file}_R1_NoIll_No4N.fastq out2=${file}_R2_NoIll_No4N.fastq; done &>bbduk_No4N.log

##only keeping reads that start with the ITS2 primer
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1_NoIll_No4N.fastq in2=${file}_R2_NoIll_No4N.fastq k=15 restrictleft=21 literal=GTGAATTGCAGAACTCCGTG,CCTCCGCTTACTTATATGCTT outm1=${file}_R1_NoIll_No4N_ITS.fastq outu1=${file}_R1_check.fastq outm2=${file}_R2_NoIll_No4N_ITS.fastq outu2=${file}_R2_check.fastq; done &>bbduk_ITS.log
##higher k = more reads removed, but can't surpass k=20 or 21

##using cutadapt to remove primer:
#module load cutadapt

# for file in $(cat samples.list)
# do
# cutadapt -g GTGAATTGCAGAACTCCGTG -a AAGCATATAAGTAAGCGGAGG -G CCTCCGCTTACTTATATGCTT -A CACGGAGTTCTGCAATTCAC -n 2 --discard-untrimmed -o ${file}_R1.fastq -p ${file}_R2.fastq ${file}_R1_NoIll_No4N_ITS.fastq ${file}_R2_NoIll_No4N_ITS.fastq
# done &> clip.log
##-g regular 5' forward primer 
##-G regular 5' reverse primer
##-o forward out
##-p reverse out
##-max-n 0 means 0 Ns allowed
##this overwrote my original renamed files 

##did sftp of *_R1.fastq & *_R2.fastq files to the folder to be used in dada2

#~########################~#
##### DADA2 BEGINS #########
#~########################~#

#installing/loading packages:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.8")
library(dada2)
#packageVersion("dada2")
#I have version 1.10.1 - tutorial says 1.8 but I think that's OK, can't find a version 1.10 walkthrough
library(ShortRead)
#packageVersion("ShortRead")
library(Biostrings)
#packageVersion("Biostrings")

path <- "/Users/nicolakriefall/fl_its2/rads" # CHANGE ME to the directory containing the fastq files after unzipping.

fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))
#see if these look right
fnFs
fnRs

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)
sample.names

#### check for primers ####
FWD <- "GTGAATTGCAGAACTCCGTG"  ## CHANGE ME to your forward primer sequence
REV <- "CCTCCGCTTACTTATATGCTT"  ## CHANGE ME...

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[6]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[6]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[6]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[6]]))
#all 0s

###### Visualizing raw data

#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs[c(1,3,4,5)])
#samples 4, 7, 8 not working here
#the ones that work look good up to ~200

#Then look at quality profile of R2 reads
plotQualityProfile(fnRs[c(1,3,4,5)])
#all worked except 7 & 8
#the ones that work look good up to ~150

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#changing a bit from default settings - maxEE=1 (1 max expected error, more conservative), truncating length at 200 bp for both forward & reverse, added "trimleft" to cut off primers [20 for forward, 21 for reverse]
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(200,180), #leave overlap
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     minLen = 50,
                     #trimLeft=c(20,21), #N nucleotides to remove from the start of each read
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchIDs=TRUE, #enforce mtching between id-line sequence identifiers of F and R reads
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)

#~############################~#
##### Learn Error Rates ########
#~############################~#

#setDadaOpt(MAX_CONSIST=30) #if necessary, increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE) 
plotErrors(errR, nominalQ=TRUE) 

#~############################~#
##### Dereplicate reads ########
#~############################~#
#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#~###############################~#
##### Infer Sequence Variants #####
#~###############################~#

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[1]]
dadaRs[[1]]

#~############################~#
##### Merge paired reads #######
#~############################~#

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

summary((mergers[[1]]))

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

#~##################################~#
##### Construct sequence table #######
#~##################################~#
#a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
rowSums(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#mostly at 300 bp, which is just what was expected [271-303 bp]

plot(table(nchar(getSequences(seqtab))))

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing

#if I wanted to remove some lengths - not recommended by dada2 for its2 data
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(268,327)] 

table(nchar(getSequences(seqtab)))
dim(seqtab)
plot(table(nchar(getSequences(seqtab))))

#~############################~#
##### Remove chimeras ##########
#~############################~#
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#Identified 965 bimeras out of 1126 input sequences.

sum(seqtab.nochim)/sum(seqtab)
#0.7769868
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. 

#~############################~#
##### Track Read Stats #########
#~############################~#

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="~/Desktop/pdam_2020/its2_reads.csv",row.names=TRUE,quote=FALSE)

#plotting for no reason

#manually added raw counts from the raw .fastq file using the following loop in Terminal:
# for file in *R1_001.fastq
# do
# echo $file >> raw_r1_names
# grep @M0 $file | wc -l >> raw_r1_counts
# done
setwd("~/moorea_holobiont/mr_ITS2/")
reads <- read.csv("its2_reads_renamed.csv")
#counts1 = raw
#counts2 = input
#counts3 = filtered
#counts4 = denoised
#counts5 = merged
#counts6 = nonchim

reread <- reshape(reads, varying = c("counts1","counts2", "counts3", "counts4", "counts5", "counts6"), timevar = "day",idvar = "sample", direction = "long", sep = "")
reread$sample <- as.factor(reread$sample)
ggplot(reread,aes(x=day,y=counts,color=sample))+
  geom_point()+
  geom_path()

library(Rmisc)
reread.se <- summarySE(data=reread,measurevar="counts",groupvars=c("day"))
ggplot(reread.se,aes(x=day,y=counts))+
  geom_point()+
  geom_path()+
  geom_errorbar(aes(ymin=counts-se,ymax=counts+se),width=0.3)

#~############################~#
##### Assign Taxonomy ##########
#~############################~#

taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/GeoSymbio_ITS2_LocalDatabase_verForPhyloseq.fasta",tryRC=TRUE,minBoot=50,verbose=TRUE)
unname(head(taxa))

#to come back to later
saveRDS(seqtab.nochim, file="~/Desktop/pdam_2020/pdam_its2_seqtab.nochim.rds")
saveRDS(taxa, file="~/Desktop/pdam_2020/pdam_its2_taxa.rds")

write.csv(seqtab.nochim, file="~/Desktop/pdam_2020/pdam_its2_seqtab.nochim.csv")
write.csv(taxa, file="~/Desktop/pdam_2020/pdam_its2_taxa.csv")

#### Reading in prior data files ####
setwd("~/Desktop/pdam_2020")
seqtab.nochim <- readRDS("pdam_its2_seqtab.nochim.rds")
taxa <- readRDS("pdam_its2_taxa.rds")

#~############################~#
##### handoff 2 phyloseq #######
#~############################~#

#BiocManager::install("phyloseq")
library('phyloseq')
library('ggplot2')
library('Rmisc')

#import dataframe holding sample information
samdf<-read.csv("samplesnames.csv")
head(samdf)
rownames(samdf) <- samdf$colony

#Construct phyloseq object (straightforward from dada2 outputs)
# ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
#                sample_data(samdf),
#                tax_table(taxa))
# 
# ps

#phyloseq object with shorter names - doing this one instead of one above
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
#making output fasta file for lulu step & maybe other things, before giving new ids to sequences
#path='~/moorea_holobiont/mr_ITS2/mrits2.fasta'
#uniquesToFasta(seqtab.nochim, path, ids = ids, mode = "w", width = 20000)
colnames(seqtab.nochim)<-ids
taxa2 <- cbind(taxa, rownames(taxa)) #retaining raw sequence info before renaming
rownames(taxa2)<-ids

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               #sample_data(samdf), 
               tax_table(taxa))

ps

#### Alpha diversity #####
library(ggplot2)
#install.packages("extrafontdb")
library(extrafontdb)
library(extrafont)
library(Rmisc)
library(cowplot)
library(ggpubr)
#font_import(paths = "/Library/Fonts/")
loadfonts()

#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)
#Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).
#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index (1 − λ).

plot_richness(ps_no87, x="site", measures=c("Shannon", "Simpson"), color="zone") + theme_bw()

df.div <- data.frame(estimate_richness(ps_no87, split=TRUE, measures =c("Shannon","InvSimpson")))
df.div

df.div$Sample <- rownames(df.div)
df.div <- merge(df.div,samdf.no87,by="Sample") #add sample data

#write.csv(df.div,file="~/moorea_holobiont/mr_ITS2/mrits_diversity.csv") #saving
#df.div <- read.csv("~/Desktop/mrits/mrits_diversity.csv") #reading back in 

diver.sh <- summarySE(data=df.div,measurevar=c("Shannon"),groupvars=c("zone","site"))
diver.si <- summarySE(data=df.div,measurevar=c("InvSimpson"),groupvars=c("zone","site"))

quartz()
gg.sh <- ggplot(diver.sh, aes(x=site, y=Shannon,color=zone,shape=zone))+
  geom_errorbar(aes(ymin=Shannon-se,ymax=Shannon+se),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  geom_point(aes(colour=zone, shape=zone),size=4,position=position_dodge(0.5))+
  xlab("Site")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  theme(text=element_text(family="Times"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))
gg.sh

gg.si <- ggplot(diver.si, aes(x=site, y=InvSimpson,color=zone,shape=zone))+
  geom_errorbar(aes(ymin=InvSimpson-se,ymax=InvSimpson+se),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  geom_point(aes(colour=zone, shape=zone),size=4,position=position_dodge(0.5))+
  xlab("Site")+
  ylab("Inv. Simpson diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  theme(text=element_text(family="Times"))
gg.si

quartz()
ggarrange(gg.sh,gg.si,common.legend=TRUE,legend="right")

mnw <- subset(df.div,site=="MNW")
mse <- subset(df.div,site=="MSE")
tah <- subset(df.div,site=="TNW")

wilcox.test(Shannon~zone,data=mnw)
#not different
wilcox.test(InvSimpson~zone,data=mnw)
#p = 0.01

wilcox.test(Shannon~zone,data=mse)
#p = 0.001
wilcox.test(InvSimpson~zone,data=mse)
#p = 0.01

wilcox.test(Shannon~zone,data=tah)
#p = 0.04
wilcox.test(InvSimpson~zone,data=tah)
#p = 0.019

##checking out how diversity works with variables - nothing interesting
#df.size <- read.csv("mr_size.csv",header=TRUE)
#df.div$coral_id <- row.names(df.div)
# 
# mergeddf <- merge(df.div, df.size, by="coral_id", sort=FALSE)
# plot(Shannon~size,data=mergeddf)
# shapiro.test(mergeddf$Shannon)
# shapiro.test(mergeddf$size)
# kruskal.test(Shannon~size,data=mergeddf)
## no significant relationship between coral size & diversity! interesting

## now by sample host heterozygosity 
# df.het <- read.table("host_het.txt") #read in data
# df.het$het <- df.het$V3/(df.het$V2+df.het$V3) #heterozygosity calculations
# df.het$V1 <- sub("TO","TNWO",df.het$V1) #just renaming some sites
# df.het$V1 <- sub("TI","TNWI",df.het$V1) #just renaming some sites
# df.het$site <- substr(df.het$V1, 0, 4)
# df.het$site <- as.factor(df.het$site)
# df.het$site <- sub("I","-B", df.het$site)
# df.het$site <- sub("O","-F", df.het$site)
# str(df.het)

## just getting sample names on the same page
# df.het$V1 <- sub(".trim.bt2.bam.out.saf.idx.ml","",df.het$V1)
# df.het$coral_id <- substr(df.het$V1, 6, 8)
# df.div$coral_id <- row.names(df.div)
# 
# mergeddf2 <- merge(df.div, df.het, by="coral_id", sort=TRUE)
# plot(Shannon~het,data=mergeddf2)
# kruskal.test(Shannon~het,data=mergeddf2)
## not significant!

#### Bar plot - raw table ####
ps.rel <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
plot_bar(ps.rel, x="Sample",fill="Class")

top <- names(sort(taxa_sums(ps_no87), decreasing=TRUE))[1:30]
ps.top <- transform_sample_counts(ps_no87, function(OTU) OTU/sum(OTU))
ps.top <- prune_taxa(top, ps.top)
plot_bar(ps.top, x="Sample",fill="Class") + facet_wrap(~zone, scales="free_x")

bot <- names(sort(taxa_sums(ps_no87), decreasing=TRUE))[5:118]
ps.bot <- transform_sample_counts(ps_no87, function(OTU) OTU/sum(OTU))
ps.bot <- prune_taxa(bot, ps.bot)
plot_bar(ps.bot, x="Sample",fill="Class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

#bar plot but without the lines between OTUs & no "NAs"
ps.acuta <- subset_samples(ps, species=="P. acuta")

ps_glom <- tax_glom(ps.acuta, "Phylum")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
plot_bar(ps0, x="colony",fill="Phylum")+
  facet_grid(clone ~ ., scales = "free", space = "free",drop=TRUE)

ps_glom <- tax_glom(ps.acuta, "Class")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
plot_bar(ps0, x="colony",fill="Class")+
  facet_grid(clone ~ ., scales = "free", space = "free",drop=TRUE)

# #ps1 <- merge_samples(ps0, "colony")
# ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
# plot_bar(ps2, fill="Class")

#getting raw average relative abundances
ps.rel <- transform_sample_counts(ps_no87, function(x) x / sum(x))
plot_bar(ps.rel,fill="Class")
ps.glom <- tax_glom(ps.rel, "Class")

c3k <- subset_taxa(ps.glom,Class==" C3k")
c3k.otu <- as.data.frame(c3k@otu_table)
mean(c3k.otu$sq1)

#all the seqs - 9 total
cspc <- subset_taxa(ps.rel,Class==" Cspc")
#rel abundance
cspc <- subset_taxa(ps.glom,Class==" Cspc")
cspc.otu <- as.data.frame(cspc@otu_table)
mean(cspc.otu$sq2)

#background ones
back <- subset_taxa(ps.rel,is.na(Class))


all.otu <- ps2@otu_table
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   Kingdom        Phylum     Class     
# sq1  "Symbiodinium" " Clade C" " C3k"  NA
# sq2  "Symbiodinium" " Clade C" " Cspc" NA
# sq7  "Symbiodinium" " Clade A" " A1"   NA
# sq33 "Symbiodinium" " Clade A" " A3"   NA
# sq66 "Symbiodinium" " Clade B" " B1"   NA
# sq81 "Symbiodinium" " Clade C" " C116" NA
all.otu2 <- as.data.frame(all.otu)
mean(all.otu2$sq1)

#### Pcoa - raw data ####

library(vegan)
#install.packages("ggforce")
#library(ggforce)
#not sure if I need this one^
library(ggpubr)
library(cowplot)

# creating a log-transfromed normalized dataset for PCoA:
df.seq <- as.data.frame(seqtab.no87)
all.log=logLin(data=df.seq)

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
all.dist=vegdist(all.log,method="bray")
all.pcoa=pcoa(all.dist)
#32.7%, 23.2%

# plotting:
scores=all.pcoa$vectors[,1:2]
scorez <- as.data.frame(scores)
scorez$Sample <- rownames(scorez)
pcoa.all <- merge(scorez,samdf.no87)

levels(pcoa.all$site) <- c("Moorea NW","Moorea SE","Tahiti NW")
ggplot(pcoa.all,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  xlab('Axis 1 (32.7%)')+
  ylab('Axis 2 (23.2%)')+
  stat_ellipse()+
  facet_wrap(~site)+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))

#now by site - unnecessary actually thanks to facet_wrap
pcoa.mnw <- subset(pcoa.all,site=="MNW")
gg.mnw <- ggplot(pcoa.mnw,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  xlab("Axis 1 ()")
gg.mnw

pcoa.mse <- subset(pcoa.all,site=="MSE")
gg.mse <- ggplot(pcoa.mse,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))
gg.mse

pcoa.tah <- subset(pcoa.all,site=="TNW")
gg.tah <- ggplot(pcoa.tah,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))
gg.tah

quartz()
ggarrange(gg.mnw,gg.mse,gg.tah,nrow=1,common.legend=TRUE,legend="right")

#relative abundance instead of absolute abundance
t.relabun <- scale(t(seqtab.no87), center=F, scale=colSums(t(seqtab.no87)))
relabun <- t(t.relabun)
all.dist=vegdist(relabun,method="manhattan")
all.pcoa=pcoa(all.dist)

scores=all.pcoa$vectors[,1:2]
scorez <- as.data.frame(scores)
scorez$Sample <- rownames(scorez)
pcoa.all <- merge(scorez,samdf.no87,by="Sample")

ggplot(pcoa.all,aes(x=Axis.1,y=Axis.2,color=site,shape=zone))+
  geom_point()+
  stat_ellipse()

#alt method
ps.rel <- transform_sample_counts(ps_no87, function(OTU) OTU/sum(OTU))
iDist <- distance(ps.rel, method="bray")
iMDS  <- ordinate(ps.rel, "MDS", distance=iDist)
plot_ordination(ps.rel, iMDS, color="site")+
  stat_ellipse()

#~##########################################~#
###### Apply LULU to cluster ASVs ############
#~##########################################~#

##Necessary pre-steps after producing ASV table and associated fasta file:
##Produce a match list using BLASTn

##IN TERMINAL#

##First produce a blastdatabase with the OTUs
#module load blast+
#makeblastdb -in mrits2.fasta -parse_seqids -dbtype nucl

##Then blast the OTUs against the database to produce the match list 
#blastn -db mrits2.fasta -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query mrits2.fasta
##HSP = high scoring pair 
##perc_identity = percent of nucleotides in the highly similar pairings that match

##transfer match_list.txt to R working directory
##BACK IN R#

#first, read in ASV table
#install.packages("remotes")
library("remotes")
#install_github("https://github.com/tobiasgf/lulu.git")
library("lulu")
#rarefying
library(vegan)

#must setting up table for lulu, saving it as a csv then reading back in 
#write.csv(seqtab.no87,"seqtab.no87.csv")
setwd("~/moorea_holobiont/mr_ITS2/")
seq.lulu <- read.csv("seqtab.no87.csv",row.names=1)

rarecurve(seq.lulu,step=100,label=FALSE)
total <- rowSums(seq.lulu)
total
subset(total, total <1994)
summary(total)

row.names.remove <- c("117","311","402","414","505","513","530","58","72","76")
lulu.rare <- seq.lulu[!(row.names(seq.lulu) %in% row.names.remove),]
samdf.rare <- samdf.no87[!(row.names(samdf.no87) %in% row.names.remove), ]
#85 samples left

counts.rare <- rrarefy(lulu.rare,sample=1994)
rarecurve(counts.rare,step=100,label=FALSE)
#save
#write.csv(counts.rare,"~/moorea_holobiont/mr_ITS2/seqtabno87.rare.csv")
counts.rare <- read.csv("~/moorea_holobiont/mr_ITS2/seqtabno87.rare.csv",row.names=1)

#back to lulu sequence
ASVs<-data.frame(t(counts.rare)) 
head(ASVs)

#just made this in terminal
matchList<-read.table("match_list.txt")

#had some zeroes to get rid of after rarefying - not the case with the raw table
ASVs <- subset(ASVs,rowSums(ASVs)!=0)
#Now, run the LULU curation
curated_result <- lulu(ASVs, matchList,minimum_match=99,minimum_relative_cooccurence=0.70)
#default: minimum_relative_cooccurence = 0.95 default, changed to 0.70 to see what happens, nothing
#default: minimum_match = 84 default, only 1 OTU different between 97 & 84
summary(curated_result)
#me: leaves 7 OTUs with 84% match - take home story is that they're all very similar
#19 otus with match of 99%
#17 otus with rarefied result

#Pull out the curated OTU list, re-transpose
lulu.out <- data.frame(t(curated_result$curated_table))

#Continue on to your favorite analysis
#write.csv(lulu.out,"~/moorea_holobiont/mr_ITS2/lulu_output.csv")
write.csv(lulu.out,"~/moorea_holobiont/mr_ITS2/lulu_output.rare.csv")
lulu.out <- read.csv("~/moorea_holobiont/mr_ITS2/lulu_output.csv",row.names=1)
lulu.out <- read.csv("~/moorea_holobiont/mr_ITS2/lulu_output.rare.csv",row.names=1)

#removing X from row names to match up with previous data frames, was only necesary from original out, is fine in my .csv file
rownames(lulu.out) <- sub("X","",rownames(lulu.out))

#phyloseq object with lulu results
ps.lulu <- phyloseq(otu_table(lulu.out, taxa_are_rows=FALSE), 
                    sample_data(samdf.no87), 
                    tax_table(taxa2))

#### Bar plot & alpha diversity - clustered results ####
library(cowplot)

ps.top <- transform_sample_counts(ps.lulu, function(OTU) OTU/sum(OTU))
plot_bar(ps.top, x="Sample",fill="Class") + facet_wrap(~zone, scales="free_x")

df.div <- data.frame(estimate_richness(ps.lulu, split=TRUE, measures =c("Shannon","InvSimpson")))
df.div

df.div$Sample <- rownames(df.div)
df.div <- merge(df.div,samdf.no87,by="Sample") #add sample data

diver.sh <- summarySE(data=df.div,measurevar=c("Shannon"),groupvars=c("zone","site"))
diver.si <- summarySE(data=df.div,measurevar=c("InvSimpson"),groupvars=c("zone","site"))

quartz()
gg.sh <- ggplot(diver.sh, aes(x=site, y=Shannon,color=zone,shape=zone))+
  geom_errorbar(aes(ymin=Shannon-se,ymax=Shannon+se),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  geom_point(aes(colour=zone, shape=zone),size=4,position=position_dodge(0.5))+
  xlab("Site")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  theme(text=element_text(family="Times"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))
gg.sh

gg.si <- ggplot(diver.si, aes(x=site, y=InvSimpson,color=zone,shape=zone))+
  geom_errorbar(aes(ymin=InvSimpson-se,ymax=InvSimpson+se),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  geom_point(aes(colour=zone, shape=zone),size=4,position=position_dodge(0.5))+
  xlab("Site")+
  ylab("Inv. Simpson diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  theme(text=element_text(family="Times"))
gg.si

#maybe a boxplot instead
df.div$zone <- gsub("Backreef","Back",df.div$zone)
df.div$zone <- gsub("Forereef","Fore",df.div$zone)
gg.si <- ggplot(df.div, aes(x=zone, y=InvSimpson,color=zone,shape=zone))+
  #geom_errorbar(aes(ymin=InvSimpson-se,ymax=InvSimpson+se),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  #geom_point(aes(colour=zone, shape=zone),size=4,position=position_dodge(0.5))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Reef zone")+
  ylab("Inv. Simpson diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  #guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  theme(text=element_text(family="Times"),legend.position="none")+
  geom_jitter(alpha=0.5)+
  facet_wrap(~site)
#geom_point()
#geom_jitter(aes(color=zone,shape=zone,group=zone))
gg.si

gg.sh <- ggplot(df.div, aes(x=zone, y=Shannon,color=zone,shape=zone))+
  #geom_errorbar(aes(ymin=InvSimpson-se,ymax=InvSimpson+se),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  #geom_point(aes(colour=zone, shape=zone),size=4,position=position_dodge(0.5))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Reef zone")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  theme(text=element_text(family="Times"),legend.position="none")+
  geom_jitter(alpha=0.5)+
  facet_wrap(~site)
gg.sh


df.div$shlog <- log(df.div$Shannon+1)
shapiro.test(df.div$shlog) #worse
library(bestNormalize)
bestNormalize(df.div$Shannon) #says to leave it 

summary(aov(Shannon~zone + Error(site),data=df.div))
wilcox.test(Shannon~zone,data=df.div)
ggplot(df.div,aes(x=site_zone,y=Shannon))+
  geom_boxplot()

mnw <- subset(df.div,site=="MNW")
mse <- subset(df.div,site=="MSE")
tah <- subset(df.div,site=="TNW")

wilcox.test(Shannon~zone,data=mnw)
#0.07459
wilcox.test(InvSimpson~zone,data=mnw)
#p-value = 0.05286

wilcox.test(Shannon~zone,data=mse)
#p = 0.004816
wilcox.test(InvSimpson~zone,data=mse)
#p-value = 0.01058

wilcox.test(Shannon~zone,data=tah)
#p-value = 0.5502
wilcox.test(InvSimpson~zone,data=tah)
#p-value = 0.6848

#install.packages("dunn.test")
library(dunn.test)
dunn.test(df.div$Shannon,df.div$site_zone,method="bh")

#### MCMC.OTU to remove underrepresented ASVs ####
library(MCMC.OTU)

#added a column with a blank name in the beginning, with 1-95 in the column, mcmc.otu likes this
#also removed the X from the beginning of sample names
lulu.mcmc <- read.csv("~/moorea_holobiont/mr_ITS2/lulu_formcmc.csv")
lulu.mcmc <- read.csv("~/moorea_holobiont/mr_ITS2/lulu_rare_formcmc.csv")
#& reading back in things

goods <- purgeOutliers(lulu.mcmc,count.columns=3:19,otu.cut=0.001,zero.cut=0.02)
#otu.cut = 0.1% of reads represented by ASV 
#zero.cut = present in more than 1 sample (2% of samples)
colnames(goods)
#sq 1, 3, 6, 7 retained with min 84% matching in lulu
#sq 1, 2, 3, 6, 7, 12, 18, 24, 32, 33 with min 99% matching in lulu
#sq 33 is only in one sample, but it has 800 reads there so I believe it

#### Bar plot - clustered & trimmed ####
rownames(goods) <- goods$sample
counts <- goods[,3:11]

#mcmc.otu removed 3 undersequenced samples: "513" "530" "76"
remove <- c("513","530","76")
samdf.mcmc <- samdf.no87[!row.names(samdf.no87)%in%remove,]

write.csv(samdf.mcmc,"~/Desktop/mr_samples.csv")
write.csv(counts,"~/Desktop/mr_counts.csv")

ps.mcmc <- phyloseq(otu_table(counts, taxa_are_rows=FALSE), 
                    sample_data(samdf.mcmc), 
                    tax_table(taxa2))

ps.rel <- transform_sample_counts(ps.mcmc, function(x) x / sum(x))
plot_bar(ps.rel,fill="Class")
ps.glom <- tax_glom(ps.rel, "Class")

#c3k.mcmc <- subset_taxa(ps.glom,Class==" C3k")
#c3k.otu <- as.data.frame(c3k.mcmc@otu_table)
#mean(c3k.otu$sq1)

#cs.mcmc <- subset_taxa(ps.glom,Class==" Cspc")
#cs.otu <- as.data.frame(cs.mcmc@otu_table)
#mean(cs.otu$sq2)
#range(cs.otu$sq2)

#bar plot
ps_glom <- tax_glom(ps.mcmc, "Class")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "site_zone")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Class")

ps.mcmc.melt <- psmelt(ps.mcmc)
ps.mcmc.melt$newclass <- paste(ps.mcmc.melt$Class,ps.mcmc.melt$OTU,sep="_")

#boxplot
ggplot(ps.mcmc.melt,aes(x=site_zone,y=Abundance,color=newclass))+
  geom_boxplot()

#individually
sq3 <- subset(ps.mcmc.melt,newclass==" C3k_sq3")
ggplot(sq3,aes(x=site_zone,y=Abundance,color=Class))+
  geom_boxplot()

library(dplyr)
ps.all <- transform_sample_counts(ps.mcmc, function(OTU) OTU/sum(OTU))
pa <- psmelt(ps.all)
tb <- psmelt(ps.all)%>%
  filter(!is.na(Abundance))%>%
  group_by(site_zone,OTU)%>%
  summarize_at("Abundance",mean)

#some more grouping variables
tb <- psmelt(ps.all)%>%
  filter(!is.na(Abundance))%>%
  group_by(site,zone,site_zone,OTU)%>%
  summarize_at("Abundance",mean)

ggplot(tb,aes(x=site_zone,y=Abundance,fill=OTU))+
  geom_bar(stat="identity", colour="black")+
  theme_cowplot()

#renaming

#ps.all@tax_table
# Taxonomy Table:     [9 taxa by 4 taxonomic ranks]:
#   Kingdom        Phylum     Class  
# sq1  "Symbiodinium" " Clade C" " C3k" - 1
# sq2  "Symbiodinium" " Clade C" " Cspc" - 1
# sq3  "Symbiodinium" " Clade C" " C3k" - 2
# sq6  "Symbiodinium" " Clade C" " C3k" - 3
# sq7  "Symbiodinium" " Clade A" " A1"  
# sq12 "Symbiodinium" " Clade C" NA - 1    
# sq18 "Symbiodinium" " Clade C" NA - 2    
# sq24 "Symbiodinium" " Clade C" " C3k" - 4
# sq32 "Symbiodinium" " Clade C" " Cspc" - 2

tb$sym <- gsub("sq12","C. - 1",tb$OTU)
tb$sym <- gsub("sq18","C. - 2",tb$sym)
tb$sym <- gsub("sq1","C3k - 1",tb$sym)
tb$sym <- gsub("sq24","C3k - 4",tb$sym)
tb$sym <- gsub("sq2","Cspc - 1",tb$sym)
tb$sym <- gsub("sq32","Cspc - 2",tb$sym)
tb$sym <- gsub("sq3","C3k - 2",tb$sym)
tb$sym <- gsub("sq6","C3k - 3",tb$sym)
tb$sym <- gsub("sq7","A1",tb$sym)

tb$zone <- gsub("Forereef","Fore",tb$zone)
tb$zone <- gsub("Backreef","Back",tb$zone)

tb$sym <- factor(tb$sym, levels=c("C3k - 1","C3k - 2","C3k - 3","C3k - 4","Cspc - 1", "Cspc - 2","C. - 1","C. - 2","A1"))
quartz()
gg.bp <- ggplot(tb,aes(x=zone,y=Abundance,fill=sym))+
  geom_bar(stat="identity")+
  theme_cowplot()+
  theme(text=element_text(family="Times"))+
  xlab('Reef zone')+
  #  scale_fill_manual(name="Sym.",values=c("seagreen1","seagreen2","seagreen3","seagreen4","blue","darkblue","orange","yellow","purple"))
  scale_fill_manual(name="Algal symbiont",values=c("#1E9C89FF","#25AB82FF","#58C765FF","#7ED34FFF","#365d8dff","#287d8eff","darkorchid4", "darkorchid1","#d4e21aff"))+
  facet_wrap(~site)
gg.bp

#"#1E9C89FF","#25AB82FF","#58C765FF","#7ED34FFF","#365d8dff","#287d8eff","#440154FF", "#48196bff","#d4e21aff"

quartz()
ggarrange(gg.bp,
          ggarrange(gg.sh,gg.si,ncol=2,labels=c("B.","C."),common.legend=T,legend="none",font.label=list(family="Times")),nrow=2,labels="A.",font.label=list(family="Times"))

#new & improved:
ggarrange(gg.bp,gg.si,ncol=1,nrow=2,labels=c("A.","B."))

#### Pcoa - clustered & trimmed ####
library(vegan)
library(cowplot)
#install.packages("ggforce")
#library(ggforce)
#not sure if I need this one^

# creating a log-transfromed normalized dataset for PCoA:
df.seq <- as.data.frame(counts)
all.log=logLin(data=df.seq)

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
all.dist=vegdist(all.log,method="manhattan")
all.pcoa=pcoa(all.dist)

# plotting:
scores=all.pcoa$vectors[,1:2]
scorez <- as.data.frame(scores)
scorez$Sample <- rownames(scorez)
pcoa.all <- merge(scorez,samdf.mcmc)

ggplot(pcoa.all,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  facet_wrap(~site)

#now by site
pcoa.mnw <- subset(pcoa.all,site=="MNW")
gg.mnw <- ggplot(pcoa.mnw,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  xlab("Axis 1 ()")
gg.mnw

pcoa.mse <- subset(pcoa.all,site=="MSE")
gg.mse <- ggplot(pcoa.mse,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  xlab("Axis 1 ()")
gg.mse

pcoa.tah <- subset(pcoa.all,site=="TNW")
gg.tah <- ggplot(pcoa.tah,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  xlab("Axis 1 ()")
gg.tah

ggarrange(gg.mnw,gg.mse,gg.tah,nrow=1,common.legend=TRUE,legend="right")

#another method
plot_ordination(ps.mcmc, ordinate(ps.mcmc, "PCoA"), color = "site") + geom_point(size = 5)

ps.mnw <- subset_samples(ps.mcmc,site=="MNW")
gg.mnw <- plot_ordination(ps.mnw, ordinate(ps.mnw, "PCoA"), color = "zone")+ 
  geom_point()+
  stat_ellipse(level=0.95)
gg.mnw

ps.mse <- subset_samples(ps.mcmc,site=="MSE")
gg.mse <- plot_ordination(ps.mse, ordinate(ps.mse, "PCoA"), color = "zone")+ 
  geom_point()+
  stat_ellipse(level=0.95)

ps.tah <- subset_samples(ps.mcmc,site=="TNW")
gg.tah <- plot_ordination(ps.tah, ordinate(ps.tah, "PCoA"), color = "zone")+ 
  geom_point()+
  stat_ellipse(level=0.95)

ggarrange(gg.mnw,gg.mse,gg.tah,nrow=1,common.legend=TRUE)

#### Deseq differentially abundant ####
library(DESeq2)

#checking if any significant 
ps.mnw = subset_samples(ps.mcmc, site=="MNW")
ds.mnw = phyloseq_to_deseq2(ps.mnw, ~ zone)
dds.mnw <- estimateSizeFactors(ds.mnw,type="poscounts")
stat.mnw = DESeq(dds.mnw, test="Wald", fitType="parametric")

res = results(stat.mnw, cooksCutoff = FALSE)
alpha = 0.05
sigtab.mnw = res[which(res$padj < alpha), ]
sigtab.mnw = cbind(as(sigtab.mnw, "data.frame"), as(tax_table(ps.mnw)[rownames(sigtab.mnw), ], "matrix"))
head(sigtab.mnw)
dim(sigtab.mnw)
#none

ps.mse = subset_samples(ps.mcmc, site=="MSE")
ds.mse = phyloseq_to_deseq2(ps.mse, ~ zone)
dds.mse <- estimateSizeFactors(ds.mse,type="poscounts")
stat.mse = DESeq(dds.mse, test="Wald", fitType="parametric")

res = results(stat.mse, cooksCutoff = FALSE)
alpha = 0.05
sigtab.mse = res[which(res$padj < alpha), ]
sigtab.mse = cbind(as(sigtab.mse, "data.frame"), as(tax_table(ps.mse)[rownames(sigtab.mse), ], "matrix"))
head(sigtab.mse)
dim(sigtab.mse)
#sq 2 & 3

ps.t = subset_samples(ps.mcmc, site=="TNW")
ds.t = phyloseq_to_deseq2(ps.t, ~ zone)
dds.t <- estimateSizeFactors(ds.t,type="poscounts")
stat.t = DESeq(dds.t, test="Wald", fitType="parametric")

res = results(stat.t, cooksCutoff = FALSE)
alpha = 0.05
sigtab.t = res[which(res$padj < alpha), ]
sigtab.t = cbind(as(sigtab.t, "data.frame"), as(tax_table(ps.t)[rownames(sigtab.t), ], "matrix"))
dim(sigtab.t)
#sq 7

#### Heat map ####
#install.packages("pheatmap")
library(pheatmap)
library(dplyr)

#transform to relative abundance rather than absolute
counts$Sample <- rownames(counts)
newnames <- merge(samdf.mcmc,counts, by="Sample")

sq1 <- newnames %>% 
  group_by(site_zone) %>% 
  summarise(sq1 = sum(sq1))
sq2 <- newnames %>% 
  group_by(site_zone) %>% 
  summarise(sq2 = sum(sq2))
sq3 <- newnames %>% 
  group_by(site_zone) %>% 
  summarise(sq3 = sum(sq3))
sq6 <- newnames %>% 
  group_by(site_zone) %>% 
  summarise(sq6 = sum(sq6))
sq7 <- newnames %>% 
  group_by(site_zone) %>% 
  summarise(sq7 = sum(sq7))
sq12 <- newnames %>% 
  group_by(site_zone) %>% 
  summarise(sq12 = sum(sq12))
sq18 <- newnames %>% 
  group_by(site_zone) %>% 
  summarise(sq18 = sum(sq18))
sq24 <- newnames %>% 
  group_by(site_zone) %>% 
  summarise(sq24 = sum(sq24))
sq32 <- newnames %>% 
  group_by(site_zone) %>% 
  summarise(sq32 = sum(sq32))
allsq <- newnames %>% 
  group_by(site_zone) %>% 
  summarise(all = sum(sq1, sq2, sq3, sq6, sq7, sq12, sq18, sq24, sq32))

df1 <- merge(sq1, sq2, by="site_zone")
df2 <- merge(df1,sq3,by="site_zone")
df3 <- merge(df2,sq6,by="site_zone")
df4 <- merge(df3,sq7,by="site_zone")
df5 <- merge(df4,sq12,by="site_zone")
df6 <- merge(df5,sq18,by="site_zone")
df7 <- merge(df6,sq24,by="site_zone")
df.all <- merge(df7,sq32,by="site_zone")
rownames(df.all) <- df.all$site_zone
df.counts <- df.all[,2:10]

df.counts.t <- t(df.counts)
#relative abundance
dat <- scale(df.counts.t, center=F, scale=colSums(df.counts.t))
quartz()
pheatmap(dat,colorRampPalette(c('white','chartreuse3','darkgreen'))(50),cluster_cols=F)

#without all the stuff I just did
counts.t <- t(counts)
dat2 <- scale(counts.t, center=F, scale=colSums(counts.t))

pheatmap(dat2,colorRampPalette(c('white','blue'))(50))

#### rarefy #####
library(vegan)

rarecurve(counts,step=100,label=FALSE) #after clustering & trimming

total <- rowSums(counts)
total
subset(total, total <1994)
summary(total)

row.names.remove <- c("117","311","402","414","505","58","72")
counts.rare <- counts[!(row.names(counts) %in% row.names.remove),]
samdf.rare <- samdf.mcmc[!(row.names(samdf.mcmc) %in% row.names.remove), ]
#85 samples left

seq.rare <- rrarefy(counts.rare,sample=1994)
rarecurve(seq.rare,step=100,label=FALSE)

rarecurve(seqtab.no87,step=100,label=FALSE)#before clustering & trimming

#save
#write.csv(seq.rare,"~/moorea_holobiont/mr_ITS2/seqtab.rare_1994.csv")
write.csv(seq.rare,"~/moorea_holobiont/mr_ITS2/seqtab.rare_1994_all.csv")
#read back in
seq.rare <- read.csv("~/moorea_holobiont/mr_ITS2/seqtab.rare_1994.csv",row.names=1)

#phyloseq object
ps.rare <- phyloseq(otu_table(seq.rare, taxa_are_rows=FALSE), 
                    sample_data(samdf.rare), 
                    tax_table(taxa2))

ps_glom <- tax_glom(ps.rare, "Class")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "site_zone")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Class")

plot_ordination(ps.rare, ordinate(ps.rare,method="PCoA"), color = "zone")+
  geom_point()+
  facet_wrap(~site)+
  theme_cowplot()+
  stat_ellipse()
#ugly

#### PCoA - rarefied, clustered, trimmed ####
df.seq <- as.data.frame(seq.rare)
all.log=logLin(data=df.seq)

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
all.dist=vegdist(all.log,method="manhattan")
all.pcoa=pcoa(all.dist)

write.csv(counts,"~/Desktop/mr_counts_rare.csv")
write.csv(samdf.rare,"~/Desktop/mr_samples_rare.csv")


# plotting:
scores=all.pcoa$vectors[,1:2]
scorez <- as.data.frame(scores)
scorez$Sample <- rownames(scorez)
pcoa.all <- merge(scorez,samdf.rare)

pcoa.all$site <- gsub("MNW","Moorea NW",pcoa.all$site)
pcoa.all$site <- gsub("MSE","Moorea SE",pcoa.all$site)
pcoa.all$site <- gsub("TNW","Tahiti NW",pcoa.all$site)

quartz()
ggplot(pcoa.all,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  facet_wrap(~site)+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  xlab("Axis 1 (54.14%)")+
  ylab("Axis 2 (33.58%)")
#looks the same as before rarefying? 

#now by site
pcoa.mnw <- subset(pcoa.all,site=="MNW")
gg.mnw <- ggplot(pcoa.mnw,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  xlab("Axis 1 ()")
gg.mnw

#3 outliers: 178, 116, 113
row.names.remove <- c("178","116","113")
pcoa.mnw.less <- pcoa.mnw[!(pcoa.mnw$Sample %in% row.names.remove),]
gg.mnw <- ggplot(pcoa.mnw.less,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  xlab("Axis 1 ()")
gg.mnw

#### Stats ####
library(vegan)
#help on adonis here:
#https://thebiobucket.blogspot.com/2011/04/assumptions-for-permanova-with-adonis.html#more

#all
#dist.seqtab <- vegdist(seqtab.no87)
#anova(betadisper(dist.seqtab,samdf.no87$zone))
adonis(seqtab.no87 ~ zone, strata=samdf.no87$site, data=samdf.no87, permutations=999)
#0.01 **

#clustered but not trimmed
adonis(lulu.out ~ zone, strata=samdf.no87$site, data=samdf.no87, permutations=999)
#p 0.009

#clustered & trimmed 
adonis(counts ~ zone, strata=samdf.mcmc$site, data=samdf.mcmc, permutations=999)
#0.05 .

#relative abundance, does this by columns so must transform
t.relabun <- scale(t(counts), center=F, scale=colSums(t(counts)))
#un-transform
relabun <- t(t.relabun)

adonis(relabun ~ zone, strata=samdf.mcmc$site, data=samdf.mcmc, permutations=999)
#0.02 *

#rarefied, clustered, trimmed
dist.rare <- vegdist(seq.rare)
anova(betadisper(dist.rare,samdf.rare$site))

adonis(counts ~ zone, strata=samdf.rare$site, data=samdf.rare, permutations=999)
adonis(seq.rare ~ zone, strata=samdf.rare$site, data=samdf.rare, permutations=999)
#0.07 . - 0.1

#Moorea NW
mnw.rare <- subset(samdf.rare,site=="MNW")
seq.rare2 <- as.data.frame(counts)
seq.rare2$sample <- rownames(seq.rare2)
mnw.seq <- subset(seq.rare2, sample %in% mnw.rare$Sample)
mnw.seq2 <- mnw.seq[,1:9]

adonis(mnw.seq2 ~ zone, data=mnw.rare, permutations=999)
#0.044 *

#Moorea SE
mse.rare <- subset(samdf.rare,site=="MSE")
mse.seq <- subset(seq.rare2, sample %in% mse.rare$Sample)
mse.seq2 <- mse.seq[,1:9]

adonis(mse.seq2 ~ zone, strata=mse.rare$site, data=mse.rare, permutations=999)
#0.018 *

#Tahiti
tah.rare <- subset(samdf.rare,site=="TNW")
tah.seq <- subset(seq.rare2, sample %in% tah.rare$Sample)
tah.seq2 <- tah.seq[,1:9]

dist.tah <- vegdist(tah.seq2)
anova(betadisper(dist.tah,tah.rare$zone))
adonis(tah.seq2 ~ zone, strata=tah.rare$site, data=tah.rare, permutations=999)
#0.055 rarefied
#0.043 *

#log-normalized
df.seq <- as.data.frame(seqtab.no87)
all.log=logLin(data=df.seq)

all.dist=vegdist(all.log,method="manhattan")
all.pcoa=pcoa(all.dist)

# plotting:
scores=all.pcoa$vectors[,1:2]
scorez <- as.data.frame(scores)
scorez$Sample <- rownames(scorez)
pcoa.all <- merge(scorez,samdf.no87)

ggplot(pcoa.all,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  facet_wrap(~site)
adonis(df.seq ~ zone, strata=samdf.no87$site, data=samdf.no87, permutations=999)
#0.012 *

#rarefied, but not clustered & trimmed
rarecurve(counts.rare,label=F) #yep def rarefied
adonis(counts.rare ~ zone, strata=samdf.rare$site, data=samdf.rare, permutations=999)
#0.01 **

#### bray-curtis ####
iDist <- distance(ps.rare, method="bray")
iMDS  <- ordinate(ps.rare, "MDS", distance=iDist)
plot_ordination(ps.rare, iMDS, color="zone", shape="site")
#ugly

#by site
ps.mnw <- subset_samples(ps.rare,site=="MNW")
iDist <- distance(ps.mnw, method="bray")
iMDS  <- ordinate(ps.mnw, "MDS", distance=iDist)
plot_ordination(ps.mnw, iMDS, color="zone")+
  stat_ellipse()

#relative abundance
t.relabun <- scale(t(counts), center=F, scale=colSums(t(counts)))
relabun <- t(t.relabun)
all.dist=vegdist(relabun,method="bray")
all.pcoa=pcoa(all.dist)

# plotting:
scores=all.pcoa$vectors[,1:2]
scorez <- as.data.frame(scores)
scorez$Sample <- rownames(scorez)
pcoa.all <- merge(scorez,samdf.no87)

ggplot(pcoa.all,aes(x=Axis.1,y=Axis.2,color=zone,shape=site))+
  geom_point()+
  stat_ellipse()

#now by site
pcoa.mnw <- subset(pcoa.all,site=="MNW")
gg.mnw <- ggplot(pcoa.mnw,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  xlab("Axis 1 ()")
gg.mnw

#### CCA ####
library(adegenet) # for transp()

pp0=capscale(seq.rare~1)
pp=capscale(seq.rare~zone%in%site,data=samdf.rare)
anova(pp, alpha=0.05)

axes2plot=c(1,2)  
quartz()
cmd=pp #change to pp for CAP, pp0 for MDS
plot(cmd,choices=axes2plot) # choices - axes to display
points(cmd,choices=axes2plot)
ordihull(cmd,choices= axes2plot,groups=samdf.rare$site_zone,draw="polygon",label=F)
#ordispider(cmd,choices= axes2plot,groups=samdf.rare$site,col="grey80")
#ordiellipse(cmd,choices= axes2plot,groups=samdf.rare$zone,draw="polygon",label=T)

#### indicator species ####
#install.packages("indicspecies")
library(indicspecies)
library(vegan)
#first testing out k means clustering
mcmc.km <- kmeans(counts,centers=2)
groupskm = mcmc.km$cluster
groupskm

all.log=logLin(data=counts)
all.dist=vegdist(all.log,method="bray")
all.pcoa=pcoa(all.dist)

km.log <- kmeans(all.log,centers=3)
groupskm = km.log$cluster

scores=all.pcoa$vectors[,1:2]
scorez <- as.data.frame(scores)
scorez$Sample <- rownames(scorez)
pcoa.all <- merge(scorez,samdf.mcmc)
grps <- as.data.frame(groupskm)
grps$Sample <- rownames(grps)
pcoa2 <- merge(grps,pcoa.all,by="Sample")

pcoa2$groupskm <- as.factor(pcoa2$groupskm)
ggplot(pcoa2,aes(x=Axis.1,y=Axis.2,color=site,shape=groupskm))+
  geom_point()+
  stat_ellipse()
#interesting - not sure what to do here

#anyway back to indicspecies
groups <- samdf.mcmc$site
groups <- samdf.mcmc$zone
indval <- multipatt(counts, groups, control = how(nperm=999))
summary(indval,alpha=1)
#doesn't really work with the clustered ASV

#unclustered
groups <- samdf.no87$zone
indval <- multipatt(seqtab.no87, groups, control = how(nperm=999))
summary(indval)
# Group Backreef  #sps.  1 
# stat p.value  
# sq22 0.478   0.017 *
#   
#   Group Forereef  #sps.  7 
# stat p.value    
# sq15 0.711   0.001 ***
#   sq9  0.526   0.007 ** 
#   sq35 0.410   0.027 *  
#   sq37 0.388   0.019 *  
#   sq19 0.380   0.019 *  
#   sq17 0.377   0.024 *  
#   sq45 0.354   0.036 *  

#now rarefied
counts.rare <- read.csv(file="~/moorea_holobiont/mr_ITS2/seqtabno87.rare.csv",row.names=1)
groups <- samdf.rare$zone
groups <- samdf.rare$site_zone
indval <- multipatt(counts.rare, groups, control = how(nperm=999))
summary(indval)
# Group Backreef  #sps.  4 
# stat p.value   
# sq10 0.620   0.004 **
#   sq29 0.594   0.005 **
#   sq22 0.556   0.003 **
#   sq66 0.420   0.038 * 
#   
#   Group Forereef  #sps.  4 
# stat p.value    
# sq15 0.729   0.001 ***
#   sq9  0.544   0.013 *  
#   sq35 0.436   0.019 *  
#   sq37 0.404   0.017 *  
#^ the same ones as unrarefied, but lost some & gained some