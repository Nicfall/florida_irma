#Nicola's FL Irma analysis - 16S
#Almost entirely based on DADA2 Pipeline 1.8 Walkthrough:
#https://benjjneb.github.io/dada2/tutorial.html
#with edits by Carly D. Kenkel and modifications for my data by Nicola Kriefall
#3/19/20

#this looks promising:
#https://github.com/RichieJu520/Co-occurrence_Network_Analysis

#~########################~#
##### PRE-PROCESSING #######
#~########################~#

#fastq files should have R1 & R2 designations for PE reads
#Also - some pre-trimming. Retain only PE reads that match amplicon primer. Remove reads containing Illumina sequencing adapters

#in Terminal home directory:
#following instructions of installing BBtools from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/
#1. download BBMap package, sftp to installation directory
#2. untar: 
#tar -xvzf BBMap_(version).tar.gz
#3. test package:
#cd bbmap
#~/bin/bbmap/stats.sh in=~/bin/bbmap/resources/phix174_ill.ref.fa.gz

# my adaptors for 16S, which I saved as "adaptors.fasta"
# >forward
# AATGATACGGCGACCAC
# >forwardrc
# GTGGTCGCCGTATCATT
# >reverse
# CAAGCAGAAGACGGCATAC
# >reverserc
# GTATGCCGTCTTCTGCTTG

#primers for 16S: 
# >forward
# GTGYCAGCMGCCGCGGTA
# >reverse
# GGACTACHVGGGTWTCTAAT

##Still in terminal - making a sample list based on the first phrase before the underscore in the .fastq name
#ls *R1_001.fastq | cut -d '_' -f 1 > samples.list
#ls *R1_001.fastq | cut -d '_' -f 2 > samples.list

##cuts off the extra words in the .fastq files
#for file in $(cat samples.list); do  mv ${file}_*R1*.fastq ${file}_R1.fastq; mv ${file}_*R2*.fastq ${file}_R2.fastq; done 
#for file in $(cat samples.list); do  mv *${file}_*R1*.fastq ${file}_R1.fastq; mv *${file}_*R2*.fastq ${file}_R2.fastq; done 

##gets rid of reads that still have the adaptor sequence, shouldn't be there, I didn't have any
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1.fastq in2=${file}_R2.fastq ref=adaptors.fasta out1=${file}_R1_NoIll.fastq out2=${file}_R2_NoIll.fastq; done &>bbduk_NoIll.log

##getting rid of first 4 bases (degenerate primers created them)
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1_NoIll.fastq in2=${file}_R2_NoIll.fastq ftl=4 out1=${file}_R1_NoIll_No4N.fastq out2=${file}_R2_NoIll_No4N.fastq; done &>bbduk_No4N.log

##only keeping reads that start with the 16S primer
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1_NoIll_No4N.fastq in2=${file}_R2_NoIll_No4N.fastq restrictleft=20 k=10 literal=GTGYCAGCMGCCGCGGTA,GGACTACHVGGGTWTCTAAT copyundefined=t outm1=${file}_R1_NoIll_No4N_16S.fastq outu1=${file}_R1_check.fastq outm2=${file}_R2_NoIll_No4N_16S.fastq outu2=${file}_R2_check.fastq; done &>bbduk_16S.log
##higher k = more reads removed, but can't surpass k=20 or 21

##using cutadapt to remove primer
# for file in $(cat samples.list)
# do
# cutadapt -g GTGYCAGCMGCCGCGGTA -a ATTAGAWACCCVHGTAGTCC -G GGACTACHVGGGTWTCTAAT -A TACCGCGGCKGCTGRCAC -n 2 --discard-untrimmed -o ${file}_R1.fastq -p ${file}_R2.fastq ${file}_R1_NoIll_No4N_16S.fastq ${file}_R2_NoIll_No4N_16S.fastq
# done &> clip.log
##-g regular 5' forward primer 
##-G regular 5' reverse primer
##-o forward out
##-p reverse out
##-max-n 0 means 0 Ns allowed
##this overwrote my original renamed files 

# did sftp of *_R1.fastq & *_R2.fastq files to the folder to be used in dada2

#~########################~#
##### DADA2 BEGINS #########
#~########################~#

#installing/loading packages:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.8")
library(dada2); packageVersion("dada2")
#I have version 1.10.0 - tutorial says 1.8 but I think that's OK, can't find a version 1.10 walkthrough
library(ShortRead)
#packageVersion("ShortRead")
library(Biostrings)
#packageVersion("Biostrings")
path <- "/Volumes/TOSHIBA EXT/flirma_apr2020/files" # CHANGE ME to the directory containing the fastq files after unzipping.

fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)
sample.names

#### check for primers ####
FWD <- "GTGYCAGCMGCCGCGGTA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACHVGGGTWTCTAAT"  ## CHANGE ME...

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
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
#no primers - amazing

#### Visualizing raw data ####

#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs.filtN[c(1,2,3,4)])
plotQualityProfile(fnFs.filtN[c(362,363,364,365)])
#looks mostly good up to 200

#Then look at quality profile of R2 reads
plotQualityProfile(fnRs.filtN[c(1,2,3,4)])
plotQualityProfile(fnRs.filtN[c(362,363,364,365)])
#180

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#changing a bit from default settings - maxEE=1 (1 max expected error, more conservative), truncating length at 200 bp for both forward & reverse [leaves ~50bp overlap], added "trimleft" to cut off primers [18 for forward, 20 for reverse]
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(200,180), #leaves ~50bp overlap
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     #trimLeft=c(18,20), #N nucleotides to remove from the start of each read
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
                     compress=TRUE, multithread=TRUE,verbose=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)

#~############################~#
##### Learn Error Rates ########
#~############################~#

#setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
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

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

plot(table(nchar(getSequences(seqtab)))) #real variants appear to be right in that 244-264 window

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(200,300)] #again, being fairly conservative wrt length

#~############################~#
##### Remove chimeras ##########
#~############################~#

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#Identified 1375 bimeras out of 80386 input sequences.

sum(seqtab.nochim)/sum(seqtab)
#0.9963629
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. 

#save
saveRDS(seqtab.nochim, file="flirma_seqtab.nochim.rds")
write.csv(seqtab.nochim, file="flirma_seqtab.nochim.csv")

seqtab.nochim <- readRDS("flirma_seqtab.nochim.rds")

#~############################~#
##### Track Read Stats #########
#~############################~#

#### need to do this :(, I done goofed #### 

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="mr16s_readstats.csv",row.names=TRUE,quote=FALSE)

#~############################~#
##### Assign Taxonomy ##########
#~############################~#

library(dada2)

# #Using package DECIPHER as an alternatie to 'assignTaxonomy'
# 
# #BiocManager::install("DECIPHER")
# library(DECIPHER); packageVersion("DECIPHER")
# #citation("DECIPHER")
# 
# #http://DECIPHER.codes/Downloads.html. Download the SILVA SSU r132 (modified) file to follow along.
# 
# dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
# load("~/Downloads/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
# ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE, threshold=50) # use all processors
# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxid <- t(sapply(ids, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa
# }))
# colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

#doing other taxonomy method:
#Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/silva_nr_v132_train_set.fa.gz",tryRC=TRUE)
unname(head(taxa))
taxa.plus <- addSpecies(taxa, "~/Downloads/silva_species_assignment_v132.fa.gz",tryRC=TRUE,verbose=TRUE)
#507 out of 79011 were assigned to the species level.
#Of which 436 had genera consistent with the input table.

setwd("~/florida_irma")                                                                    
saveRDS(taxa.plus, file="mr16s_taxaplus.rds")
saveRDS(taxa, file="mr16s_taxa.rds")
write.csv(taxa.plus, file="mr16s_taxaplus.csv")
write.csv(taxa, file="mr16s_taxa.csv")

saveRDS(seqtab.nochim, file="mr16s_seqtab.nochim.rds")
write.csv(seqtab.nochim, file="mr16s_seqtab.nochim.csv")
write.csv(seqtab.nochim, file="mr16s_seqtab.nochim_renamed.csv")

#### Read in previously saved datafiles ####
setwd("~/florida_irma/fl16S")
seqtab.nochim <- readRDS("flirma_seqtab.nochim.rds")
taxa <- readRDS("fl16s_taxa.rds")
taxa.plus <- readRDS("fl16s_taxaplus.rds")

#~############################~#
##### handoff 2 phyloseq #######
#~############################~#

#BiocManager::install("phyloseq")
library('phyloseq')
library('ggplot2')
library('Rmisc')
library('cowplot')

#import dataframe holding sample information
samdf<-read.csv("fl16s_sampledata.csv")
head(samdf)
rownames(samdf) <- samdf$sample_16s

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa.plus))

ps

#first look at data
ps_glom <- tax_glom(ps, "Family")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "transect_zone")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, x="transect_zone", fill="Family")+
  theme(legend.position="none")

#phyloseq object with shorter names - doing this one instead of one above
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
#making output fasta file for lulu step & maybe other things, before giving new ids to sequences
#path='~/moorea_holobiont/mr_16s/mr16s.fasta'
#uniquesToFasta(seqtab.nochim, path, ids = ids, mode = "w", width = 20000)

colnames(seqtab.nochim)<-ids
taxa2 <- cbind(taxa.plus, rownames(taxa.plus)) #retaining raw sequence info before renaming
rownames(taxa2)<-ids

#phyloseq object with new taxa ids
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa2))

ps #79011 taxa

#### remove mitochondria, chloroplasts, non-bacteria #### 
ps.mito <- subset_taxa(ps, (Family=="Mitochondria"))
ps.mito #944 taxa to remove
ps.chlor <- subset_taxa(ps, (Order=="Chloroplast"))
ps.chlor #1142 taxa to remove
ps.notbact <- subset_taxa(ps, (Kingdom!="Bacteria") | is.na(Kingdom))
ps.notbact #3912 taxa to remove

ps.nomito <- subset_taxa(ps, (Family!="Mitochondria") | is.na(Family))
ps.nomito #78067 taxa
ps.nochlor <- subset_taxa(ps.nomito, (Order!="Chloroplast") | is.na(Order))
ps.nochlor #76925 taxa
ps.clean <- subset_taxa(ps.nochlor, (Kingdom=="Bacteria"))
ps.clean #73013 taxa

#just archaea
ps.arch <- subset_taxa(ps.nomito, (Kingdom=="Archaea"))
ps.arch #2608 taxa

seqtab.clean <- data.frame(otu_table(ps.clean))
write.csv(seqtab.clean,file="fl16s_seqtab.cleaned.csv")

saveRDS(seqtab.clean, file="fl16s_seqtab.cleaned.rds")
saveRDS(taxa2,file="fl16s_taxa2.rds")

#### Read in previously saved datafiles ####
setwd("~/florida_irma/fl16S")
seqtab.clean <- readRDS("fl16s_seqtab.cleaned.rds")
taxa2 <- readRDS("fl16s_taxa2.rds")
samdf <- read.csv("fl16s_sampledata.csv")
rownames(samdf) <- samdf$sample_16s

library("phyloseq")

#re-make phyloseq object if necessary
ps.clean <- phyloseq(otu_table(seqtab.clean, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa2))
ps.clean #should be 73013 taxa

#remove some pdam & core samples
ps.clean2 <- subset_samples(ps.clean, species!="pdam" & species!="cores")

#just sids or rads
ps.sid <- subset_samples(ps.clean2, species=="ssid")
ps.rad <- subset_samples(ps.clean2, species=="srad")

#save sids & rads tables
seqtab.sid <- data.frame(otu_table(ps.sid))
seqtab.rad <- data.frame(otu_table(ps.rad))

#many are 0 total counts
sid.zero = seqtab.sid[,colSums(seqtab.sid) == 0]
ncol(sid.zero) #44389
removecols.sid <- c(colnames(sid.zero))
seqtab.sid.no0 <- seqtab.sid[,!(colnames(seqtab.sid) %in% removecols.sid)]

#many are 0 total counts
rad.zero = seqtab.rad[,colSums(seqtab.rad) == 0]
ncol(rad.zero) #14008
removecols.rad <- c(colnames(rad.zero))
seqtab.rad.no0 <- seqtab.rad[,!(colnames(seqtab.rad) %in% removecols.rad)]

saveRDS(seqtab.sid.no0, file="fl16s_seqtab.sid.rds")
saveRDS(seqtab.rad.no0, file="fl16s_seqtab.rad.rds")

#### read in files here too ####
seqtab.sid.no0 <- readRDS(file="fl16s_seqtab.sid.rds")
seqtab.rad.no0 <- readRDS(file="fl16s_seqtab.rad.rds")

samdf.sid <- data.frame(sample_data(ps.sid))
samdf.rad <- data.frame(sample_data(ps.rad))

write.csv(samdf.sid, file="fl16s_samdf.sid.csv")
write.csv(samdf.rad, file="fl16s_samdf.rad.csv")

samdf.sid <- read.csv("fl16s_samdf.sid.csv",row.names=1)
samdf.rad <- read.csv("fl16s_samdf.rad.csv",row.names=1)

# #### put somewhere ####
# #just1516 <- subset_samples(ps.clean2, year==15 | year==16)
# #just18 <- subset_samples(ps.clean2, year==18)
# 
# ps_glom <- tax_glom(ps.sid, "Family")
# ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
# ps1 <- merge_samples(ps0, "transect_zone_year")
# ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
# plot_bar(ps2, x="transect_zone_year", fill="Family")+
#   theme(legend.position="none")+
#   facet_wrap(~transect)

#### rarefy ####
library(vegan)

rarecurve(seqtab.sid.no0,step=100,label=FALSE) #after removing contaminants

total <- rowSums(seqtab.sid.no0)
toofew <- subset(total, total <9100)
#20 samples to remove, ouch

row.names.remove <- c(names(toofew))
seqtab.sid.less <- seqtab.sid.no0[!(row.names(seqtab.sid.no0) %in% row.names.remove),]
samdf.sid.rare <- samdf.sid[!(row.names(samdf.sid) %in% row.names.remove), ]
#173 samples left (168 first round through, rarefying to 11627)

seqtab.sid.rare <- rrarefy(seqtab.sid.less,sample=9100)
#rarecurve(seqtab.rare,step=100,label=FALSE)
saveRDS(seqtab.sid.rare, file="seqtab.sid.rare9100.rds")
seqtab.sid.rare <- readRDS(file="seqtab.sid.rare9100.rds")

sid.rare.zero = seqtab.sid.rare[,colSums(seqtab.sid.rare) == 0]
ncol(sid.rare.zero) #3661, 4098 for 9100
removecols.sid <- c(colnames(sid.rare.zero))
seqtab.sid.rare.no0 <- seqtab.sid.rare[,!(colnames(seqtab.sid.rare) %in% removecols.sid)]

samdf.sid.rare$year <- as.factor(samdf.sid.rare$year)
samdf.sid.rare$year <- gsub("16","15",samdf.sid.rare$year)
samdf.sid.rare$year <- as.factor(samdf.sid.rare$year)
samdf.sid.rare.nobis <- subset(samdf.sid.rare,transect!="Biscayne")
ggplot(data.frame(samdf.sid.rare.nobis), aes(x=year)) +
  geom_bar()

#phyloseq object but rarefied
ps.sid.rare <- phyloseq(otu_table(seqtab.sid.rare.no0, taxa_are_rows=FALSE), 
                    sample_data(samdf.sid.rare.nobis), 
                    tax_table(taxa2))
ps.sid.rare #24526 taxa after removing 0s, 153 samples after rarefying

#### alpha diversity ####
#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)
#Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).
#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index (1 − λ).
plot_richness(ps.sid, x="transect_zone_year", measures=c("Shannon", "Simpson","Observed"), color="zone") + theme_bw()

df <- data.frame(estimate_richness(ps.sid, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))
df
df <- data.frame(estimate_richness(ps.sid.rare, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))

df$sample_16s <- rownames(df)
df$sample_16s <- gsub("X","",df$sample_16s)
df.div <- merge(df,samdf.sid.rare,by="sample_16s") #add sample data

write.csv(df,file="fl16s_diversity_sids.csv") #saving
df.div <- read.csv("fl16s_diversity_sids.csv") #reading back in 

#concatenating 15 with 16
df.div$year <- gsub("16","15",df.div$year)
df.div$year <- as.factor(df.div$year)
df.div.nobis <- subset(df.div,transect!="Biscayne")

quartz()
#zone
gg.sh <- ggplot(df.div.nobis, aes(x=zone, y=Shannon,color=zone,shape=zone))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Reef zone")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  facet_wrap(~transect*year)
gg.sh
  # scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  # scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  # #guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  # theme(text=element_text(family="Times"),legend.position="none")+
  # geom_jitter(alpha=0.5)+

gg.sh <- ggplot(df.div.nobis, aes(x=year,y=Shannon,color=year))+
  geom_boxplot()+
  xlab("Year")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  #facet_wrap(~transect)+
  scale_color_manual(values=c("#781C6D","#F8850F","#AE305C"),name="Time point",labels=c("Baseline","Irma","Post-Irma"))+
  scale_x_discrete(labels=c("2015","2017","2018"))+
  theme(axis.text.x=element_text(angle=45,hjust=1))
quartz()
gg.sh

gg.si <- ggplot(df.div.nobis, aes(x=transect, y=InvSimpson,color=year,shape=year))+
  geom_boxplot()+
  xlab("Reef zone")+
  ylab("Inv. Simpson diversity")+
  theme_cowplot()+
  ylim(0,100)#+
  #scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  #scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  #guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  #theme(text=element_text(family="Times"),legend.position="none")+
  #geom_jitter(alpha=0.5)#+
  #facet_wrap(~site)
gg.si

gg.obs <- ggplot(df.div, aes(x=zone, y=Observed,color=zone,shape=zone))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Reef zone")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  #guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  theme(text=element_text(family="Times"),legend.position="none")+
  geom_jitter(alpha=0.5)+
  facet_wrap(~site)
gg.obs

#stats
shapiro.test(df.div.nobis$Observed)
df.div.nobis$obs.log <- log(df.div.nobis$Observed)
shapiro.test(df.div.nobis$obs.log)

a.all <- aov(obs.log~year*transect,data=df.div.nobis)
summary(a.all)
TukeyHSD(a.all) #NOTHING

#stats
shapiro.test(df.div.nobis$Shannon)
df.div.nobis$sh.log <- log(df.div.nobis$Shannon)
shapiro.test(df.div.nobis$sh.log)

a.all <- aov(Shannon~year*transect,data=df.div.nobis)
summary(a.all)
TukeyHSD(a.all) #nothing

a.all <- aov(InvSimpson~year*transect,data=df.div.nobis)
summary(a.all)
TukeyHSD(a.all) #nothing

div.low <- subset(df.div.nobis,transect=="Lower")
div.mid <- subset(df.div.nobis,transect=="Middle")
div.upp <- subset(df.div.nobis,transect=="Upper")

#lower transect - OTU number
shapiro.test(div.low$Observed) #nope
div.low$obs.log <- log(div.low$Observed)
shapiro.test(div.low$obs.log) #yessss

a.div.low <- aov(obs.log~year,data=div.low)
summary(a.div.low)
TukeyHSD(a.div.low) #nothing

#lower transect - Shannon
shapiro.test(div.low$Shannon) #nope
div.low$sh.log <- logb(div.low$Shannon, base=10)
shapiro.test(div.low$sh.log) 

a.div.low <- aov(sh.log~year*zone,data=div.low)
summary(a.div.low)
TukeyHSD(a.div.low) #nothing

#lower transect - simpson
shapiro.test(div.low$InvSimpson) #nope
div.low$si.log <- log(div.low$InvSimpson)
shapiro.test(div.low$si.log) #no

a.div.low <- aov(InvSimpson~year*zone,data=div.low)
summary(a.div.low)
TukeyHSD(a.div.low) #nothing

div.low15 <- subset(div.low,year=="15")
div.low17 <- subset(div.low,year=="17")
div.low18 <- subset(div.low,year=="18")

wilcox.test(InvSimpson~zone,data=div.low15) #yes
wilcox.test(InvSimpson~zone,data=div.low17) #nope
wilcox.test(InvSimpson~zone,data=div.low18) #nope

#mid transect - OTU number
shapiro.test(div.mid$Observed) #nope
div.mid$obs.log <- log(div.mid$Observed)
shapiro.test(div.mid$obs.log) #yessss

a.div.mid <- aov(obs.log~year*zone,data=div.mid)
summary(a.div.mid)
TukeyHSD(a.div.mid) #nothing

#mid transect - Shannon
shapiro.test(div.mid$Shannon) #good

a.div.mid <- aov(Shannon~year*zone,data=div.mid)
summary(a.div.mid)
TukeyHSD(a.div.mid) #nothing

#mid transect - simpson
shapiro.test(div.mid$InvSimpson) #nope
div.mid$si.log <- log(div.mid$InvSimpson)
shapiro.test(div.mid$si.log) #yes

a.div.mid <- aov(InvSimpson~year*zone,data=div.mid)
summary(a.div.mid)
TukeyHSD(a.div.mid) #nothing

#upp transect - OTU number
shapiro.test(div.upp$Observed) #nope
div.upp$obs.log <- log(div.upp$Observed)
shapiro.test(div.upp$obs.log) #yessss

a.div.upp <- aov(obs.log~year*zone,data=div.upp)
summary(a.div.upp)
TukeyHSD(a.div.upp) #nothing

#upp transect - Shannon
shapiro.test(div.upp$Shannon) #good

a.div.upp <- aov(Shannon~year*zone,data=div.upp)
summary(a.div.upp)
TukeyHSD(a.div.upp) #nothing

#upp transect - simpson
shapiro.test(div.upp$InvSimpson) #nope
div.upp$si.log <- log(div.upp$InvSimpson)
shapiro.test(div.upp$si.log) #yes

a.div.upp <- aov(InvSimpson~year*zone,data=div.upp)
summary(a.div.upp)
TukeyHSD(a.div.upp) #nothing

#### trim underrepresented ASVs ####
library(MCMC.OTU)

#formatting the table for mcmc.otu - requires one first column that's 1 through whatever
#& has "X" as column name
nums.sid <- 1:nrow(seqtab.sid.rare.no0) #did it another time with untrimmed, unrarefied seqtab.sid.no0
samples.sid <- rownames(seqtab.sid.rare.no0)

sid.int <- cbind(sample = 0, seqtab.sid.rare.no0)
sid.formcmc <- cbind(X = 0, sid.int)

sid.formcmc$X <- nums.sid
sid.formcmc$sample <- samples.sid

#now for radians
nums.rad <- 1:nrow(seqtab.rad.no0)
samples.rad <- rownames(seqtab.rad.no0)

rad.int <- cbind(sample = 0, seqtab.rad.no0)
rad.formcmc <- cbind(X = 0, rad.int)

rad.formcmc$X <- nums.rad
rad.formcmc$sample <- samples.rad

#regular
seq.trim.sid <- purgeOutliers(sid.formcmc,count.columns=3:28626,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)
#3 bad samples - "C10" "E8"  "H6"
#673 ASVs passing cutoff for reads
#584 ASVs show up in 2% of samples
#remove bad samples from sample data frame
row.names.remove <- c("C10","E8","H6")
samdf.trim.sid <- samdf.sid[!(row.names(samdf.sid) %in% row.names.remove), ]

#### radians - leaving for now ####
seq.trim.rad <- purgeOutliers(rad.formcmc,count.columns=3:28626,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.015)

#remove sample info
seq.trim.sid <- seq.trim.sid[,3:584] #regular (not rarefied)

write.csv(seq.trim.sid,file="seq.trim.sid.1e4.csv")
seq.trim.sid <- read.csv("seq.trim.sid.1e4.csv",row.names=1)

#remake phyloseq object
ps.sid.trim <- phyloseq(otu_table(seq.trim.sid, taxa_are_rows=FALSE), 
                        sample_data(samdf.sid), 
                        tax_table(taxa2))
ps.sid.trim #611 taxa left, 582 after more strict 2% of samples

#### trying alpha diversity again post-trim ####
df <- data.frame(estimate_richness(ps.sid.trim, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))

df$sample_16s <- rownames(df)
df$sample_16s <- gsub("X","",df$sample_16s)
df.div <- merge(df,samdf.sid,by="sample_16s") #add sample data

write.csv(df,file="fl16s_diversity_sids.csv") #saving
df.div <- read.csv("fl16s_diversity_sids.csv") #reading back in 

#concatenating 15 with 16
df.div$year <- gsub("16","15",df.div$year)
df.div.nobis <- subset(df.div,transect!="Biscayne")

quartz()
#zone
gg.sh <- ggplot(df.div.nobis, aes(x=zone, y=Shannon,color=zone,shape=zone))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Reef zone")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  facet_wrap(~transect*year)
gg.sh

#year
gg.sh <- ggplot(df.div.nobis, aes(x=year, y=Observed))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Reef zone")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  facet_wrap(~transect)
gg.sh

#### rarefy trimmed ####
library(vegan)

rarecurve(seq.trim.sid,step=100,label=FALSE) #after removing contaminants

total <- rowSums(seq.trim.sid)
toofew <- subset(total, total <2175)
toofew
#19 samples to remove ~10%)

row.names.remove <- c(names(toofew))
seq.trim.sid.less <- seq.trim.sid[!(row.names(seq.trim.sid) %in% row.names.remove),]
samdf.trim.sid.rare <- samdf.sid[!(row.names(samdf.sid) %in% row.names.remove), ]
#174 samples left being strict

seq.trim.sid.rare <- rrarefy(seq.trim.sid.less,sample=2175)
rarecurve(seq.trim.sid.rare,step=100,label=FALSE)

#phyloseq object but rarefied
samdf.trim.sid.rare$year <- as.factor(samdf.trim.sid.rare$year)
samdf.trim.sid.rare$year <- gsub("16","15",samdf.trim.sid.rare$year)
samdf.trim.sid.rare.nobis <- subset(samdf.trim.sid.rare,transect!="Biscayne")
ps.trim.sid.rare <- phyloseq(otu_table(seq.trim.sid.rare, taxa_are_rows=FALSE), 
                        sample_data(samdf.trim.sid.rare.nobis), 
                        tax_table(taxa2))
ps.trim.sid.rare #582 taxa, 171 samples

#### pcoa plots - batch effect :( ####
library(stats)
library(MCMC.OTU)

#transform to relative abundance
ps.trim.sid.rel <- transform_sample_counts(ps.trim.sid.rare, function(x) x / sum(x))
ord <- ordinate(ps.trim.sid.rel, "PCoA", "bray")
plot_ordination(ps.trim.sid.rel, ord,color="year", shape="year",axes=1:2)+
  stat_ellipse()#+
  facet_wrap(~transect)

#saving
seq.trim.rel <- data.frame(otu_table(ps.trim.rel))

ord <- ordinate(ps.trim.sid.rare, "PCoA", "bray")
plot_ordination(ps.trim.sid.rare, ord,color="transect", shape="year")+
  stat_ellipse()

#now by site
ps.low <- subset_samples(ps.trim.sid.rare,transect=="Lower")
#ps.mnw <- subset_samples(ps.rare.trim,site=="MNW")
ord.low <- ordinate(ps.low, "PCoA", "bray")
gg.low <- plot_ordination(ps.low, ord.low,color="year",axes=1:2)+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()+  
  scale_color_manual(name="Time point",values=c("#781C6D","#F8850F","#AE305C"),labels=c("Baseline","Irma","Post-Irma"))

gg.low
gg.mid
gg.upp

library(ggpubr)
quartz()
ggarrange(gg.low,gg.mid,gg.upp,nrow=1,common.legend=TRUE,legend="right",labels="AUTO")

#now by site
ps.low <- subset_samples(ps.trim.sid.rare,transect=="Upper")
#ps.mnw <- subset_samples(ps.rare.trim,site=="MNW")
ord.low <- ordinate(ps.low, "PCoA", "bray")
plot_ordination(ps.low, ord.low,color="year",axes=1:2)+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()+
  scale_color_manual(name="Time point",values=c("#781C6D","#F8850F","#AE305C"),labels=c("Baseline","Irma","Post-Irma"))

#### pcoa stats ####
library("funfuns")
sid.low.otu <- data.frame(ps.low@otu_table)
sid.low.dat <- data.frame(ps.low@sam_data)
pairwise.adonis(sid.low.otu, factors = sid.low.dat$year, permutations = 999)

dist.seqtab <- vegdist(sid.low.otu)
bet.all <- betadisper(dist.seqtab,sid.low.dat$year)
anova(bet.all)
#not sig - rel or not
plot(bet.all)
adonis(sid.low.otu ~ year, data=sid.low.dat, permutations=999)
#low - all sig different
#mid - all sig different
#upp - all sig different

#overall
row.names.keep <- c(row.names(seq.trim.sid.rare))
samdf.trim.sid.rare.nobis2 <- samdf.sid[(row.names(samdf.sid) %in% row.names.keep), ]
adonis(seq.trim.sid.rare~year*transect, data=samdf.trim.sid.rare.nobis2,permutations=999)
#all sig

ps.mse <- subset_samples(ps.trim.rel,site=="MSE")
#ps.mse <- subset_samples(ps.trim,site=="MSE")
ord.mse <- ordinate(ps.mse, "PCoA", "bray")
gg.mse <- plot_ordination(ps.mse, ord.mse,color="zone", shape="zone")+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  ggtitle("Moorea SE")+
  annotate(geom="text", x=-0.35, y=0.55, label="p < 0.01**",size=4)+#change for non-rarefied
  #xlab("Axis 1 (33.1%)")+ #rarefied
  #ylab("Axis 2 (22.2%)")+ #rarefied
  xlab("Axis 1 (27.9%)")+#non-rarefied
  ylab("Axis 2 (19%)")+#non-rarefied
  theme(axis.text=element_text(size=10))
gg.mse

ps.tnw <- subset_samples(ps.trim.rel,site=="TNW")
#ps.tnw <- subset_samples(ps.trim,site=="TNW")
ord.tnw <- ordinate(ps.tnw, "PCoA", "bray")
gg.tnw <- plot_ordination(ps.tnw, ord.tnw,color="zone", shape="zone")+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  ggtitle("Tahiti NW")+
  #annotate(geom="text", x=-0.4, y=0.6, label="p < 0.01**",size=4)+ #rarefied
  annotate(geom="text", x=-0.25, y=0.6, label="p < 0.01**",size=4)+ #not rarefied
  #xlab("Axis 1 (42.6%)")+ #rarefied
  #ylab("Axis 2 (27.9%)")+ #rarefied
  xlab("Axis 1 (34.3%)")+#non-rarefied
  ylab("Axis 2 (26.5%)")+#non-rarefied
  theme(axis.text=element_text(size=10))
gg.tnw

#### pcoa plot ####
library(ggpubr)

quartz()
ggarrange(gg.mnw,gg.mse,gg.tnw,nrow=1,common.legend=TRUE,legend="right",labels="AUTO")

#### some subsampling for CS506 project: ####
samdf.sid.1516 <- subset(samdf.sid, year==15 | year==16 & transect!="Biscayne")
samdf.sid.17 <- subset(samdf.sid,year==17 & transect!="Biscayne")
samdf.sid.18 <- subset(samdf.sid,year==18 & transect!="Biscayne")

samples.sid.1516 <- rownames(samdf.sid.1516)
samples.sid.17 <- rownames(samdf.sid.17)
samples.sid.18 <- rownames(samdf.sid.18)

seq.trim.sid.1516 <- seq.trim.sid[(row.names(seq.trim.sid) %in% samples.sid.1516), ]
zeros.1516 = seq.trim.sid.1516[,colSums(seq.trim.sid.1516) == 0]
removecols.1516 <- c(colnames(zeros.1516))
seq.trim.sid.1516.no0 <- seq.trim.sid.1516[,!(colnames(seq.trim.sid.1516) %in% removecols.1516)]

seq.trim.sid.17 <- seq.trim.sid[(row.names(seq.trim.sid) %in% samples.sid.17), ]
zeros.17 = seq.trim.sid.17[,colSums(seq.trim.sid.17) == 0]
removecols.17 <- c(colnames(zeros.17))
seq.trim.sid.17.no0 <- seq.trim.sid.17[,!(colnames(seq.trim.sid.17) %in% removecols.17)]

seq.trim.sid.18 <- seq.trim.sid[(row.names(seq.trim.sid) %in% samples.sid.18), ]
zeros.18 = seq.trim.sid.18[,colSums(seq.trim.sid.18) == 0]
removecols.18 <- c(colnames(zeros.18))
seq.trim.sid.18.no0 <- seq.trim.sid.18[,!(colnames(seq.trim.sid.18) %in% removecols.18)]

write.csv(seq.trim.sid.1516.no0,file="seq.trim.sid.1516.no0.csv")
write.csv(seq.trim.sid.17.no0,file="seq.trim.sid.17.no0.csv")
write.csv(seq.trim.sid.18.no0,file="seq.trim.sid.18.no0.csv")

#### round 2 - conglomerate by family ####
#### rename ASVs ####
library(rlang)
library(stringr)

tax <- as.data.frame(ps.trim.sid.rare@tax_table@.Data)

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "D_0__",""),
                        Phylum = str_replace(tax[,2], "D_1__",""),
                        Class = str_replace(tax[,3], "D_2__",""),
                        Order = str_replace(tax[,4], "D_3__",""),
                        Family = str_replace(tax[,5], "D_4__",""),
                        Genus = str_replace(tax[,6], "D_5__",""),
                        Species = str_replace(tax[,7], "D_6__",""),
                        stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
####### Fill holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

tax_table(ps.trim.sid.rare) <- as.matrix(tax.clean)

ps_glom <- tax_glom(ps.trim.sid.rare, "Family")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
plot_bar(ps0, x="year", fill="Family")+
  theme(legend.position="none")+
  facet_wrap(~transect)

seq.glom <- data.frame(otu_table(ps_glom))
sqs <- rownames(tax_table(ps_glom))
tax.glom <- as.data.frame(ps_glom@tax_table@.Data)

rownames(tax.glom) == sqs
colnames(seq.glom) <- tax.glom$Family

samdf.sid.15 <- subset(samdf.sid.rare, year==15 | year==16 & transect!="Biscayne") #57
samdf.sid.17 <- subset(samdf.sid.rare,year==17 & transect!="Biscayne") #56
samdf.sid.18 <- subset(samdf.sid.rare,year==18 & transect!="Biscayne") #35

samples.sid.15 <- rownames(samdf.sid.15)
samples.sid.17 <- rownames(samdf.sid.17)
samples.sid.18 <- rownames(samdf.sid.18)

seq.glom.15 <- seq.glom[(row.names(seq.glom) %in% samples.sid.15), ]
zeros.15 = seq.glom.15[,colSums(seq.glom.15) == 0]
removecols.15 <- c(colnames(zeros.15))
seq.glom.15.no0 <- seq.glom.15[,!(colnames(seq.glom.15) %in% removecols.15)]

seq.glom.17 <- seq.glom[(row.names(seq.glom) %in% samples.sid.17), ]
zeros.17 = seq.glom.17[,colSums(seq.glom.17) == 0]
removecols.17 <- c(colnames(zeros.17))
seq.glom.17.no0 <- seq.glom.17[,!(colnames(seq.glom.17) %in% removecols.17)]

seq.glom.18 <- seq.glom[(row.names(seq.glom) %in% samples.sid.18), ]
zeros.18 = seq.glom.18[,colSums(seq.glom.18) == 0]
removecols.18 <- c(colnames(zeros.18))
seq.glom.18.no0 <- seq.glom.18[,!(colnames(seq.glom.18) %in% removecols.18)]

write.csv(seq.glom.15.no0,file="seq.glom.15.no0.csv")
write.csv(seq.glom.17.no0,file="seq.glom.17.no0.csv")
write.csv(seq.glom.18.no0,file="seq.glom.18.no0.csv")



