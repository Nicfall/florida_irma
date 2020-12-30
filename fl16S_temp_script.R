#fl16s temp
#overview: just sids data frame - already had chloroplasts etc removed
#remove underrepresented shits
#rarefy
#alpha diversity
#PCOA, unifrac shits
#community bar graph

#only sids - 0 counts removed
setwd("~/florida_irma/fl16s")
seqtab.sid.no0 <- readRDS(file="fl16s_seqtab.sid.rds") #193 samples, 28.6k taxa
samdf.sid <- read.csv("fl16s_samdf.sid.csv",row.names=1) #193 samples
taxa2 <- readRDS("fl16s_taxa2.rds") #79k taxa

#### rarefy ####
library(vegan)

rarecurve(seq.trim.sid,step=100,label=FALSE) #after removing contaminants & underrepresented
#did seq.trim.sid a different time - very batchy

total <- rowSums(seqtab.sid.no0)
total
toofew <- subset(total, total <9100)
#20 samples to remove, ouch

row.names.remove <- c(names(toofew))
seqtab.sid.less <- seqtab.sid.no0[!(row.names(seqtab.sid.no0) %in% row.names.remove),]
samdf.sid.rare <- samdf.sid[!(row.names(samdf.sid) %in% row.names.remove), ]
#170 samples left (rarefying to 6k)
#173 left rarefying to 9100

seqtab.sid.rare <- rrarefy(seqtab.sid.less,sample=9100)
#rarecurve(seqtab.sid.rare,step=100,label=FALSE)
saveRDS(seqtab.sid.rare, file="fltemp_seqtab.sid.rare9100.rds")
seqtab.sid.rare <- readRDS(file="fltemp_seqtab.sid.rare9100.rds")

sid.rare.zero = seqtab.sid.rare[,colSums(seqtab.sid.rare) == 0]
ncol(sid.rare.zero) #4105 rarefiying to 9100 - none when trimming before
removecols.sid <- c(colnames(sid.rare.zero))
seqtab.sid.rare.no0 <- seqtab.sid.rare[,!(colnames(seqtab.sid.rare) %in% removecols.sid)]
#173 samples, 24.5k taxa
seqtab.sid.rare.no0 <- data.frame(seqtab.sid.rare.no0)

samdf.sid.rare$year <- as.factor(samdf.sid.rare$year)
samdf.sid.rare$year <- gsub("16","15",samdf.sid.rare$year)
samdf.sid.rare$year <- as.factor(samdf.sid.rare$year)
samdf.sid.rare <- subset(samdf.sid.rare,transect!="Biscayne")

#phyloseq object but rarefied
library('phyloseq')

samdf.sid.rare$year <- as.factor(samdf.sid.rare$year)
ps.sid.rare <- phyloseq(otu_table(seqtab.sid.rare.no0, taxa_are_rows=FALSE), 
                        sample_data(samdf.sid.rare), 
                        tax_table(taxa2))
ps.sid.rare #582 taxa & 170 samples when trimming first
#24.5k taxa & 153 samples rarefying to 9k & removing biscayne

#### trim underrepresented ASVs ####
library(MCMC.OTU)

#formatting the table for mcmc.otu - requires one first column that's 1 through whatever
#& has "X" as column name
nums.sid <- 1:nrow(seqtab.sid.rare.no0)
samples.sid <- rownames(seqtab.sid.rare.no0)

sid.int <- cbind(sample = 0, seqtab.sid.rare.no0)
sid.formcmc <- cbind(X = 0, sid.int)

sid.formcmc$X <- nums.sid
sid.formcmc$sample <- samples.sid

#regular
library(MCMC.OTU)
seq.trim.sid <- purgeOutliers(sid.formcmc,count.columns=3:24521,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)
#3 bad samples - "C10" "E8"  "H6"
#673 ASVs passing cutoff for reads
#584 ASVs show up in 2% of samples
#remove bad samples from sample data frame
#row.names.remove <- c("C10","E8","H6")
#samdf.trim.sid <- samdf.sid[!(row.names(samdf.sid) %in% row.names.remove), ]
##9100 rare, then trim: no bad samples, 549 ASVs left

#remove sample info
seq.trim.sid <- seq.trim.sid[,3:584] #regular (not rarefied)
seq.trim.sid <- seq.trim.sid[,3:551] #rarefied

#write.csv(seq.trim.sid,file="seq.trim.sid.1e4.csv")
#seq.trim.sid <- read.csv("seq.trim.sid.1e4.csv",row.names=1)

#remake phyloseq object
library(phyloseq)
ps.sid.trim <- phyloseq(otu_table(seq.trim.sid, taxa_are_rows=FALSE), 
                        sample_data(samdf.sid), 
                        tax_table(taxa2))
ps.sid.trim #611 taxa left, 582 after more strict 2% of samples
#9100 rare then trim = 549 taxa & 173 samples left

#### rarefy rarefied trimmed ####
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


#### alpha diversity ####
#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)
#Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).
#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index (1 − λ).
plot_richness(ps.sid.rare, x="transect_zone_year", measures=c("Shannon", "Simpson","Observed"), color="zone") + theme_bw()

#df <- data.frame(estimate_richness(ps.sid, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))
#df
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
gg.sh <- ggplot(df.div.nobis, aes(x=year,y=Shannon,color=year))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Year")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  #facet_wrap(~transect)+
  scale_color_manual(values=c("#781C6D","#F8850F","#AE305C"),name="Time point",labels=c("Baseline","Irma","Post-Irma"))+
  scale_x_discrete(labels=c("2015","2017","2018"))+
  theme(axis.text.x=element_text(angle=45,hjust=1))
gg.sh

gg.ob <- ggplot(df.div.nobis, aes(x=year,y=Observed,color=year))+
  geom_boxplot()+
  xlab("Year")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  #facet_wrap(~transect)+
  scale_color_manual(values=c("#781C6D","#F8850F","#AE305C"),name="Time point",labels=c("Baseline","Irma","Post-Irma"))+
  scale_x_discrete(labels=c("2015","2017","2018"))+
  theme(axis.text.x=element_text(angle=45,hjust=1))
gg.ob

gg.si <- ggplot(df.div.nobis, aes(x=year,y=InvSimpson,color=year))+
  geom_boxplot()+
  xlab("Year")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  #facet_wrap(~transect)+
  scale_color_manual(values=c("#781C6D","#F8850F","#AE305C"),name="Time point",labels=c("Baseline","Irma","Post-Irma"))+
  scale_x_discrete(labels=c("2015","2017","2018"))+
  theme(axis.text.x=element_text(angle=45,hjust=1))
gg.si

#shannon
shapiro.test(df.div.nobis$Shannon)

a.all <- aov(Shannon~year,data=df.div.nobis)
summary(a.all)
TukeyHSD(a.all) #nothing

#observed
shapiro.test(df.div.nobis$Observed)
df.div.nobis$ob.log <- log(df.div.nobis$Observed)
shapiro.test(df.div.nobis$ob.log)

k.ob <- kruskal.test(Observed~year,data=df.div.nobis)
k.ob

w.ob <- pairwise.wilcox.test(df.div.nobis$Observed, df.div.nobis$year,p.adjust.method="BH")
w.ob
# 15      17   
# 17 0.069   -    
#   18 4.6e-05 0.017

#anova just to check
a.ob <- aov(Observed~year,data=df.div.nobis)
summary(a.ob)
TukeyHSD(a.ob) #same thing - marginally significant between '15 & '17, super with the '18

#simpson
shapiro.test(log(df.div.nobis$InvSimpson))

a.si <- aov(log(InvSimpson)~year,data=df.div.nobis)
summary(a.si)
TukeyHSD(a.si) #17-18 not different but the rest are..

#### pcoa ####
library('ggplot2')
library('cowplot')

ord.sid.rare <- ordinate(ps.trim.sid.rare, "PCoA", "bray")
gg.pc.all <- plot_ordination(ps.trim.sid.rare, ord.sid.rare,color="year", shape="year")+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()
gg.pc.all

ps.low <- subset_samples(ps.sid.rare,transect=="Upper")
#ps.mnw <- subset_samples(ps.rare.trim,site=="MNW")
ord.low <- ordinate(ps.low, "PCoA", "bray")
gg.low <- plot_ordination(ps.low, ord.low,color="year",axes=1:2)+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()+  
  scale_color_manual(name="Time point",values=c("#781C6D","#F8850F","#AE305C"),labels=c("Baseline","Irma","Post-Irma"))

gg.low
