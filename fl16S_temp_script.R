#fl16s temp
#overview: just sids data frame
#remove underrepresented shits
#rarefy
#alpha diversity
#PCOA, unifrac shits
#community bar graph

#only sids - 0 counts removed
seqtab.sid.no0 <- readRDS(file="fl16s_seqtab.sid.rds")
samdf.sid <- read.csv("fl16s_samdf.sid.csv",row.names=1)
taxa2 <- readRDS("fl16s_taxa2.rds")



