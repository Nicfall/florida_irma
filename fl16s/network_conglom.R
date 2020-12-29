# Install Required packages 
#install.packages("igraph") 
library("igraph")
#install.packages("qgraph") 
library("qgraph")
#install.packages("MCL")
library("MCL")
#install.packages("vegan")
library("vegan")
# Install SpiecEasi package 
#install.packages("devtools") 
library("devtools") 
#install_github("zdk123/SpiecEasi")
library("SpiecEasi")
library("ggplot2")
library("ggpubr")

### Reading in & renaming data

data_15 <- read.csv("seq.glom.15.no0.csv",row.names=1) #from 2015&2016
data_17 <- read.csv("seq.glom.17.no0.csv",row.names=1) #from 2017
data_18 <- read.csv("seq.glom.18.no0.csv",row.names=1) #from 2018

#removing samples with less than an average of 2 reads per OTU - recommended by authors
data_15_2 <- data_15[,-which(colMeans(data_15) <= 2)]
data_17_2 <- data_17[,-which(colMeans(data_17) <= 2)]
data_18_2 <- data_18[,-which(colMeans(data_18) <= 2)]

### Building the network

This is a network of correlations of the abundances of the types of bacteria across samples. Following instructions [here](https://www-ncbi-nlm-nih-gov.ezproxy.bu.edu/pubmed/30298259).

# SparCC network
sparcc.matrix_15 <- sparcc(data_15_2)
sparcc.matrix_17 <- sparcc(data_17_2)
sparcc.matrix_18 <- sparcc(data_18_2)

sparcc.cutoff <- 0.3
sparcc.adj_15 <- ifelse(abs(sparcc.matrix_15$Cor) >= sparcc.cutoff, 1, 0)
sparcc.adj_17 <- ifelse(abs(sparcc.matrix_17$Cor) >= sparcc.cutoff, 1, 0)
sparcc.adj_18 <- ifelse(abs(sparcc.matrix_18$Cor) >= sparcc.cutoff, 1, 0)

# Add OTU names to rows and columns
rownames(sparcc.adj_15) <- colnames(data_15_2) 
colnames(sparcc.adj_15) <- colnames(data_15_2) 

rownames(sparcc.adj_17) <- colnames(data_17_2) 
colnames(sparcc.adj_17) <- colnames(data_17_2) 

rownames(sparcc.adj_18) <- colnames(data_18_2) 
colnames(sparcc.adj_18) <- colnames(data_18_2) 

# Build network from adjacency
sparcc.net_15 <- graph.adjacency(sparcc.adj_15,
                                   mode = "undirected",
                                   diag = FALSE)
sparcc.net_17 <- graph.adjacency(sparcc.adj_17,
                                 mode = "undirected",
                                 diag = FALSE)
sparcc.net_18 <- graph.adjacency(sparcc.adj_18,
                                 mode = "undirected",
                                 diag = FALSE)

# Use sparcc.net for the rest of the method 
net_15 <- sparcc.net_15
net_17 <- sparcc.net_17
net_18 <- sparcc.net_18
# Hub detection
net_15.cn <- closeness(net_15)
net_15.bn <- betweenness(net_15) 
net_15.pr <- page_rank(net_15)$vector 
net_15.hs <- hub_score(net_15)$vector
net_15.hs.sort <- sort(net_15.hs, decreasing = TRUE)
net_15.hs.top5 <- head(net_15.hs.sort, n = 5)

net_17.cn <- closeness(net_17)
net_17.bn <- betweenness(net_17) 
net_17.pr <- page_rank(net_17)$vector 
net_17.hs <- hub_score(net_17)$vector
net_17.hs.sort <- sort(net_17.hs, decreasing = TRUE)
net_17.hs.top5 <- head(net_17.hs.sort, n = 5)

net_18.cn <- closeness(net_18)
net_18.bn <- betweenness(net_18) 
net_18.pr <- page_rank(net_18)$vector 
net_18.hs <- hub_score(net_18)$vector
net_18.hs.sort <- sort(net_18.hs, decreasing = TRUE)
net_18.hs.top5 <- head(net_18.hs.sort, n = 5)

### Plotting the network

# Function 2: Plot network with node size scaled to hubbiness 
plot.net <- function(net, scores, outfile, title) {
  # Convert node label from names to numerical IDs. 
  features <- V(net)$name
  col_ids <- seq(1, length(features))
  V(net)$name <- col_ids
  node.names <- features[V(net)] # Nodes’ color.
  V(net)$color <- "white"
  # Define output image file.
  outfile <- paste(outfile, "jpg", sep=".") # Image properties.
  jpeg(outfile, width = 4800, height = 9200, res = 300, quality = 100)
  par(oma = c(4, 1, 1, 1))
  # Main plot function.
  plot(net, vertex.size = (scores*5)+4, vertex.label.cex = 1) 
  title(title, cex.main = 4)
  # Plot legend containing OTU names.
  #labels = paste(as.character(V(net)), node.names, sep = ") ") 
  #legend("bottom", legend = labels, xpd = TRUE, ncol = 5, cex = 1.2) 
  dev.off()
}

plot.net(net_15, net_15.hs, outfile = "network_family_15", title = "Network - 2015 & 2016")
plot.net(net_17, net_17.hs, outfile = "network_family_17", title = "Network - 2017")
plot.net(net_18, net_18.hs, outfile = "network_family_18", title = "Network - 2018")
```

### Looking at results

eigen.15 <- eigen_centrality(net_15)
eigen.15.tosort <- data.frame(eigen.15[["vector"]])
eigen.15.tosort$vector <- eigen.15.tosort$eigen.15...vector...
eigen.15.sorted <- eigen.15.tosort[order(-eigen.15.tosort$vector),]
eigen.15.sorted$sqs <- rownames(eigen.15.sorted)
eigen.15.top10 <- eigen.15.sorted[1:10,]
sqs_order <- eigen.15.top10$sqs
eigen.15.top10$sqs <- factor(eigen.15.top10$sqs,levels=sqs_order)
gg.eig.15 <- ggplot(eigen.15.top10,aes(x=sqs,y=vector))+
  geom_bar(stat="identity")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  ylab("Eigenvector centrality")+
  xlab("Bacterial family")+
  ggtitle("Pre-Irma (2015/16)")+
  geom_hline(yintercept=0.7885872,color="#781C6D")
gg.eig.15

eigen.17 <- eigen_centrality(net_17)
eigen.17.tosort <- data.frame(eigen.17[["vector"]])
eigen.17.tosort$vector <- eigen.17.tosort$eigen.17...vector...
eigen.17.sorted <- eigen.17.tosort[order(-eigen.17.tosort$vector),]
eigen.17.sorted$sqs <- rownames(eigen.17.sorted)
eigen.17.top10 <- eigen.17.sorted[1:10,]
sqs_order <- eigen.17.top10$sqs
eigen.17.top10$sqs <- factor(eigen.17.top10$sqs,levels=sqs_order)
gg.eig.17 <- ggplot(eigen.17.top10,aes(x=sqs,y=vector))+
  geom_bar(stat="identity")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  ylab("Eigenvector centrality")+
  xlab("Bacterial family")+
  ggtitle("Irma (2017)")+
  geom_hline(yintercept=0.7885872,color="#781C6D")+
  geom_hline(yintercept=0.7140983,color="#AE305C")
gg.eig.17

eigen.18 <- eigen_centrality(net_18)
eigen.18.tosort <- data.frame(eigen.18[["vector"]])
eigen.18.tosort$vector <- eigen.18.tosort$eigen.18...vector...
eigen.18.sorted <- eigen.18.tosort[order(-eigen.18.tosort$vector),]
eigen.18.sorted$sqs <- rownames(eigen.18.sorted)
eigen.18.top10 <- eigen.18.sorted[1:10,]
sqs_order <- eigen.18.top10$sqs
eigen.18.top10$sqs <- factor(eigen.18.top10$sqs,levels=sqs_order)
gg.eig.18 <- ggplot(eigen.18.top10,aes(x=sqs,y=vector))+
  geom_bar(stat="identity")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  ylab("Eigenvector centrality")+
  xlab("Bacterial family")+
  ggtitle("Post-Irma (2018)")+
  geom_hline(yintercept=0.7140983,color="#AE305C")
gg.eig.18

#colors:
#baseline, irma, recovery
#"#781C6D","#F8850F","#AE305C"

quartz()
gg.eig.15
quartz()
gg.eig.17
quartz()
gg.eig.18

ggarrange(gg.eig.15,gg.eig.17,gg.eig.18,ncol=3,heights=c(2,1,2))

#are any shared?
merge(eigen.18.top10,eigen.17.top10,by="sqs")
#Burkh
#Sandara
#Woesei
merge(eigen.17.top10,eigen.15.top10,by="sqs")
#Rhizo
#Sand
#Woe
wow <- merge(eigen.18.top10,eigen.15.top10,by="sqs")
#Caulobacteraceae        
#Kiloniellaceae          
#Order_Thalassobaculales
#Sandaracinaceae         
#Woeseiaceae 

## Eigenvector centrality

```{r network networks 2}
btwn.1516.tosort <- data.frame(betweenness(net_1516))
btwn.1516.tosort$vector <- btwn.1516.tosort$betweenness.net_1516.
btwn.1516.sorted <- btwn.1516.tosort[order(-btwn.1516.tosort$vector),]
btwn.1516.sorted$sqs <- rownames(btwn.1516.sorted)
btwn.1516.top10 <- btwn.1516.sorted[1:10,]
sqs_order <- btwn.1516.top10$sqs
btwn.1516.top10$sqs <- factor(btwn.1516.top10$sqs,levels=sqs_order)
gg.btwn.1516 <- ggplot(btwn.1516.top10,aes(x=sqs,y=vector))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  ylab("Betweenness centrality")+
  xlab("Bacterial sequence")+
  ggtitle("2015 & 2016 (pre-Irma)")

btwn.17.tosort <- data.frame(betweenness(net_17))
btwn.17.tosort$vector <- btwn.17.tosort$betweenness.net_17.
btwn.17.sorted <- btwn.17.tosort[order(-btwn.17.tosort$vector),]
btwn.17.sorted$sqs <- rownames(btwn.17.sorted)
btwn.17.top10 <- btwn.17.sorted[1:10,]
sqs_order <- btwn.17.top10$sqs
btwn.17.top10$sqs <- factor(btwn.17.top10$sqs,levels=sqs_order)
gg.btwn.17 <- ggplot(btwn.17.top10,aes(x=sqs,y=vector))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  ylab("Betweenness centrality")+
  xlab("Bacterial sequence")+
  ggtitle("2017 (Irma)")

btwn.18.tosort <- data.frame(betweenness(net_18))
btwn.18.tosort$vector <- btwn.18.tosort$betweenness.net_18.
btwn.18.sorted <- btwn.18.tosort[order(-btwn.18.tosort$vector),]
btwn.18.sorted$sqs <- rownames(btwn.18.sorted)
btwn.18.top10 <- btwn.18.sorted[1:10,]
sqs_order <- btwn.18.top10$sqs
btwn.18.top10$sqs <- factor(btwn.18.top10$sqs,levels=sqs_order)
gg.btwn.18 <- ggplot(btwn.18.top10,aes(x=sqs,y=vector))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  ylab("Betweenness centrality")+
  xlab("Bacterial sequence")+
  ggtitle("2018 (post-Irma)")

ggarrange(gg.btwn.1516,gg.btwn.17,gg.btwn.18)
```

Similar pattern as in the eigenvector centrality plots

Left to do:
- More metric analysis & all the statistics!


