library(Biostrings)
library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)
library(reshape2)
library(ggbeeswarm)
library(umap)
library(Rtsne)
library(factoextra)
library(org.Hs.eg.db)

donor78 <- read.table("Infos_BCOR-Mut_BCOR_WT_FinalList.txt", sep = "\t", header = T, as.is = T)
rownames(donor78) <- donor78$Sample
donor78 <- donor78[,-1]
donor78$Genotype <- factor(donor78$Genotype)
donor78$Gender <- factor(donor78$Gender)
donor78$Group <- factor(donor78$Group)
donor78$MDS <- factor(donor78$MDS)
donor78$BCOR_vaf <- factor(donor78$BCOR_vaf)

expression_ips <- read.table("INSIGNIA_78_Samples_Study.txt", sep = "\t", header = T, as.is = T)
donor78 <- donor78[names(expression_ips[,-1]),]

rownames(expression_ips) <- expression_ips$GeneID
dds <- DESeqDataSetFromMatrix(countData = expression_ips[,-1],
                              colData = donor78,
                             design = ~ BCOR_vaf)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#dds <- DESeq(dds) # no effect to dds

# vst transformation
vstd <- vst(dds, blind=FALSE)
pcadata<-plotPCA(vstd,intgroup=c("BCOR_vaf","Gender"), returnData = TRUE)
#print(plot_pca)
percentVar <- round(100 * attr(pcadata, "percentVar")) 
pdf(file="PCA_expression_BCOR.pdf", onefile=TRUE,height=3.5,width=5, useDingbats=FALSE,encoding = "ISOLatin2.enc")
p <- ggplot(pcadata, aes(x=PC1, y=PC2,colour=BCOR_vaf, shape=Gender)) +  geom_point(size=3)+
  scale_shape_manual(values = c(Male = 17, Female = 16))+
  scale_colour_manual(values=c("#D55E00","#56B4E9", "#E69F00", "#999999"))+
#  geom_text_repel(aes_string(label = "name"),size=3.5)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))
p <- p + theme(axis.text.x=element_text(size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
print(p)
dev.off()

require("ggrepel")
pdf(file="PCA_expression_BCOR_withlabel.pdf", onefile=TRUE,height=10,width=15, useDingbats=FALSE,encoding = "ISOLatin2.enc")
p <- ggplot(pcadata, aes(x=PC1, y=PC2,colour=BCOR_vaf, shape=Gender)) +  geom_point(size=3)+
  scale_shape_manual(values = c(Male = 17, Female = 16))+
  scale_colour_manual(values=c("#D55E00","#56B4E9", "#E69F00", "#999999"))+
  geom_text_repel(aes_string(label = "name"),size=3.5)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))
p <- p + theme(axis.text.x=element_text(size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
print(p)
dev.off()


################################
# Differentiation_D0_D6_D12_D27
# 34,40
################################
expression_2ips <- read.table("Counts_Directed_Differentiation_D0_D6_D12_D27.txt", sep = "\t", header = T, as.is = T)

# Time course of PCA of all gene expression for iPSCs MSH34i and MSH40i
donor2 <- read.table("Info_AllCells_BCOR_Low_vs_High_34_40_Samples.txt", sep = "\t", header = T, as.is = T)
rownames(donor2) <- donor2$Sample
donor2 <- donor2[,2:4]
donor2$BCOR_Expression <- factor(donor2$BCOR_Expression)
donor2$CellLine <- factor(donor2$CellLine)
donor2$CellType <- factor(donor2$CellType)

rownames(expression_2ips) <- expression_2ips$GeneID
dds <- DESeqDataSetFromMatrix(countData = expression_2ips[,-1],
                              colData = donor2,
                              design = ~ BCOR_Expression)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#dds <- DESeq(dds) # no effect to dds

# vst transformation
vstd <- vst(dds, blind=FALSE)
pcadata<-plotPCA(vstd,intgroup=c("BCOR_Expression","CellType"), returnData = TRUE)
#print(plot_pca)
percentVar <- round(100 * attr(pcadata, "percentVar")) 
pdf(file="PCA_expression_BCOR_34_40.pdf", onefile=TRUE,height=3.5,width=5, useDingbats=FALSE,encoding = "ISOLatin2.enc")
p <- ggplot(pcadata, aes(x=PC1, y=PC2,colour=BCOR_Expression, shape=CellType)) +  geom_point(size=3)+
  scale_shape_manual(values = c(0,3,4,8))+
  scale_colour_manual(values=c("#D55E00","#56B4E9", "#E69F00", "#999999"))+
  #  geom_text_repel(aes_string(label = "name"),size=3.5)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))
p <- p + theme(axis.text.x=element_text(size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
print(p)
dev.off()


# Time course of marker genes
expression_2ips <- read.table("Counts_Directed_Differentiation_D0_D6_D12_D27.txt", sep = "\t", header = T, as.is = T)
markers <- read.table("Heatmap_Markers_Ordered.txt", sep = "\t", header = T, as.is = T)
marker_expression <- merge(expression_2ips,markers,by="GeneID")
rownames(marker_expression) <- marker_expression$GeneName
#marker_expression <- marker_expression[order(marker_expression$Group, marker_expression$GeneName),]

marker_expression <- marker_expression[match(markers$GeneName,marker_expression$GeneName),]
donor2 <- read.table("Info_AllCells_BCOR_Low_vs_High_34_40_Samples.txt", sep = "\t", header = T, as.is = T)
rownames(donor2) <- donor2$Sample
donor2 <- donor2[,2:4]
donor2$BCOR_Expression <- factor(donor2$BCOR_Expression)
donor2$CellLine <- factor(donor2$CellLine)
donor2$CellType <- factor(donor2$CellType)


dds2 <- DESeqDataSetFromMatrix(countData = marker_expression[,2:25],
                              colData = donor2,
                              design = ~ BCOR_Expression)

# Pre-filtering
keep <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep,]

dds2 <- estimateSizeFactors(dds2)
dds2_nor <- counts(dds2, normalized=TRUE)
dds2_nor_relative <- log2(dds2_nor+1)-rowMeans(log2(dds2_nor+1))
msh34_nor_relative <- dds2_nor_relative[,c(1:3,7:9,13:15,19:21)]
msh40_nor_relative <- dds2_nor_relative[,c(4:6,10:12,16:18,22:24)]


sample_nor_relatvie <- msh40_nor_relative

# heatmap
#colnames(sample_nor_relatvie) <- sig_total$MutationType
sample_nor_relatvie <- round(sample_nor_relatvie,3)
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("blue","white", "red"))(n = 128)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(-14,  # for blue
               0,           # for white
               1,max(sample_nor_relatvie))             # for green

pdf(file="msh40_nor_relative_heatmap_log2_v32.pdf", h=8, w=10, onefile=TRUE)
gplots::heatmap.2(sample_nor_relatvie,
                 # hclustfun=function(x) hclust(x,method = 'median'),
                  reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean), # Reorder dendrogram by branch mean rather than sum
                #  main = "Heatmap of substitution signature", # heat map title
                  notecol="black",      # change font color of cell labels to black
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",         # turns off trace lines inside the heat map
                  margins =c(12,9),     # widens margins around plot
                  keysize=1,
                  col=my_palette,       # use on color palette defined earlier
                  key = TRUE, 
                  ColSideColors=rep(c("#CCCCCC","#999999","#666666","#333333"), each=3),
                  RowSideColors=rep(c("#E4BFF7","#EB7CF7","#F200FF","#820793"), c(9,10,12,7)),
    #              breaks=col_breaks,    # enable color transition at specified limits
                 # dendrogram="row",     # only draw a row dendrogram
                  Rowv = "NA",           # turn off row clustering
                  Colv="NA")            # turn off column clustering
dev.off()


####################### END ####################### 





