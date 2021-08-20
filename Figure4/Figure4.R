source("/rds/project/sn206/rds-sn206-nik-zainal/users/xz388/cancer_archive04/b_1176/15_subs_20180907/00_data/sub_common.R")

##############
#  No SNP
##############
# HipSci Skin drivers
hipsci_disub <- read.table("hipsci_disubs_FATHMM_uniq_pathogenic_petr_final.txt", sep="\t", header = T, as.is = T)
hdisub_s <- hipsci_disub[hipsci_disub$SNP=="n" & hipsci_disub$ips_nALT>0,c("ips","Fibro","Fibro_nALT","ips_nALT","gene","mutationAA","FATHMM","FATHMM_score")] # 58
hdisub_s$type <- "disub"
  
hipsci_sub <- read.table("hipsci_snvs_FATHMM_uniq_pathogenic_petr_final.txt", sep="\t", header = T, as.is = T)
hsub_s <- hipsci_sub[hipsci_sub$SNP=="n" & hipsci_sub$ips_nALT>0,c("ips","Fibro","Fibro_nALT","ips_nALT","gene","mutationAA","FATHMM","FATHMM_score")] # 247
hsub_s$type <- "sub"

hipsci_indel <- read.table("hipsci_indel_cancergenes_info_petr_frameshift_final.validated.txt", sep="\t", header = T, as.is = T)
hindel_s <- hipsci_indel[hipsci_indel$Validate%in%c("tp","exp"),c("ips","Fibro","Fibro_nALT","ips_nALT","gene")]  # 45
hindel_s$mutationAA <- "frameshift"
hindel_s$FATHMM <- "PATHOGENIC"
hindel_s$FATHMM_score <- 1
hindel_s$type <- "indel"

h_muts <- rbind(hdisub_s,hsub_s)
h_muts <- rbind(h_muts,hindel_s)
write.table(h_muts,"hipsci_driver_noSNP.txt",sep = "\t", col.names = T, row.names = F, quote = F)



h_muts <- read.table("hipsci_driver_noSNP.txt", sep="\t", header = T, as.is = T)
h_muts$driver <- 1
h_muts_freq <- data.frame(table(h_muts$gene))
h_muts_freq <- h_muts_freq[order(h_muts_freq$Freq,decreasing = F),]
names(h_muts_freq) <- c("gene","Freq")
h_muts_2 <- h_muts[h_muts$gene%in% as.character(h_muts_freq[h_muts_freq$Freq>2,1]),]
h_muts_3 <- data.frame(table(h_muts_2$ips,h_muts_2$gene))
names(h_muts_3) <- c("ips","gene","driver")

# plot
pdf(file="hipsci_driver.pdf", onefile=TRUE,width = 10,height = 4)
g <-ggplot(h_muts_3, aes(y=factor(gene), x=factor(ips))) + geom_tile(aes(fill=driver),colour="black")
g <- g+xlab("F-hiPSCs")+ylab("Gene")
g <- g+scale_y_discrete(limits = as.character(h_muts_freq[h_muts_freq$Freq>2,1]))+scale_x_discrete(position = "top") 
g <- g+scale_fill_gradient2(high="#CA9F92", low="white",limits=c(0, 1))
g <- g+theme(axis.text.x=element_text(size=8,angle=90, hjust=0.9,vjust=0.9,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.border = element_rect(colour = "black", fill=NA)
)
print(g)
#ggMarginal(g, type = "histogram", fill="transparent")
dev.off()


pdf(file="hipsci_driver_freq.pdf", onefile=TRUE,width = 2,height = 4)
g <-ggplot(h_muts_freq[h_muts_freq$Freq>2,], aes(x=factor(gene), y=Freq)) + geom_bar(fill="#CA9F92",stat="identity")
g <- g+xlab("Freq")+ylab("Gene")
g <- g+scale_x_discrete(limits = as.character(h_muts_freq[h_muts_freq$Freq>2,1]))+ coord_flip()
#g <- g+scale_fill_gradient2(high="#CA9F92", low="white",limits=c(0, 1))
g <- g+theme(axis.text.x=element_text(size=8,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA)
)
print(g)
dev.off()


# mutation burden
hipsci_sub_burden <- read.table("sample_summary_3experiments.txt",sep="\t", header = T, as.is = T)
h_sub_burden <- hipsci_sub_burden[hipsci_sub_burden$ips %in% h_muts_3$ips,]
h_sub_burden <- h_sub_burden[order(h_sub_burden$wgs,decreasing = T),]

# Insignia blood drivers
insignia_muts <- read.table("insignia_codingmuts_score.txt", sep="\t", header = T, as.is = T)
insignia_muts <- insignia_muts[,c("Sample","patient","celltype","gene_name")]





##############
#  Wwith SNP
##############
hipsci_disub <- read.table("hipsci_disubs_FATHMM_uniq_pathogenic_petr_final.txt", sep="\t", header = T, as.is = T)
hdisub_s <- hipsci_disub[hipsci_disub$ips_nALT>0,c("ips","Fibro","Fibro_nALT","ips_nALT","gene","mutationAA","FATHMM","FATHMM_score")] # 58
hdisub_s$type <- "disub"

hipsci_sub <- read.table("hipsci_snvs_FATHMM_uniq_pathogenic_petr_final.txt", sep="\t", header = T, as.is = T)
hsub_s <- hipsci_sub[hipsci_sub$ips_nALT>0,c("ips","Fibro","Fibro_nALT","ips_nALT","gene","mutationAA","FATHMM","FATHMM_score")] # 247
hsub_s$type <- "sub"

hipsci_indel <- read.table("hipsci_indel_cancergenes_info_petr_frameshift_final.validated.txt", sep="\t", header = T, as.is = T)
hindel_s <- hipsci_indel[hipsci_indel$Validate%in%c("tp","exp"),c("ips","Fibro","Fibro_nALT","ips_nALT","gene")]  # 45
hindel_s$mutationAA <- "frameshift"
hindel_s$FATHMM <- "PATHOGENIC"
hindel_s$FATHMM_score <- 1
hindel_s$type <- "indel"

h_muts <- rbind(hdisub_s,hsub_s)
h_muts <- rbind(h_muts,hindel_s)
write.table(h_muts,"hipsci_driver_withSNP.txt",sep = "\t", col.names = T, row.names = F, quote = F)


h_muts <- read.table("hipsci_driver_withSNP.txt", sep="\t", header = T, as.is = T)
h_muts$driver <- 1
h_muts_freq <- data.frame(table(h_muts$gene))
h_muts_freq <- h_muts_freq[order(h_muts_freq$Freq,decreasing = F),]
names(h_muts_freq) <- c("gene","Freq")
h_muts_2 <- h_muts[h_muts$gene%in% as.character(h_muts_freq[h_muts_freq$Freq>2,1]),]
write.table(h_muts_2,"hipsci_driver_withSNP_freq_g2.txt",sep = "\t", col.names = T, row.names = F, quote = F)

h_muts_2_freq <- data.frame(table(h_muts_2$gene))
h_muts_2_freq <- h_muts_2_freq[order(h_muts_2_freq$Freq,decreasing = F),]
names(h_muts_2_freq) <- c("gene","geneFreq")
h_muts_2_freq$rank <- seq_len(nrow(h_muts_2_freq))

# order ips according to gene frequence, use the max rank for fibro
h_muts_2_order <- merge(h_muts_2,h_muts_2_freq,by="gene")
h_muts_2_ipsorder <- unique(h_muts_2[,c("ips","Fibro")])

h_muts_2_fibroorder <- unique(h_muts_2_order[,c("gene","Fibro","rank")])
h_muts_2_fibroorder <- dcast(h_muts_2_fibroorder,Fibro~gene)
h_muts_2_fibroorder[is.na(h_muts_2_fibroorder)] <- 0
h_muts_2_fibroorder$maxrank <- apply(h_muts_2_fibroorder[,2:dim(h_muts_2_fibroorder)[2]], 1, max)

h_muts_2_ipsorder <- merge(h_muts_2_ipsorder, h_muts_2_fibroorder,by="Fibro")
h_muts_2_ipsorder <- h_muts_2_ipsorder[order(h_muts_2_ipsorder$maxrank, h_muts_2_ipsorder$ips, decreasing = T),]


h_muts_3 <- data.frame(table(h_muts_2$ips,h_muts_2$gene))
names(h_muts_3) <- c("ips","gene","driver")
# plot
pdf(file="hipsci_driver_SNP.pdf", onefile=TRUE,width = 10,height = 4)
g <-ggplot(h_muts_3, aes(y=factor(gene), x=factor(ips))) + geom_tile(aes(fill=driver),colour="black")
g <- g+xlab("F-hiPSCs")+ylab("Gene")+scale_x_discrete(limits = as.character(h_muts_2_ipsorder$ips), position = "top")
g <- g+scale_y_discrete(limits = as.character(h_muts_2_freq$gene))
g <- g+scale_fill_gradient2(high="#59ADD0", low="white",limits=c(0, 1))
g <- g+theme(axis.text.x=element_text(size=8,angle=90, hjust=0.9,vjust=0.9,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.border = element_rect(colour = "black", fill=NA)
)
print(g)
#ggMarginal(g, type = "histogram", fill="transparent")
dev.off()

# genes with freq >2
h_muts_g2 <- h_muts[h_muts$gene %in% h_muts_freq[h_muts_freq$Freq>2,1],]
h_muts_g2_freq <- data.frame(table(h_muts_g2$gene,h_muts_g2$type))
names(h_muts_g2_freq) <- c("gene","type","Freq")
h_muts_g2_freq$type <- as.character(h_muts_g2_freq$type)
h_muts_g2_freq[which(h_muts_g2_freq$type=="indel"),"type"] <- "aindel"
pdf(file="hipsci_driver_SNP_freq.pdf", onefile=TRUE,width = 3,height = 4)
g <-ggplot(h_muts_g2_freq, aes(x=factor(gene), y=Freq,fill=type)) + geom_bar(stat="identity",width=0.8)+ylim(0,35)
g <- g+xlab("Freq")+ylab("Gene")
g <- g+scale_x_discrete(limits = as.character(h_muts_freq[h_muts_freq$Freq>2,1]))+ coord_flip()
g <- g+scale_fill_manual(values=c("#CA9F92","#F9CD97","#E3D9B0"))
g <- g+theme(axis.text.x=element_text(size=8,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA)
)
print(g)
dev.off()


# mutation burden
hipsci_sub_burden <- read.table("sample_summary_3experiments.txt",sep="\t", header = T, as.is = T)
h_sub_burden <- hipsci_sub_burden[hipsci_sub_burden$ips %in% h_muts_3$ips,]
h_sub_burden <- h_sub_burden[order(h_sub_burden$wgs,decreasing = T),]



# assign signature to snv drivers



#########################
# Insignia blood drivers
#########################
insignia_muts <- read.table("insignia_muts_pathogenic_validated.txt", sep="\t", header = T, as.is = T)
insignia_muts <- insignia_muts[,c("Sample","patient","celltype","gene_name")]
names(insignia_muts) <- c("Sample","patient","celltype","gene")
insignia_muts$driver <- 1
insignia_muts_freq <- data.frame(table(insignia_muts$gene))
insignia_muts_freq <- insignia_muts_freq[order(insignia_muts_freq$Freq,decreasing = F),]
names(insignia_muts_freq) <- c("gene","Freq")
insignia_muts_2 <- insignia_muts[insignia_muts$gene%in% as.character(insignia_muts_freq[insignia_muts_freq$Freq>=2,1]),]
write.table(insignia_muts_2,"insignia_muts_g2.txt",sep = "\t", col.names = T, row.names = F, quote = F)

#cave_1000gfilter <- read.table("/rds/project/sn206/rds-sn206-nik-zainal/users/xz388/cancer_archive04/g_1133/12_subs/00_data/cave_11142018.tab.1000gFiltered", sep = "\t", header = T, as.is = T)
#samples <- data.frame(table(cave_1000gfilter$Sample))
#names(samples) <- c("Sample","Freq")
#donor21 <- read.table("donor21_included.txt", sep = "\t", header = F, as.is = T)
#samples$patient <- sub("\\i.*","",samples$Sample)
#samples$patient <- sub("e","",samples$patient)
#samples_donor21 <- samples[samples$patient%in% donor21[,1],]
#write.table(samples_donor21,"samples_donor21.txt", sep = "\t", col.names = T, row.names = F, quote = F)

samples_donor21 <- read.table("samples_donor21_final.txt", sep="\t", header = T, as.is = T)
samples_donor21$celltype <- ""
samples_donor21[grepl("i", samples_donor21$Sample, fixed=TRUE),"celltype"] <- "ips_parent"
samples_donor21[grepl("_", samples_donor21$Sample, fixed=TRUE),"celltype"] <- "ips_subclone"


sample_gene <- data.frame("Sample"=samples_donor21$Sample, "BCOR"=0,"DNMT3A"=0, "CBL"=0, "ACVR2A"=0,"FIP1L1"=0, "AKAP9"=0)
sample_gene <- melt(sample_gene,c("Sample"))
names(sample_gene) <- c("Sample","gene","Freq")
sample_gene <- sample_gene[,c("Sample","gene")]
insignia_muts_3 <- merge(sample_gene, insignia_muts_2[,c("Sample","gene","driver")],all.x=T)
insignia_muts_3[is.na(insignia_muts_3)] <- 0



insignia_muts_2_freq <- data.frame(table(insignia_muts_2$gene))
insignia_muts_2_freq <- insignia_muts_2_freq[order(insignia_muts_2_freq$Freq,decreasing = F),]
names(insignia_muts_2_freq) <- c("gene","geneFreq")
insignia_muts_2_freq$rank <- seq_len(nrow(insignia_muts_2_freq))

# order ips according to gene frequence, use the max rank for fibro
insignia_muts_2_order <- merge(insignia_muts_2,insignia_muts_2_freq,by="gene")
insignia_muts_2_ipsorder <- unique(insignia_muts_2[,c("Sample","patient")])

insignia_muts_2_patientorder <- unique(insignia_muts_2_order[,c("gene","patient","rank")])
insignia_muts_2_patientorder <- dcast(insignia_muts_2_patientorder,patient~gene)
insignia_muts_2_patientorder[is.na(insignia_muts_2_patientorder)] <- 0
insignia_muts_2_patientorder$maxrank <- apply(insignia_muts_2_patientorder[,2:dim(insignia_muts_2_patientorder)[2]], 1, max)

insignia_muts_3_ipsorder <- unique(samples_donor21[,c("Sample","patient")])
insignia_muts_3_ipsorder <- merge(insignia_muts_3_ipsorder, insignia_muts_2_patientorder,by="patient", all.x=T)
insignia_muts_3_ipsorder[is.na(insignia_muts_3_ipsorder)] <- 0
insignia_muts_3_ipsorder <- insignia_muts_3_ipsorder[order(insignia_muts_3_ipsorder$maxrank, insignia_muts_3_ipsorder$Sample, decreasing = c(TRUE, FALSE), method="radix"),]


pdf(file="insignia_driver_SNP.pdf", onefile=TRUE,width = 10,height = 2)
g <-ggplot(insignia_muts_3, aes(y=factor(gene), x=factor(Sample))) + geom_tile(aes(fill=driver),colour="black")
g <- g+xlab("B-hiPSCs")+ylab("Gene")
g <- g+scale_y_discrete(limits = as.character(insignia_muts_2_freq$gene))+scale_x_discrete(limits = as.character(insignia_muts_3_ipsorder$Sample), position = "top") 
g <- g+scale_fill_gradient2(high="#DB8DB2", low="white",limits=c(0, 1))
g <- g+theme(axis.text.x=element_text(size=8,angle=90, hjust=0.9,vjust=0.9,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA)
)
print(g)
#ggMarginal(g, type = "histogram", fill="transparent")
dev.off()

insignia_muts <- read.table("insignia_muts_pathogenic_validated.txt", sep="\t", header = T, as.is = T)
insignia_muts_g2 <- insignia_muts[insignia_muts$gene_name %in% insignia_muts_freq[insignia_muts_freq$Freq>=2,1],]
insignia_muts_g2$type <- as.character(insignia_muts_g2$type)
insignia_muts_g2[which(insignia_muts_g2$type!="Sub"),"type"] <- "aindel"
insignia_muts_g2_freq <- data.frame(table(insignia_muts_g2$gene_name,insignia_muts_g2$type))
names(insignia_muts_g2_freq) <- c("gene","type","Freq")

pdf(file="insignia_driver_SNP_freq.pdf", onefile=TRUE,width = 3,height = 1.5)
g <-ggplot(insignia_muts_g2_freq, aes(x=factor(gene), y=Freq,fill=type)) + geom_bar(stat="identity",width=0.8) +ylim(0,35)
g <- g+xlab("Freq")+ylab("Gene")
g <- g+scale_x_discrete(limits = as.character(insignia_muts_freq[insignia_muts_freq$Freq>=2,1]))+ coord_flip()
g <- g+scale_fill_manual(values=c("#CA9F92","#E3D9B0"))
g <- g+theme(axis.text.x=element_text(size=8,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA)
)
print(g)
dev.off()



##################################
# combine two and bar plot of hits in each gene
##################################
i <- insignia_muts_freq[insignia_muts_freq$Freq>=2,]
names(i) <- c("gene","freq")
i$dataset <- "Insignia_B_hiPSCs"
h <- h_muts_freq[h_muts_freq$Freq>2,]
names(h) <- c("gene","freq")
h$dataset <- "Hipsci_F_hiPSCs"
ih <- rbind(i,h)
ih_freq <- data.frame(table(as.character(ih$gene)))

ih_gene_hits <- data.frame("gene"=ih_freq[,1],"Insignia_B_hiPSCs"=0,"Hipsci_F_hiPSCs"=0) 
ih_gene_hits <- melt(ih_gene_hits,c("gene"))
names(ih_gene_hits) <- c("gene","dataset","freq1")

ih_gene_hits <- merge(ih_gene_hits[,c("gene","dataset")], ih,by=c("gene","dataset"),all.x=T)
ih_gene_hits[is.na(ih_gene_hits)] <- 0
ih_gene_hits_dcast <- dcast(ih_gene_hits,gene~dataset)
ih_gene_hits_dcast$total <- ih_gene_hits_dcast$Insignia_B_hiPSCs + ih_gene_hits_dcast$Hipsci_F_hiPSCs
ih_gene_hits_dcast <- ih_gene_hits_dcast[order(ih_gene_hits_dcast$total,decreasing = F),]

pdf(file="ih_driver_SNP.pdf", onefile=TRUE,width = 2.5,height = 5)
g <-ggplot(ih_gene_hits, aes(y=factor(gene), x=factor(dataset))) + geom_tile(aes(fill=freq),colour="black")
g <- g+ylab("Gene")
g <- g+scale_y_discrete(limits = as.character(ih_gene_hits_dcast$gene))+scale_x_discrete(position = "top") 
g <- g+scale_fill_gradient2(high="#B1C27A", low="white",limits=c(0, 45))
g <- g+theme(axis.text.x=element_text(size=8,angle=90, hjust=0.9,vjust=0.9,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA)
)
print(g)
#ggMarginal(g, type = "histogram", fill="transparent")
dev.off()

pdf(file="ih_driver_SNP_barplot.pdf", onefile=TRUE,width = 5,height = 4)
g <-ggplot(ih_gene_hits, aes(x=factor(gene), y=freq,fill=dataset)) +geom_bar(stat="identity",width=0.8)
g <- g+xlab("Gene")
g <- g+scale_x_discrete(limits = as.character(ih_gene_hits_dcast$gene))+ coord_flip()
g <- g+scale_fill_manual(values=c("#DB8DB2","#59ADD0"))
g <- g+theme(axis.text.x=element_text(size=10,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA)
)
print(g)
#ggMarginal(g, type = "histogram", fill="transparent")
dev.off()


##################################
# WGS sub and indel burden distribution of two 
##################################

hipsci_muts_burden <- read.table("hipsci_ips_summary_muts.txt", sep = "\t", header = T, as.is = T)
hipsci_muts_burden$dataset <- "hipsci"
insignia_muts_burden <- read.table("insignia_muts_summary_subclone.txt", sep = "\t", header = T, as.is = T)
insignia_muts_burden <- insignia_muts_burden[,c("Samples","sub_num","disub_num","indel_num")]
names(insignia_muts_burden) <- c("ips","a_sub_num","b_disub_num","c_indel_num")
insignia_muts_burden$dataset <- "insignia"
both_muts <- rbind(hipsci_muts_burden,insignia_muts_burden)
both_muts_melt <- melt(both_muts,c("ips","dataset"))
names(both_muts_melt) <- c("ips","dataset","type","Freq")
#both_muts_melt$freq_log <- log10(both_muts_melt)
pdf(file="ih_muts_burden.pdf", onefile=TRUE,width = 8,height = 4)
g <-ggplot(both_muts_melt, aes(x=dataset, y=Freq,colour=dataset)) + geom_quasirandom()#+geom_boxplot(outlier.size=0.5)
g <- g+  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom="pointrange", color = "black")
#g <- g+scale_x_discrete(limits = as.character(ih_gene_hits_dcast$gene))
g <- g+scale_colour_manual(values=c("#59ADD0","#DB8DB2"))#+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
g <- g+theme(axis.text.x=element_text(size=10,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA)
)
g <- g+facet_wrap(~type,ncol=3,scales = "free")
print(g)
dev.off()

##################################
# WGS sub and indel burden distribution of samples with drivers 
##################################

hipsci_muts_burden <- read.table("hipsci_ips_summary_muts.txt", sep = "\t", header = T, as.is = T)
hipsci_muts_burden$dataset <- "hipsci"
insignia_muts_burden <- read.table("insignia_muts_summary_subclone.txt", sep = "\t", header = T, as.is = T)
insignia_muts_burden <- insignia_muts_burden[,c("Samples","sub_num","disub_num","indel_num")]
names(insignia_muts_burden) <- c("ips","a_sub_num","b_disub_num","c_indel_num")
insignia_muts_burden$dataset <- "insignia"
both_muts <- rbind(hipsci_muts_burden,insignia_muts_burden)

insignia_muts_2 <- read.table("insignia_muts_g2.txt", sep = "\t", header = T, as.is = T)
hipsci_muts_2 <- read.table("hipsci_driver_withSNP_freq_g2.txt", sep = "\t", header = T, as.is = T)

sample_bcor <- c(insignia_muts_2[insignia_muts_2$gene=="BCOR","Sample"], hipsci_muts_2[hipsci_muts_2$gene=="BCOR","ips"])

both_muts_bcor <- both_muts[both_muts$ips %in%sample_bcor, ]
both_muts_bcor_melt <- melt(both_muts_bcor,c("ips","dataset"))
names(both_muts_bcor_melt) <- c("ips","dataset","type","Freq")

pdf(file="ih_muts_burden_BCOR.pdf", onefile=TRUE,width = 8,height = 4)
g <-ggplot(both_muts_bcor_melt, aes(x=dataset, y=Freq,colour=dataset)) + geom_quasirandom()#+geom_boxplot(outlier.size=0.5)
g <- g+  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom="pointrange", color = "black")
#g <- g+scale_x_discrete(limits = as.character(ih_gene_hits_dcast$gene))
g <- g+scale_colour_manual(values=c("#DB8DB2","#59ADD0"))#+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
g <- g+theme(axis.text.x=element_text(size=10,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA)
)
g <- g+facet_wrap(~type,ncol=3,scales = "free")
print(g)
dev.off()

###
both_muts2 <- both_muts
names(both_muts2) <- c("ips","csub","bdisub","aindel","dataset")
both_muts2 <- both_muts2[,c("ips","aindel","bdisub","csub","dataset")]
both_muts2$bcor <- 0.7
both_muts2[both_muts2$ips %in%sample_bcor, "bcor"] <- 1
both_muts2$total <- rowSums(both_muts2[,2:4])
both_muts2_order <- both_muts2[order(both_muts2$dataset, both_muts2$total, decreasing = T),]

both_muts2_melt <- melt(both_muts2,c("ips","dataset","bcor","total"))
names(both_muts2_melt) <- c("ips","dataset","bcor","total","type","Freq")
pdf(file="ih_muts_burden_BCOR_barplot.pdf", onefile=TRUE,width = 12,height = 2)
g <-ggplot(both_muts2_melt, aes(x=ips, y=Freq,fill=factor(type),alpha=bcor)) +geom_bar(stat="identity",width=0.8)
g <- g+scale_x_discrete(limits = as.character(both_muts2_order$ips))
g <- g+scale_fill_manual(values=c("#CA9F92","#F9CD97","#E3D9B0"))+scale_alpha(range = c(0.2, 1))
g <- g+theme(axis.text.x=element_blank(),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA)
)
print(g)
dev.off()

pdf(file="ih_muts_burden_BCOR_barplot_black.pdf", onefile=TRUE,width = 12,height = 2)
g <-ggplot(both_muts2, aes(x=ips, y=total,fill=dataset, alpha=bcor)) +geom_bar(stat="identity",width=0.8)
g <- g+scale_x_discrete(limits = as.character(both_muts2_order$ips))
g <- g+scale_fill_manual(values=c("#59ADD0","#DB8DB2"))+scale_alpha(range = c(0.2, 1))
g <- g+theme(axis.text.x=element_blank(),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA)
)
print(g)
dev.off()
