source("/rds/project/sn206/rds-sn206-nik-zainal/users/xz388/cancer_archive04/b_1176/15_subs_20180907/00_data/sub_common.R")


# dNdScv analysis of the HIPSCI data

setwd("/rds/project/rds-acvgC5qEG2w/users/xz388/cancer_archive04/g_1133/54_foad_78donors/dNdScv/")
library(dndscv)

## 1. INSIGNIA (blood) data
mutations = read.table("../cave_pindel_Insignia_1000gFiltered_muts_donor78_parent.txt", header=1, sep="\t", stringsAsFactors=F) # Loading the mutation file
mutations = mutations[ ,c(1,4,5,6,7)] # 5-column mutation format
mutations$Donor = sapply(strsplit(mutations$Sample, split = "_"), function(x) x[[1]]) # Extracting the DonorID to work with unique mutations per donor
mutations = unique(mutations[ ,c(6,2,3,4,5)]) # Unique mutations per donor (to avoid double counting)

# Default dNdScv
dndsout = dndscv(mutations, outmats = T) # Running dNdScv with default parameters
#dndsout = dndscv(mutations, outmats = T,cv=NULL) # Running dNdScv with default parameters; without covariates
print(dndsout$globaldnds) # The global dN/dS ratios suggest considerable germline SNP contamination (can they be filtered out with a SNP database?)
#name       mle     cilow   cihigh
#wmis wmis 0.9860407 0.8644260 1.124765
#wnon wnon 0.9991679 0.7305408 1.366572
#wspl wspl 0.9217558 0.6532505 1.300625
#wtru wtru 0.9632501 0.7546569 1.229500
#wall wall 0.9843728 0.8642538 1.121187

# Adding RHT q-values
known_cancergenes = read.table("apriori_list_known_cancer_genes.txt", header=0, sep="\t", stringsAsFactors = F)[,1]
knowndrivers = dndsout$sel_cv$gene_name %in% known_cancergenes # Index of the known cancer genes in the sel_cv output table
dndsout$sel_cv$qglobal_RHT[knowndrivers] = p.adjust(dndsout$sel_cv$pglobal_cv[knowndrivers], method="BH")
selcv_insignia = dndsout$sel_cv[which(dndsout$sel_cv$qglobal_cv<0.5 | dndsout$sel_cv$qglobal_RHT<0.5), ]
print(selcv_insignia)
write.table(selcv_insignia, file="dNdScv_Insignia_blood.txt", col.names=T, row.names=F, sep="\t", quote=F)

# Attempting to run RHT on known hotspots: No known substitution hotspots observed in the data
data("knownhotspots_hg19", package = "dndscv")
hotspots_siteRHT = sitednds(dndsout, site_list = known_hotspots, method = "LNP", min_recurr = 1)
print(hotspots_siteRHT)


## 2. HIPSCI (fibroblasts) data

m1 = read.table("wgs.SNVs.530bams.2019-09-02.txt", header=0, sep="\t", stringsAsFactors=F)
m2 = read.table("wgs.indels.mpileup.530bams.2019-09-02.txt", header=0, sep="\t", stringsAsFactors=F)
m3 = read.table("wes.SNVs.637bams.chrX.2021-01-06.txt", header=0, sep="\t", stringsAsFactors=F)
m4 = read.table("wes.indels.mpileup.637bams.2019-09-02.txt", header=0, sep="\t", stringsAsFactors=F)
m5 = read.table("wes.SNVs.highcov-137bams.2019-09-02.txt", header=0, sep="\t", stringsAsFactors=F)
m6 = read.table("wes.SNVs.637bams.2019-09-02.txt", header=0, sep="\t", stringsAsFactors=F)

mutations = rbind(m1, m2, m3, m4, m5, m6) # Merging all substitution and indel calls available (WES and WGS)
mutations = mutations[,c(4,1,2,6,7,11,21)] # Selecting the columns of interest
colnames(mutations) = c("sampleID","chr","pos","ref","mut","ips_nALT","FILTER") # Adding column names
mutations = mutations[order(mutations$sampleID, mutations$chr, mutations$pos), ]
mutations = mutations[mutations$ips_nALT>0, ] # Filtering out mutations without read support in column 11 (ips:nALT)
#mutations = mutations[mutations$FILTER=="PASS", ] # Only PASS variants

# Collapsing duplicate mutations within a donor
mutations$donor = sapply(strsplit(mutations$sampleID, split = "_"), function(x) x[[1]]) # Extracting the donorID to work with unique mutations per donor
mutsfordnds = unique(mutations[ ,c("donor","chr","pos","ref","mut")])

# Default dNdScv
dndsout = dndscv(mutsfordnds, outmats = T) # Running dNdScv with default parameters
print(dndsout$globaldnds) # The global dN/dS ratios are OK (they suggest no considerable SNP contamination)

# Adding RHT q-values
known_cancergenes = read.table("apriori_list_known_cancer_genes.txt", header=0, sep="\t", stringsAsFactors = F)[,1]
knowndrivers = dndsout$sel_cv$gene_name %in% known_cancergenes # Index of the known cancer genes in the sel_cv output table
dndsout$sel_cv$qglobal_RHT[knowndrivers] = p.adjust(dndsout$sel_cv$pglobal_cv[knowndrivers], method="BH")
selcv_hipsci = dndsout$sel_cv[which(dndsout$sel_cv$qglobal_cv<0.5 | dndsout$sel_cv$qglobal_RHT<0.5), ]
write.table(selcv_hipsci, file="dNdScv_HipSci_fibroblasts.txt", col.names=T, row.names=F, sep="\t", quote=F)
print(selcv_hipsci)

# Attempting to run RHT on known hotspots: No known substitution hotspots observed in the data
data("knownhotspots_hg19", package = "dndscv")
hotspots_siteRHT = sitednds(dndsout, site_list = known_hotspots, method = "LNP", min_recurr = 1)
print(hotspots_siteRHT)



##############
#  HipSci 
##############
hipsci_disub <- read.table("hipsci_disubs_FATHMM_uniq_pathogenic_petr_final.txt", sep="\t", header = T, as.is = T)
hdisub_s <- hipsci_disub[hipsci_disub$ips_nALT>0,c("ips","Fibro","Fibro_nREF","Fibro_nALT","ips_nREF","ips_nALT","gene","mutationAA","FATHMM","FATHMM_score")] # 58
hdisub_s$type <- "disub"

hipsci_sub <- read.table("hipsci_snvs_FATHMM_uniq_pathogenic_petr_final.txt", sep="\t", header = T, as.is = T)
hsub_s <- hipsci_sub[hipsci_sub$ips_nALT>0,c("ips","Fibro","Fibro_nREF","Fibro_nALT","ips_nREF","ips_nALT","gene","mutationAA","FATHMM","FATHMM_score")] # 247
hsub_s$type <- "sub"

hipsci_indel <- read.table("hipsci_indel_cancergenes_info_petr_frameshift_final.validated.txt", sep="\t", header = T, as.is = T)
hindel_s <- hipsci_indel[hipsci_indel$Validate%in%c("tp","exp"),c("ips","Fibro","Fibro_nREF","Fibro_nALT","ips_nREF","ips_nALT","gene")]  # 45
hindel_s$mutationAA <- "frameshift"
hindel_s$FATHMM <- "PATHOGENIC"
hindel_s$FATHMM_score <- 1
hindel_s$type <- "indel"

h_muts <- rbind(hdisub_s,hsub_s)
h_muts <- rbind(h_muts,hindel_s)
write.table(h_muts,"hipsci_driver_withSNP.txt",sep = "\t", col.names = T, row.names = F, quote = F)


h_muts <- read.table("hipsci_driver_withSNP.txt", sep="\t", header = T, as.is = T)

# BCOR iPSCs
bcor_hipsci <- h_muts[h_muts$gene=="BCOR",]
bcor_hipsci <- bcor_hipsci[order(bcor_hipsci$type,bcor_hipsci$Fibro_nREF),]

# plot fibroblast reads
pdf(file="hipsci_bcor_fibro_reads.pdf", onefile=TRUE,width = 5,height = 4)
g <-ggplot(bcor_hipsci, aes(y=Fibro_nREF, x=factor(ips))) + geom_bar(stat="identity", width = .6)
g <- g+scale_x_discrete(limits = as.character(bcor_hipsci$ips))+coord_flip()
g <- g+scale_y_continuous(breaks = seq(0,300,50))
#g <- g+scale_y_discrete(limits = as.character(h_muts_2_freq$gene))
#g <- g+scale_fill_gradient2(high="#59ADD0", low="white",limits=c(0, 1))
g <- g+theme(axis.text.x=element_text(size=8,angle=90, hjust=0.9,vjust=0.9,colour="black"),
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
#ggMarginal(g, type = "histogram", fill="transparent")
dev.off()

# plot ips reads
bcor_hipsci_trans <- melt(bcor_hipsci[,c("ips","ips_nREF","ips_nALT","type")],c("ips","type"))

pdf(file="hipsci_bcor_ips_reads.pdf", onefile=TRUE,width = 5,height = 4)
g <-ggplot(bcor_hipsci_trans, aes(y=value, x=factor(ips),fill=variable)) + geom_bar(stat="identity", width = .6,position = position_stack(reverse = TRUE))
g <- g+scale_x_discrete(limits = as.character(bcor_hipsci$ips))+coord_flip()
g <- g+scale_y_continuous(breaks = seq(0,300,50))
#g <- g+scale_y_discrete(limits = as.character(h_muts_2_freq$gene))
#g <- g+scale_fill_gradient2(high="#59ADD0", low="white",limits=c(0, 1))
g <- g+theme(axis.text.x=element_text(size=8,angle=90, hjust=0.9,vjust=0.9,colour="black"),
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
#ggMarginal(g, type = "histogram", fill="transparent")
dev.off()



##############
#  Insignia/Hipsci mutant mutation burden
##############
ips_bcorinfo <- read.table("sample_bcorinfo_78donors.txt", sep="\t", header = T, as.is = T)

##################################
# WGS sub and indel burden distribution of samples with drivers 
##################################

hipsci_muts_burden <- read.table("hipsci_ips_summary_muts.txt", sep = "\t", header = T, as.is = T)
hipsci_muts_burden$dataset <- "hipsci"
hipsci_muts_burden$mut_num <- hipsci_muts_burden$a_sub_num+hipsci_muts_burden$b_disub_num+hipsci_muts_burden$c_indel_num
insignia_muts_burden <- read.table("Insignia_BhiPSc_mut_burden.txt", sep = "\t", header = T, as.is = T)
names(insignia_muts_burden) <- c("ips","mut_num")
insignia_muts_burden$dataset <- "insignia"
both_muts <- rbind(hipsci_muts_burden[,c("ips","mut_num","dataset")],insignia_muts_burden)

insignia_muts_2 <- read.table("ips_subclones_donor78_bcormut.txt", sep = "\t", header = T, as.is = T)
hipsci_muts_2 <- read.table("hipsci_driver_withSNP_freq_g2.txt", sep = "\t", header = T, as.is = T)

sample_bcor <- c(insignia_muts_2[insignia_muts_2$BCOR_staus_ips=="BCOR_MT","ips"], hipsci_muts_2[hipsci_muts_2$gene=="BCOR","ips"])

both_muts_bcor <- both_muts[both_muts$ips %in%sample_bcor, ]
both_muts_bcor_melt <- melt(both_muts_bcor,c("ips","dataset"))
names(both_muts_bcor_melt) <- c("ips","dataset","type","Freq")

###
both_muts2 <- both_muts
both_muts2$bcor <- 0.7
both_muts2[both_muts2$ips %in%sample_bcor, "bcor"] <- 1
both_muts2_order <- both_muts2[order(both_muts2$dataset, both_muts2$mut_num, decreasing = T),]

pdf(file="ih_muts_burden_BCOR_barplot.pdf", onefile=TRUE,width = 6,height = 2)
g <-ggplot(both_muts2_order, aes(x=ips, y=mut_num,fill=factor(dataset),alpha=bcor)) +geom_bar(stat="identity",width=0.8)
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

### pvalue between bcor wt and mt


pdf(file="ih_muts_burden_BCOR_quasirandom.pdf", onefile=TRUE,width = 8,height = 4)
g <-ggplot(both_muts2_order, aes(x=bcor, y=mut_num,fill=factor(dataset),alpha=bcor)) + geom_quasirandom()#+geom_boxplot(outlier.size=0.5)
#g <- g+  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom="pointrange", color = "black")
#g <- g+scale_x_discrete(limits = as.character(ih_gene_hits_dcast$gene))
g <- g+scale_fill_manual(values=c("#DB8DB2","#59ADD0"))+scale_alpha(range = c(0.2, 1))
g <- g+theme(axis.text.x=element_text(size=10,colour="black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),                                                               
             axis.title.y = element_text(size=15),                                                               
             plot.title = element_text(size=10),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA)
)
g <- g+facet_wrap(~dataset,ncol=2)
print(g)
dev.off()

#p=0.5504784
p <- wilcox.test(both_muts2_order[both_muts2_order$dataset=="insignia"&both_muts2_order$bcor==1,"mut_num"], both_muts2_order[both_muts2_order$dataset=="insignia"&both_muts2_order$bcor==0.7,"mut_num"])$p.value
#p=0.77
p <- wilcox.test(both_muts2_order[both_muts2_order$dataset=="hipsci"&both_muts2_order$bcor==1,"mut_num"], both_muts2_order[both_muts2_order$dataset=="hipsci"&both_muts2_order$bcor==0.7,"mut_num"])$p.value



