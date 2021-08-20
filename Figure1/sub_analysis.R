source("/rds/project/sn206/rds-sn206-nik-zainal/users/xz388/cancer_archive04/b_1176/15_subs_20180907/00_data/sub_common.R")
source("/rds/project/sn206/rds-sn206-nik-zainal/users/xz388/cancer_archive04/b_1176/29_indels_fullscale_final/indel_common.R")

subs <- read.table("./cave_asmd.tab",sep = "\t",header = T, as.is = T)
###############################################
#
#    VAF for each subclone
#
###############################################
filename <- paste0("sub_VAF_distribution", ".pdf")
pdf(file=filename, onefile=TRUE,width=10,height=8)
p <- ggplot(data=subs, aes(x=PM.Tum))+ geom_histogram(binwidth=0.05)
p <- p+theme(axis.text.x=element_text(size=10,colour = "black"),
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
p <- p+facet_wrap(~Sample,ncol=3, scales="free")
print(p)
dev.off()

# ASMD
filename <- paste0("sub_ASMD_distribution", ".pdf")
pdf(file=filename, onefile=TRUE,width=10,height=8)
p <- ggplot(data=subs, aes(x=ASMD))+ geom_histogram(binwidth=0.05)
p <- p+theme(axis.text.x=element_text(size=10,colour = "black"),
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
p <- p+facet_wrap(~Sample,ncol=3, scales="free")
print(p)
dev.off()


subs <- subs[subs$ASMD>=90,]
#######################################
# Check mutation number in each sample
#######################################
SampleMutNum(subs,6,8,"muts_num")

########################################################################
# Plot mutational profiles for each sample using ALL mutations
########################################################################
mut_catalogue <- GenCatalogue(subs,"Sample")
write.table(mut_catalogue,"all_muts_catalogue.txt",sep = "\t",col.names = T, row.names = F, quote = F)
plotCountbasis(mut_catalogue[,-2],2,12,15,"all_muts_catalogue.pdf")


###############################################
#
#    Shared mutation
#
###############################################
subs2 <- subs
subs2$Treatment <- "S7"
subs2[subs2$Sample %in% c("S2_RE19_P2","S2_RE5_P1","S2_SF2_P2","S2_SF3_P2"),"Treatment"] <- "S2"

ShareMuts_summary(subs2,"Treatment","Sample")
ShareMuts_UpSetR(subs2,10,15,"ShareMuts_UpSetR")

############################################
# Remove shared mutations
############################################
denovo_muts <- Remove_sharedmuts(subs2,"Treatment")
ShareMuts_UpSetR(denovo_muts,10,10,"ShareMuts_UpSetR_denovo")
write.table(denovo_muts,"denovo_muts.txt", sep = "\t", col.names = T, row.names = F, quote = F)

########################################################################
# Plot mutational profiles for each sample using denovo mutations
########################################################################
mut_catalogue <- GenCatalogue(denovo_muts,"Sample")
write.table(mut_catalogue,"denovo_muts_catalogue.txt",sep = "\t",col.names = T, row.names = F, quote = F)
plotCountbasis(mut_catalogue[,-2],2,12,15,"denovo_muts_catalogue.pdf")

#######################################
# Check de novo mutation number in each sample
#######################################
SampleMutNum(denovo_muts,6,8,"denovo_muts_num")



########################################################################
# Plot mutational profiles for shared mutations
########################################################################
shared_muts <- Sharedmuts(subs2,"Treatment")
mut_catalogue <- GenCatalogue(shared_muts,"Sample")
write.table(mut_catalogue,"shared_muts_catalogue.txt",sep = "\t",col.names = T, row.names = F, quote = F)
plotCountbasis(mut_catalogue[,-2],2,12,15,"shared_muts_catalogue.pdf")


###############################################
#
#   2131 + 2195
#
###############################################
subs2131 <- read.table("/nfs/cancer_archive04/xz3/o_foad/11_subs/cave_asmd.tab",sep = "\t",header = T, as.is = T)
subs2131 <- subs2131[subs2131$ASMD>=90,]
subs2131 <- subs2131[subs2131$PM.Tum>0.1,]
subs2131$mut_ID <- paste0(subs2131$Sample,"_",subs2131$Chrom,"_",subs2131$Pos)
subs2131 <- subs2131[subs2131$Sample %in% c("S2_RE19_P2","S2_RE5_P1","S2_SF2_P2","S2_SF3_P2","S7_RE11_P2","S7_RE14_P2","S7_RE17_P2","S7_RE2_P2"), ]


subs2195 <- read.table("/nfs/cancer_archive04/xz3/o_foad/2195/01_subs/cave_asmd.tab",sep = "\t",header = T, as.is = T)
sample_info <- read.table("/nfs/cancer_archive04/xz3/o_foad/2195/01_subs/sample_info.txt",sep = "\t",header = T, as.is = T)
subs2195 <- merge(subs2195, sample_info,by="Sample")
subs2195$Sample <- subs2195$Sample2
subs2195 <- subs2195[,-dim(subs2195)[2]]
subs2195 <- subs2195[subs2195$ASMD>=90,]
subs2195 <- subs2195[subs2195$PM.Tum>0.1,]
subs2195$mut_ID <- paste0(subs2195$Sample,"_",subs2195$Chrom,"_",subs2195$Pos)

subs <- rbind(subs2131,subs2195)
write.table(subs,"subs_2131_2195.txt",sep = "\t", col.names = T, row.names = F, quote = F)
###############################################
# Find double substitutions
###############################################
subs <- read.table("subs_2131_2195.txt", sep = "\t", header = T, as.is = T)
doublesubs <- FindDinucleotides(subs)

doublesubs_list <- c(paste0(doublesubs$Sample,"_",doublesubs$Chrom,"_",doublesubs$Pos))

dbsubs <- subs[subs$mut_ID %in% doublesubs_list,]

write.table(dbsubs,"doublesubs_2131_2195.txt",sep = "\t", col.names = T, row.names = F, quote = F)

# remove double substitutions from subs

subs <- subs[!subs$mut_ID %in% doublesubs_list,]

doublesub_num <- data.frame(table(doublesubs$Sample))
sub_num <- data.frame(table(subs$Sample))
names(doublesub_num) <- c("Sample","doublesubs")
names(sub_num) <- c("Sample","subs")
ds_subs <- merge(doublesub_num,sub_num,by="Sample")
ds_subs$tissue <- "skin"
ds_subs[ds_subs$Sample%in%c("S2_RE19_P2","S2_RE5_P1","S7_RE11_P2","S7_RE14_P2","S7_RE17_P2","S7_RE2_P2"),"tissue"] <- "blood"
write.table(ds_subs,"skin_blood_subs_burden.txt",sep = "\t", col.names = T, row.names = F, quote = F)

ds_subs_melt <- melt(ds_subs, c("Sample","tissue"))
names(ds_subs_melt) <- c("Sample","tissue","type","freq")
ds_subs_melt <- ds_subs_melt[order(ds_subs_melt$tissue,ds_subs_melt$freq,decreasing = T),]
pdf(file="skin_blood_subs_burden.pdf", onefile=TRUE,height=3,width=5, useDingbats=FALSE)
p <- ggplot(data=ds_subs_melt, aes(x=Sample, y=freq, fill=type))+ geom_bar(stat="identity",width=.8)+xlab("Sample")+ylab("Count")
#p <- p+scale_x_discrete(limits = ds_subs$Sample,labels = ds_subs$Sample)
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=8,colour = "black",hjust = 1),
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
p <- p+facet_grid(.~tissue, scales="free", space="free")
print(p)
dev.off()


# 96 profile of double subs
doublesubs <- read.table("doublesubs_2131_2195.txt",sep = "\t", header = T, as.is = T)
sub_catalouge <- GenCatalogue(doublesubs,"Sample")
plotCountbasis(sub_catalouge[,-2],2,12,15,"doublesub_catalogue2131_2195.pdf")



#######################################
#
# Six mutation types
#  
########################################
subs <- read.table("subs_2131_2195.txt", sep = "\t", header = T, as.is = T)
doublesubs <- FindDinucleotides(subs)
doublesubs_list <- c(paste0(doublesubs$Sample,"_",doublesubs$Chrom,"_",doublesubs$Pos), paste0(doublesubs$Sample,"_",doublesubs$Chrom,"_",doublesubs$Pos_neighbor))

# remove double substitutions from subs

subs <- subs[!subs$mut_ID %in% doublesubs_list,]


sub_catalouge <- Gen6Catalogue(subs,"Sample")
subs_melt <- melt(sub_catalouge, c("Mutation"))
names(subs_melt) <- c("Mutation","Sample","freq")
subs_melt$tissue <- "skin"
subs_melt[subs_melt$Sample%in%c("S2_RE19_P2","S2_RE5_P1","S7_RE11_P2","S7_RE14_P2","S7_RE17_P2","S7_RE2_P2"),"tissue"] <- "blood"
write.table(subs_melt,"skin_blood_subs_6Mutation.txt",sep = "\t", col.names = T, row.names = F, quote = F)
mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")

pdf(file="skin_blood_subs_6Mutation.pdf", onefile=TRUE,height=3,width=5, useDingbats=FALSE)
p <- ggplot(data=subs_melt, aes(x=Sample, y=freq, fill=Mutation))+ geom_bar(stat="identity",width=.8)+xlab("Sample")+ylab("Count")+scale_fill_manual(values=mypalette)
#p <- p+scale_x_discrete(limits = ds_subs$Sample,labels = ds_subs$Sample)
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=8,colour = "black",hjust = 1),
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
p <- p+facet_grid(.~tissue, scales="free", space="free")
print(p)
dev.off()





#######################################
#
# Mutational signatures
#  
########################################
# Fit Skin-spefic signatures
PancanSig <- read.table("/nfs/cancer_archive04/xz3/b_1176/70_MSI/00_data/Pancan_signatures_subs_final.txt",sep = "\t",header = T, as.is = T)
names(PancanSig)[1] <- "MutationType"
Skin_tissueSig <- PancanSig[,c(1,which(sub("_[^_]+$","",names(PancanSig))=="Skin"))]

subs <- read.table("subs_2131_2195.txt", sep = "\t", header = T, as.is = T)
doublesubs <- FindDinucleotides(subs)
doublesubs_list <- c(paste0(doublesubs$Sample,"_",doublesubs$Chrom,"_",doublesubs$Pos), paste0(doublesubs$Sample,"_",doublesubs$Chrom,"_",doublesubs$Pos_neighbor))

# remove double substitutions from subs

subs <- subs[!subs$mut_ID %in% doublesubs_list,]


sub_catalouge <- GenCatalogue(subs,"Sample")
plotCountbasis(sub_catalouge[,-2],2,12,15,"sub_catalogue2131_2195.pdf")

mut_sig <- merge(Skin_tissueSig,sub_catalouge[,-2],by="MutationType")
mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
sig_cat <- mut_sig[,2:dim(Skin_tissueSig)[2]]
mut_cat <- mut_sig[,(dim(Skin_tissueSig)[2]+1):dim(mut_sig)[2]]

a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
write.table(a$E_median_filtered,paste0("Skinsig_exposure_2131_2195",".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
Exposure <-a$E_median_filtered
Exposure <- read.table("Skinsig_exposure_2131_2195.txt", sep = "\t", header = T, as.is = T)
sample_exposure <- as.data.frame(Exposure)
sample_exposure$Sig <- rownames(sample_exposure)
sample_exposure_melt <- melt(sample_exposure,c("Sig"))
names(sample_exposure_melt) <- c("Sig","Sample","exposure")
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_A","Sig"] <- "Culture"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_D","Sig"] <- "UV_Skin_D"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_J","Sig"] <- "UV_Skin_J"
sample_exposure_melt$tissue <- "skin"
sample_exposure_melt[sample_exposure_melt$Sample%in%c("S2_RE19_P2","S2_RE5_P1","S7_RE11_P2","S7_RE14_P2","S7_RE17_P2","S7_RE2_P2"),"tissue"] <- "blood"
cbbPalette <- c("#000000", "#999999","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#FF66FF", "#CC79A7")
pdf(file="Skinsig_exposure_2131_2195_SigExposur.pdf", onefile=TRUE,height=3,width=5, useDingbats=FALSE)
p <- ggplot(sample_exposure_melt,aes(x=Sample,y=exposure,fill=Sig))+geom_bar(stat="identity", width=.8)+scale_fill_manual(values=cbbPalette)
#p <- p+scale_x_discrete(limits = as.character(summary_subs[,"ips"]))+xlab("ips")
#p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=8,colour = "black",hjust = 1),
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
p <- p+facet_grid(.~tissue, scales="free", space="free")

print(p)
dev.off()



#######################################
#
#  COSMIC SIGNATURE
#  Compare with COSMIC signatures
########################################
cosmic_sig <- read.table("/nfs/cancer_archive04/xz3/a_1242/00_common/signatures_probabilities.txt",sep = "\t", header = T, as.is = T)
cosmic_sig <- cosmic_sig[,3:33]
colnames(cosmic_sig)[1] <- "MutationType"
mut_catalogue <- read.table("denovo_muts_catalogue.txt", sep = "\t", header = T, as.is = T)
mut_catalogue <- mut_catalogue[,-2] 
mut_catalogue[,-which(colnames(mut_catalogue)=="MutationType")] <- mut_catalogue[,-which(colnames(mut_catalogue)=="MutationType")]/colSums(mut_catalogue[,-which(colnames(mut_catalogue)=="MutationType")])[col(mut_catalogue[,-which(colnames(mut_catalogue)=="MutationType")])]

Cossimi_CompareSig(mut_catalogue,cosmic_sig,2,10,2,"COSMIC Signatures","Cossimi_Cosmic")

#######################################
#
#  New figure 1
#  
########################################
subs <- read.table("subs_2131_2195.txt", sep = "\t", header = T, as.is = T)
doublesubs <- FindDinucleotides(subs)

# double substitutions
doublesubs_list <- c(paste0(doublesubs$Sample,"_",doublesubs$Chrom,"_",doublesubs$Pos))
dbsubs <- subs[subs$mut_ID %in% doublesubs_list,]

# remove double substitutions from subs
subs <- subs[!subs$mut_ID %in% doublesubs_list,]

# A. mutation burden of double subs and single subs
doublesub_num <- data.frame(table(doublesubs$Sample))
sub_num <- data.frame(table(subs$Sample))
names(doublesub_num) <- c("Sample","doublesubs")
names(sub_num) <- c("Sample","subs")
ds_subs <- merge(doublesub_num,sub_num,by="Sample")
#ds_subs$tissue <- "skin"
#ds_subs[ds_subs$Sample%in%c("S2_RE19_P2","S2_RE5_P1","S7_RE11_P2","S7_RE14_P2","S7_RE17_P2","S7_RE2_P2"),"tissue"] <- "blood"
ds_subs$donor <- sub("\\_.*",'', ds_subs$Sample)
ds_subs[!ds_subs$donor %in% c("S2","S7"),"donor"] <- "dHPSI_Skin"
ds_subs[ds_subs$donor =="S7","donor"] <- "aS7_Blood"
ds_subs[ds_subs$donor =="S2","donor"] <- "cS2_Skin"
ds_subs[ds_subs$Sample %in%c("S2_RE19_P2","S2_RE5_P1"),"donor"] <- "bS2_Blood"
write.table(ds_subs,"skin_blood_subs_burden.txt",sep = "\t", col.names = T, row.names = F, quote = F)

ds_subs_melt <- melt(ds_subs, c("Sample","donor"))
names(ds_subs_melt) <- c("Sample","donor","type","freq")
ds_subs <- ds_subs[order(ds_subs$donor,ds_subs$subs,decreasing = F),]
pdf(file="skin_blood_subs_burden_donor.pdf", onefile=TRUE,height=3,width=5, useDingbats=FALSE)
p <- ggplot(data=ds_subs_melt, aes(x=Sample, y=freq, fill=type))+ geom_bar(stat="identity",width=.8)+xlab("Sample")+ylab("Count")
#p <- p+scale_x_discrete(limits = ds_subs$Sample,labels = ds_subs$Sample)
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=8,colour = "black",hjust = 1),
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
p <- p+facet_grid(.~donor, scales="free", space="free")
print(p)
dev.off()

# B. 6 mutation types of single subs
sub_catalouge <- Gen6Catalogue(subs,"Sample")
subs_melt <- melt(sub_catalouge, c("Mutation"))
names(subs_melt) <- c("Mutation","Sample","freq")
subs_melt$donor <- sub("\\_.*",'', subs_melt$Sample)
subs_melt[!subs_melt$donor %in% c("S2","S7"),"donor"] <- "dHPSI_Skin"
subs_melt[subs_melt$donor =="S7","donor"] <- "aS7_Blood"
subs_melt[subs_melt$donor =="S2","donor"] <- "cS2_Skin"
subs_melt[subs_melt$Sample %in%c("S2_RE19_P2","S2_RE5_P1"),"donor"] <- "bS2_Blood"
write.table(subs_melt,"skin_blood_subs_6Mutation.txt",sep = "\t", col.names = T, row.names = F, quote = F)
mypalette <- c("deepskyblue","black", "red", "grey", "seagreen2", "lightpink2")

pdf(file="skin_blood_subs_6Mutation_donor.pdf", onefile=TRUE,height=3,width=5, useDingbats=FALSE)
p <- ggplot(data=subs_melt, aes(x=Sample, y=freq, fill=Mutation))+ geom_bar(stat="identity",width=.8)+xlab("Sample")+ylab("Count")+scale_fill_manual(values=mypalette)
#p <- p+scale_x_discrete(limits = ds_subs$Sample,labels = ds_subs$Sample)
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=8,colour = "black",hjust = 1),
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
p <- p+facet_grid(.~donor, scales="free", space="free")
print(p)
dev.off()

# C. Skin signature exposure
Exposure <- read.table("Skinsig_exposure_2131_2195.txt", sep = "\t", header = T, as.is = T)
sample_exposure <- as.data.frame(Exposure)
sample_exposure$Sig <- rownames(sample_exposure)
sample_exposure_melt <- melt(sample_exposure,c("Sig"))
names(sample_exposure_melt) <- c("Sig","Sample","exposure")
sample_exposure_melt[,"Sig1"] <- "tOther signatures"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_A","Sig1"] <- "Signature 18"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_D","Sig1"] <- "UV_Skin_D"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_J","Sig1"] <- "UV_Skin_J"

sample_exposure_melt$donor <- sub("\\_.*",'', sample_exposure_melt$Sample)
sample_exposure_melt[!sample_exposure_melt$donor %in% c("S2","S7"),"donor"] <- "dHPSI_Skin"
sample_exposure_melt[sample_exposure_melt$donor =="S7","donor"] <- "aS7_Blood"
sample_exposure_melt[sample_exposure_melt$donor =="S2","donor"] <- "cS2_Skin"
sample_exposure_melt[sample_exposure_melt$Sample %in%c("S2_RE19_P2","S2_RE5_P1"),"donor"] <- "bS2_Blood"

cbbPalette <- c("#000000","#F0E442", "#CC79A7","#CC79A7")
pdf(file="Skinsig_exposure_2131_2195_SigExposur_donor.pdf", onefile=TRUE,height=3,width=5, useDingbats=FALSE)
p <- ggplot(sample_exposure_melt,aes(x=Sample,y=exposure,fill=Sig1))+geom_bar(stat="identity", width=.8)+scale_fill_manual(values=cbbPalette)
#p <- p+scale_x_discrete(limits = as.character(summary_subs[,"ips"]))+xlab("ips")
#p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=8,colour = "black",hjust = 1),
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
p <- p+facet_grid(.~donor, scales="free", space="free")

print(p)
dev.off()

# D. Indel types
indels <- read.table("indels_2131_2195.txt", sep = "\t", header = T, as.is = T)

a <- data.frame(table(indels$Sample,indels$classification))
names(a) <- c("Sample","type","Freq")
a$donor <- sub("\\_.*",'', a$Sample)
a[!a$donor %in% c("S2","S7"),"donor"] <- "dHPSI_Skin"
a[a$donor =="S7","donor"] <- "aS7_Blood"
a[a$donor =="S2","donor"] <- "cS2_Skin"
a[a$Sample %in%c("S2_RE19_P2","S2_RE5_P1"),"donor"] <- "bS2_Blood"

indel_mypalette_fill <- c("orchid","pink", "lightgreen","skyblue","skyblue","orange","orange",
                          "purple","deeppink","darkgreen","blue","blue","tomato","tomato",
                          "grey")
indel_positions <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                     "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                     "Complex") 

indel_mypalette_fill <- c("orchid","pink","lightgreen","skyblue","orange","purple","deeppink","darkgreen","blue","tomato","grey")
#cbbPalette <- c("#000000", "#999999","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#FF66FF", "#CC79A7")
pdf(file="indels_2131_2195_catalogue_donor.pdf", onefile=TRUE,height=3,width=5, useDingbats=FALSE)
p <- ggplot(a,aes(x=Sample,y=Freq,fill=type))+geom_bar(stat="identity", width=.8)+scale_fill_manual(values=indel_mypalette_fill)
#p <- p+scale_x_discrete(limits = as.character(summary_subs[,"ips"]))+xlab("ips")
#p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=8,colour = "black",hjust = 1),
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
p <- p+facet_grid(.~donor, scales="free", space="free")

print(p)
dev.off()

# E. Rg number
#####################
rgs_2131 <- read.table("/rds/project/sn206/rds-sn206-nik-zainal/users/xz388/cancer_archive04/o_foad/30_rgs/Rgs_BrassI_2131.txt", sep = "\t", header = T, as.is = T)
names(rgs_2131)[11] <- "Sample"
rgs_2195 <- read.table("/rds/project/sn206/rds-sn206-nik-zainal/users/xz388/cancer_archive04/o_foad/2195/03_rgs/Rgs_BrassI.txt", sep = "\t", header = T, as.is = T)
info_2195 <- read.table("/rds/project/sn206/rds-sn206-nik-zainal/users/xz388/cancer_archive04/o_foad/2195/01_subs/sample_info.txt", sep = "\t", header = T, as.is = T)
names(info_2195) <- c("sample","Sample")
rgs_2195 <- merge(rgs_2195,info_2195, by="sample")
rgs_2195 <- rgs_2195[,-1]
rgs <- rbind(rgs_2131,rgs_2195)
rgs_brassII <- rgs[rgs$assembly_score !="_",]
rgs_brassII <- rgs_brassII[!grepl(",",rgs_brassII$Sample),]
rgs_brassII <- rgs_brassII[!grepl("S7RE14",rgs_brassII$Sample),]

write.table(rgs_brassII, "rgs_brassII_2131_2195.txt", sep = "\t", col.names = T, row.names = F, quote = F)
#########################

rgs_brassII <- read.table("rgs_brassII_2131_2195.txt", sep="\t", header = T, as.is = T)
a <- data.frame(table(rgs_brassII$Sample))
names(a) <- c("Sample","Freq")
a$donor <- sub("\\_.*",'', a$Sample)
a[!a$donor %in% c("S2","S7"),"donor"] <- "dHPSI_Skin"
a[a$donor =="S7","donor"] <- "aS7_Blood"
a[a$donor =="S2","donor"] <- "cS2_Skin"
a[a$Sample %in%c("S2_RE19_P2","S2_RE5_P1"),"donor"] <- "bS2_Blood"

pdf(file="rgs_2131_2195_donor.pdf", onefile=TRUE,height=3,width=5, useDingbats=FALSE)
p <- ggplot(a,aes(x=Sample,y=Freq))+geom_bar(stat="identity", width=.8,fill="#00909E")
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=8,colour = "black",hjust = 1),
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
p <- p+facet_grid(.~donor, scales="free", space="free")

print(p)
dev.off()


