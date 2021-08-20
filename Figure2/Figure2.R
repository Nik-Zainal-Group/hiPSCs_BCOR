source("/rds/project/sn206/rds-sn206-nik-zainal/users/xz388/cancer_archive04/b_1176/15_subs_20180907/00_data/sub_common.R")

#########################
# New figures 2
#########################
# Fig2B
# sub, disub, indel
ips_denovo_final <- read.table("ips_final.txt", sep = "\t", header = T, as.is = T)
summary_subs <- data.frame(table(ips_denovo_final$ips))
names(summary_subs) <- c("ips","Freq")
summary_subs <- summary_subs[order(summary_subs$Freq, decreasing = T),]
write.table(summary_subs,"ips_summary_subs.txt",sep = "\t", col.names = T, row.names = F, quote = F)

indels <- read.table("WGS_indels_final.txt", sep = "\t", header = T, as.is = T)
ips_indels <- indels[indels$VAF>0,] 
ips_denovo_final <- ips_indels[ips_indels$VAF>0.1,]

summary_subs <- data.frame(table(ips_denovo_final$ips))
names(summary_subs) <- c("ips","Freq")
summary_subs <- summary_subs[order(summary_subs$Freq, decreasing = T),]
write.table(summary_subs,"ips_summary_indels.txt",sep = "\t", col.names = T, row.names = F, quote = F)


disubs <- read.table("WGS_disubs_final.txt", sep = "\t", header = T, as.is = T)
ips_disubs <- disubs[disubs$VAF>0,] 
ips_denovo_final <- ips_disubs[ips_disubs$VAF>0.1,]

summary_subs <- data.frame(table(ips_denovo_final$ips))
names(summary_subs) <- c("ips","Freq")
summary_subs <- summary_subs[order(summary_subs$Freq, decreasing = T),]
write.table(summary_subs,"ips_summary_disubs.txt",sep = "\t", col.names = T, row.names = F, quote = F)


summary_subs <- read.table("ips_summary_subs.txt", sep = "\t", header = T, as.is = T)
names(summary_subs) <- c("ips","a_sub_num")
summary_indels <- read.table("ips_summary_indels.txt", sep = "\t", header = T, as.is = T)
names(summary_indels) <- c("ips","c_indel_num")
summary_disubs <- read.table("ips_summary_disubs.txt", sep = "\t", header = T, as.is = T)
names(summary_disubs) <- c("ips","b_disub_num")

summary_muts <- merge(summary_subs,summary_disubs,by="ips", all.x=T)
summary_muts <- merge(summary_muts,summary_indels, by="ips", all.x=T)
summary_muts[is.na(summary_muts)] <- 0
write.table(summary_muts,"ips_summary_muts.txt", sep = "\t", col.names = T, row.names = F, quote = F)

summary_muts_melt <- melt(summary_muts,c("ips"))

summary_muts_melt[summary_muts_melt$value==0,"value"] <- 1 # just for plotting on log scale
pdf(file="ips_muts_distribution_quassirandom.pdf", onefile=TRUE,height=3,width=5.5, useDingbats=FALSE)
#p <- ggplot(summary_muts_melt,aes(x=variable,y=value))+geom_jitter(aes(colour=variable),position=position_jitter(0.2))+  #+geom_dotplot(aes(fill=variable),binaxis = "y", stackdir = "center")+
stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom="pointrange", color = "black")
p <- ggplot(summary_muts_melt,aes(x=variable,y=value))+geom_quasirandom(aes(colour=variable))+  
  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom="pointrange", color = "black")
p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_text(size=10,colour = "black"),
             axis.text.y=element_text(size=10,colour = "black"),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
#p <- p+facet_grid(.~ Sig)
print(p)
dev.off()

# Figure2C
# Fit Skin-spefic signatures
PancanSig <- read.table("/nfs/cancer_archive04/xz3/b_1176/70_MSI/00_data/Pancan_signatures_subs_final.txt",sep = "\t",header = T, as.is = T)
names(PancanSig)[1] <- "MutationType"
Skin_tissueSig <- PancanSig[,c(1,which(sub("_[^_]+$","",names(PancanSig))=="Skin"))]

mut_catalogue <- GenCatalogue(ips_denovo_final,"ips")
write.table(mut_catalogue,"ips_catalogue_final.txt",sep = "\t",col.names = T, row.names = F, quote = F)

sub_catalouge <- read.table("ips_catalogue_final.txt", sep = "\t", header = T, as.is = T)
mut_sig <- merge(Skin_tissueSig,sub_catalouge[,-2],by="MutationType")
mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
sig_cat <- mut_sig[,2:dim(Skin_tissueSig)[2]]
mut_cat <- mut_sig[,(dim(Skin_tissueSig)[2]+1):dim(mut_sig)[2]]

a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
write.table(a$E_median_filtered,paste0("ips_exposure_","Skin",".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
Exposure <-a$E_median_filtered
sample_exposure <- as.data.frame(Exposure)
sample_exposure$Sig <- rownames(sample_exposure)
sample_exposure_melt <- melt(sample_exposure,c("Sig"))
names(sample_exposure_melt) <- c("Sig","ips","exposure")
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_A","Sig"] <- "Culture"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_D","Sig"] <- "UV_Skin_D"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_J","Sig"] <- "UV_Skin_J"
sample_exposure_melt$ips <- chartr(".", "-", sample_exposure_melt$ips)
summary_subs <- read.table("../ips_summary_subs.txt",sep = "\t", header = T, as.is = T)
#cbbPalette <- c("#000000", "#999999","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#FF66FF", "#CC79A7")
cbbPalette <- c("#000000", "#F0E442","#F0E442", "#F0E442", "#F0E442", "#F0E442", "#F0E442", "#F0E442","#CC79A7", "#CC79A7")

pdf(file="ips_sub_SigExposure2.pdf", onefile=TRUE,height=3,width=8, useDingbats=FALSE)
p <- ggplot(sample_exposure_melt,aes(x=ips,y=exposure,fill=Sig))+geom_bar(stat="identity", width=.8, position="fill")+scale_fill_manual(values=cbbPalette)
p <- p+scale_x_discrete(limits = as.character(summary_subs[,"ips"]))+xlab("ips")
#p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_blank(),
             axis.text.y=element_text(size=10,colour = "black"),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
print(p)
dev.off()

pdf(file="ips_sub_SigExposure.pdf", onefile=TRUE,height=3,width=8, useDingbats=FALSE)
p <- ggplot(sample_exposure_melt,aes(x=ips,y=exposure,fill=Sig))+geom_bar(stat="identity", width=.8)+scale_fill_manual(values=cbbPalette)
p <- p+scale_x_discrete(limits = as.character(summary_subs[,"ips"]))+xlab("ips")
#p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_blank(),
             axis.text.y=element_text(size=10,colour = "black"),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
print(p)
dev.off()



# Figure2D
# UV exposure vs double subs
Exposure_ips <- read.table("ips_exposure_Skin.txt", sep = "\t", header = T, as.is = T)
Exposure_ips <- rbind(Exposure_ips,colSums(Exposure_ips))
rownames(Exposure_ips)[11] <- "Skin_UV"
Exposure_ips_t <- as.data.frame(t(Exposure_ips))
Exposure_ips_t$ips <- gsub("\\.","-",rownames(Exposure_ips_t))

summary_disubs <- read.table("ips_summary_disubs.txt", sep = "\t", header = T, as.is = T)
names(summary_disubs) <- c("ips","b_disub_num")

summary_disubs_subexposure <- merge(summary_disubs, Exposure_ips_t[,c("ips", "Skin_UV")], by="ips", all.y=T)
summary_disubs_subexposure[is.na(summary_disubs_subexposure)] <- 0
res <- cor.test(summary_disubs_subexposure$b_disub_num, summary_disubs_subexposure$Skin_UV, method = "pearson") # cor: 0.9731135 
pdf(file="ips_disub_UVexposure_correlation.pdf", onefile=TRUE,height=3,width=4, useDingbats=FALSE)
p <- ggplot(summary_disubs_subexposure,aes(x=Skin_UV,y=b_disub_num))+geom_point()+geom_smooth(method="lm")
p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_text(size=10,colour = "black"),
             axis.text.y=element_text(size=10,colour = "black"),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
#p <- p+facet_grid(.~types)
print(p)
dev.off()


# Figure 2G
mut_burden_ips <- read.table("../ips_summary_subs.txt",sep = "\t", header = T, as.is = T)
mut_burden_ips$Fibro <- sub('.*\\-','',sub('\\_.*', '', mut_burden_ips$ips))

sex_age <- read.table("../ips_sex_age.txt", sep = "\t", header = T, as.is = T)
names(sex_age) <- c("ips", "Sex","Age")
sex_age_mutburden <- merge(sex_age, mut_burden_ips, by="ips")

# ips 
pdf(file="ips_sub_age_correlation.pdf", onefile=TRUE,height=3,width=5, useDingbats=FALSE)
p <- ggplot(sex_age_mutburden,aes(x=Age,y=Freq))+geom_violin(trim = FALSE)+geom_jitter(position=position_jitter(0.1))#+scale_color_manual(values=cbbPalette)
#p <- p+scale_x_discrete(limits = as.character(summary_subs[,"Fibro"]))+xlab("Fibro")
p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_text(size=8,colour = "black"),
             axis.text.y=element_text(size=10,colour = "black"),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
#p <- p+facet_grid(.~ Sig)
print(p)
dev.off()

# Figure 2H
pdf(file="ips_sub_sex_correlation.pdf", onefile=TRUE,height=3,width=3, useDingbats=FALSE)
p <- ggplot(sex_age_mutburden,aes(x=Sex,y=Freq))+geom_violin(trim = FALSE)+geom_jitter(position=position_jitter(0.1))#+scale_color_manual(values=cbbPalette)
#p <- p+scale_x_discrete(limits = as.character(summary_subs[,"Fibro"]))+xlab("Fibro")
p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_text(size=8,colour = "black"),
             axis.text.y=element_text(size=10,colour = "black"),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
#p <- p+facet_grid(.~ Sig)
print(p)
dev.off()

wilcox.test(sex_age_mutburden[sex_age_mutburden$Sex=="Female","Freq"], sex_age_mutburden[sex_age_mutburden$Sex=="Male","Freq"]) # p-value = 0.03114
wilcox.test(Freq~Sex, data=sex_age_mutburden) # p-value = 0.03114


###### Supplementary Figures
# Figure S3 
# Correlation between UV exposure vs  15 indel subtypes
Exposure_ips <- read.table("ips_exposure_Skin.txt", sep = "\t", header = T, as.is = T)
Exposure_ips_2 <- rbind(Exposure_ips[1,],colSums(Exposure_ips[c("Skin_D","Skin_J"),]))
Exposure_ips_2 <- rbind(Exposure_ips_2,colSums(Exposure_ips[c("Skin_B","Skin_C","Skin_E","Skin_F","Skin_G","Skin_H","Skin_I"),]))
rownames(Exposure_ips_2) <- c("sig18","sig7","sig.others")

subclones_catalogue <- read.table("ips_indel_catalogue15.txt",sep = "\t", header = T, as.is = T)
rownames(subclones_catalogue) <- subclones_catalogue$indelsubtype

SubSig_IndelSig_ips <- rbind(Exposure_ips_2, subclones_catalogue[,-1])

SubSig_IndelSig_ips_t <- t(SubSig_IndelSig_ips)
cor.test(SubSig_IndelSig_ips_t[,2], SubSig_IndelSig_ips_t[,10])

corr_matrix <- matrix(rep(0,45),nrow=3)
for(i in 1:3){
  a <- SubSig_IndelSig_ips[i,]
  
  for(j in 4:18){
    b <- SubSig_IndelSig_ips[j,]
    corr_matrix[i,j-3] <- cor.test(t(a), t(b), method = "pearson")$estimate
    
  }
}
corr_matrix <- data.frame(corr_matrix)
corr_matrix[is.na(corr_matrix)] <- 0
names(corr_matrix) <- rownames(SubSig_IndelSig_ips)[4:18]
rownames(corr_matrix) <- rownames(SubSig_IndelSig_ips)[1:3]
corr_matrix$Sub_sig <- rownames(corr_matrix)
corr_matrix_melt <- melt(corr_matrix,c("Sub_sig"))
names(corr_matrix_melt) <- c("Sub_sig","Indel_type","corr")

pvalue_matrix <- matrix(rep(0,45),nrow=3)
for(i in 1:3){
  a <- SubSig_IndelSig_ips[i,]
  
  for(j in 4:18){
    b <- SubSig_IndelSig_ips[j,]
    pvalue_matrix[i,j-3] <- cor.test(t(a), t(b), method = "pearson")$p.value
    
  }
}
pvalue_matrix <- data.frame(pvalue_matrix)
pvalue_matrix[is.na(pvalue_matrix)] <- 0
names(pvalue_matrix) <- rownames(SubSig_IndelSig_ips)[4:18]
rownames(pvalue_matrix) <- rownames(SubSig_IndelSig_ips)[1:3]
pvalue_matrix$Sub_sig <- rownames(pvalue_matrix)
pvalue_matrix_melt <- melt(pvalue_matrix,c("Sub_sig"))
names(pvalue_matrix_melt) <- c("Sub_sig","Indel_type","pvalue")

corr_pvalue_matrix_melt <- merge(corr_matrix_melt,pvalue_matrix_melt,by=c("Sub_sig","Indel_type"))
corr_pvalue_matrix_melt[corr_pvalue_matrix_melt$pvalue<0.0001,"pvalue"] <- 0.0001
corr_pvalue_matrix_melt$pvalue_log10 <- -log10(corr_pvalue_matrix_melt$pvalue)
indel_positions <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                     "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                     "Complex") 

pdf(file="corr_pvalue_matrix_melt2.pdf", onefile=TRUE,height=3,width=10, useDingbats=FALSE)
p <- ggplot(corr_pvalue_matrix_melt,aes(x=Indel_type,y=Sub_sig))+geom_point(aes(size=pvalue_log10,colour=corr))+ scale_radius()
#p <- p + scale_colour_gradient2(name = "Correlation", trans = "log10",low = "blue",midpoint=-2, high = "red",breaks = c(1,0.1,0.01,0.001), labels = c(1,0.1,0.01,0.001))
p <- p + scale_colour_gradient2(name = "Correlation",low = "blue",midpoint=0, high = "red",breaks = c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1),limits=c(-1,1))
p <- p+scale_x_discrete(limits = indel_positions,labels = indel_positions) + ylab("Signature")+ xlab("Indel type")
p <- p+theme(axis.text.x=element_text(angle=45, vjust=1, size=10,colour = "black",hjust=1),
             axis.text.y=element_text(size=10,colour = "black"),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
#p <- p+facet_grid(.~types)
print(p)
dev.off()

# Figure S4
# Signature fitting for fibroblasts
# Fibroblast substitution burden
#######################################
# mutations in fibro Subs
#######################################
subs <- read.table("../WGS_subs_final.txt", sep = "\t", header = T, as.is = T)
fibro_subs <- subs[subs$fibro_VAF>0,] 

fibro_subs$ID <- paste0(fibro_subs$Chrom,"_",fibro_subs$Pos)
fibro_subsID_summary <- data.frame(table(fibro_subs$ID))
names(fibro_subsID_summary) <- c("ID","Freq")
fibro_subs <- merge(fibro_subs,fibro_subsID_summary,by="ID")
fibro_subs_uniq <- unique(fibro_subs[,c("Chrom","Pos","Fibro","Ref","Alt","Fibro_nREF","Fibro_nALT","fibro_VAF")])
write.table(fibro_subs_uniq,"fibro_subs_uniq.txt",sep = "\t", col.names = T, row.names = F, quote = F)
fibro_subs_uniq <- read.table("fibro_subs_uniq.txt", sep = "\t", header = T, as.is = T)
# Profiles
mut_catalogue <- GenCatalogue(fibro_subs_uniq,"Fibro")
write.table(mut_catalogue,"fibro_subs_catalogue.txt",sep = "\t",col.names = T, row.names = F, quote = F)

summary_subs <- data.frame(table(fibro_subs_uniq$Fibro))
names(summary_subs) <- c("Fibro","Freq")
summary_subs <- summary_subs[order(summary_subs$Freq, decreasing = T),]
write.table(summary_subs,"fibro_summary_subs.txt",sep = "\t", col.names = T, row.names = F, quote = F)



PancanSig <- read.table("/nfs/cancer_archive04/xz3/b_1176/70_MSI/00_data/Pancan_signatures_subs_final.txt",sep = "\t",header = T, as.is = T)
names(PancanSig)[1] <- "MutationType"
Skin_tissueSig <- PancanSig[,c(1,which(sub("_[^_]+$","",names(PancanSig))=="Skin"))]
sub_catalouge <- read.table("fibro_subs_catalogue.txt", sep = "\t", header = T, as.is = T)

mut_sig <- merge(Skin_tissueSig,sub_catalouge[,-2],by="MutationType")
mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
sig_cat <- mut_sig[,2:dim(Skin_tissueSig)[2]]
mut_cat <- mut_sig[,(dim(Skin_tissueSig)[2]+1):dim(mut_sig)[2]]

a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
write.table(a$E_median_filtered,paste0("fibro_exposure_","Skin",".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
Exposure <-a$E_median_filtered
sample_exposure <- as.data.frame(Exposure)
sample_exposure$Sig <- rownames(sample_exposure)
sample_exposure_melt <- melt(sample_exposure,c("Sig"))
names(sample_exposure_melt) <- c("Sig","Fibro","exposure")
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_A","Sig"] <- "Culture"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_D","Sig"] <- "UV_Skin_D"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_J","Sig"] <- "UV_Skin_J"
sample_exposure_melt$Fibro <- chartr(".", "-", sample_exposure_melt$Fibro)
summary_subs <- read.table("fibro_summary_subs.txt",sep = "\t", header = T, as.is = T)
#cbbPalette <- c("#000000", "#999999","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#FF66FF", "#CC79A7")
cbbPalette <- c("#000000", "#F0E442","#F0E442", "#F0E442", "#F0E442", "#F0E442", "#F0E442", "#F0E442","#CC79A7", "#CC79A7")

pdf(file="fibro_sub_SigExposure2.pdf", onefile=TRUE,height=3,width=8, useDingbats=FALSE)
p <- ggplot(sample_exposure_melt,aes(x=Fibro,y=exposure,fill=Sig))+geom_bar(stat="identity", width=.8, position="fill")+scale_fill_manual(values=cbbPalette)
p <- p+scale_x_discrete(limits = as.character(summary_subs[,"Fibro"]))+xlab("Fibro")
#p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_blank(),
             axis.text.y=element_text(size=10,colour = "black"),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
print(p)
dev.off()

pdf(file="fibro_sub_SigExposure.pdf", onefile=TRUE,height=3,width=8, useDingbats=FALSE)
p <- ggplot(sample_exposure_melt,aes(x=Fibro,y=exposure,fill=Sig))+geom_bar(stat="identity", width=.8)+scale_fill_manual(values=cbbPalette)
p <- p+scale_x_discrete(limits = as.character(summary_subs[,"Fibro"]))+xlab("Fibro")
#p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_blank(),
             axis.text.y=element_text(size=10,colour = "black"),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
print(p)
dev.off()







