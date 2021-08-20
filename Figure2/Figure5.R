source("/rds/project/sn206/rds-sn206-nik-zainal/users/xz388/cancer_archive04/b_1176/15_subs_20180907/00_data/sub_common.R")

#########################
# New figures 5
#########################
PancanSig <- read.table("./Pancan_signatures_subs_final.txt",sep = "\t",header = T, as.is = T)
names(PancanSig)[1] <- "MutationType"
Skin_tissueSig <- PancanSig[,c(1,which(sub("_[^_]+$","",names(PancanSig))=="Skin"))]

ips_denovo_final <- read.table("ips_final.txt", sep = "\t", header = T, as.is = T)
ips_realdenovo_final <- ips_denovo_final
ips_realdenovo_final$denovo <- "shared"
ips_realdenovo_final[ips_realdenovo_final$Fibro_nALT==0,]$denovo <- "denovo"

# shared
sub_catalouge <- GenCatalogue(ips_realdenovo_final[ips_realdenovo_final$denovo=="shared",],"ips")
mut_sig <- merge(Skin_tissueSig,sub_catalouge[,-2],by="MutationType")
mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
sig_cat <- mut_sig[,2:dim(Skin_tissueSig)[2]]
mut_cat <- mut_sig[,(dim(Skin_tissueSig)[2]+1):dim(mut_sig)[2]]

a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
write.table(a$E_median_filtered,paste0("ips_exposure_","Skin_shared",".txt"),sep = "\t",row.names = T, col.names = T, quote = F)

# denovo
sub_catalouge <- GenCatalogue(ips_realdenovo_final[ips_realdenovo_final$denovo=="denovo",],"ips")
mut_sig <- merge(Skin_tissueSig,sub_catalouge[,-2],by="MutationType")
mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
sig_cat <- mut_sig[,2:dim(Skin_tissueSig)[2]]
mut_cat <- mut_sig[,(dim(Skin_tissueSig)[2]+1):dim(mut_sig)[2]]

a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
write.table(a$E_median_filtered,paste0("ips_exposure_","Skin_denovo",".txt"),sep = "\t",row.names = T, col.names = T, quote = F)

summary_subs <- read.table("ips_summary_subs.txt",sep = "\t", header = T, as.is = T)

# Fig2A
Exposure_shared <- read.table("ips_exposure_Skin_shared.txt", sep = "\t", header = T, as.is = T)
sample_exposure <- as.data.frame(Exposure_shared)
sample_exposure$Sig <- rownames(sample_exposure)
sample_exposure_melt <- melt(sample_exposure,c("Sig"))
names(sample_exposure_melt) <- c("Sig","ips","exposure")
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_A","Sig"] <- "Culture"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_D","Sig"] <- "UV_Skin_D"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_J","Sig"] <- "UV_Skin_J"
sample_exposure_melt$ips <- chartr(".", "-", sample_exposure_melt$ips)
sample_exposure_melt$denovo <- "shared"
sample_exposure_melt_shared <- sample_exposure_melt

if(TRUE){
  cbbPalette <- c("#000000", "#F0E442","#F0E442", "#F0E442", "#F0E442", "#F0E442", "#F0E442", "#F0E442","#CC79A7", "#CC79A7")
  pdf(file="ips_sub_SigExposure_sharedwithFibro.pdf", onefile=TRUE,height=3,width=8, useDingbats=FALSE)
  p <- ggplot(sample_exposure_melt_shared,aes(x=ips,y=exposure,fill=Sig))+geom_bar(stat="identity", width=.8, position="fill")+scale_fill_manual(values=cbbPalette)
  p <- p+scale_x_discrete(limits = as.character(summary_subs[,"ips"]))+xlab("ips")
  p <- p+ scale_y_continuous(labels=percent)
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
  
}

Exposure_denovo <- read.table("ips_exposure_Skin_denovo.txt", sep = "\t", header = T, as.is = T)
sample_exposure <- as.data.frame(Exposure_denovo)
sample_exposure$Sig <- rownames(sample_exposure)
sample_exposure_melt <- melt(sample_exposure,c("Sig"))
names(sample_exposure_melt) <- c("Sig","ips","exposure")
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_A","Sig"] <- "Culture"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_D","Sig"] <- "UV_Skin_D"
sample_exposure_melt[sample_exposure_melt$Sig=="Skin_J","Sig"] <- "UV_Skin_J"
sample_exposure_melt$ips <- chartr(".", "-", sample_exposure_melt$ips)
sample_exposure_melt$denovo <- "denovo"
sample_exposure_melt_denovo <- sample_exposure_melt

if(TRUE){
  cbbPalette <- c("#000000", "#F0E442","#F0E442", "#F0E442", "#F0E442", "#F0E442", "#F0E442", "#F0E442","#CC79A7", "#CC79A7")
  pdf(file="ips_sub_SigExposure_denovo.pdf", onefile=TRUE,height=3,width=8, useDingbats=FALSE)
  p <- ggplot(sample_exposure_melt_denovo,aes(x=ips,y=exposure,fill=Sig))+geom_bar(stat="identity", width=.8, position="fill")+scale_fill_manual(values=cbbPalette)
  p <- p+scale_x_discrete(limits = as.character(summary_subs[,"ips"]))+xlab("ips")
  p <- p+ scale_y_continuous(labels=percent)
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
  
}
sample_exposure_melt_all <- rbind(sample_exposure_melt_denovo,sample_exposure_melt_shared)
write.table(sample_exposure_melt_all, "ips_sample_sig_exposure_denovo_shared.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# Figure5B
# all sig exposure
sample_exposure_melt_all <- read.table("ips_sample_sig_exposure_denovo_shared.txt", sep = "\t", header = T, as.is = T)
sample_exposure <- dcast(sample_exposure_melt_all,ips+denovo~Sig,value.var="exposure")
sample_exposure$UV_Skin_D_J <- sample_exposure$UV_Skin_D + sample_exposure$UV_Skin_J
sample_exposure$Others <- sample_exposure$Skin_B + sample_exposure$Skin_C+ sample_exposure$Skin_E+ sample_exposure$Skin_F+ sample_exposure$Skin_G+ sample_exposure$Skin_H+sample_exposure$Skin_I

sample_exposure_melt2 <- melt(sample_exposure,c("ips","denovo"))
sample_exposure_melt2 <- sample_exposure_melt2[sample_exposure_melt2$variable %in%c("Culture","UV_Skin_D_J","Others"),]
names(sample_exposure_melt2) <- c("Sample","denovo","type","Freq")
a2 <- dcast(sample_exposure_melt2,Sample+denovo~type,value.var="Freq")
a2[,3:dim(a2)[2]] <-  a2[,3:dim(a2)[2]]/rowSums(a2[,3:dim(a2)[2]])[row(a2[,3:dim(a2)[2]])]
a3 <- melt(a2,c("Sample","denovo"))
names(a3) <- c("Sample","denovo","type","expo")
cbbPalette <- c("#E69F00", "#56B4E9")
pdf(file="SubSigs_ips_subclonecatalogue_exposure_share_uniqu2.pdf", onefile=TRUE,height=2.3,width=4, useDingbats=FALSE)
p <- ggplot(a3,aes(x=type, y=expo, fill=denovo))+geom_boxplot(outlier.size=0.5, size=0.5)+xlab("Number of samples")+scale_fill_manual(values=cbbPalette)+scale_y_continuous(labels=percent)
#p <- p+stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), geom="pointrange", color = denovo)
#p <- p+scale_y_continuous(breaks = seq(-350,350,50), labels = c(seq(350,0,-50),seq(50,350,50))) 
p <- p+theme(axis.text.x=element_text(size=5,angle=90,vjust=0.5,colour = "black"),
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

wt <- wilcox.test(sample_exposure_melt2[sample_exposure_melt2$variable=="Culture" & sample_exposure_melt2$denovo=="denovo","value"], sample_exposure_melt2[sample_exposure_melt2$variable=="Culture"& sample_exposure_melt2$denovo=="shared","value"]) #p.value:1.277813e-105
wt <- wilcox.test(sample_exposure_melt2[sample_exposure_melt2$variable=="UV_Skin_D_J" & sample_exposure_melt2$denovo=="denovo","value"], sample_exposure_melt2[sample_exposure_melt2$variable=="UV_Skin_D_J"& sample_exposure_melt2$denovo=="shared","value"]) #p.value:0.08493461
wt <- wilcox.test(sample_exposure_melt2[sample_exposure_melt2$variable=="Others" & sample_exposure_melt2$denovo=="denovo","value"], sample_exposure_melt2[sample_exposure_melt2$variable=="Others"& sample_exposure_melt2$denovo=="shared","value"]) #p.value:1.277813e-105


# signature emerge order
sample_exposure_melt_all <- read.table("ips_sample_sig_exposure_denovo_shared.txt", sep = "\t", header = T, as.is = T)
sample_exposure <- dcast(sample_exposure_melt_all,ips+denovo~Sig,value.var="exposure")
sample_exposure$UV_Skin_D_J <- sample_exposure$UV_Skin_D + sample_exposure$UV_Skin_J
sample_exposure$Others <- sample_exposure$Skin_B + sample_exposure$Skin_C+ sample_exposure$Skin_E+ sample_exposure$Skin_F+ sample_exposure$Skin_G+ sample_exposure$Skin_H+sample_exposure$Skin_I

sample_exposure_melt2 <- melt(sample_exposure,c("ips","denovo"))
sample_exposure_melt2 <- sample_exposure_melt2[sample_exposure_melt2$variable %in%c("Culture","UV_Skin_D_J","Others"),]
sample_exposure_melt2$sig_2 <- 1
sample_exposure_melt2[sample_exposure_melt2$value==0,"sig_2"] <- 0
sample_exposure_melt2_dcast <- dcast(sample_exposure_melt2,variable~denovo,value.var="sig_2",fun.aggregate=sum)
sample_exposure_melt3 <- melt(sample_exposure_melt2_dcast,c("variable"))
names(sample_exposure_melt3) <- c("Sig","type","Freq")
sample_exposure_melt3[sample_exposure_melt3$type=="shared","Freq"] <- -1*sample_exposure_melt3[sample_exposure_melt3$type=="shared","Freq"] 

cbbPalette <- c("#E69F00", "#56B4E9")
pdf(file="Subsigs_ips_subclonecatalogue_emerge_share_uniqu2.pdf", onefile=TRUE,height=2.3,width=5, useDingbats=FALSE)
p <- ggplot(sample_exposure_melt3,aes(x=Sig, y=Freq, fill=type))+geom_bar(stat="identity",width=0.6)+xlab("Number of samples")+scale_fill_manual(values=cbbPalette)
#p <- p+scale_y_continuous(breaks = seq(-350,350,50), labels = c(seq(350,0,-50),seq(50,350,50))) 
p <- p+scale_y_continuous(limits=c(-350,350))
p <- p+ coord_flip()+theme(axis.text.x=element_text(size=10,colour = "black"),
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



# Figure5C




# Figure5D
pass <- read.table("./hipsci_files.tsv", sep = "\t", header = T, as.is = T,check.names = F)
names(pass) <- c("ips","celltype","assay","passage","disease","sex")
Exposure <- read.table("../ips_exposure_Skin.txt", sep = "\t", header = T, as.is = T)
sample_exposure <- as.data.frame(Exposure)
sample_exposure$Sig <- rownames(sample_exposure)
sample_exposure_melt <- melt(sample_exposure,c("Sig"))
names(sample_exposure_melt) <- c("Sig","ips","exposure")
sample_exposure_melt$ips <- chartr(".", "-", sample_exposure_melt$ips)
sample_exposure_sig <- dcast(sample_exposure_melt,ips~Sig,value.var="exposure")
sample_exposure_sig_pass <- merge(sample_exposure_sig,pass[,c(1,4)],by="ips")
sample_exposure_sig_pass <- sample_exposure_sig_pass[sample_exposure_sig_pass$Skin_A>0,]
sample_exposure_sig_pass$allmutation <- rowSums(sample_exposure_sig_pass[,-1]) 
sample_exposure_sig_pass$UV <- sample_exposure_sig_pass$Skin_D+sample_exposure_sig_pass$Skin_J

# ips ~ culture
pdf(file="ips_sub_passage_correlation_culture.pdf", onefile=TRUE,height=5,width=5, useDingbats=FALSE)
p <- ggplot(sample_exposure_sig_pass,aes(x=passage,y=Skin_A))+geom_point()+geom_smooth(method="lm", colour="black")
#p <- p+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
p <- p+theme(axis.text.x=element_text(size=10,colour = "black"),
             axis.text.y=element_text(size=10,colour = "black"),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
print(p)
dev.off()

