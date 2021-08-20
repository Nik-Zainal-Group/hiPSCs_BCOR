source("/rds/project/sn206/rds-sn206-nik-zainal/users/xz388/cancer_archive04/b_1176/15_subs_20180907/00_data/sub_common.R")

#########################
# New figures 3
#########################

# shared ips with cosine similarity
Samedonor_ips <- read.table("Samedonor_ips_number_cossim.txt", sep="\t",header = T, as.is = T)
b <- read.table("Samedonor_ips_sharedmuts.txt",sep = "\t", header = T, as.is = T)
b_dcast <- dcast(b,ips~Flag,value.var="Freq")
b_dcast$shared <- b_dcast$d_shared_inFibro+b_dcast$c_shared_notinFibro
b_dcast <- b_dcast[order(b_dcast$shared, b_dcast$ips, decreasing = T),]
b_dcast$Fibro <- sub('\\_.*', '', b_dcast$ips)
Samedonor_ips_shared <- merge(Samedonor_ips,unique(b_dcast[,c("Fibro","shared")]),by="Fibro")
Samedonor_ips_shared$share_group <- "1"
Samedonor_ips_shared[Samedonor_ips_shared$shared<=10,]$share_group <- "0"

pdf(file="Samedonor_ips_number_cossim_sharedinfo.pdf", onefile=TRUE,height=2.8,width=4.2, useDingbats=FALSE)
p <- ggplot(Samedonor_ips_shared,aes(x=min,y=max,colour=cossim, shape=share_group))+geom_point()+scale_colour_gradientn(colours=rainbow(4))+scale_shape_manual(values=c(1, 16))
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
print(p)
dev.off()


# Signature of shared mutations





