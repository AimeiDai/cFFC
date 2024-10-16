setwd("C:/Users/Administrator/OneDrive/1.cFF/manuscript/2024-3-18-first_revision/conservation")
setwd("/Users/aimeidai/OneDrive/1.cFF/manuscript/2024-3-18-first_revision/conservation")
options(stringsAsFactors=F)
library(reshape)
library(ggplot2)
require(RColorBrewer)
# cols <- c("darkgrey",rev(brewer.pal(9, "Set1")[c(1,3,2)]),"black")
# names(cols) =c("miR->non-TF","non-redundant","low-redundant","high-redundant","miR->TF")
cols2 <- c("darkgrey","black")
names(cols2) =c("non.cFFC.sites","cFFC.sites")

cols <- rev(brewer.pal(9, "Set1")[c(1,3,2)])
names(cols) = c("non-redundant","low-redundant","high-redundant")

mir_expr = read.table("../TF_miRNA_expression/miRFami.expr.txt")
expr_m = cbind("Embryonic"=rowMeans(mir_expr[,1:7]),"L3"=mir_expr[,9])

############ larva 3
mir_tg = read.table("l3_consmiR.lowThreshold.redundancy",header = T,quote="",sep="\t")
mir_tg$type2=ifelse(mir_tg$type%in%c("non-redundant","low-redundant","high-redundant"),"cFFC.sites","non.cFFC.sites")

cff_sits = subset(mir_tg,type%in%c("non-redundant","low-redundant","high-redundant"))
cff_sits$type2=cff_sits$type
#### site type distribution in cFFC redundancy
count = cbind(table(cff_sits[,c("Seed.match","type2")]),background=table(cff_sits$Seed.match))
mosaicplot(count, shade = T, las=2, main = "Frequency of targeting sites",color=T)

for (i in 1:3){
  tb = chisq.test(cbind(count[,i],count[,4]-count[,i]))$p.value
  print(ifelse(tb<=0.001,"***",ifelse(tb<=0.01,"**",ifelse(tb<=0.05,"*","ns"))))
}

prop = melt(sweep(count,2,colSums(count),"/"))

p <- ggplot(prop) + geom_bar(aes(x=X2,y=value,fill=X1),stat="identity")+ theme_bw(12) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",x="", y="Proportion of targeting sites") + scale_fill_manual(values=rev(brewer.pal(9, "Set1")[c(1,3,2)]))
ggsave("l3_consmiR.lowThreshold.SiteType.distribution.pdf",width = 2.8,height = 3.5)
#######################
all = rbind(data.frame(mir_tg,group="non.cFFC.sites"),data.frame(cff_sits,group="cFFC.sites"))
all$type2=factor(all$type2,levels=c("non.cFFC.sites","cFFC.sites","non-redundant","low-redundant","high-redundant"))
all$group = factor(all$group,levels=c("non.cFFC.sites","cFFC.sites"))
p <- ggplot(all) + geom_boxplot(aes(x=type2,y=as.numeric(PCT),color=type2))+ facet_wrap(~group,scales = "free_x") + scale_color_manual(values=c(cols,cols2))+ theme_classic(12) + theme(legend.background=element_rect(fill="transparent"), panel.grid=element_blank(),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), legend.position="topright", axis.text.x = element_text(angle=45,hjust=1))+ xlab("") + ylab("PCT")+scale_y_log10(limits=c(0.008,5))
ggsave("l3_consmiR.lowThreshold.redundancy-2.pdf",width = 2.5,height = 3.2)

data = all
data$log10PCT = log10(as.numeric(data$PCT))
data = subset(data,!is.na(log10PCT)&is.finite(log10PCT))
10^by(data$log10PCT,data[,c("type2")],median)

x1 = as.numeric(subset(all,type2=="non.cFFC.sites")$PCT)
x2 = as.numeric(subset(all,type2=="cFFC.sites")$PCT)
wilcox.test(x1,x2)$p.value
ix=combn(c("non-redundant","low-redundant","high-redundant"),2)
p_values=vector(length = 3)
for(i in 1:ncol(ix)){
  x1 = as.numeric(subset(all,type2==ix[1,i])$PCT)
  x2 = as.numeric(subset(all,type2==ix[2,i])$PCT)
  #print(ix[,i])
  p_values[i]=wilcox.test(x1,x2)$p.value
  #print(ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","NS"))))
}
p.adjust(p_values, method = "BH")

summary(aov(as.numeric(PCT)~type2,cff_sits))
######## targeting efficacy by site type
p <- ggplot(data=subset(all,group=="cFFC.sites")) + geom_boxplot(aes(x=type2,y=as.numeric(PCT),color=type2))+ facet_wrap(~Seed.match,scales = "free_x") + scale_color_manual(values=c(cols,cols2))+ theme_classic(12) + theme(legend.background=element_rect(fill="transparent"), panel.grid=element_blank(),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), legend.position="topright", axis.text.x = element_text(angle=45,hjust=1))+ xlab("") + ylab("PCT")+scale_y_log10(limits=c(0.01,5))
ggsave("l3_consmiR.lowThreshold.redundancy-siteType.pdf",width = 3.4,height = 2.9)
data=subset(all,group=="cFFC.sites")
data$log10PCT = log10(as.numeric(data$PCT))
data = subset(data,!is.na(log10PCT)&is.finite(log10PCT))
summary(aov(log10PCT~Seed.match*type,data))
# Df Sum Sq Mean Sq F value Pr(>F)    
# Seed.match          2   2130  1065.2 666.153 <2e-16 ***
#   type                2    463   231.6 144.864 <2e-16 ***
#   Seed.match:type     4      1     0.3   0.183  0.947    
# Residuals       37449  59882     1.6                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
pa = combn(c("non-redundant","low-redundant","high-redundant"),2)
for(i in unique(data$Seed.match)){
  print(i)
  for(j in 1:ncol(pa)){
    x1 = subset(data,type==pa[1,j]&Seed.match==i)
    x2 = subset(data,type==pa[2,j]&Seed.match==i)
    p=wilcox.test(x1$log10PCT,x2$log10PCT)$p.value
    cat(ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<=0.05,"*","ns")))," ")
  }
  cat("\n")
}

#####################
mir_red = unique(all[,c("miRs","type2","group")])
mir_red = mir_red[mir_red[,1]!="miR-2a-2/2c",]
mir_red[,4] = expr_m[mir_red[,1],2]

p <- ggplot(mir_red) + geom_boxplot(aes(x=type2,y=as.numeric(V4),color=type2)) + facet_wrap(~group,scales="free_x") + scale_color_manual(values=c(cols2,cols))+ theme_bw(10) + theme(legend.background=element_rect(fill="transparent"), panel.grid=element_blank(),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), legend.position="topright", axis.text.x = element_text(angle=45,hjust=1))+ xlab("") + ylab("miRNA expression \n(RPM)")+scale_y_log10(limits=c(0.001,10^7))
ggsave("l3.lowThreshold.redundancyvsExpresion.pdf",width = 2.4,height = 2.5)

for(i in 1:ncol(ix)){
  x1 = as.numeric(subset(mir_red,type2==ix[1,i])$V4)
  x2 = as.numeric(subset(mir_red,type2==ix[2,i])$V4)
  print(ix[,i])
  p=wilcox.test(x1,x2)$p.value
  print(ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","NS"))))
}
######### by expression level
# expr = read.table("../TF_miRNA_expression/miR_expression_clusters_heatmap.txt",header = T,row.names = 1)
# 
# ix = ifelse(expr$ifLarva==0,"ns",ifelse(expr$cluster==5,"high",ifelse(expr$cluster==4|expr$cluster==3,"median","low")))
# names(ix)=rownames(expr)
# all$ix = ix[all$miRs]
# data=unique(subset(all,group=="cFFC.sites"&ix!="ns")[,c("miRs","type","ix")])
# data$ix=factor(data$ix,levels=c("high","median","low"))
# p <- ggplot(data=subset(all,group=="cFFC.sites")) + geom_boxplot(aes(x=type2,y=as.numeric(PCT),color=type2))+ facet_wrap(~ix,scales = "free_x") + scale_color_manual(values=c(cols,cols2))+ theme_classic(12) + theme(legend.background=element_rect(fill="transparent"), panel.grid=element_blank(),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), legend.position="topright", axis.text.x = element_text(angle=45,hjust=1))+ xlab("") + ylab("PCT")+scale_y_log10(limits=c(0.01,5))
# ggsave("l3_consmiR.lowThreshold.redundancy-miRExpr.pdf",width = 3.4,height = 2.9)
# data=subset(all,group=="cFFC.sites")
# data$log10PCT = log10(as.numeric(data$PCT))
# data = subset(data,!is.na(log10PCT)&is.finite(log10PCT))
# for(i in c("low","median","high")){
#   data2=subset(data,ix==i)
#   print(summary(aov(log10PCT~type,data2)))
#}



############## embryonic stage

mir_tg = read.table("Embryonic_consmiR.lowThreshold.redundancy",header = T,quote="",sep="\t")
mir_tg$type2=ifelse(mir_tg$type%in%c("non-redundant","low-redundant","high-redundant"),"cFFC.sites","non.cFFC.sites")

cff_sits = subset(mir_tg,type%in%c("non-redundant","low-redundant","high-redundant"))
cff_sits$type2=cff_sits$type
#### site type distribution in cFFC redundancy
count = cbind(table(cff_sits[,c("Seed.match","type2")]),background=table(cff_sits$Seed.match))
for (i in 1:3){
  tb = chisq.test(cbind(count[,i],count[,4]-count[,i]))$p.value
  print(ifelse(tb<=0.001,"***",ifelse(tb<=0.01,"**",ifelse(tb<=0.05,"*","ns"))))
}
prop = melt(sweep(count,2,colSums(count),"/"))

p <- ggplot(prop) + geom_bar(aes(x=X2,y=value,fill=X1),stat="identity")+ theme_bw(12) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",x="", y="Proportion of targeting sites") + scale_fill_manual(values=rev(brewer.pal(9, "Set1")[c(1,3,2)]))
ggsave("Embryo_consmiR.lowThreshold.SiteType.distribution.pdf",width = 2.8,height = 3.5)
##########################
all = rbind(data.frame(mir_tg,group="non.cFFC.sites"),data.frame(cff_sits,group="cFFC.sites"))
all$type2=factor(all$type2,levels=c("non.cFFC.sites","cFFC.sites","non-redundant","low-redundant","high-redundant"))
all$group = factor(all$group,levels=c("non.cFFC.sites","cFFC.sites"))
p <- ggplot(all) + geom_boxplot(aes(x=type2,y=as.numeric(PCT),color=type2))+ facet_wrap(~group,scales = "free_x") + scale_color_manual(values=c(cols,cols2))+ theme_classic(12) + theme(legend.background=element_rect(fill="transparent"), panel.grid=element_blank(),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), legend.position="topright", axis.text.x = element_text(angle=45,hjust=1))+ xlab("") + ylab("PCT")+scale_y_log10(limits=c(0.01,5))
ggsave("Embryonic_consmiR.lowThreshold.redundancy.pdf",width = 2.5,height = 3.2)

data = all
data$log10PCT = log10(as.numeric(data$PCT))
data = subset(data,!is.na(log10PCT)&is.finite(log10PCT))
10^by(data$log10PCT,data[,c("type2")],median)

ix=combn(c("non.cFFC.sites","cFFC.sites","non-redundant","low-redundant","high-redundant"),2)
for(i in 1:ncol(ix)){
  x1 = as.numeric(subset(all,type2==ix[1,i])$PCT)
  x2 = as.numeric(subset(all,type2==ix[2,i])$PCT)
  print(ix[,i])
  p=wilcox.test(x1,x2)$p.value
  print(ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","NS"))))
}
summary(aov(as.numeric(PCT)~type,cff_sits))

######## targeting efficacy by site type
p <- ggplot(data=subset(all,group=="cFFC.sites")) + geom_boxplot(aes(x=type2,y=as.numeric(PCT),color=type2))+ facet_wrap(~Seed.match,scales = "free_x") + scale_color_manual(values=c(cols,cols2))+ theme_classic(12) + theme(legend.background=element_rect(fill="transparent"), panel.grid=element_blank(),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), legend.position="topright", axis.text.x = element_text(angle=45,hjust=1))+ xlab("") + ylab("PCT")+scale_y_log10(limits=c(0.01,5))
ggsave("embryonic_consmiR.lowThreshold.redundancy-siteType.pdf",width = 3.4,height = 2.9)
data=subset(all,group=="cFFC.sites")
data$log10PCT = log10(as.numeric(data$PCT))
data = subset(data,!is.na(log10PCT)&is.finite(log10PCT))
summary(aov(log10PCT~Seed.match*type,data))
# Df Sum Sq Mean Sq F value Pr(>F)    
# Seed.match          2    379  189.69 626.718 <2e-16 ***
#   type                2     71   35.69 117.933 <2e-16 ***
#   Seed.match:type     4      0    0.11   0.375  0.827    
# Residuals       34282  10376    0.30                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
pa = combn(c("non-redundant","low-redundant","high-redundant"),2)
for(i in unique(data$Seed.match)){
  print(i)
  for(j in 1:ncol(pa)){
    x1 = subset(data,type==pa[1,j]&Seed.match==i)
    x2 = subset(data,type==pa[2,j]&Seed.match==i)
    p=wilcox.test(x1$log10PCT,x2$log10PCT)$p.value
    cat(ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<=0.05,"*","ns")))," ")
  }
  cat("\n")
}

#########################
mir_red = unique(all[,c("miRs","type2","group")])
mir_red = mir_red[mir_red[,1]!="miR-2a-2/2c",]
mir_red[,4] = expr_m[mir_red[,1],1]

p <- ggplot(mir_red) + geom_boxplot(aes(x=type2,y=as.numeric(V4),color=type2)) + facet_wrap(~group,scales="free_x") + scale_color_manual(values=c(cols2,cols))+ theme_bw(10) + theme(legend.background=element_rect(fill="transparent"), panel.grid=element_blank(),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), legend.position="topright", axis.text.x = element_text(angle=45,hjust=1))+ xlab("") + ylab("miRNA expression \n(RPM)")+scale_y_log10(limits=c(0.001,10^7))
ggsave("embryonic.lowThreshold.redundancyvsExpresion.pdf",width = 2.4,height = 2.5)

for(i in 1:ncol(ix)){
  x1 = as.numeric(subset(mir_red,type2==ix[1,i])$V4)
  x2 = as.numeric(subset(mir_red,type2==ix[2,i])$V4)
  print(ix[,i])
  p=wilcox.test(x1,x2)$p.value
  print(ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","NS"))))
}







