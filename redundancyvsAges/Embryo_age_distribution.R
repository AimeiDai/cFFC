library(ggplot2)
library(reshape)
options(stringsAsFactors=F)
setwd("C:/Users/Administrator/OneDrive/1.cFF/manuscript/2024-3-18-first_revision/redundancyvsAges")
require(RColorBrewer)
cols <- c("darkgrey",rev(brewer.pal(9, "Set1")[c(1,3,2)]))
names(cols) =c("non-cFFC","non-redundant","low-redundant","high-redundant")

kondoages = read.table("kondoages.csv",sep=",",header=T)
redundancy=read.table("../expression_noise/embryo_consmiR.lowThreshold.allTargetType",sep="\t",header = T)
redundancy$type= ifelse(redundancy$type%in%c("miR->non-TF","miR->TF"),"non-cFFC",redundancy$type)
rownames(redundancy)=redundancy[,1]
utr_len = read.table("../expression_noise/UTR_length.txt",sep="\t")
utr_len_s = data.frame(len=tapply(utr_len[,3],utr_len[,2],max))
utr_len_s$geneId = rownames(utr_len_s)
utr_len_s$logUTR = log2(utr_len_s$len)
utr_len_s$Cut = as.character(cut(utr_len_s$logUTR, seq(0,16,1), include.lowest = T,right=T,ordered_result=T))
utr_len_s$type = redundancy[rownames(utr_len_s),2]
#red_s = utr_len_s
age_groups = read.table("age_groups.txt",sep="\t")

length(intersect(rownames(utr_len_s),kondoages[,1]))
length(rownames(utr_len_s))
length(unique(kondoages[,1]))

aa = unique(merge(utr_len_s[,c(2,5)],kondoages[,c(1,6,7)],by.x=1,by.y=1))
aa = subset(aa,type%in%c("non-cFFC","non-redundant","low-redundant","high-redundant"))
aa$type = factor(aa$type,levels=c("non-cFFC","non-redundant","low-redundant","high-redundant"))

p<-ggplot() + geom_boxplot(data=aa,aes(x=type,y=tau,color=type))+ theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="none", legend.title=element_blank()) + labs(title="",y=expression(paste("Tissue specificity (",tau,")",sep="")), x="") +scale_color_manual(values=cols) + ylim(0,1.35)
ggsave("embryo_redundancyvsTau.pdf",width=1.5,height=3)

by(aa[,3],aa[,2],median)
# aa[, 2]: non-cFFC
# [1] 0.798
# ------------------------------------------------------------------------------ 
#   aa[, 2]: non-redundant
# [1] 0.747
# ------------------------------------------------------------------------------ 
#   aa[, 2]: low-redundant
# [1] 0.539
# ------------------------------------------------------------------------------ 
#   aa[, 2]: high-redundant
# [1] 0.492
median(subset(aa,type!="non-cFFC")$tau)
# 0.545
median(subset(aa,type!="non-cFFC"&type!="non-redundant")$tau)
# 0.517
pairs=combn(c("non-cFFC","non-redundant","low-redundant","high-redundant"),2)
p_vals = vector(length=ncol(pairs))
for(i in 1:ncol(pairs)){
  x1 = subset(aa,type==pairs[1,i])
  x2 = subset(aa,type==pairs[2,i])
  p_vals[i]=wilcox.test(x1[,3],x2[,3])$p.value
  #print(ifelse(P<=0.001,"***",ifelse(P<=0.01,"**",ifelse(P<=0.05,"*","NS"))))
}
P = p.adjust(p_vals,method = "BH")
print(ifelse(P<=0.001,"***",ifelse(P<=0.01,"**",ifelse(P<=0.05,"*","NS"))))
#[1] "*"   "***" "***" "***" "***" "***"

bb = merge(aa,age_groups,by.x=4,by.y=2)
r=rbind(table(bb[,c(3,5)]))
r_gg = melt(r/rowSums(r))
r_gg$X2 = factor(r_gg$X2,levels=c("pre-Bilateria","pre-Drosophilid","Drosophilid"))
r_gg$X1 = factor(r_gg$X1,levels=c("non-cFFC","non-redundant","low-redundant","high-redundant"))

p<-ggplot() + geom_bar(data=r_gg,aes(x=X1,y=value,fill=X2),stat="identity",color="white")+ theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of genes", x="") +scale_fill_manual(values=c("darkred","lightseagreen","darkgrey"))
ggsave("Embryo_redundancyvsAges.pdf",width=3,height=3)

for(i in 2:4){
	for (j in 1:3){
		m = matrix(c(sum(r[i,-j]),r[i,j],sum(r[i-1,-j]),r[i-1,j]),nrow=2,byrow=T)
		p = fisher.test(m)$p.value
		cat(i,j,": ")
		print(ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<=0.05,"*","NS"))))
	}
}
# 2 1 : [1] "*"
# 2 2 : [1] "*"
# 2 3 : [1] "NS"
# 3 1 : [1] "***"
# 3 2 : [1] "***"
# 3 3 : [1] "NS"
# 4 1 : [1] "***"
# 4 2 : [1] "***"
# 4 3 : [1] "NS"

bb$V1 = factor(bb$V1,levels=c("pre-Bilateria","pre-Drosophilid","Drosophilid"))
p <-ggplot() + geom_boxplot(data=bb,aes(x=type,y=tau,color=V1)) + facet_wrap(~V1)+ theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y=expression(paste("Tissue specificity (",tau,")",sep="")), x="") +scale_color_manual(values=c("darkred","lightseagreen","darkgrey"))
ggsave("Embryo_redundancyvsAgesTau.pdf",width=4.5,height=3)

fit=aov(tau~type*V1,data=bb)
summary(fit)
posthoc <- TukeyHSD(fit)
print(posthoc)

bb$tau_rg = cut(bb$tau,c(0,0.8,0.9,1))
x = rbind(t(table(bb[,c("type","tau_rg")])))
#,background=table(bb[,5])
x_gg = melt(x/rowSums(x))
x_gg$X2 = factor(x_gg$X2,levels=c("non-cFFC","non-redundant","low-redundant","high-redundant"))
x_gg$X1=factor(x_gg$X1,levels=c("(0,0.8]","(0.8,0.9]","(0.9,1]"))
p <-ggplot() + geom_bar(data=x_gg,aes(x=X1,y=value,fill=X2),stat="identity",color="white")+ theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",x=expression(paste("Tissue specificity (",tau,")",sep="")), y="Proportion of genes") +scale_fill_manual(values=cols)
ggsave("Embryo_tauGroupsRedundancy.pdf",width=2.8,height=3)

#non-motif
for(i in 2:3){
	for(j in 1:4){
		m = matrix(c(sum(x[i,-j]),x[i,j],sum(x[i-1,-j]),x[i-1,j]),nrow=2,byrow=T)
		p = fisher.test(m)$p.value
		cat(i,j,": ")
		print(ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<=0.05,"*","NS"))))
	}
}

# 2 1 : [1] "***"
# 2 2 : [1] "***"
# 2 3 : [1] "NS"
# 2 4 : [1] "***"
# 3 1 : [1] "***"
# 3 2 : [1] "NS"
# 3 3 : [1] "NS"
# 3 4 : [1] "***"

bb$redundancy = ifelse(bb$type=="non-cFFC",0,ifelse(bb$type=="non-redundant",1,ifelse(bb$type=="low-redundant",2,3)))
pdf("Embryo_TukeyHSD.pdf",width=4,height=7)
par(mfrow = c(3, 1),mar=c(5, 14, 2, 2) + 0.1)
for(i in unique(bb$V1)){
	y = subset(bb,V1==i)
	fit = aov(tau~type,data=y)
	print(summary(fit))
	posthoc <- TukeyHSD(fit)
	print(posthoc)
	plot(TukeyHSD(fit, conf.level=.95), las = 2)
}
dev.off()
# Df Sum Sq Mean Sq F value   Pr(>F)    
# type           3    0.8 0.26682   6.304 0.000294 ***
#   Residuals   2486  105.2 0.04233                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Df Sum Sq Mean Sq F value   Pr(>F)    
# type          3  0.882 0.29416   7.039 0.000121 ***
#   Residuals   508 21.229 0.04179                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Df Sum Sq Mean Sq F value Pr(>F)    
# type           3   4.94  1.6468   47.78 <2e-16 ***
#   Residuals   1262  43.50  0.0345                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Tukey multiple comparisons of means
# 95% family-wise confidence level