library(ggplot2)
library(reshape)
options(stringsAsFactors=F)
require(RColorBrewer)
cols <- c("darkgrey",rev(brewer.pal(9, "Set1")[c(1,3,2)]))
names(cols) =c("miR->non-TF","non-redundant","low-redundant","high-redundant")

mir_tg = read.table("/public1/users/aimei/cFF/2.miR_targets/Conserved_Family_Info_converted_filtered.txt",header=T)

utr_len = read.table("UTR_length.txt")
utr_len_s = utr_len[utr_len[,1]%in%mir_tg[,5],]
rownames(utr_len_s) = utr_len_s[,2]

redundancy=read.table("../4.conservation/motif_with_redundancy.txt")
### number of redundant and non-redundant cFFs
red_s = melt(tapply(redundancy[,4],redundancy[,1],function(x){
	if(length(x)>1){
		return(ifelse(any(x=="high-redundant"),"high-redundant",ifelse(any(x=="low-redundant"),"low-redundant",ifelse(any(x=="non-redundant"),"non-redundant",ifelse(any(x=="miR->TF"),"miR->TF","miR->non-TF")))))
	}
	return(x)
	}))
rownames(red_s)=red_s[,1]

utr_len_s$type = red_s[utr_len_s[,1],2]

utr_len_s$logUTR = log2(utr_len_s$V3)
utr_len_s$Cut = as.character(cut(utr_len_s$logUTR, seq(4,16,1)))
########### expression variance
genes = subset(utr_len_s,(Cut=="(8,9]"|Cut=="(9,10]"|Cut=="(10,11]")&type!="miR->TF")
expr = read.table("/public1/users/aimei/cFF/expression_data/lanwq_stress/hs_stress.conf.tpm",header=T)

expr_gg = melt(as.matrix(expr))
# ############## Fold changes
expr_gg[,4:6]=matrix(unlist(strsplit(as.character(expr_gg[,2]),"[.]")),ncol=3,byrow=T)
data_summary <- function(x) {
   m <- mean(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}
##########################
# fc = by(expr_gg[,c(1,3,5)],expr_gg[,4],function(x){
# 	#browser()
# 	m=tapply(x[,2],x[,c(1,3)],function(y){mean(y,na.rm=T)})
# 	return(melt(log2(m[,2:5]/m[,1])))
# 	})
# fc_gg = rbind(data.frame(fc$Ctrl25,temp="C25"),data.frame(fc$Ctrl30,temp="C30"))
# fc_gg$type=factor(genes[as.character(fc_gg[,1]),5],levels=c("miR->TG","miR->TF","non-redundant","redundant"))
# fc_gg$V5 = factor(fc_gg$V5, levels=c("3h","6h","9h","12h"))

# p <- ggplot(subset(fc_gg,!is.na(type)), aes(x=V5,y=value, color=type))+ geom_hline(yintercept=0,color="red",linetype="dashed")+ geom_boxplot(outlier.shape=NA,width=0.5,position=position_dodge()) + facet_grid(type~temp)+ theme_bw(10) + ylab("log2 Fold Change") + theme(legend.background=element_rect(fill="transparent"),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), axis.text.x = element_text(angle=45,hjust=1), panel.grid=element_blank()) + coord_cartesian(ylim = c(-2,3))+scale_color_manual(values=c("darkgrey","lightseagreen","darkred","black"))
# ggsave("FCTimevs0h_lwq_stress.pdf",width=4,height=4)

#########################

# fc = by(expr_gg[,3:4],expr_gg[,c(1,5)],function(x){
# 	#browser()
# 	x1 = subset(x,V4=="Ctrl25")
# 	x2 = subset(x, V4=="Ctrl30")
# 	return(log2(mean(x2[,1])/mean(x1[,1])))
# 	})

# fcx = do.call("rbind",list(fc))
# fcx_sub = melt(subset(fcx,rownames(fcx)%in%genes[,2]))
# fcx_sub$type=factor(genes[as.character(fcx_sub[,1]),5],levels=c("miR->TG","miR->TF","non-redundant","redundant"))
# fcx_sub$X2 = factor(fcx_sub$X2, levels=c("0h","3h","6h","9h","12h"))


# p <- ggplot(fcx_sub, aes(x=X2,y=value, color=type))+ geom_hline(yintercept=0,color="red",linetype="dashed")+ stat_summary(fun.data=data_summary,position=position_dodge(0.5)) + theme_bw(10) + ylab("log2 Fold Change") + theme(legend.background=element_rect(fill="transparent"),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), axis.text.x = element_text(angle=0,hjust=1), panel.grid=element_blank()) + coord_cartesian(ylim = c(-2,2))+scale_color_manual(values=c("darkgrey","lightseagreen","darkred","black"))
# ggsave("FC_lwq_stress.pdf",width=5,height=3)


fc = by(expr_gg[,3:4],expr_gg[,c(1,5)],function(x){
	#browser()
	x1 = subset(x,V4=="Ctrl25")
	x2 = subset(x, V4=="Ctrl30")
	return((mean(x2[,1])-mean(x1[,1]))*2/(mean(x1[,1])+mean(x2[,1])))
	})

fcx = do.call("rbind",list(fc))
fcx_sub = melt(subset(fcx,rownames(fcx)%in%genes[,2]))
fcx_sub$Cut=genes[as.character(fcx_sub[,1]),"Cut"]
fcx_sub$X2 = factor(fcx_sub$X2, levels=c("0h","3h","6h","9h","12h"))
fcx_sub$type=factor(genes[as.character(fcx_sub[,1]),4],levels=c("non-motif","non-redundant","low-redundant","high-redundant"))
fcx_sub$Cut=factor(fcx_sub$Cut,levels=c("(8,9]","(9,10]","(10,11]"))

tmp = subset(fcx_sub,X2!="0h")
x=melt(tapply(tmp[,3],tmp[,c(2,4,5)],function(x){mean(x,na.rm=T)}))
#apply(x,2,function(x){sd(x,na.rm=T)/mean(x,na.rm=T)})
#by(fcx_sub[,3],fcx_sub[,4:5],function(x) {sd(x,na.rm=T)/mean(x,na.rm=T)})
x_gg = melt(tapply(x[,4],x[,2:3],function(x){mean(x,na.rm=T)}))
x_gg$type = factor(x_gg$type, levels=c("non-motif","non-redundant","low-redundant","high-redundant"))
x_gg$Cut=factor(x_gg$Cut,levels=c("(8,9]","(9,10]","(10,11]"))
#        Cut           type     value
# 1  (10,11] high-redundant 0.1722526
# 2    (8,9] high-redundant 0.2142443
# 3   (9,10] high-redundant 0.1790839
# 4  (10,11]  low-redundant 0.1869047
# 5    (8,9]  low-redundant 0.1763595
# 6   (9,10]  low-redundant 0.1762853
# 7  (10,11]    miR->non-TF 0.2553138
# 8    (8,9]    miR->non-TF 0.3534658
# 9   (9,10]    miR->non-TF 0.3410195
# 10 (10,11]  non-redundant 0.1129459
# 11   (8,9]  non-redundant 0.2049688
# 12  (9,10]  non-redundant 0.2785666

p <- ggplot()+ geom_hline(data=x_gg,aes(yintercept=value,color=type),linetype="dashed")+ stat_summary(data=fcx_sub, aes(x=X2,y=value, color=type),fun.data=data_summary,position=position_dodge(0.7)) + facet_grid(Cut~.) + theme_bw(12) + ylab("Expression difference") + theme(legend.background=element_rect(fill="transparent"),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), axis.text.x = element_text(angle=0,hjust=0.5), panel.grid=element_blank()) +scale_color_manual(values=cols)
ggsave("ED_lwq_stress.pdf",width=4.5,height=5)

by(fcx_sub[,c(2,3,5)],fcx_sub[,4],function(x){
	by(x[,c(2,3)],x[,1],function(y){
		#browser()
		by(y[,1],y[,2],function(z){mean(z,na.rm=T)})
		})
	})

for (c in c("(8,9]","(9,10]","(10,11]")){
	print(c)
	for (t in c("0h","3h","6h","9h","12h")){
		print(t)
		for(i in 1:ncol(pairs)){
			x1 = subset(fcx_sub,type==pairs[1,i]&X2==t&Cut==c)
			x2 = subset(fcx_sub,type==pairs[2,i]&X2==t&Cut==c)
			print(wilcox.test(x1[,3],x2[,3])$p.value)
		}
	}
}

write.table(subset(genes,Cut=="(10,11]"&type=="non-redundant")[,2],"non-redundant.10_11.txt",row.names=F,col.names=F, sep="\t",quote=F)
write.table(subset(genes,Cut=="(8,9]"&type=="high-redundant")[,2],"high-redundant.8_9.txt",row.names=F,col.names=F, sep="\t",quote=F)

# ########## expression variance
# cvx = melt(tapply(expr_gg[,3],expr_gg[,c(1,4,5)],function(x){mean(x)/sd(x)}))

# cvx_sub = subset(cvx,X1%in%genes[,2])


# cvx_sub$type=factor(genes[as.character(cvx_sub[,1]),5],levels=c("miR->TG","miR->TF","non-redundant","redundant"))
# cvx_sub$V5 = factor(cvx_sub$V5, levels=c("0h","3h","6h","9h","12h"))

# p <- ggplot(cvx_sub, aes(x=V5,y=value,color=V4)) + geom_boxplot(outlier.shape = NA,position=position_dodge(0.9)) + facet_grid(type~.) + theme_bw(10) + ylab("Coefficient of variance") +scale_color_manual(values=c("black","darkred")) + theme(legend.background=element_rect(fill="transparent"),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), axis.text.x = element_text(angle=45,hjust=1), panel.grid=element_blank()) + coord_cartesian(ylim = c(0,40)) 
# ggsave("CV_lwq_stress.pdf",width=3.5,height=5)

# # by(cvx_sub[,3:5],cvx_sub[,2],function(x){
# # 	by(x[,2:3],x[,1],function(y){
# # 		x1=subset(y,type=="miR->TG")[,1]
# # 		x2=subset(y,type=="miR->TF")[,1]
# # 		x3=subset(y,type=="non-redundant")[,1]
# # 		x4=subset(y,type=="redundant")[,1]
# # 		return(c(wilcox.test(x1,x2)$p.value,
# # 		wilcox.test(x2,x3)$p.value,
# # 		wilcox.test(x3,x4)$p.value,
# # 		wilcox.test(x2,x4)$p.value))
# # 	})
# # 	})

# by(cvx_sub[,2:4],cvx_sub[,5],function(x){
# 	by(x[,c(1,3)],x[,2],function(y){
# 		wilcox.test(y[,2]~y[,1])$p.value
# 		})
# 	})

