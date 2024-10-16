
options(stringsAsFactors=F)

gtf = read.table("/public1/users/aimei/genomes/Drosophila_melanogaster/BDGP6_r35/protein_coding.gtf",stringsAsFactors=F,sep="\t",quote="")
tmp = regexpr("gene_id \"\\S+?\";",gtf[,9])
gtf$gene_id = substr(gtf[,9],tmp+9,tmp+attr(tmp,"match.length")-3)
proteins=unique(gtf$gene_id)

tgs=read.table("../5.polymorphism/TSFly_r6.19_3UTRs.gff")
tgs$gene_id = do.call("rbind",strsplit(tgs[,9],":"))[,2]
tgs = intersect(proteins,tgs$gene_id)

##############
all = subset(read.table("Conserved_Family_Info.txt",stringsAsFactors=F, header=T,sep="\t",quote=""),Species.ID==7227)
all = subset(all,Gene.ID%in%tgs)
miRfami = unique(subset(read.table("../2.miR_TF_expression/miR_expression/miR_Family_Info_modified.txt",stringsAsFactors=F, header=T,sep="\t",quote=""),Family.Conservation.!=-1)[,1:2])
rownames(miRfami)=miRfami[,2]
all$miRs = miRfami[all[,1],1]


motifs = read.table("embryonic_motifs.consmiR.lowThreshold",stringsAsFactors=F)
motifs = subset(motifs,V3%in%tgs)

type = do.call("rbind",by(motifs,motifs[,3],function(x){
	#browser()
	ntfs = length(unique(x[,2]))
	nmirs = length(unique(x[,1]))
	ncFFCs = nrow(x)
	if(ntfs>1){
		if (nmirs == 1){return(c(ncFFCs,"B"))}
		if (nmirs > 1){return(c(ncFFCs,"D"))}
	}
	if (ntfs == 1){
		if (nmirs == 1){return(c(ncFFCs,"A"))}
		if (nmirs > 1){return(c(ncFFCs,"C"))}
	}
}))

type = data.frame(ncFFCs=as.numeric(type[,1]),type=type[,2],stringsAsFactors=F)

nmirs=tapply(motifs[,1],motifs[,3],function(x){
	return(length(unique(x)))
	})
ntfs=tapply(motifs[,2],motifs[,3],function(x){
	return(length(unique(x)))
	})

red_type=data.frame(nmirs,ntfs,type[names(nmirs),])
red_type$redundancy = ifelse(red_type$type=="A"|red_type$type=="B","non-redundant","redundant")
table(red_type$type)
by(red_type[,1],red_type[,4],sum) ## number of miRNAs 
#by(red_type[,2],red_type[,4],sum) ## number of tfs
by(red_type[,3],red_type[,4],sum) ## number of cFFCs

table(red_type[,5]) ## number of cFFC TGs
by(red_type[,1],red_type[,5],sum) ## number of miRNAs 
#by(red_type[,2],red_type[,5],sum) ## number of tfs
by(red_type[,3],red_type[,5],sum) ## number of cFFCs


th=7.24
red_type$type2 = ifelse(red_type[,1]==1,"non-redundant",ifelse(red_type[,1]>=th,"high-redundant","low-redundant"))

tf_tg = read.table("../1.TF_bindings/Embryo_tf_bindings.txt",stringsAsFactors=F,sep="\t",quote="")

all_tgs = data.frame(tgs,row.names=tgs)
# expr_mirs = read.table("../2.miR_TF_expression/miR_expression/embryo_miRs.cons.low.txt",stringsAsFactors=F)
# expr_all = subset(all, miRs%in%expr_mirs[,1])
expr_all = all
all_tgs$type = ifelse(all_tgs$tgs%in%expr_all[,2]&all_tgs$tgs%in%tf_tg[,2],"miR->TF",ifelse(all_tgs$tgs%in%expr_all[,2],"miR->non-TF","other"))

all_tgs[rownames(red_type),"type"] = red_type[,"type2"]
table(all_tgs[,2])
# high-redundant  low-redundant        miR->TF    miR->non-TF  non-redundant 
#           2549           5209              2            222           2324 
#          other 
#           3120 

write.table(all_tgs,"Embryo_consmiR.lowThreshold.allTargetType",row.names=F,col.names=T,sep="\t",quote=F)

tf_tg = tf_tg[tf_tg[,3] %in% tgs,]
expressed_tfs = read.table("../2.miR_TF_expression/TF_expression/filtered_embryonic_TFs.txt",stringsAsFactors=F)
tf_tg = tf_tg[tf_tg[,2]%in%expressed_tfs[,1],2:3]

mir_tg = expr_all


type=apply(mir_tg,1,function(x){
		#browser()
		y=motifs[motifs[,1]==x[12]&motifs[,3]==x[2],]
		if(nrow(y)==0){ ## non-motif sites
			if(x[2]%in%tf_tg[,1]){return("miR->TF")
			}else{return("miR->non-TF")}
		}
		if(nrow(y)>0){
			#browser()
			z = red_type[unique(y[,3]),]
			if(nrow(z)>1){
				return("error")
			}
			if(z[1,1]>th){return("high-redundant")}
			if(z[1,1]>1&z[1,1]<=th){return("low-redundant")}
			if(z[1,1]==1){return("non-redundant")}
		}
		})

save(type,file="Embryonic_consmiR.lowThreshold.redundancy_type.RData")

#load("Embryonic_redundancy_type.RData")
mir_tg$type=type
# x=subset(mir_tg,PCT!="NULL")
# x1=subset(x,type=="high-redundant")
# by(x$PCT,x$type,function(x){median(as.numeric(x))})

write.table(mir_tg,"Embryonic_consmiR.lowThreshold.redundancy",row.names=F,col.names=T,quote=F,sep="\t")



# bitwise = t(apply(mir_tg,1,function(x){
# 	#browser()
# 	as.numeric(unlist(strsplit(x[8],"")))
# 	}))
# mir_tg[,14:19] = bitwise

# conservation=do.call("rbind",by(mir_tg[,14:19],mir_tg$type,colSums))
# cons_freq = conservation/conservation[,1]
# colnames(cons_freq)=c("0","4","10","14","27","40")
# cons_freq_gg = melt(cons_freq)
# cons_freq_gg$X1=factor(cons_freq_gg$X1,levels=c("non-redundant","miR->non-TF","low-redundant","high-redundant","miR->TF"))

# require(RColorBrewer)
# cols <- c("darkgrey",rev(brewer.pal(9, "Set1")[c(1,3,2)]),"black")
# names(cols) =c("miR->non-TF","non-redundant","low-redundant","high-redundant","miR->TF")

# p <- ggplot() + geom_line(data=subset(cons_freq_gg,X1!="mix"),aes(x=X2,y=value,color=X1,group=X1),size=0.3,alpha=0.6,position=position_dodge(2)) + geom_point(data=subset(cons_freq_gg,X1!="mix"),aes(x=X2,y=value,color=X1,group=X1),size=2.5,position=position_dodge(2))+ theme_bw(10) + ylab("Conservation") + theme(legend.background=element_rect(fill="transparent"),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), legend.position=c(0.7,0.7), axis.text.x = element_text(angle=45,hjust=1)) + scale_x_continuous(name="Divergence Time",breaks=c(0,4,10,14,27,40),labels=c("D.mel","D.sim","D.yak","D.ana","D.pse","D.vir")) + ylim(0,1) + scale_color_manual(values=cols)
# ggsave("conservation_fourtype.pdf",width=3,height=2.5)

# conservation=conservation[c("miR->TF","high-redundant","low-redundant","miR->non-TF","non-redundant"),]
# cons_freq=cons_freq[c("miR->TF","high-redundant","low-redundant","miR->non-TF","non-redundant"),]
# for(i in 2:ncol(conservation)){
# 	print(i)
# 	ix = order(cons_freq[,i],decreasing=T)
# 	aa=combn(ix,2)
# 	for( j in 1:ncol(aa)){
# 		m = cbind(conservation[aa[,j],i],conservation[aa[,j],1] - conservation[aa[,j],i])
# 		cat(paste(aa[,j],collapse="v"),fisher.test(m)$p.value,"\n")
# 	}
# }

# [1] 2
# 1v2 0.0003653462 
# 1v3 1.164781e-06 
# 1v4 5.922668e-08 
# 1v5 4.531387e-09 
# 2v3 8.146926e-18 
# 2v4 6.455901e-59 
# 2v5 1.46564e-11 
# 3v4 0.0001227415 
# 3v5 0.0007632131 
# 4v5 0.04126612 
# [1] 3
# 1v2 0.003753959 
# 1v3 3.791772e-07 
# 1v4 1.745059e-08 
# 1v5 1.119199e-11 
# 2v3 1.027826e-42 
# 2v4 1.171395e-112 
# 2v5 6.859879e-26 
# 3v4 0.0003225673 
# 3v5 1.907445e-07 
# 4v5 5.857712e-05 
# [1] 4
# 1v2 0.01158116 
# 1v3 5.92487e-06 
# 1v4 4.837931e-07 
# 1v5 2.610239e-11 
# 2v3 2.111539e-34 
# 2v4 2.505776e-84 
# 2v5 2.345152e-24 
# 3v4 0.00975849 
# 3v5 1.489809e-08 
# 4v5 9.618649e-07 
# [1] 5
# 1v2 0.1006917 
# 1v3 1.080948e-05 
# 1v4 5.410524e-06 
# 1v5 4.167074e-12 
# 2v3 2.496442e-61 
# 2v4 9.046792e-111 
# 2v5 6.056398e-35 
# 3v4 0.6959295 
# 3v5 2.187593e-10 
# 4v5 2.167796e-10 
# [1] 6
# 1v2 0.08469925 
# 1v4 0.0001780582 
# 1v3 1.891899e-05 
# 1v5 5.433973e-12 
# 2v4 2.279572e-56 
# 2v3 5.579504e-49 
# 2v5 6.718944e-31 
# 4v3 0.00425038 
# 4v5 6.287265e-14 
# 3v5 2.240681e-10

#x=merge(mir_tg[,c(1:2,10)],all[,c(1,3,5,9,10)],by.x=c(1,2),by.y=c(1,2))
#write.table(mir_tg[,c(5:7,13)],"motif_with_redundancy.txt",sep="\t",row.names=F,col.names=F,quote=F)
##########################

# conservation=do.call("rbind",by(mir_tg[,3:8],mir_tg$ifTF,colSums))
# cons_freq = conservation/conservation[,1]
# colnames(cons_freq)=c("0","4","10","14","27","40")
# cons_freq_gg = melt(cons_freq)
# cons_freq_gg$X1=factor(cons_freq_gg$X1,levels=c("miR->TF","miR->non-TF"))
# p <- ggplot(cons_freq_gg,aes(x=X2,y=value,color=X1)) + geom_line(size=0.5,alpha=0.6) + geom_point(size=3.5)+ theme_bw(10) + ylab("Conservation") + theme(legend.background=element_rect(fill="transparent"),legend.text = element_text(size = 10), legend.key = element_rect(fill = "transparent", colour = "transparent") ,legend.title=element_blank(), legend.position=c(0.7,0.7), axis.text.x = element_text(angle=45,hjust=1)) + scale_x_continuous(name="Divergence Time",breaks=c(0,4,10,14,27,40),labels=c("D.mel","D.sim","D.yak","D.ana","D.pse","D.vir")) + ylim(0,1) + scale_color_manual(values=c("darkred","black"))
# ggsave("conservation_TF_TG.pdf",width=3.5,height=2.5)

# for(i in 2:ncol(conservation)){
# 	print(i)
# 	hl = cbind(conservation[c(1,2),i],conservation[c(1,2),1] - conservation[c(1,2),i])
# 	print(fisher.test(hl)$p.value)
# }

# [1] 2
# [1] 3.862648e-18
# [1] 3
# [1] 1.896994e-28
# [1] 4
# [1] 1.007331e-21
# [1] 5
# [1] 6.113535e-24
# [1] 6
# [1] 6.917358e-31


# x=merge(mir_tg[,c(1:2,9)],all[,c(1,3,5,9,10)],by.x=c(1,2),by.y=c(1,2))

# write.table(x[,c(4:6,3)],"TF_TG_type.txt",sep="\t",row.names=F,col.names=F,quote=F)

