options(stringsAsFactor=F)

est = read.table("sample_table.txt",stringsAsFactor=F)
fbgn2symbol = unique(read.table(paste("stringtie_out/",est[1,1],"/gene.abundance",sep=""),sep="\t",header=T,stringsAsFactor=F,quote="")[,c(1,2)])
rownames(fbgn2symbol)=fbgn2symbol[,1]
write.table(fbgn2symbol,"fbgn2symbol.txt",row.names=F,col.names=T,sep="\t",quote=F)

etpm = data.frame(row.names=fbgn2symbol[,1])
for (i in 1:nrow(est)){
	tmp = read.table(paste("stringtie_out/",est[i,1],"/gene.abundance",sep=""),sep="\t",header=T,row.names=1,stringsAsFactor=F,quote="")
	etpm[,i]=tmp[rownames(etpm),"TPM"]
}
colnames(etpm)=est[,2]

#etfs = read.table("/public1/users/aimei/cFF_revised/1.TF_bindings/Embryo_tf_bindings.txt",stringsAsFactor=F,sep="\t",quote="")
etfs = read.table("/public1/users/aimei/cFF_revised/1.TF_bindings/SupplementaryTable1-EmbryoChIP.tsv",stringsAsFactor=F,sep=",",quote="",header=T)
etfs=unique(etfs[,2])
#etfs = etfs[etfs!=""]
etpm_s = etpm[etfs,]
filtered_tfs = rownames(etpm_s)[rowSums(etpm_s[,1:6]>=10)>=3]
write.table(filtered_tfs,"filtered_embryonic_TFs.txt",row.names=F,col.names=F,sep="\t",quote=F)
rownames(etpm_s)=fbgn2symbol[rownames(etpm_s),2]
write.table(etpm_s,"Embryo_TFs.TPM",sep="\t",row.names=T,col.names=T,quote=F)


#ltfs = read.table("/public1/users/aimei/cFF_revised/1.TF_bindings/W3L_tf_bindings.txt",stringsAsFactor=F)
ltfs=read.table("/public1/users/aimei/cFF_revised/1.TF_bindings/SupplementaryTable2-W3LChIP.tsv",stringsAsFactor=F,sep=",",quote="",header=T)
ltfs=unique(ltfs[,2])
ltpm_s = etpm[ltfs,]
rownames(ltpm_s)=fbgn2symbol[rownames(ltpm_s),2]
write.table(ltpm_s,"W3L_TFs.TPM",sep="\t",row.names=T,col.names=T,quote=F)

############## LWQ' heatshock data
expr = read.table("~/cFF/expression_data/lanwq_stress/hs_stress.conf.tpm",stringsAsFactor=F)
lexpr_s = expr[ltfs,]
filtered_tfs = rownames(lexpr_s)[rowSums(lexpr_s>=10)>=3]
write.table(filtered_tfs,"filtered_L3_TFs.txt",row.names=F,col.names=F,sep="\t",quote=F)

rownames(lexpr_s)=fbgn2symbol[rownames(lexpr_s),2]
write.table(lexpr_s,"lwq_L3_larvaTFs.TPM",sep="\t",row.names=T,col.names=T,quote=F)


eexpr_s = expr[etfs,]
rownames(eexpr_s)=fbgn2symbol[rownames(eexpr_s),2]
write.table(eexpr_s,"lwq_L3_embryoTFs.TPM",sep="\t",row.names=T,col.names=T,quote=F)


############3

