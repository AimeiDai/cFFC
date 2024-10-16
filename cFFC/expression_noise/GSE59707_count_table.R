
ortholog = read.table("dmel_orthologs_in_drosophila_species_fb_2022_01.tsv",sep="\t",quote="",stringsAsFactors=F)
dsi_ort = ortholog[grep("Dsim", ortholog$V7),]
x1 = table(dsi_ort[,1])
x2 = table(dsi_ort[,6])
one2one = dsi_ort[dsi_ort[,1] %in% names(x1[x1==1]) & dsi_ort[,6] %in% names(x2[x2==1]),c(1,6)]
rownames(one2one) = one2one[,1]
colnames(one2one) = c("dm_id","ds_id")

st = read.table("GSE59707_sample_table.txt",stringsAsFactors=F,sep="\t",header=T)

tpm = data.frame(row.names=one2one[,1])
for ( i in 1:nrow(st)){
	tmp = read.table(paste(st[i,"path"],"/",st[i,2],"/gene.abundance",sep=""),sep="\t",quote="",header=T)
	x = tapply(tmp[,"TPM"],tmp[,"Gene.ID"],sum)
	if (st[i,3] == "dme"){
		tpm[,st[i,2]] = x[rownames(tpm)]
	}else{
		tpm[,st[i,2]] = x[one2one[rownames(tpm),2]]
	}

}

tpm = tpm[-unique(which(is.na(tpm),arr.ind=T)[,1]),]

write.table(tpm,"GSE59707_dm_ds.tpm",sep="\t",row.names=T,col.names=T,quote=F)

###################
st = read.table("AMO_ctrl.sample_table.txt",stringsAsFactors=F,sep="\t",header=T)

tpm = data.frame(row.names=one2one[,1])
for ( i in 1:nrow(st)){
	tmp = read.table(paste(st[i,"path"],st[i,1],"/gene.abundance",sep=""),sep="\t",quote="",header=T)
	x = tapply(tmp[,"TPM"],tmp[,"Gene.ID"],sum)
	if (st[i,3] == "dme"){
		tpm[,st[i,1]] = x[rownames(tpm)]
	}else{
		tpm[,st[i,1]] = x[one2one[rownames(tpm),2]]
	}

}

tpm = tpm[-unique(which(is.na(tpm),arr.ind=T)[,1]),]
write.table(tpm,"AMO_ctrl_dm_ds.tpm",sep="\t",row.names=T,col.names=T,quote=F)

####################


