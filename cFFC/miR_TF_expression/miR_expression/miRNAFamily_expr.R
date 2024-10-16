setwd("/public1/users/aimei/cFF_revised/2.miR_TF_expression/miR_expression")
st = read.table("sample_table.txt",stringsAsFactor=F)

expr = read.table("mirdeep2_out/miRNAs_expressed_all_samples_now.csv",stringsAsFactor=F,header=T)
expr_t = do.call("rbind",by(expr[,40:74],expr[,1],colSums))

familyInfo = read.table("miR_Family_Info_modified.txt",stringsAsFactor=F,sep="\t",quote="",header=T)
miRs=intersect(expr[,1],familyInfo[,3])
familyInfo = familyInfo[familyInfo[,3]%in%miRs,]

#familyInfo = read.table("/public1/users/aimei/cFF/2.miR_targets/miR_Family_Info.txt",stringsAsFactor=F, sep="\t",quote="",header=T)
setdiff(familyInfo[,3],expr[,1])
setdiff(expr[,1],familyInfo[,3])

expr_fami=do.call("rbind",tapply(familyInfo[,3],familyInfo[,1],function(x){

	if(length(x)==1){return(expr_t[x,])}
		#browser()
	return(colSums(expr_t[x,]))
	}))
colnames(expr_fami)=substr(colnames(expr_fami),1,3)
#expr_fami = as.data.frame(expr_fami)

expr_m=do.call("cbind",by(st[,2],st[,3],function(x){
	if (length(x)==1){return(expr_fami[,x])}
		#browser()
	return(rowMeans(expr_fami[,x]))
	}))


expr_m = as.data.frame(expr_m)
familyInfo$ifConserved = ifelse(familyInfo[,5]==1,1,0)
consInfo = unique(familyInfo[,c(1,7)])
rownames(consInfo) = consInfo[,1]
table(consInfo[,2])

expr_m$consInfo = consInfo[rownames(expr_m),2]

expr_m$embLow=rowSums(expr_m[,1:2]>=10)>=2
expr_m$l3Low=expr_m[,9]>=10
expr_m$embMedian=rowSums(expr_m[,1:2]>=100)>=2
expr_m$l3Median=expr_m[,9]>=100
# expr_m$embHigh=rowSums(expr_m[,1:7]>=500)>=3
# expr_m$l3High=expr_m[,9]>=500

table(expr_m[,c(11,12)])
table(expr_m[,c(11,13)])
table(expr_m[,c(11,14)])
table(expr_m[,c(11,15)])
# table(expr_m[,c(11,16)])
# table(expr_m[,c(11,17)])

write.table(subset(expr_m,!is.na(expr_m$consInfo)),"miRFami.expr.txt",row.names=T,col.names=T,sep="\t",quote=F)

ix = (expr_m$consInfo==1)&!is.na(expr_m$consInfo)&(expr_m$embLow)
write.table(rownames(expr_m)[ix],"embryo_miRs.cons.low.txt",row.names=F, col.names=F, sep="\t",quote=F)
ix = (expr_m$consInfo==0)&!is.na(expr_m$consInfo)&(expr_m$embLow)
write.table(rownames(expr_m)[ix],"embryo_miRs.noncons.low.txt",row.names=F, col.names=F, sep="\t",quote=F)
ix = (expr_m$consInfo==1)&!is.na(expr_m$consInfo)&(expr_m$embMedian)
write.table(rownames(expr_m)[ix],"embryo_miRs.cons.high.txt",row.names=F, col.names=F, sep="\t",quote=F)
ix = (expr_m$consInfo==0)&!is.na(expr_m$consInfo)&(expr_m$embMedian)
write.table(rownames(expr_m)[ix],"embryo_miRs.noncons.high.txt",row.names=F, col.names=F, sep="\t",quote=F)

ix = (expr_m$consInfo==1)&!is.na(expr_m$consInfo)&(expr_m$l3Low)
write.table(rownames(expr_m)[ix],"l3_miRs.cons.low.txt",row.names=F, col.names=F, sep="\t",quote=F)
ix = (expr_m$consInfo==0)&!is.na(expr_m$consInfo)&(expr_m$l3Low)
write.table(rownames(expr_m)[ix],"l3_miRs.noncons.low.txt",row.names=F, col.names=F, sep="\t",quote=F)
ix = (expr_m$consInfo==1)&!is.na(expr_m$consInfo)&(expr_m$l3Median)
write.table(rownames(expr_m)[ix],"l3_miRs.cons.high.txt",row.names=F, col.names=F, sep="\t",quote=F)
ix = (expr_m$consInfo==0)&!is.na(expr_m$consInfo)&(expr_m$l3Median)
write.table(rownames(expr_m)[ix],"l3_miRs.noncons.high.txt",row.names=F, col.names=F, sep="\t",quote=F)





