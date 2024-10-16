gtf = read.table("/public1/users/aimei/genomes/Drosophila_melanogaster/BDGP6_r35/protein_coding.gtf",stringsAsFactors=F,sep="\t",quote="")
tmp = regexpr("gene_id \"\\S+?\";",gtf[,9])
gtf$gene_id = substr(gtf[,9],tmp+9,tmp+attr(tmp,"match.length")-3)
proteins=unique(gtf$gene_id)

tfb = read.table("Embryo_tf_bindings.txt",stringsAsFactors=F,sep="\t",quote="")
etfs = read.table("../2.miR_TF_expression/TF_expression/filtered_embryonic_TFs.txt", stringsAsFactors = F)

tfb_filtered = subset(tfb, V2 %in% etfs[,1]&V3 %in% proteins)
length(unique(tfb_filtered[,2]))
#341
length(unique(tfb_filtered[,3]))
#13660
dim(unique(tfb_filtered))
#3937911

#setdiff(etfs[,1],tfb[,2])



tfb = read.table("W3L_tf_bindings.txt",stringsAsFactors=F,sep="\t",quote="")
etfs = read.table("../2.miR_TF_expression/TF_expression/filtered_L3_TFs.txt", stringsAsFactors = F)

tfb_filtered = subset(tfb, V2 %in% etfs[,1]&V3 %in% proteins)
length(unique(tfb_filtered[,2]))
#31
length(unique(tfb_filtered[,3]))
#13579
dim(unique(tfb_filtered))
#324352
