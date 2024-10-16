library(ggplot2)
library(reshape)
options(stringsAsFactors=F)
require(RColorBrewer)

####################################
theta_pi=function(x,n){
	# n is number of samples, x is the frequency specturm with format:
	# size   count
	result=sum(2*x[,1]*(n-x[,1])/(n*(n-1))*x[,2]);
	return(result);
}
	
tajima_D=function(x,n){
	# n is number of samples, x is the frequency specturm with format:
	# size   count
	a1=sum(1/(1:(n-1)));
	a2=sum(1/((1:(n-1))^2));
	e1=1/a1*((n+1)/(3*(n-1))-1/a1);
	e2=1/(a1^2+a2)*(2*(n^2+n+3)/(9*n*(n-1))-(n+2)/(n*a1)+a2/a1^2);
	s=sum(x[,2]);
	tpi=theta_pi(x,n);
	tw=s/a1;
	varpi=e1*s+e2*s*(s-1);
	result=(tpi-tw)/sqrt(varpi);
	return(result);
}

theta_W=function(x,n){
	result=sum(x[,2])/sum(1/(1:(n-1)));
	return(result);
}


boostrapping <- function(x){
	#browser()
	n = length(x)*0.85
	ix = sample(1:length(x),n)
	return(x[ix])
}

###############################
ff = subset(read.table("l3_consmiR.lowThreshold.UTRcoor_sorted",header=F),V1!="chr4")
ff = read.table("wholeGenomeDivPolDmel/dme.ff-r6.29.site",header=T)
nrow(subset(ff,chr!="4"))
ff = read.table("wholeGenomeDivPolDmel/dme_ff.div",header=T)
nrow(subset(ff,chr!="chr4"))
ff = read.table("wholeGenomeDivPolDmel/dme_ff.pol",header=T)
nrow(ff)
ff1 = subset(ff, AF!=0&AF!=1&AC/AF>=30&CHR!="chr4"&ANC!="u")
nrow(ff1)
nrow(subset(ff1,MLAC>2))

#x = matrix(c(324004,291633,190,219,10294,8888,1306,1009,4046,3472,8030,7753),nrow=2,byrow=F)
# W3L
x = matrix(c(324004,291633,6,5,3319,2554,887,697,8095,6809,12493,12116),nrow=2,byrow=F)

# Embyonic
x = matrix(c(324004,291633,306,286,9253,8046,988,756,6030,5093,9030,8696),nrow=2,byrow=F)
for(i in 2:ncol(x)){
	m = cbind(x[,1],x[,i])
	print(fisher.test(m)$p.val)
}

#########################
ff = read.table("wholeGenomeDivPolDmel/dme_ff.pol",header=T)
ff1 = subset(ff, AF!=0&AF!=1&AC/AF>=30&CHR!="chr4"&ANC!="u")
td_ff = vector(length=1000)
for (i in 1:length(td_ff)){
	if(i%%100==1){print(i)}
	ix = sample(1:nrow(ff1),nrow(ff1)*0.85)
	ff_x = ff1[ix,]
	freq=melt(table(ff_x$AC))
	td_ff[i]=tajima_D(freq,197*2)
}

daf = read.table("l3_consmiR.lowThreshold.redundancy_DAF.tab")

freq = tapply(daf[,4],daf[,3],function(x){
	#browser()
	tds = vector(length=1000)
	for (i in 1:length(td_ff)){
		ix = sample(1:length(x),length(x)*0.85)
		ff_x = x[ix]
		freq=melt(table(ff_x))
		tds[i]=tajima_D(freq,197*2)
	}
	return(tds)
	})

cols <- c("lightslategrey","grey",rev(brewer.pal(9, "Set1")[c(1,3,2)]),"black")
names(cols) =c("fourfold","miR..non.TF","non.redundant","low.redundant","high.redundant","miR..TF")
freq_gg=melt(data.frame(do.call("cbind",freq),fourfold=td_ff))
freq_gg$variable=factor(freq_gg$variable,levels=c("fourfold","miR..TF","miR..non.TF","non.redundant","low.redundant","high.redundant"))
p <- ggplot(freq_gg) + geom_boxplot(aes(x=variable,y=value,color=variable))+ theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="none", legend.title=element_blank()) + labs(title="",y="Tajima's D", x="") +scale_color_manual(values=cols) + ylim(-2,0)
ggsave("boxplot_tajimaD_byAG.pdf",width=1.8,height=2.9)


comp = combn(c("fourfold","miR..non.TF","non.redundant","low.redundant","high.redundant","miR..TF"),2)
for(i in 1:ncol(comp)){
	x1 = subset(freq_gg,variable==comp[1,i])[,2]
	x2 = subset(freq_gg,variable==comp[2,i])[,2]
	print(wilcox.test(x1,x2)$p.value)
}

data = subset(freq_gg,variable=="non.redundant"|variable=="low.redundant"|variable=="high.redundant")

summary(aov(value~variable,data))
fit = glm(value~variable,data=data)
summary(fit)
# daf = read.table("TF_TG_DAF.tab_by_AG")

# freq = tapply(daf[,4],daf[,3],function(x){
# 	#browser()
# 	tds = vector(length=1000)
# 	for (i in 1:length(td_ff)){
# 		ix = sample(1:length(x),length(x)*0.85)
# 		ff_x = x[ix]
# 		freq=melt(table(ff_x))
# 		tds[i]=tajima_D(freq,197*2)
# 	}
# 	return(tds)
# 	})
# freq_gg=melt(data.frame(do.call("cbind",freq),fourfold=td_ff))
# freq_gg$variable=factor(freq_gg$variable,levels=c("fourfold","miR..TF","miR..non.TF"))
# p <- ggplot(freq_gg) + geom_boxplot(aes(x=variable,y=value,color=variable))+ theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="none", legend.title=element_blank()) + labs(title="",y="Tajima's D", x="") +scale_color_manual(values=c("darkgrey","black","darkred")) + ylim(-2,0)
# ggsave("boxplot_TF_TG_tajimaD_byAG.pdf",width=1.5,height=3)

