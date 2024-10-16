library(ggplot2)
library(reshape)
options(stringsAsFactors=F)
require(RColorBrewer)

cols2 <- c("grey","darkgrey","black")
names(cols2) =c("fourfold","non.cFFC.sites","cFFC.sites")

cols <- rev(brewer.pal(9, "Set1")[c(1,3,2)])
names(cols) = c("non.redundant","low.redundant","high.redundant")

setwd("C:/Users/Administrator/OneDrive/1.cFF/manuscript/2024-3-18-first_revision/polymorphism")
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
#########################
ff = read.table("dme_ff.pol",header=T)
ff1 = subset(ff, AF!=0&AF!=1&AC/AF>=30&CHR!="chr4"&ANC!="u")
td_ff = vector(length=1000)
for (i in 1:length(td_ff)){
	if(i%%100==1){print(i)}
	ix = sample(1:nrow(ff1),nrow(ff1)*0.85)
	ff_x = ff1[ix,]
	freq=melt(table(ff_x$AC))
	td_ff[i]=tajima_D(freq,197*2)
}

########## W3L, conserved miRNAs, low threshold
daf = read.table("l3_nonconsmiR.lowThreshold.redundancy_DAF.tab")
daf$type = ifelse(daf$V3=="high-redundant"|daf$V3=="low-redundant"|daf$V3=="non-redundant","cFFC sites","non-cFFC sites")
freq = tapply(daf[,5],daf[,6],function(x){
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

freq_gg1=melt(data.frame(do.call("cbind",freq),fourfold=td_ff))
freq_gg1$variable=factor(freq_gg1$variable,levels=c("fourfold","non.cFFC.sites","cFFC.sites"))

freq = tapply(daf[,5],daf[,3],function(x){
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

freq_gg2=melt(data.frame(do.call("cbind",freq)))
freq_gg2 = subset(freq_gg2,variable%in%c("non.redundant","low.redundant","high.redundant"))
freq_gg2$variable=factor(freq_gg2$variable,levels=c("non.redundant","low.redundant","high.redundant"))

freq_gg = rbind(data.frame(freq_gg1,type="non.cFFC.sites"),data.frame(freq_gg2,type="cFFC.sites"))
freq_gg$type = factor(freq_gg$type,levels=c("non.cFFC.sites","cFFC.sites"))
p <- ggplot(freq_gg) + geom_boxplot(aes(x=variable,y=value,color=variable))+ facet_wrap(~type,scales="free_x") + theme_classic(12) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="none", legend.title=element_blank()) + labs(title="",y="Tajima's D", x="") +scale_color_manual(values=c(cols2,cols)) + ylim(-2,0)
ggsave("tajimaD_l3_nonconsmiR.lowTh_2.pdf",width=2.8,height=3.4)

by(freq_gg[,2],freq_gg[,c(1,3)],median)

ix = combn(c("fourfold","non.cFFC.sites","cFFC.sites","non.redundant","low.redundant","high.redundant"),2)
for(i in 1:ncol(ix)){
  x1 = as.numeric(subset(freq_gg,variable==ix[1,i])$value)
  x2 = as.numeric(subset(freq_gg,variable==ix[2,i])$value)
  print(ix[,i])
  p=wilcox.test(x1,x2)$p.value
  print(ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","NS"))))
}


data = subset(freq_gg,variable=="non.redundant"|variable=="low.redundant"|variable=="high.redundant")

summary(aov(value~variable,data))
# fit = glm(value~variable,data=data)
# summary(fit)
#################### site type
cff_sits = subset(daf,type=="cFFC sites")
freq = tapply(cff_sits[,5],cff_sits[,3:4],function(x){
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
freq_gg3=data.frame()
for(i in rownames(freq)){
  for (j in colnames(freq)){
    freq_gg3=rbind(freq_gg3,data.frame(value=unlist(freq[i,j]),cFFType=i,siteType=j,stringsAsFactors = F))
  }
}

freq_gg3$cFFType=factor(freq_gg3$cFFType,levels=c("non-redundant","low-redundant","high-redundant"))
cols3=cols
names(cols3) = c("non-redundant","low-redundant","high-redundant")
p <- ggplot(freq_gg3) + geom_boxplot(aes(x=cFFType,y=value,color=cFFType))+ facet_wrap(~siteType,scales="free_x") + theme_bw(12) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="none", legend.title=element_blank()) + labs(title="",y="Tajima's D", x="") +scale_color_manual(values=c(cols3)) + ylim(-2,0)
ggsave("tajimaD_l3_nonconsmiR.lowTh_siteType.pdf",width=3.4,height=3.4)
pa = combn(c("non-redundant","low-redundant","high-redundant"),2)
for (i in c("7mer-a1","7mer-m8","8mer")){
  print(i)
  for (j in 1:ncol(pa)){
    x1 = subset(freq_gg3,cFFType==pa[1,j]&siteType==i)
    x2 = subset(freq_gg3,cFFType==pa[2,j]&siteType==i)
    p = wilcox.test(x1$value,x2$value)$p.value
    print(ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<0.05,"*","ns"))))
  }
}
summary(aov(value~cFFType*siteType,freq_gg3))


########## embryonic, conserved miRNAs, low threshold
daf = read.table("Embryonic_nonconsmiR.lowThreshold.redundancy_DAF.tab")
daf$type = ifelse(daf$V3=="high-redundant"|daf$V3=="low-redundant"|daf$V3=="non-redundant","cFFC sites","non-cFFC sites")

freq = tapply(daf[,5],daf[,6],function(x){
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

freq_gg1=melt(data.frame(do.call("cbind",freq),fourfold=td_ff))
freq_gg1$variable=factor(freq_gg1$variable,levels=c("fourfold","non.cFFC.sites","cFFC.sites"))


freq = tapply(daf[,5],daf[,3],function(x){
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

freq_gg2=melt(data.frame(do.call("cbind",freq)))
freq_gg2 = subset(freq_gg2,variable%in%c("non.redundant","low.redundant","high.redundant"))
freq_gg2$variable=factor(freq_gg2$variable,levels=c("non.redundant","low.redundant","high.redundant"))

freq_gg = rbind(data.frame(freq_gg1,type="non-cFFC sites"),data.frame(freq_gg2,type="cFFC sites"))
freq_gg$type = factor(freq_gg$type,levels=c("non-cFFC sites","cFFC sites"))
p <- ggplot(freq_gg) + geom_boxplot(aes(x=variable,y=value,color=variable))+facet_wrap(~type,scales = "free_x")+ theme_classic(12) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="none", legend.title=element_blank()) + labs(title="",y="Tajima's D", x="") +scale_color_manual(values=c(cols2,cols)) + ylim(-2,0)
ggsave("tajimaD_embryonic_nonconsmiR.lowTh_2.pdf",width=2.8,height=3.4)

by(freq_gg[,2],freq_gg[,c(1,3)],median)

for(i in 1:ncol(ix)){
  x1 = as.numeric(subset(freq_gg,variable==ix[1,i])$value)
  x2 = as.numeric(subset(freq_gg,variable==ix[2,i])$value)
  print(ix[,i])
  p=wilcox.test(x1,x2)$p.value
  print(ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","NS"))))
}

data = subset(freq_gg,variable=="non.redundant"|variable=="low.redundant"|variable=="high.redundant")

summary(aov(value~variable,data))
#################### site type
cff_sits = subset(daf,type=="cFFC sites")
freq = tapply(cff_sits[,5],cff_sits[,3:4],function(x){
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
freq_gg3=data.frame()
for(i in rownames(freq)){
  for (j in colnames(freq)){
    freq_gg3=rbind(freq_gg3,data.frame(value=unlist(freq[i,j]),cFFType=i,siteType=j,stringsAsFactors = F))
  }
}

freq_gg3$cFFType=factor(freq_gg3$cFFType,levels=c("non-redundant","low-redundant","high-redundant"))
cols3=cols
names(cols3) = c("non-redundant","low-redundant","high-redundant")
p <- ggplot(freq_gg3) + geom_boxplot(aes(x=cFFType,y=value,color=cFFType))+ facet_wrap(~siteType,scales="free_x") + theme_bw(12) + theme(axis.text.x=element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="none", legend.title=element_blank()) + labs(title="",y="Tajima's D", x="") +scale_color_manual(values=c(cols3)) + ylim(-2,0)
ggsave("tajimaD_embryo_nonconsmiR.lowTh_siteType.pdf",width=3.4,height=3.4)
pa = combn(c("non-redundant","low-redundant","high-redundant"),2)
for (i in c("7mer-a1","7mer-m8","8mer")){
  print(i)
  for (j in 1:ncol(pa)){
    x1 = subset(freq_gg3,cFFType==pa[1,j]&siteType==i)
    x2 = subset(freq_gg3,cFFType==pa[2,j]&siteType==i)
    p = wilcox.test(x1$value,x2$value)$p.value
    print(ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<0.05,"*","ns"))))
  }
}
summary(aov(value~cFFType*siteType,freq_gg3))

