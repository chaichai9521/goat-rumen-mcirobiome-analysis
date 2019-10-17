## stacked barplot
##generate pivot table in R
##combine OTU levels to genus level
##total 151 genera
###########################
otu=read.delim("otu_table_new.csv", header=T, sep=",")
otu1=otu[,c(2,24:ncol(otu))]

rn=colnames(otu1)[2:ncol(otu)]
rnx= as.character(rn)
rnx=rn ##rename taxnomy in standard level
for(i in 1:860){
   ww <- strsplit(rnx, "[_]")[[i]]
if(length(ww)==15) rn[i] <- paste(ww[1],ww[13],ww[14],ww[15],sep="_")
if(length(ww)==14) rn[i] <- paste(ww[1],ww[14],sep="_")
if(length(ww)==13) rn[i] <- paste(ww[1],ww[13],sep="_")
if(length(ww)==12) rn[i] <- paste(ww[1],ww[9],ww[10],ww[11],ww[12],sep="_")
if(length(ww)==11) rn[i] <- paste(ww[1],ww[11],'unclassified',sep="_")
if(length(ww)==9) rn[i] <- paste(ww[1],ww[9],'unclassified',sep="_")
if(length(ww)==7) rn[i] <- paste(ww[1],ww[7],'unclassified',sep="_")
if(length(ww)==5 )rn[i] <- paste(ww[1],ww[5],'unclassified',sep="_")
    if(length(ww)==3) rn[i] <- paste(ww[1],ww[3],'unclassified',sep="_")

 if(length(ww)==1) rn[i] <- paste(ww[1])
 }
write.csv(rnx,"rnx.csv")

otu=read.delim("otu_table_new2.csv", header=T, sep=",")
otu1=otu[,c(2,24:ncol(otu))]
write.csv(otu1,"otu.csv")

#########generate pivot table
gg=read.delim("otu.csv", header=T, sep=",")
ccc=  aggregate(gg[,-1],by=list(gg$group),FUN=sum)
write.csv(ccc,'genus2.csv')

hgg=read.delim("phy1.csv", header=T, skip=1,sep=",")
ccc=aggregate(hgg[,2:ncol(hgg)], by = list(hgg$group),FUN=sum)
write.csv(ccc,'phy2.csv')


##############
b=read.delim("genus2.csv", header=T, sep=",")
a=b[1:31,]
library(RColorBrewer)
a=as.matrix(a)
barplot(a[,-1], main = " ", ylim=c(0,1),xlab = "",ylab = "Relative abundace %",col=brewer.pal(12,"Paired"))
plot.new()
legend("center", ncol=2,bty = "n",legend=a[1:31,1],cex=0.8,fill=brewer.pal(12,"Paired"))

b=read.delim("genus2-avg.csv", header=T, sep=",")
a=b[1:31,]
library(RColorBrewer)
a=as.matrix(a)
barplot(a[,-1], main = " ", ylim=c(0,1),xlab = "",ylab = "Relative abundace %",col=brewer.pal(12,"Paired"))
##
a=read.delim("genus-litter.csv", skip=1,header=T, sep=",")
library(RColorBrewer)
a=as.matrix(a)
barplot(a[,-1], main = " ", ylim=c(0,1),xlab = "",ylab = "Relative abundace %",col=brewer.pal(12,"Paired"))
plot.new()
legend("center", ncol=2,bty = "n",legend=a[,1],cex=0.8,fill=brewer.pal(12,"Paired"))

a=read.delim("genus-litter-avg.csv", header=T, sep=",")
library(RColorBrewer)
a=as.matrix(a)
barplot(a[,-1], main = " ", ylim=c(0,1),xlab = "",ylab = "Relative abundace %",col=brewer.pal(12,"Paired"))

##
a=read.delim("phylum-ceca.csv", header=T,  skip=1,sep=",")
library(RColorBrewer)
a=as.matrix(a)
barplot(a[,-1], main = " ", ylim=c(0,1),xlab = "",ylab = "Relative abundace %",col=brewer.pal(12,"Paired"))
plot.new()
legend("center", ncol=1,bty = "n",legend=a[,1],cex=0.8,fill=brewer.pal(12,"Paired"))

a=read.delim("phylum-ceca-avg.csv", header=T, sep=",")
library(RColorBrewer)
a=as.matrix(a)
barplot(a[,-1], main = " ", ylim=c(0,1),xlab = "",ylab = "Relative abundace %",col=brewer.pal(12,"Paired"))


a=read.delim("phylum-litter.csv", header=T,  skip=1,sep=",")
library(RColorBrewer)
a=as.matrix(a)
barplot(a[,-1], main = " ", ylim=c(0,1),xlab = "",ylab = "Relative abundace %",col=brewer.pal(12,"Paired"))
plot.new()
legend("center", ncol=1,bty = "n",legend=a[,1],cex=0.8,fill=brewer.pal(12,"Paired"))

a=read.delim("phylum-litter-avg.csv", header=T, sep=",")
library(RColorBrewer)
a=as.matrix(a)
barplot(a[,-1], main = " ", ylim=c(0,1),xlab = "",ylab = "Relative abundace %",col=brewer.pal(12,"Paired"))

##
a=read.delim("otu_barplot-avg.csv", header=T, sep=",")
library(RColorBrewer)
a=as.matrix(a)
barplot(a[,-1], main = " ", ylim=c(0,1),xlab = "",ylab = "Relative abundace %",col=brewer.pal(12,"Paired"))
plot.new()
legend("center", ncol=1,bty = "n",legend=a[,1],cex=0.8,fill=brewer.pal(12,"Paired"))



