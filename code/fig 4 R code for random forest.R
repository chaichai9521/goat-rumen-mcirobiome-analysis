###################################################################################################################
##re-analysis
###3 groups together RF
a2=read.delim("otu_table_new3.csv", header=T, sep=",")
data4=a2[,c(1,24:861)]  
library(randomForest);library(pROC)
library(dplyr);library(AUCRF)
data4$trt=factor(data4$trt)
levels(data4$trt)=c(1:length(levels(data4$trt))-1)
set.seed(1)
  rf_testset_trt4 <- AUCRF(trt~., data=data4, ntree=10000, pdel=0.05, ranking="MDA")
   aucrf_test_trt4 <- AUCRFcv(rf_testset_trt4, nCV=10, M=20)
  test_held_trt4 <- predict(aucrf_test_trt4$RFopt, type='prob')[,2]
  trt_roc4 <- roc(data4$trt ~ test_held_trt4)


pdf(file='AUC-3goups.pdf', width=4, height=3)
layout(matrix(c(1,
                1), 
              nrow=1, byrow = TRUE))
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
plot(trt_roc4, col='blue', add=T, lty=1,print.thres.col="blue",print.thres=T)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
legend('bottomright', legend=c(
  sprintf('MRO-MRC-MCA, AUC = 1.00')
  
),lty=1, lwd = 2, cex=0.7, col=c( 'blue'), bty='n')
dev.off()

pdf(file='AUC-3goups-33.pdf', width=4, height=3)
layout(matrix(c(1,
                1), 
              nrow=1, byrow = TRUE))
par(mar=c(4,4,1,1))

##my code=above function
  p1=max(aucrf_test_trt4$AUCcurve$AUC);p2=(aucrf_test_trt4$AUCcurve$AUC)
  r=p2-p1;ylim=c(max(0,p2-r),min(1,p1+r))
  plot(aucrf_test_trt4$AUCcurve$AUC~aucrf_test_trt4$AUCcurve$k,type='o',col='4',pch=20,ylim=c(0,1.1),
       ylab='OOB-AUC',xlab='Number of selected variables')
  text(aucrf_test_trt4$Kopt,aucrf_test_trt4$"OOB-AUCopt","|", pos=3,col='4')
  text(aucrf_test_trt4$Kopt,aucrf_test_trt4$"OOB-AUCopt"+1.8*(ylim[1]-ylim[2])/10, paste("OOB-AUCopt = ",round(aucrf_test_trt4$"OOB-AUCopt",digits=3)," (Kopt = ",aucrf_test_trt4$Kopt,")",sep=""),
       pos=4,col='4', cex=0.7, offset=0)
  text(aucrf_test_trt4$Kopt,aucrf_test_trt4$"OOB-AUCopt"+2.8*(ylim[1]-ylim[2])/10, paste("cvAUC = ",round(aucrf_test_trt4$cvAUC,digits=3),sep=""),  col='4',pos=4, cex=0.7, offset=0)
  dev.off()

plot(aucrf_test_trt4,which=c("ranking"),maxvars=25)

imp <- sort(aucrf_test_trt4$ranking,decreasing=T)[1:25]
imp <- sort(aucrf_test_trt4$ranking,decreasing=T)[26:50]
specific=subset(data4,(colnames(data4) %in% names(imp)))
 Oname <- array(1:25,c(25,4))
 for(i in 1:25){
  Oname[i,2] <- which(colnames(data4)==names(imp[i]))
 Oname[i,3] <- names(imp[i])
 Oname[i,4] <- imp[i]
 }
 d1 <-data4[,c(1,as.integer(Oname[,2]))]
levels(d1$trt)=c('MCA','MRC','MRO')
rn <- colnames(d1)
rn= as.character(rn)
for(i in 1:26){
   ww <- strsplit(rn, "[.]")[[i]]
if(length(ww)==7) rn[i] <- paste(ww[7],sep="_")
if(length(ww)==6) rn[i] <- paste(ww[6],sep="_")
if(length(ww)==5 )rn[i] <- paste(ww[5],sep="_")
   if(length(ww)==4) rn[i] <- paste(ww[4],sep="_")
   if(length(ww)==3) rn[i] <- paste(ww[3],sep="_")
 if(length(ww)==2) rn[i] <- paste(ww[2])
 if(length(ww)==1) rn[i] <- paste(ww[1])
 }
rn1=rn;rn1= as.character(rn1)
for(i in 1:26){
   ww <- strsplit(rn, "[__]")[[i]]
if(length(ww)==7) rn[i] <- paste(ww[7],sep="_")
if(length(ww)==6) rn[i] <- paste(ww[6],sep="_")
if(length(ww)==5 )rn[i] <- paste(ww[2],ww[3],ww[4],ww[5],sep="_")
   if(length(ww)==4) rn[i] <- paste(ww[4],sep="_")
   if(length(ww)==3) rn[i] <- paste(ww[3],sep="_")
 if(length(ww)==2) rn[i] <- paste(ww[2])
 if(length(ww)==1) rn[i] <- paste(ww[1])
 }
#rn1[10]='Subdivision5_genera_incertae_sedis'

colnames(d1)=rn1

d1$trt=factor(d1$trt, levels=c('MRO','MRC','MCA'))

pdf("boxpot-trt-AE-others-auc22223.pdf",width=8, height=10.5, paper="letter")

par(mfrow=c(5,3));par(mar=c(2,2,2,2))
for(i in 2:100){
boxplot(d1[,i]~d1$trt,main=rn[i],col=c('4','red'))
}
dev.off()


###ggplot generate mutiple figures with pvalue
library(reshape2) 
library(ggplot2);library(magrittr)
library(ggpubr)
d <- melt(d1, id.var=c("trt"))

ggplot(d,aes(x=factor(trt),y=value))+geom_boxplot(aes(fill=trt))+  labs(x="", y="")+
ylab("")+theme(axis.title.y = element_text(size = rel(1.2)),plot.title = element_text(hjust = .5))+border(color = "black", size = 0.8)+
facet_wrap(. ~ variable,scale="free",nrow=5)+ guides(fill=FALSE)+
  stat_compare_means(label = "p.format",label.y.npc='top')+ theme_bw()

############
#change color manual 
ggplot(d,aes(x=factor(trt),y=value))+geom_boxplot(aes(fill=trt))+ labs(x="", y="")+
ylab("")+theme(axis.title.y = element_text(size = rel(1.2)),plot.title = element_text(hjust = .5))+border(color = "black", size = 0.8)+
facet_wrap(. ~ variable,scale="free",nrow=5)+ guides(fill=FALSE)+
  stat_compare_means(label = "p.format",label.y.npc='top')+ theme_bw()+scale_colour_manual(values=c("chartreuse4",'orangered', "dodgerblue4"), aesthetics = "fill")


#compare_means(shannon~trt, data=d1)
#my_comparisons <- list(c("others", "AECOPD"))
#comparisons=my_comparisons,
# Add pairwise comparisons p-value from Wilcoxon test
  

#######################################
##multiple figures into multiple pages
gg1=ggplot(d,aes(x=factor(trt),y=value))+geom_boxplot(aes(fill=trt))+ labs(x="", y="")+
ylab("")+theme(axis.title.y = element_text(size = rel(1.2)),plot.title = element_text(hjust = .5))+border(color = "black", size = 0.8)+
#facet_wrap(. ~ variable,scale="free",nrow=5)+ 
guides(fill=FALSE)+
  stat_compare_means(label = "p.format",label.y.npc='top')+ theme_bw()
devtools::install_github("guiastrennec/ggplus")
library(ggplus)
pdf("need10.pdf")
gg10 <- facet_multiple(plot=gg1, facets="variable", ncol = 2, nrow = 5)
dev.off()

write.csv(imp,'imp.csv')
imp <- sort(aucrf_test_trt4$ranking,decreasing=T)[1:1398]
imp <- sort(aucrf_test_trt4$ranking,decreasing=T)[1:838]
#########################################################################################
##AUCRF-pair-wise

##MRO vs MRC
a2=read.delim("otu_table_new3.csv", header=T, sep=",")
a=subset(a2,trt=='MRO'|trt=='MRC')
data=a[,c(1,24:861)] 
data$trt=droplevels(data$trt)
##> levels(data$trt)
#"MRC" "MRO"
library(randomForest);library(pROC)
library(dplyr);library(AUCRF)
data$trt=factor(data$trt)
levels(data$trt)=c(1:length(levels(data$trt))-1)
set.seed(1)
  rf_testset_trt <- AUCRF(trt~., data=data, ntree=10000, pdel=0.05, ranking="MDA")
   aucrf_test_trt <- AUCRFcv(rf_testset_trt, nCV=10, M=20)
  test_held_trt <- predict(aucrf_test_trt$RFopt, type='prob')[,2]
  trt_roc <- roc(data$trt ~ test_held_trt)
pdf(file='AUC-MRO-MRC.pdf', width=4, height=3)
layout(matrix(c(1,
                1), 
              nrow=1, byrow = TRUE))
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
plot(trt_roc, col='blue', add=T, lty=1,print.thres.col="blue",print.thres=T)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
legend('bottomright', legend=c(
  sprintf('MRO-MRC, AUC = 1.00')
  
),lty=1, lwd = 2, cex=0.7, col=c( 'blue'), bty='n')
dev.off()

pdf(file='AUC-MRO-MRC-33.pdf', width=4, height=3)
layout(matrix(c(1,
                1), 
              nrow=1, byrow = TRUE))
par(mar=c(4,4,1,1))

##my code=above function
  p1=max(aucrf_test_trt$AUCcurve$AUC);p2=(aucrf_test_trt$AUCcurve$AUC)
  r=p2-p1;ylim=c(max(0,p2-r),min(1,p1+r))
  plot(aucrf_test_trt$AUCcurve$AUC~aucrf_test_trt$AUCcurve$k,type='o',col='4',pch=20,ylim=c(0,1.1),
       ylab='OOB-AUC',xlab='Number of selected variables')
  text(aucrf_test_trt$Kopt,aucrf_test_trt$"OOB-AUCopt","|", pos=3,col='4')
  text(aucrf_test_trt$Kopt,aucrf_test_trt$"OOB-AUCopt"+1.8*(ylim[1]-ylim[2])/10, paste("OOB-AUCopt = ",round(aucrf_test_trt$"OOB-AUCopt",digits=3)," (Kopt = ",aucrf_test_trt$Kopt,")",sep=""),
       pos=4,col='4', cex=0.7, offset=0)
  text(aucrf_test_trt$Kopt,aucrf_test_trt$"OOB-AUCopt"+2.8*(ylim[1]-ylim[2])/10, paste("cvAUC = ",round(aucrf_test_trt$cvAUC,digits=3),sep=""),  col='4',pos=4, cex=0.7, offset=0)
  dev.off()

plot(aucrf_test_trt,which=c("ranking"),maxvars=25)

imp <- sort(aucrf_test_trt$ranking,decreasing=T)[1:25]
specific=subset(data,(colnames(data) %in% names(imp)))
 Oname <- array(1:25,c(25,4))
 for(i in 1:25){
  Oname[i,2] <- which(colnames(data)==names(imp[i]))
 Oname[i,3] <- names(imp[i])
 Oname[i,4] <- imp[i]
 }
 d1 <-data[,c(1,as.integer(Oname[,2]))]
levels(d1$trt)=c('MRC','MRO')

d1$trt=factor(d1$trt, levels=c('MRO','MRC'))
library(reshape2) 
library(ggplot2);library(magrittr)
library(ggpubr)
d <- melt(d1, id.var=c("trt"))
ggplot(d,aes(x=factor(trt),y=value))+geom_boxplot(aes(fill=trt))+ labs(x="", y="")+
ylab("")+theme(axis.title.y = element_text(size = rel(1.2)),plot.title = element_text(hjust = .5))+border(color = "black", size = 0.8)+
facet_wrap(. ~ variable,scale="free",nrow=5)+ guides(fill=FALSE)+
  stat_compare_means(label = "p.format",label.y.npc='top')+ theme_bw()+scale_colour_manual(values=c("chartreuse4",'orangered'), aesthetics = "fill")


##MRO vs MCA
a2=read.delim("otu_table_new3.csv", header=T, sep=",")
a=subset(a2,trt=='MRO'|trt=='MCA')
data=a[,c(1,24:861)] 
data$trt=droplevels(data$trt)
##> levels(data$trt)
#"MCA" "MRO"
library(randomForest);library(pROC)
library(dplyr);library(AUCRF)
data$trt=factor(data$trt)
levels(data$trt)=c(1:length(levels(data$trt))-1)
set.seed(1)
  rf_testset_trt <- AUCRF(trt~., data=data, ntree=10000, pdel=0.05, ranking="MDA")
   aucrf_test_trt <- AUCRFcv(rf_testset_trt, nCV=10, M=20)
  test_held_trt <- predict(aucrf_test_trt$RFopt, type='prob')[,2]
  trt_roc <- roc(data$trt ~ test_held_trt)
pdf(file='AUC-MRO-MCA.pdf', width=4, height=3)
layout(matrix(c(1,
                1), 
              nrow=1, byrow = TRUE))
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
plot(trt_roc, col='blue', add=T, lty=1,print.thres.col="blue",print.thres=T)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
legend('bottomright', legend=c(
  sprintf('MRO-MCA, AUC = 1.00')
  
),lty=1, lwd = 2, cex=0.7, col=c( 'blue'), bty='n')
dev.off()

pdf(file='AUC-MRO-MCA-33.pdf', width=4, height=3)
layout(matrix(c(1,
                1), 
              nrow=1, byrow = TRUE))
par(mar=c(4,4,1,1))

##my code=above function
  p1=max(aucrf_test_trt$AUCcurve$AUC);p2=(aucrf_test_trt$AUCcurve$AUC)
  r=p2-p1;ylim=c(max(0,p2-r),min(1,p1+r))
  plot(aucrf_test_trt$AUCcurve$AUC~aucrf_test_trt$AUCcurve$k,type='o',col='4',pch=20,ylim=c(0,1.1),
       ylab='OOB-AUC',xlab='Number of selected variables')
  text(aucrf_test_trt$Kopt,aucrf_test_trt$"OOB-AUCopt","|", pos=3,col='4')
  text(aucrf_test_trt$Kopt,aucrf_test_trt$"OOB-AUCopt"+1.8*(ylim[1]-ylim[2])/10, paste("OOB-AUCopt = ",round(aucrf_test_trt$"OOB-AUCopt",digits=3)," (Kopt = ",aucrf_test_trt$Kopt,")",sep=""),
       pos=4,col='4', cex=0.7, offset=0)
  text(aucrf_test_trt$Kopt,aucrf_test_trt$"OOB-AUCopt"+2.8*(ylim[1]-ylim[2])/10, paste("cvAUC = ",round(aucrf_test_trt$cvAUC,digits=3),sep=""),  col='4',pos=4, cex=0.7, offset=0)
  dev.off()

plot(aucrf_test_trt,which=c("ranking"),maxvars=25)

imp <- sort(aucrf_test_trt$ranking,decreasing=T)[1:25]
specific=subset(data,(colnames(data) %in% names(imp)))
 Oname <- array(1:25,c(25,4))
 for(i in 1:25){
  Oname[i,2] <- which(colnames(data)==names(imp[i]))
 Oname[i,3] <- names(imp[i])
 Oname[i,4] <- imp[i]
 }
 d1 <-data[,c(1,as.integer(Oname[,2]))]
levels(d1$trt)=c('MCA','MRO')
d1$trt=factor(d1$trt, levels=c('MRO','MCA'))
library(reshape2) 
library(ggplot2);library(magrittr)
library(ggpubr)
d <- melt(d1, id.var=c("trt"))
ggplot(d,aes(x=factor(trt),y=value))+geom_boxplot(aes(fill=trt))+ labs(x="", y="")+
ylab("")+theme(axis.title.y = element_text(size = rel(1.2)),plot.title = element_text(hjust = .5))+border(color = "black", size = 0.8)+
facet_wrap(. ~ variable,scale="free",nrow=5)+ guides(fill=FALSE)+
  stat_compare_means(label = "p.format",label.y.npc='top')+ theme_bw()+scale_colour_manual(values=c("chartreuse4",'dodgerblue4'), aesthetics = "fill")

##MRC vs MCA
a2=read.delim("otu_table_new3.csv", header=T, sep=",")
a=subset(a2,trt=='MRC'|trt=='MCA')
data=a[,c(1,24:861)] 
data$trt=droplevels(data$trt)
##> levels(data$trt)
#"MCA" "MRC"
library(randomForest);library(pROC)
library(dplyr);library(AUCRF)
data$trt=factor(data$trt)
levels(data$trt)=c(1:length(levels(data$trt))-1)
set.seed(1)
  rf_testset_trt <- AUCRF(trt~., data=data, ntree=10000, pdel=0.05, ranking="MDA")
   aucrf_test_trt <- AUCRFcv(rf_testset_trt, nCV=10, M=20)
  test_held_trt <- predict(aucrf_test_trt$RFopt, type='prob')[,2]
  trt_roc <- roc(data$trt ~ test_held_trt)
pdf(file='AUC-MRC-MCA.pdf', width=4, height=3)
layout(matrix(c(1,
                1), 
              nrow=1, byrow = TRUE))
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
plot(trt_roc, col='blue', add=T, lty=1,print.thres.col="blue",print.thres=T)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
legend('bottomright', legend=c(
  sprintf('MRC-MCA, AUC = 1.00')
  
),lty=1, lwd = 2, cex=0.7, col=c( 'blue'), bty='n')
dev.off()

pdf(file='KOPT-MRC-MCA.pdf', width=4, height=3)
layout(matrix(c(1,
                1), 
              nrow=1, byrow = TRUE))
par(mar=c(4,4,1,1))

##my code=above function
  p1=max(aucrf_test_trt$AUCcurve$AUC);p2=(aucrf_test_trt$AUCcurve$AUC)
  r=p2-p1;ylim=c(max(0,p2-r),min(1,p1+r))
  plot(aucrf_test_trt$AUCcurve$AUC~aucrf_test_trt$AUCcurve$k,type='o',col='4',pch=20,ylim=c(0,1.1),
       ylab='OOB-AUC',xlab='Number of selected variables')
  text(aucrf_test_trt$Kopt,aucrf_test_trt$"OOB-AUCopt","|", pos=3,col='4')
  text(aucrf_test_trt$Kopt,aucrf_test_trt$"OOB-AUCopt"+1.8*(ylim[1]-ylim[2])/10, paste("OOB-AUCopt = ",round(aucrf_test_trt$"OOB-AUCopt",digits=3)," (Kopt = ",aucrf_test_trt$Kopt,")",sep=""),
       pos=4,col='4', cex=0.7, offset=0)
  text(aucrf_test_trt$Kopt,aucrf_test_trt$"OOB-AUCopt"+2.8*(ylim[1]-ylim[2])/10, paste("cvAUC = ",round(aucrf_test_trt$cvAUC,digits=3),sep=""),  col='4',pos=4, cex=0.7, offset=0)
  dev.off()

plot(aucrf_test_trt,which=c("ranking"),maxvars=25)

imp <- sort(aucrf_test_trt$ranking,decreasing=T)[1:25]
specific=subset(data,(colnames(data) %in% names(imp)))
 Oname <- array(1:25,c(25,4))
 for(i in 1:25){
  Oname[i,2] <- which(colnames(data)==names(imp[i]))
 Oname[i,3] <- names(imp[i])
 Oname[i,4] <- imp[i]
 }
 d1 <-data[,c(1,as.integer(Oname[,2]))]
levels(d1$trt)=c('MCA','MRC')

d1$trt=factor(d1$trt, levels=c('MRC','MCA'))
library(reshape2) 
library(ggplot2);library(magrittr)
library(ggpubr)
d <- melt(d1, id.var=c("trt"))
ggplot(d,aes(x=factor(trt),y=value))+geom_boxplot(aes(fill=trt))+ labs(x="", y="")+
ylab("")+theme(axis.title.y = element_text(size = rel(1.2)),plot.title = element_text(hjust = .5))+border(color = "black", size = 0.8)+
facet_wrap(. ~ variable,scale="free",nrow=5)+ guides(fill=FALSE)+
  stat_compare_means(label = "p.format",label.y.npc='top')+ theme_bw()+scale_colour_manual(values=c("orangered",'dodgerblue4'), aesthetics = "fill")