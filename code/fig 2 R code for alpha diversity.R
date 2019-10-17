a2=read.delim("alpha.csv", header=T, sep=",")
a2$trt=factor(a2$trt, levels=c('MRO','MRC','MCA'))
library(ggplot2);library(magrittr)
library(ggpubr)
compare_means(shannon~trt, data=a2)
my_comparisons <- list(c('MRO','MRC'), c('MRC','MCA'), c('MRO','MCA'))

ggboxplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Shannon Index",legend="none",y="shannon", add ="boxplot", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   #stat_compare_means(comparisons=my_comparisons, label.y = c(5.5, 6, 6.5,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y =7.5)+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))


my_comparisons <- list(c('MRO','MRC'),  c('MRO','MCA'))
ggboxplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Observed_species",legend="none",y="observed_species", add ="boxplot", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   stat_compare_means(comparisons=my_comparisons, label.y = c(375,  395,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means(label.y =75)+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))
#############
library(dunn.test);library(FSA);library("tibble")
kruskal.test(a2$observed_species~a2$trt)
sss=dunnTest(a2$observed_species~a2$trt,method="bonferroni",digits = 3)
ppp=sss$res
foo <- data.frame(do.call('rbind', strsplit(as.character(ppp$Comparison),' - ',fixed=TRUE)))
colnames(foo)=c('group2','group1')
ppp[,5:6]=foo[,1:2]
ppp[,4]=round(ppp[,4],3)
as_data_frame(ppp)

ggboxplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Observed_species",legend="none",y="observed_species", add ="boxplot", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
  stat_compare_means(label.y =75)+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))+
stat_pvalue_manual(    data = ppp, label = "P.adj",y.position = c(375,395,385))






