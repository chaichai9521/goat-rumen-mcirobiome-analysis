a2=read.delim("otu_table_for_VFA.csv", header=T, sep=",")
library(ggplot2);library(magrittr)
library(ggpubr)
a2$trt=factor(a2$trt, levels=c('MRO','MRC','MCA'))

my_comparisons <- list(c('MRO','MRC'), c('MRC','MCA'), c('MRO','MCA'))

ammonia=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Ammonia Nitrogen",legend="none",y="Ammonia_nitrogen", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   stat_compare_means(comparisons=my_comparisons, label.y = c(22, 15, 23,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means()+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))

MCP=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Microbial Protein",legend="none",y="MCP", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   stat_compare_means(comparisons=my_comparisons, label.y = c(1.4, 1.5, 1.6,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means()+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))

TVFA=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Total VFA",legend="none",y="TVFA", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   stat_compare_means(comparisons=my_comparisons, label.y = c(90, 101, 105,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means()+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))

acetate=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Acetate",legend="none",y="acetate", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   stat_compare_means(comparisons=my_comparisons, label.y = c(42, 50, 52,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means()+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))

propionate=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Propionate",legend="none",y="propionate", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   stat_compare_means(comparisons=my_comparisons, label.y = c(20, 27, 29,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means()+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))

butyrate=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Butyrate",legend="none",y="butyrate", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   stat_compare_means(comparisons=my_comparisons, label.y = c(15, 22, 23,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = c(18))+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))

ggarrange(ammonia,MCP,TVFA,acetate,propionate,butyrate,
          labels = c("A", "B",'C','D','E','F'),
          ncol = 2, nrow = 3) %>%
  ggexport(filename = "fig 1 VFA.pdf",width = 130, height = 960)

########################################################################
##using anova test
a2=read.delim("otu_table_for_VFA.csv", header=T, sep=",")
library(ggplot2);library(magrittr)
library(ggpubr)
a2$trt=factor(a2$trt, levels=c('MRO','MRC','MCA'))

library(rstatix)  # https://github.com/kassambara/rstatix
stat.test <- aov(Ammonia_nitrogen~ trt, data = a2) %>%
  tukey_hsd()
library(car)
leveneTest(butyrate~ trt, data = a2,center=median)

ammonia=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Ammonia Nitrogen",legend="none",y="Ammonia_nitrogen", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   #stat_compare_means(comparisons=my_comparisons, label.y = c(22, 15, 23,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova")+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))+stat_pvalue_manual(
    stat.test, label = "p.adj", 
    y.position = c(22, 20, 15))

stat.test <- aov(MCP~ trt, data = a2) %>%
  tukey_hsd()
MCP=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Microbial Protein",legend="none",y="MCP", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   #stat_compare_means(comparisons=my_comparisons, label.y = c(1.4, 1.5, 1.6,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova",label.y = 1.5)+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))+
stat_pvalue_manual(stat.test, label = "p.adj",y.position = c(1.1, 1.4, 1.25))

stat.test <- aov(TVFA~ trt, data = a2) %>%
  tukey_hsd()
TVFA=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Total VFA",legend="none",y="TVFA", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   #stat_compare_means(comparisons=my_comparisons, label.y = c(90, 101, 105,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova")+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))+
stat_pvalue_manual(stat.test, label = "p.adj",y.position = c(85, 105, 100))

stat.test <- aov(acetate~ trt, data = a2) %>%
  tukey_hsd()
acetate=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Acetate",legend="none",y="acetate", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   #stat_compare_means(comparisons=my_comparisons, label.y = c(42, 50, 52,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova")+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))+
stat_pvalue_manual(stat.test, label = "p.adj",y.position = c(42, 51, 48))

stat.test <- aov(propionate~ trt, data = a2) %>%
  tukey_hsd()
propionate=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Propionate",legend="none",y="propionate", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   #stat_compare_means(comparisons=my_comparisons, label.y = c(20, 27, 29,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova")+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))+
stat_pvalue_manual(stat.test, label = "p.adj",y.position = c(22, 28, 26))

stat.test <- aov(butyrate~ trt, data = a2) %>%
  tukey_hsd()
butyrate=ggbarplot(a2, cex=0.25,x="trt",xlab="",ylab="",title="Butyrate",legend="none",y="butyrate", add ="mean_se", error.plot=  "errorbar",fill="trt", color ="black",palette = c("chartreuse4",'orangered', "dodgerblue4"))+
   #stat_compare_means(comparisons=my_comparisons, label.y = c(22, 23, 15,4, 4.5, 5,6, 6.5, 7))+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova",label.y = c(18))+
 border(color = "black", size = 0.8)+
 theme(plot.title = element_text(hjust = .5))+
stat_pvalue_manual(stat.test, label = "p.adj",y.position = c(15, 20, 22))
require("magrittr")
ggarrange(ammonia,MCP,TVFA,acetate,propionate,butyrate,
          labels = c("A", "B",'C','D','E','F'),
          ncol = 2, nrow = 3) %>%
  ggexport(filename = "fig 1 VFA-new.pdf",width = 800, height = 1960)

