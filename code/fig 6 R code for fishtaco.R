

require(FishTacoPlot); require(ggplot2); require(scales); require(grid)
library(devtools)
install_github("borenstein-lab/fishtaco-plot") 
##########################################################33

fishtaco_out_main_output_SCORE_wilcoxon_ASSESSMENT_multi_taxa

####L3 MRO--MCA
p <- MultiFunctionTaxaContributionPlots(input_dir="/data/MRO-MCA", input_prefix="fishtaco_out",
input_taxa_taxonomy="/data/MRO-MCA/taxnomy333.txt", input_permutation="single_taxa",sort_by="list", plot_type="bars",
input_function_filter_list=c("ABC_transporters",'Biosynthesis_of_unsaturated_fatty_acids','Carbohydrate_metabolism','Insulin_signaling_pathway','Lipid_metabolism','Selenocompound_metabolism'), 
min_cont_as_separate=0.025,add_predicted_da_markers=TRUE)

p <- p +guides(fill = guide_legend(ncol=4)) + ylab("Wilcoxon test statistic (W)") +
theme(plot.title=element_blank(), axis.title.x=element_text(size=12,colour="black",face="plain"),
axis.text.x=element_text(size=10,colour="black",face="plain"), axis.title.y=element_blank(),
axis.text.y=element_text(size=10,colour="black",face="plain"), axis.ticks.y=element_blank(),
axis.ticks.x=element_blank(), panel.grid.major.x = element_line(colour="light gray"), panel.grid.major.y = element_line(colour="light gray"),
panel.grid.minor.x = element_line(colour="light gray"), panel.grid.minor.y = element_line(colour="light gray"),
panel.background = element_rect(fill="transparent",colour=NA), panel.border = element_rect(fill="transparent",colour="black"),
legend.background=element_rect(colour="black"), legend.title=element_text(size=10), legend.text=element_text(size=8,face="plain"),
legend.key.size=unit(0.8,"line"), legend.margin=unit(0.1,"line"), legend.position="bottom")
p

ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]

ggsave("FishTaco_T2D.png", p, width=11.8, height=6.07, units="in")
####################################
#make plot of last pathway
p <- MultiFunctionTaxaContributionPlots(input_dir="/data/MRO-MCA", input_prefix="fishtaco_out",
input_taxa_taxonomy="taxonomy33.txt", sort_by="list", plot_type="bars",
input_function_filter_list=c("Arachidonic_acid_metabolism",'Carbohydrate_metabolism','Flavone_and_flavonol_biosynthesis','Lipid_metabolism','Renal_cell_carcinoma',
'Signal_transduction_mechanisms','Vibrio_cholerae_pathogenic_cycle','Meiosis___yeast','Transporters','Ethylbenzene_degradation','Bacterial_motility_proteins','Glycosyltransferases',
'Pentose_and_glucuronate_interconversions','Transcription_factors','Bisphenol_degradation','Glutathione_metabolism','Other_ion_coupled_transporters','Pathways_in_cancer','Plant_pathogen_interaction',
'Selenocompound_metabolism','Sulfur_metabolism','DNA_repair_and_recombination_proteins','Transcription_machinery'), add_predicted_da_markers=TRUE)

p <- p +guides(fill = guide_legend(ncol=4)) + ylab("Wilcoxon test statistic (W)") +
theme(plot.title=element_blank(), axis.title.x=element_text(size=12,colour="black",face="plain"),
axis.text.x=element_text(size=10,colour="black",face="plain"), axis.title.y=element_blank(),
axis.text.y=element_text(size=10,colour="black",face="plain"), axis.ticks.y=element_blank(),
axis.ticks.x=element_blank(), panel.grid.major.x = element_line(colour="light gray"), panel.grid.major.y = element_line(colour="light gray"),
panel.grid.minor.x = element_line(colour="light gray"), panel.grid.minor.y = element_line(colour="light gray"),
panel.background = element_rect(fill="transparent",colour=NA), panel.border = element_rect(fill="transparent",colour="black"),
legend.background=element_rect(colour="black"), legend.title=element_text(size=10), legend.text=element_text(size=8,face="plain"),
legend.key.size=unit(0.8,"line"), legend.margin=unit(0.1,"line"), legend.position="bottom")
p


#######################################################
####L3 MRO--Mrc
p <- MultiFunctionTaxaContributionPlots(input_dir="/data/MRO-MRC", input_prefix="fishtaco_out",
input_taxa_taxonomy="/data/MRO-MRC/taxnomy333.txt", input_permutation="single_taxa",sort_by="list", plot_type="bars",
input_function_filter_list=c("Transporters",'Biosynthesis_of_unsaturated_fatty_acids','Carbohydrate_digestion_and_absorption','Carbohydrate_metabolism','Insulin_signaling_pathway','Lipid_metabolism','Selenocompound_metabolism'), 
min_cont_as_separate=0.025,add_predicted_da_markers=TRUE)

p <- p +guides(fill = guide_legend(ncol=4)) + ylab("Wilcoxon test statistic (W)") +
theme(plot.title=element_blank(), axis.title.x=element_text(size=12,colour="black",face="plain"),
axis.text.x=element_text(size=10,colour="black",face="plain"), axis.title.y=element_blank(),
axis.text.y=element_text(size=10,colour="black",face="plain"), axis.ticks.y=element_blank(),
axis.ticks.x=element_blank(), panel.grid.major.x = element_line(colour="light gray"), panel.grid.major.y = element_line(colour="light gray"),
panel.grid.minor.x = element_line(colour="light gray"), panel.grid.minor.y = element_line(colour="light gray"),
panel.background = element_rect(fill="transparent",colour=NA), panel.border = element_rect(fill="transparent",colour="black"),
legend.background=element_rect(colour="black"), legend.title=element_text(size=10), legend.text=element_text(size=8,face="plain"),
legend.key.size=unit(0.8,"line"), legend.margin=unit(0.1,"line"), legend.position="bottom")
p
##########################
p <- MultiFunctionTaxaContributionPlots(input_dir="/data/MRO-MRC", input_prefix="fishtaco_out",
input_taxa_taxonomy="/data/MRO-MRC/taxnomy333.txt",input_permutation="single_taxa",
min_cont_as_separate=0.025,sort_by="predicted_da", plot_type="bars", add_predicted_da_markers=TRUE, add_original_da_markers=TRUE)
p = p + scale_fill_manual(values = c("yellow2",'red','pink','purple','steelblue3','darkgreen',
"cyan",'coral','darkred','chartreuse','darkolivegreen3','lightslateblue',"lightblue1",'mediumseagreen',
'mediumpurple1','yellow','maroon4','lightgoldenrodyellow','blue','gray21','darkslategray','gray74'))
p <- p + guides(fill = guide_legend(ncol=5)) + ylab("Wilcoxon test statistic (W)") +
theme(plot.title=element_blank(), axis.title.x=element_text(size=12,colour="black",face="plain"),
axis.text.x=element_text(size=10,colour="black",face="plain"), axis.title.y=element_blank(),
axis.text.y=element_text(size=10,colour="black",face="plain"), axis.ticks.y=element_blank(),
axis.ticks.x=element_blank(), panel.grid.major.x = element_line(colour="light gray"), panel.grid.major.y = element_line(colour="light gray"),
panel.grid.minor.x = element_line(colour="light gray"), panel.grid.minor.y = element_line(colour="light gray"),
panel.background = element_rect(fill="transparent",colour=NA), panel.border = element_rect(fill="transparent",colour="black"),
legend.background=element_rect(colour="black"), legend.title=element_text(size=10), legend.text=element_text(size=8,face="plain"),
legend.key.size=unit(0.9,"line"), legend.margin=unit(0.2,"line"), legend.position="bottom")
p
ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
ggsave("MRO-MRC-f.pDF", p, width=11.8, height=15, units="in")
show_only_diff_abun_taxa=TRUE,
##
p <- MultiFunctionTaxaContributionPlots(input_dir="/data/MRO-MCA", input_prefix="fishtaco_out",
input_taxa_taxonomy="/data/MRO-MCA/taxnomy333.txt",input_permutation="single_taxa",
min_cont_as_separate=0.025,sort_by="predicted_da", plot_type="bars", add_predicted_da_markers=TRUE, add_original_da_markers=TRUE)
p = p + scale_fill_manual(values = c("yellow2",'red','pink','purple','steelblue3','darkgreen',
"cyan",'coral','darkred','chartreuse','darkolivegreen3','lightslateblue',"lightblue1",'mediumseagreen',
'mediumpurple1','yellow','maroon4','lightgoldenrodyellow','blue','plum1','darkslategray','honeydew','gray0','gray47','gray100','gray74'))

p <- p + guides(fill = guide_legend(ncol=5)) + ylab("Wilcoxon test statistic (W)") +
theme(plot.title=element_blank(), axis.title.x=element_text(size=12,colour="black",face="plain"),
axis.text.x=element_text(size=10,colour="black",face="plain"), axis.title.y=element_blank(),
axis.text.y=element_text(size=10,colour="black",face="plain"), axis.ticks.y=element_blank(),
axis.ticks.x=element_blank(), panel.grid.major.x = element_line(colour="light gray"), panel.grid.major.y = element_line(colour="light gray"),
panel.grid.minor.x = element_line(colour="light gray"), panel.grid.minor.y = element_line(colour="light gray"),
panel.background = element_rect(fill="transparent",colour=NA), panel.border = element_rect(fill="transparent",colour="black"),
legend.background=element_rect(colour="black"), legend.title=element_text(size=10), legend.text=element_text(size=8,face="plain"),
legend.key.size=unit(0.9,"line"), legend.margin=unit(0.2,"line"), legend.position="bottom")
p
ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
ggsave("MRO-MRC-f.pDF", p, width=11.8, height=15, units="in")
#######################################################################################################################################

##########new  change color manually
p <- MultiFunctionTaxaContributionPlots(input_dir="/data/MRO-MCA", input_prefix="fishtaco_out",
input_taxa_taxonomy="/data/MRO-MCA/taxnomy333.txt", input_permutation="single_taxa",sort_by="list", plot_type="bars",
input_function_filter_list=c("ABC_transporters",'Biosynthesis_of_unsaturated_fatty_acids','Carbohydrate_metabolism','Insulin_signaling_pathway','Lipid_metabolism','Selenocompound_metabolism'), 
min_cont_as_separate=0.025,add_predicted_da_markers=TRUE)


##p = p + scale_fill_manual(values = fishtaco_palette)
p = p + scale_fill_manual(values = c("yellow2",'red','pink','purple','steelblue3','darkgreen',
"cyan",'coral','darkred','chartreuse','deeppink','khaki1',"lightblue1",'mediumseagreen',
'mediumpurple1','lightslateblue','navy','gray74'))


p <- p +guides(fill = guide_legend(ncol=4)) + ylab("Wilcoxon test statistic (W)") +
theme(plot.title=element_blank(), axis.title.x=element_text(size=12,colour="black",face="plain"),
axis.text.x=element_text(size=10,colour="black",face="plain"), axis.title.y=element_blank(),
axis.text.y=element_text(size=10,colour="black",face="plain"), axis.ticks.y=element_blank(),
axis.ticks.x=element_blank(), panel.grid.major.x = element_line(colour="light gray"), panel.grid.major.y = element_line(colour="light gray"),
panel.grid.minor.x = element_line(colour="light gray"), panel.grid.minor.y = element_line(colour="light gray"),
panel.background = element_rect(fill="transparent",colour=NA), panel.border = element_rect(fill="transparent",colour="black"),
legend.background=element_rect(colour="black"), legend.title=element_text(size=10), legend.text=element_text(size=8,face="plain"),
legend.key.size=unit(0.8,"line"), legend.margin=unit(0.1,"line"), legend.position="bottom")
p



#########
p <- MultiFunctionTaxaContributionPlots(input_dir="/data/MRO-MRC", input_prefix="fishtaco_out",
input_taxa_taxonomy="/data/MRO-MRC/taxnomy333.txt", input_permutation="single_taxa",sort_by="list", plot_type="bars",
input_function_filter_list=c("Transporters",'Biosynthesis_of_unsaturated_fatty_acids','Carbohydrate_digestion_and_absorption','Carbohydrate_metabolism','Insulin_signaling_pathway','Lipid_metabolism','Selenocompound_metabolism'), 
min_cont_as_separate=0.025,add_predicted_da_markers=TRUE)


##p = p + scale_fill_manual(values = fishtaco_palette)
p = p + scale_fill_manual(values = c("yellow2",'red','pink','purple','steelblue3','darkgreen',
"cyan",'coral','darkred','chartreuse','khaki1','lightslateblue',"lightblue1",'gray74','mediumseagreen',
'mediumpurple1','deeppink4','maroon4'))


p <- p +guides(fill = guide_legend(ncol=4)) + ylab("Wilcoxon test statistic (W)") +
theme(plot.title=element_blank(), axis.title.x=element_text(size=12,colour="black",face="plain"),
axis.text.x=element_text(size=10,colour="black",face="plain"), axis.title.y=element_blank(),
axis.text.y=element_text(size=10,colour="black",face="plain"), axis.ticks.y=element_blank(),
axis.ticks.x=element_blank(), panel.grid.major.x = element_line(colour="light gray"), panel.grid.major.y = element_line(colour="light gray"),
panel.grid.minor.x = element_line(colour="light gray"), panel.grid.minor.y = element_line(colour="light gray"),
panel.background = element_rect(fill="transparent",colour=NA), panel.border = element_rect(fill="transparent",colour="black"),
legend.background=element_rect(colour="black"), legend.title=element_text(size=10), legend.text=element_text(size=8,face="plain"),
legend.key.size=unit(0.8,"line"), legend.margin=unit(0.1,"line"), legend.position="bottom")
p
