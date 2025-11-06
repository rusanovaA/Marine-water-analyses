library(phyloseq)
library(ggplot2)
library(ggrepel)
library(dplyr)
library("readxl")
library(pairwiseAdonis)
library(RColorBrewer)
library(vegan) #PERNANOVA
library(ggpubr)
library(microbiome) #
library(ranacapa) #ggrare
library(MicrobeR) #PCoA 3D
library(umap)
library(dendextend) #dendrogram

results_dir_ad <- "C:/Your/path/Additional/"
###Creating Phyloseq
seqtab.nochim_OTU <- read.table(sep = "\t", file = "C:/Your/path/all_OTU_frequency.tsv", header = TRUE, check.names=FALSE, row.names=1)
taxa_OTU <- read.table(sep = "\t", file = "C:/Your/path/all_OTUs_phylogeny.tsv", header = TRUE, check.names=FALSE, row.names=1)
taxa_OTU <- as.matrix(taxa_OTU)
sample.data <- read_excel("C:/Your/path/Supplementary_Table_1_St_Anna_Trough_metadata_0m_50m.xlsx")
sample.data <- sample.data %>% tibble::column_to_rownames("Sample_ID")
water <- phyloseq(otu_table(seqtab.nochim_OTU, taxa_are_rows = TRUE), sample_data(sample.data), tax_table(taxa_OTU))
water

###Find and samples with Fam Pseudoalteromonadaceae more than 10% (filtration)
water_fam <- tax_glom(water, "Family")
water_fam_rel <- transform_sample_counts(water_fam, function(x) 100 * x / sum(x))
tbl <- psmelt(water_fam_rel)
pseudo_percent <- tbl %>%
  filter(Family == "Pseudoalteromonadaceae") %>%
  select(Sample, Abundance)
pseudo_label <- ifelse(pseudo_percent$Abundance > 10, "High", "Low")
names(pseudo_label) <- pseudo_percent$Sample
sample_data(water)$Pseudoaltero_Status <- NA
samples_in_ph <- sample_names(water)
sample_data(water)[samples_in_ph, "Pseudoaltero_Status"] <- pseudo_label[samples_in_ph]
water

###Rarefaction curves
rare_water <- ggrare(water, step = 100, color = "Pseudoaltero_Status", se = FALSE, plot = TRUE) +
  theme_minimal()
rare_water <- rare_water + scale_color_manual(values = c("indianred2", "deepskyblue3"))
print(rare_water)
ggsave(file.path(results_dir_ad, "rare_water_wo_filtering.svg"), plot = rare_water, device = 'svg',dpi = 300,
       width = 5000, height = 1000, units = 'px')
ggsave(file.path(results_dir_ad, "rare_water_wo_filtering.pdf"), plot = rare_water, device = 'pdf',dpi = 300,
       width = 1500, height = 1000, units = 'px')
#ggsave(file.path(resuils_dir, "rare_water.png"), plot = rare_water, device = 'png',dpi = 300,
       # width = 1500, height = 1000, units = 'px')

###Filtering taxa
water <- subset_taxa(water, Kingdom %in% c("Bacteria", "Archaea"))
water <- subset_taxa(water, !(Family %in% c("Mitochondria")))
water

###Normalization and rare taxa removal
total = median(sample_sums(water))
standf = function(x, t=total) round(t * (x / sum(x)))
water = transform_sample_counts(water, standf)
water
total_ab_001 = 1e-3
reads_per_OTU_water <- taxa_sums(water)
keepTaxa_001 = ((reads_per_OTU_water / sum(reads_per_OTU_water))*100) > total_ab_001
water = prune_taxa(keepTaxa_001, water)
water

###Rarefaction curves
rare_water <- ggrare(water, step = 100, color = "Pseudoaltero_Status", se = FALSE, plot = TRUE) +
  theme_minimal()
rare_water <- rare_water + scale_color_manual(values = c("indianred2", "deepskyblue3"))
print(rare_water)
ggsave(file.path(results_dir_ad, "rare_water.svg"), plot = rare_water, device = 'svg',dpi = 300,
       width = 5000, height = 1000, units = 'px')
ggsave(file.path(results_dir_ad, "rare_water.pdf"), plot = rare_water, device = 'pdf',dpi = 300,
       width = 1500, height = 1000, units = 'px')

###Count fraction of contaminating OTUs - OTU_133 and OTU_196 -
###(Escherichia-Shigella and Pseudomonas, respectively) fractions in all samples
water_clean <- subset_samples(water, Sample_common_name != '6Z Kit control')
otu_mat <- as(otu_table(water_clean), "matrix")
if (!taxa_are_rows(water_clean)) {
  otu_mat <- t(otu_mat)
}
cont_otu_name <- c("OTU_133", "OTU_195")
cont_otu_abund <- colSums(otu_mat[cont_otu_name, ])
sample_sums <- colSums(otu_mat)
cont_otu_fraction <- (cont_otu_abund / sample_sums) * 100
sample_metadata <- as.data.frame(sample_data(water_clean))
sample_metadata_ordered <- sample_metadata[match(colnames(otu_mat), rownames(sample_metadata)), ]
result_df <- data.frame(
  SampleID = rownames(sample_metadata),
  Station = sample_metadata$Station,
  Cont_OTU_fraction = cont_otu_fraction,
  Sample_type = sample_metadata_ordered$Sample_type
)
write.table(result_df, file = "C:/Your/path/Additional/Contamination_fraction_per_sample_exp.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cont_otu_plot <- ggplot(result_df, aes(x = Sample_type, y = cont_otu_fraction)) + 
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5, aes(color = Sample_type)) +
  ylab("Contamination fraction, %") + xlab("Sample_type type") +
  theme_minimal() +
  ggtitle("Fraction of contaminating OTUs across samples")
cont_otu_plot
ggsave(file.path(results_dir_ad, "Contamination_plot.pdf"), plot = cont_otu_plot, device = 'pdf',dpi = 300,
       width = 1500, height = 1000, units = 'px')

###Barplot for all samples
#Aggregates data at the Genus level
GP.genus <- tax_glom(water, 'Genus')

#Selects top 40 most abundant genera
N <- 40
top_40 <- names(sort(taxa_sums(GP.genus), decreasing = TRUE))[1:N]

#Calculates relative abundance
GP.genus.prop <- transform_sample_counts(GP.genus, function(x) x / sum(x) )
GP.genus.prop.top <- prune_taxa(top_40, GP.genus.prop)

otu_mat <- as(otu_table(GP.genus.prop.top), "matrix")
if(taxa_are_rows(GP.genus.prop.top)) {
  otu_mat <- t(otu_mat)
}
sample_sums <- rowSums(otu_mat)
print(sample_sums)
min(sample_sums)
max(sample_sums)
median(sample_sums)

GP.genus.prop.top <- microbiome::transform(GP.genus.prop.top, "compositional")
GP.genus.prop.top <- aggregate_taxa(GP.genus.prop.top, level = 'Genus')

#concatinate taxonomy names
tax_df <- as.data.frame(tax_table(GP.genus.prop.top))
tax_df$full_taxonomy <- apply(tax_df[, c("Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) {
  x <- x[!is.na(x)]
  paste(x, collapse = "/")
})
taxa_names(GP.genus.prop.top) <- tax_df$full_taxonomy

#Creates bar plots showing relative abundance per group
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 8,
                                                                       face = "italic", 
                                                                       colour = "Black", 
                                                                       angle = 0)))


ra_water <- plot_composition(GP.genus.prop.top, sample.sort = 'Water_mass', 
                             group_by = 'Water_mass', x.label ="Station") +
  scale_fill_manual("Genus", values = c("indianred2", "darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", 
                                        "darkolivegreen1","darkorange","royalblue2", "darksalmon", 
                                        "palegreen1", "darkslategrey", "pink2", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                                        "firebrick4", "darkolivegreen4","deepskyblue3", "darkgoldenrod3","tomato3", "peachpuff",
                                        "skyblue4","khaki2","orange","royalblue4","goldenrod1", "slateblue", "olivedrab","thistle",
                                        "steelblue1","palevioletred", "darkgreen", "gold", "purple2", "sienna",
                                        "slategray2", "maroon", "mediumaquamarine")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = 10),
        panel.background = element_blank(), 
        legend.position = "bottom") +
  ggtitle("Relative abundance") +
  guide_italics
print(ra_water)
ggsave(file.path(results_dir_ad,"RA_water_all_samples_with_PA.svg"), plot = ra_water, 
       device = 'svg', width = 8000, height = 1750, units = 'px')
ggsave(file.path(results_dir_ad,"RA_water_all_samples_with_PA.pdf"), plot = ra_water, 
       device = 'pdf', width = 8000, height = 1750, units = 'px')
ggsave(file.path(results_dir_ad,"RA_water_all_samples_with_PA.png"), plot = ra_water, 
       device = 'png', width = 8000, height = 1750, units = 'px')

###Generates PCoA plots all with Pseudoalteromonadacea axes 1-2
mass_pallet <- c("darkkhaki","pink2", "slateblue","paleturquoise3","goldenrod1")
ord_bc_0m_12 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_12 <- plot_ordination(water, ord_bc_0m_12, color="Water_mass", shape = "Pseudoaltero_Status",  axes=c(1, 2)) +
  scale_colour_manual(values = mass_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_12 <- PCoA_water_bc_0m_12 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_12)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_Pa_12.pdf"), plot = PCoA_water_bc_0m_12, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots all axes 1-3
ord_bc_0m_13 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_13 <- plot_ordination(water, ord_bc_0m_13, color="Water_mass", shape = "Pseudoaltero_Status", axes=c(1, 3)) +
  scale_colour_manual(values = mass_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_13 <- PCoA_water_bc_0m_13 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_13)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_Pa_13.pdf"), plot = PCoA_water_bc_0m_13, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots all axes 2-3
ord_bc_0m_23 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_23 <- plot_ordination(water, ord_bc_0m_13, color="Water_mass", shape = "Pseudoaltero_Status", axes=c(2, 3)) +
  scale_colour_manual(values = mass_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_23 <- PCoA_water_bc_0m_23 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_23)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_Pa_23.pdf"), plot = PCoA_water_bc_0m_23, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

### Delete samples with Fam Pseudoalteromonadaceae more than 10% (filtration)
samples_to_keep <- pseudo_percent %>%
  filter(Abundance <= 10) %>%
  pull(Sample)
water <- prune_samples(samples_to_keep, water)
water

###Barplot for Lab
#Aggregates data at the Genus level
water_filtr <- subset_samples(water, Sample_type == 'Lab control')
water_filtr
GP.genus <- tax_glom(water_filtr, 'Genus')

#Selects top 40 most abundant genera
N <- 8
top_8 <- names(sort(taxa_sums(GP.genus), decreasing = TRUE))[1:N]

#Calculates relative abundance
GP.genus_top <- prune_taxa(top_8, GP.genus)
GP.genus_top_prop <- microbiome::transform(GP.genus_top, "compositional")
tax_df <- as.data.frame(tax_table(GP.genus_top_prop))
tax_df$full_taxonomy <- apply(tax_df[, c("Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) {
  x <- x[!is.na(x)]
  paste(x, collapse = "_")
})
if(any(duplicated(tax_df$full_taxonomy))) {
  stop("Duplicates")
}
taxa_names(GP.genus_top_prop) <- tax_df$full_taxonomy
GP.genus.prop.top <- GP.genus_top_prop

#Creates bar plots showing relative abundance per group
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 4,
                                                                       face = "italic", 
                                                                       colour = "Black", 
                                                                       angle = 0)))


ra_water <- plot_composition(GP.genus.prop.top, sample.sort = 'Sample_type', 
                             group_by = 'Sample_type', x.label ="Sample_type") +
  scale_fill_manual("Genus", values = c("indianred2", "darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", 
                                        "darkolivegreen1","darkorange","royalblue2", "darksalmon", 
                                        "palegreen1", "darkslategrey", "pink2", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                                        "firebrick4", "darkolivegreen4","deepskyblue3", "darkgoldenrod3","tomato3", "peachpuff",
                                        "skyblue4","khaki2","orange","royalblue4","goldenrod1", "slateblue", "olivedrab","thistle",
                                        "steelblue1","palevioletred", "darkgreen", "gold", "purple2", "sienna",
                                        "slategray2", "maroon", "mediumaquamarine")) +
  theme_bw() + 
  geom_bar(stat="identity", width=0.01) +
  theme(axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = 10),
        panel.background = element_blank(), 
        legend.position = "bottom") +
  ggtitle("Relative abundance") +
  guide_italics
print(ra_water)
ggsave(file.path(results_dir_ad,"RA_LAB.svg"), plot = ra_water, 
       device = 'svg', width = 5000, height = 1000, units = 'px')
ggsave(file.path(results_dir_ad,"RA_LAB.pdf"), plot = ra_water, 
       device = 'pdf', width = 5000, height = 1000, units = 'px')
ggsave(file.path(results_dir_ad,"RA_LAB.png"), plot = ra_water, 
       device = 'png', width = 6000, height = 1750, units = 'px')


###
sample_data(water)$Temperature_C <- as.numeric(as.character(sample_data(water)$Temperature_C))
sample_data(water)$Salinity_ppt <- as.numeric(as.character(sample_data(water)$Salinity_ppt))
sample_data(water)$Depth_m <- as.numeric(as.character(sample_data(water)$Depth_m))

###Generates PCoA plots all axes 1-2 temperature
ord_bc_0m_12 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_12 <- plot_ordination(water, ord_bc_0m_12, color="Temperature_C", axes=c(1, 2)) +
  scale_color_gradient(low="paleturquoise1", high="indianred2") +
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_12 <- PCoA_water_bc_0m_12 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_12)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_vs_temperature_12.pdf"), plot = PCoA_water_bc_0m_12, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots all axes 1-3
ord_bc_0m_13 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_13 <- plot_ordination(water, ord_bc_0m_13, color="Temperature_C", axes=c(1, 3)) +
  scale_color_gradient(low="paleturquoise1", high="indianred2") +
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_13 <- PCoA_water_bc_0m_13 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_13)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_vs_temperature_13.pdf"), plot = PCoA_water_bc_0m_13, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots all axes 2-3
ord_bc_0m_23 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_23 <- plot_ordination(water, ord_bc_0m_13, color="Temperature_C", axes=c(2, 3)) +
  scale_color_gradient(low="paleturquoise1", high="indianred2") + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_23 <- PCoA_water_bc_0m_23 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_23)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_vs_temperature_23.pdf"), plot = PCoA_water_bc_0m_23, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

###Generates PCoA plots all axes 1-2 salinity
ord_bc_0m_12 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_12 <- plot_ordination(water, ord_bc_0m_12, color="Salinity_ppt", axes=c(1, 2)) +
  scale_color_gradient(low="darkblue", high="darkgoldenrod2") +
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_12 <- PCoA_water_bc_0m_12 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_12)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_vs_salinity_12.pdf"), plot = PCoA_water_bc_0m_12, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots all axes 1-3
ord_bc_0m_13 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_13 <- plot_ordination(water, ord_bc_0m_13, color="Salinity_ppt", axes=c(1, 3)) +
  scale_color_gradient(low="darkblue", high="darkgoldenrod2") +
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_13 <- PCoA_water_bc_0m_13 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_13)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_vs_salinity_13.pdf"), plot = PCoA_water_bc_0m_13, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots all axes 2-3
ord_bc_0m_23 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_23 <- plot_ordination(water, ord_bc_0m_13, color="Salinity_ppt", axes=c(2, 3)) +
  scale_color_gradient(low="darkblue", high="darkgoldenrod2") + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_23 <- PCoA_water_bc_0m_23 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_23)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_vs_salinity_23.pdf"), plot = PCoA_water_bc_0m_23, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

###Generates PCoA plots all axes 1-2 depth
ord_bc_0m_12 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_12 <- plot_ordination(water, ord_bc_0m_12, color="Depth_m", axes=c(1, 2)) +
  scale_color_gradient(low="paleturquoise3", high="darkblue") +
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_12 <- PCoA_water_bc_0m_12 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_12)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_vs_depth_12.pdf"), plot = PCoA_water_bc_0m_12, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots all axes 1-3
ord_bc_0m_13 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_13 <- plot_ordination(water, ord_bc_0m_13, color="Depth_m", axes=c(1, 3)) +
  scale_color_gradient(low="paleturquoise3", high="darkblue") +
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_13 <- PCoA_water_bc_0m_13 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_13)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_vs_depth_13.pdf"), plot = PCoA_water_bc_0m_13, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots all axes 2-3
ord_bc_0m_23 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_23 <- plot_ordination(water, ord_bc_0m_13, color="Depth_m", axes=c(2, 3)) +
  scale_color_gradient(low="paleturquoise3", high="darkblue") + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_23 <- PCoA_water_bc_0m_23 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_23)
ggsave(file.path(results_dir_ad,"PCoA_water_bc_vs_depth_23.pdf"), plot = PCoA_water_bc_0m_23, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

###Alpha
mass_pallet <- c("darkkhaki","pink2", "slateblue","paleturquoise3","goldenrod1")
mass_clusters_pallet <- c("#f2bac6", "#f389a1", "#887bd7", "#abd7d7")

#Alpha mass

alpha_water_Chao <- plot_richness(water, 'Water_mass', color = 'Water_mass', 
                                  measures = "Chao1")
alpha_water_Shan <- plot_richness(water, 'Water_mass', color = 'Water_mass', 
                                  measures = "Shannon")


alpha_water_Chao <- alpha_water_Chao +
  theme_minimal() +
  geom_boxplot(aes(fill = Water_mass), color = "black", coef = 0) +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = mass_pallet) + 
  stat_compare_means(comparisons = list(c("BSBW","HAW", "KC", "MW", 
                                          "SAW")), label = "p.signif", 
                     size = 2.9, method = "t.test", exact=FALSE, vjust = 0.1, 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*","NS"))) +
  facet_wrap(~variable, scales = "free_x")
print(alpha_water_Chao)


alpha_water_Shan <- alpha_water_Shan +
  theme_minimal() +
  geom_boxplot(aes(fill = Water_mass), color = "black", coef = 0) +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = mass_pallet) +
  facet_wrap(~variable, scales = "free_x")
print(alpha_water_Shan)

alpha_water_Chao <- alpha_water_Chao + theme(legend.position = "none", aspect.ratio=1)
alpha_water_Shan <- alpha_water_Shan + theme(legend.position = "none", aspect.ratio=1)

ggsave(file.path(results_dir_ad,"Alpha_water_mass_Chao.pdf"), 
       plot = alpha_water_Chao, device = 'pdf',
       width = 1350, height = 1000, units = 'px')
ggsave(file.path(results_dir_ad,"Alpha_water_mass_Shan.pdf"), 
       plot = alpha_water_Shan, device = 'pdf',
       width = 1350, height = 1000, units = 'px')

#Alpha clusters

alpha_water_Chao <- plot_richness(water, 'Microbiome_cluster', color = 'Microbiome_cluster', 
                                  measures = "Chao1")
alpha_water_Shan <- plot_richness(water, 'Microbiome_cluster', color = 'Microbiome_cluster', 
                                  measures = "Shannon")

alpha_water_Chao <- alpha_water_Chao +
  theme_minimal() +
  geom_boxplot(aes(fill = Microbiome_cluster), color = "black", coef = 0) +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = mass_clusters_pallet) + 
  facet_wrap(~variable, scales = "free_x")
print(alpha_water_Chao)


alpha_water_Shan <- alpha_water_Shan +
  theme_minimal() +
  geom_boxplot(aes(fill = Microbiome_cluster), color = "black", coef = 0) +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = mass_clusters_pallet) +
  facet_wrap(~variable, scales = "free_x")
print(alpha_water_Shan)

alpha_water_Chao <- alpha_water_Chao + theme(legend.position = "none", aspect.ratio=1)
alpha_water_Shan <- alpha_water_Shan + theme(legend.position = "none", aspect.ratio=1)

ggsave(file.path(results_dir_ad,"Alpha_water_clusters_Chao.pdf"), 
       plot = alpha_water_Chao, device = 'pdf',
       width = 1350, height = 1000, units = 'px')
ggsave(file.path(results_dir_ad,"Alpha_water_clusters_Shan.pdf"), 
       plot = alpha_water_Shan, device = 'pdf',
       width = 1350, height = 1000, units = 'px')
