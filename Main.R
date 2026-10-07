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
library('microbiomeMarker') #LefSe
library(pheatmap)
library(gridExtra)
library(VennDiagram) # For Venn diagrams
#SpeSpeNet graph
library(pals)
library(tidygraph)
library(ggraph)
library(ggdark)
library(igraph)


#packageVersion("package_name")

results_dir_m <- "C:/Your/path/Main/"
###Creating Phyloseq
seqtab.nochim_OTU <- read.table(sep = "\t", file = "C:/Your/path/all_OTU_frequency.tsv", header = TRUE, check.names=FALSE, row.names=1)
taxa_OTU <- read.table(sep = "\t", file = "C:/Your/path/all_OTUs_phylogeny.tsv", header = TRUE, check.names=FALSE, row.names=1)
taxa_OTU <- as.matrix(taxa_OTU)
sample.data <- read_excel("C:/Your/path/Supplementary_Table_1_St_Anna_Trough_metadata_0m_50m.xlsx")
sample.data <- sample.data %>% tibble::column_to_rownames("Sample_ID")
water <- phyloseq(otu_table(seqtab.nochim_OTU, taxa_are_rows = TRUE), sample_data(sample.data), tax_table(taxa_OTU))
water

###Filtering taxa
water <- subset_taxa(water, Kingdom %in% c("Bacteria", "Archaea"))
water <- subset_taxa(water, !(Family %in% c("Mitochondria")))
water

###Normalization and rare taxa removal
total = median(sample_sums(water))
standf = function(x, t=total) round(t * (x / sum(x)))
water = transform_sample_counts(water, standf)
total_ab_001 = 1e-3
reads_per_OTU_water <- taxa_sums(water)
keepTaxa_001 = ((reads_per_OTU_water / sum(reads_per_OTU_water))*100) > total_ab_001
water = prune_taxa(keepTaxa_001, water)
water

###Find and delete samples with Fam Pseudoalteromonadaceae more than 10% (filtration)
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
#delete
samples_to_keep <- pseudo_percent %>%
  filter(Abundance <= 10) %>%
  pull(Sample)
water <- prune_samples(samples_to_keep, water)
water

###Palettes
big_pallet <- c( "steelblue1","gold","palevioletred", "darkgreen","skyblue4","khaki2","orange","royalblue4",
                    "goldenrod1", "slateblue", "olivedrab","thistle", "purple2", "sienna", "plum1",
                    "slategray2", "maroon", "mediumaquamarine", "darkgoldenrod3","tomato3", "peachpuff",
                    "firebrick4", "deepskyblue3","darkolivegreen4","indianred2", "darkblue", "darkgoldenrod2", "darkorchid",
                    "paleturquoise3", "darkolivegreen1","darkorange","royalblue2", "darksalmon",
                    "palegreen1", "pink2", "darkkhaki", "darkslategrey", "seagreen", "burlywood","paleturquoise1")

mass_0_pallet <- c("goldenrod1", "darkkhaki", "slateblue","paleturquoise3")
mass_50_pallet <- c("darkkhaki","pink2","tomato3", "slateblue")
mass_pallet <- c("goldenrod1", "darkkhaki", "pink2", "tomato3", "slateblue", "paleturquoise3")

###Group samples
water <- subset_samples(water, Station != '6Z Kit')
water
water_0m <- subset_samples(water, Depth_m == '0')
water_0m
water_50m <- subset_samples(water, Depth_m != '0')
water_50m
<<<<<<< HEAD
water_core <- subset_samples(water, Microbiome_core %in% c("MW core", "KC core", 
                                                                   "CDW core", 
                                                                   "CDW+BSBW core", 
                                                           "BSBW core"))
=======
water_core <- subset_samples(water, Microbiome_cluster_core %in% c("MW core", "KC core", 
                                                                   "CDW core", 
                                                                   "CDW+BSBW core"))

water_core_0m <- subset_samples(water_core, Depth_m == '0')

water_core_50m <- subset_samples(water_core, Depth_m != '0')

################################################################################


###3D for 0 m
meta <- as(sample_data(water_0m), "data.frame")
otu <- otu_table(water_0m)
pcoa3d_water <- PCoA3D(METRIC="braycurtis", METADATA=meta, FEATURES=otu, COLOR="Water_mass", AXIS=c(1,2,3), PALETTE=mass_0_pallet)
pcoa3d_water
pcoa3d_water <- PCoA3D(METRIC="braycurtis", METADATA=meta, FEATURES=otu, COLOR="Station_ID", AXIS=c(1,2,3), PALETTE=big_pallet)
pcoa3d_water

###3D for 50 m
meta <- as(sample_data(water_50m), "data.frame")
otu <- otu_table(water_50m)
pcoa3d_water <- PCoA3D(METRIC="braycurtis", METADATA=meta, FEATURES=otu, COLOR="Water_mass", AXIS=c(1,2,3), PALETTE=mass_50_pallet)
pcoa3d_water
pcoa3d_water <- PCoA3D(METRIC="braycurtis", METADATA=meta, FEATURES=otu, COLOR="Station_ID", AXIS=c(1,2,3), PALETTE=big_pallet)
pcoa3d_water

###3D for 0 m and 50 m
meta <- as(sample_data(water), "data.frame")
otu <- otu_table(water)
pcoa3d_water <- PCoA3D(METRIC="braycurtis", METADATA=meta, FEATURES=otu, COLOR="Water_mass", AXIS=c(1,2,3), PALETTE=mass_pallet)
pcoa3d_water
pcoa3d_water <- PCoA3D(METRIC="braycurtis", METADATA=meta, FEATURES=otu, COLOR="Station_ID", AXIS=c(1,2,3), PALETTE=big_pallet)
pcoa3d_water

###Bray tree
dist_bc <- phyloseq::distance(water, method = "bray")
hc <- hclust(dist_bc, method = "ward.D2")
dend <- as.dendrogram(hc)
meta <- data.frame(sample_data(water))
labels(dend) <- meta$Station[match(labels(dend), rownames(meta))]
pdf(file.path(results_dir_m, "Dendrogram_all_stations.pdf"), width = 15, height = 3)
plot(dend)
dev.off()

#Get labels order and create new order for samples
dend_order <- labels(dend)
dend_order_rev <- rev(dend_order)

###Barplot for all samples
#Aggregates data at the Genus level
GP.genus <- tax_glom(water, 'Genus')

#Selects top 40 most abundant genera
N <- 40
top_40 <- names(sort(taxa_sums(GP.genus), decreasing = TRUE))[1:N]

#Calculates relative abundance
GP.genus.prop <- transform_sample_counts(GP.genus, function(x) x / sum(x) )
GP.genus.prop.top <- prune_taxa(top_40, GP.genus.prop)
GP.genus.prop.top <- microbiome::transform(GP.genus.prop.top, "compositional")
GP.genus.prop.top <- aggregate_taxa(GP.genus.prop.top, level = 'Genus')

#concatinate taxonomy names
tax_df <- as.data.frame(tax_table(GP.genus.prop.top))
tax_df$full_taxonomy <- apply(tax_df[, c("Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) {
  x <- x[!is.na(x)]
  paste(x, collapse = "_")
})
taxa_names(GP.genus.prop.top) <- tax_df$full_taxonomy

#New order for samples
meta <- sample_data(GP.genus.prop.top)
meta_df <- as.data.frame(meta)
meta_df$Station <- factor(meta_df$Station, levels = dend_order_rev)
sample_data(GP.genus.prop.top) <- sample_data(meta_df)

#Creates bar plots showing relative abundance per group
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 8,
                                                                       face = "italic", 
                                                                       colour = "Black", 
                                                                       angle = 0)))

#group_by = 'Water_mass'
ra_water <- plot_composition(GP.genus.prop.top, sample.sort = 'Station'
                             , x.label ="Station") +
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
ggsave(file.path(results_dir_m,"RA_water_all_samples.svg"), plot = ra_water, 
       device = 'svg', width = 8000, height = 1750, units = 'px')
ggsave(file.path(results_dir_m,"RA_water_all_samples.pdf"), plot = ra_water, 
       device = 'pdf', width = 9000, height = 1950, units = 'px')
ggsave(file.path(results_dir_m,"RA_water_all_samples.png"), plot = ra_water, 
       device = 'png', width = 8000, height = 1750, units = 'px')

###Temperature and Salinity heatmap
meta_df <- as.data.frame(sample_data(water))
meta_df$Station <- factor(meta_df$Station, levels = dend_order_rev)
meta_df <- meta_df[order(meta_df$Station), ]
meta_df$Temperature_C <- as.numeric(as.character(meta_df$Temperature_C))
meta_df$Salinity_PSU <- as.numeric(as.character(meta_df$Salinity_PSU))
mat_temp <- as.matrix(meta_df[, "Temperature_C", drop=FALSE])
rownames(mat_temp) <- meta_df$Station
mat_sal <- as.matrix(meta_df[, "Salinity_PSU", drop=FALSE])
rownames(mat_sal) <- meta_df$Station

p1 <- pheatmap(mat_temp,
               color = colorRampPalette(c("paleturquoise1","indianred2"))(50),
               cluster_rows = FALSE, cluster_cols = FALSE,
               fontsize_row = 8, main = "Temperature_C")

p2 <- pheatmap(mat_sal,
               color = colorRampPalette(c("darkblue","darkgoldenrod2"))(50),
               cluster_rows = FALSE, cluster_cols = FALSE,
               fontsize_row = 8, main = "Salinity_PSU")

grid.arrange(p1$gtable, p2$gtable, ncol = 2)

pdf(file.path(results_dir_m,"Heatmaps_water_temp_salinity.pdf"), width = 5, height = 8)
grid.arrange(p1$gtable, p2$gtable, ncol = 2)
dev.off()

####Generates PCoA plots 0m axes 1-2
ord_bc_0m_12 <- ordinate(water_0m, "PCoA", "bray")
PCoA_water_bc_0m_12 <- plot_ordination(water_0m, ord_bc_0m_12, color="Water_mass", axes=c(1, 2)) +
  scale_colour_manual(values = mass_0_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_12 <- PCoA_water_bc_0m_12 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_12)
ggsave(file.path(results_dir_m,"PCoA_water_bc_0m_12.pdf"), plot = PCoA_water_bc_0m_12, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#legend
ord_bc_0m_12 <- ordinate(water_0m, "PCoA", "bray")
PCoA_water_bc_0m_12 <- plot_ordination(water_0m, ord_bc_0m_12, color="Water_mass", axes=c(1, 2)) +
  scale_colour_manual(values = mass_0_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_12 <- PCoA_water_bc_0m_12 + theme(legend.position = "top", aspect.ratio=1)
print(PCoA_water_bc_0m_12)
ggsave(file.path(results_dir_m,"PCoA_water_bc_0m_12_legend.pdf"), plot = PCoA_water_bc_0m_12, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots 0m axes 1-3
ord_bc_0m_13 <- ordinate(water_0m, "PCoA", "bray")
PCoA_water_bc_0m_13 <- plot_ordination(water_0m, ord_bc_0m_13, color="Water_mass", axes=c(1, 3)) +
  scale_colour_manual(values = mass_0_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_13 <- PCoA_water_bc_0m_13 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_13)
ggsave(file.path(results_dir_m,"PCoA_water_bc_0m_13.pdf"), plot = PCoA_water_bc_0m_13, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots 0m axes 2-3
ord_bc_0m_23 <- ordinate(water_0m, "PCoA", "bray")
PCoA_water_bc_0m_23 <- plot_ordination(water_0m, ord_bc_0m_13, color="Water_mass", axes=c(2, 3)) +
  scale_colour_manual(values = mass_0_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_23 <- PCoA_water_bc_0m_23 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_23)
ggsave(file.path(results_dir_m,"PCoA_water_bc_0m_23.pdf"), plot = PCoA_water_bc_0m_23, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

###Generates PCoA plots 50m axes 1-2
ord_bc_0m_12 <- ordinate(water_50m, "PCoA", "bray")
PCoA_water_bc_0m_12 <- plot_ordination(water_50m, ord_bc_0m_12, color="Water_mass", axes=c(1, 2)) +
  scale_colour_manual(values = mass_50_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_12 <- PCoA_water_bc_0m_12 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_12)
ggsave(file.path(results_dir_m,"PCoA_water_bc_50m_12.pdf"), plot = PCoA_water_bc_0m_12, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#legend
ord_bc_0m_12 <- ordinate(water_50m, "PCoA", "bray")
PCoA_water_bc_0m_12 <- plot_ordination(water_50m, ord_bc_0m_12, color="Water_mass", axes=c(1, 2)) +
  scale_colour_manual(values = mass_50_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_12 <- PCoA_water_bc_0m_12 + theme(legend.position = "top", aspect.ratio=1)
print(PCoA_water_bc_0m_12)
ggsave(file.path(results_dir_m,"PCoA_water_bc_50m_12_legend.pdf"), plot = PCoA_water_bc_0m_12, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots 50m axes 1-3
ord_bc_0m_13 <- ordinate(water_50m, "PCoA", "bray")
PCoA_water_bc_0m_13 <- plot_ordination(water_50m, ord_bc_0m_13, color="Water_mass", axes=c(1, 3)) +
  scale_colour_manual(values = mass_50_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
#PCoA_water_bc_0m_13 <- PCoA_water_bc_0m_13 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_13)
ggsave(file.path(results_dir_m,"PCoA_water_bc_50m_13.pdf"), plot = PCoA_water_bc_0m_13, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots 50m axes 2-3
ord_bc_0m_23 <- ordinate(water_50m, "PCoA", "bray")
PCoA_water_bc_0m_23 <- plot_ordination(water_50m, ord_bc_0m_13, color="Water_mass", axes=c(2, 3)) +
  scale_colour_manual(values = mass_50_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_23 <- PCoA_water_bc_0m_23 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_23)
ggsave(file.path(results_dir_m,"PCoA_water_bc_50m_23.pdf"), plot = PCoA_water_bc_0m_23, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

###Generates PCoA plots all axes 1-2
ord_bc_0m_12 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_12 <- plot_ordination(water, ord_bc_0m_12, color="Water_mass", axes=c(1, 2)) +
  scale_colour_manual(values = mass_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_12 <- PCoA_water_bc_0m_12 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_12)
ggsave(file.path(results_dir_m,"PCoA_water_bc_12.pdf"), plot = PCoA_water_bc_0m_12, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#legend
ord_bc_0m_12 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_12 <- plot_ordination(water, ord_bc_0m_12, color="Water_mass", axes=c(1, 2)) +
  scale_colour_manual(values = mass_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_12 <- PCoA_water_bc_0m_12 + theme(legend.position = "top", aspect.ratio=1)
print(PCoA_water_bc_0m_12)
ggsave(file.path(results_dir_m,"PCoA_water_bc_12_legend.pdf"), plot = PCoA_water_bc_0m_12, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots all axes 1-3
ord_bc_0m_13 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_13 <- plot_ordination(water, ord_bc_0m_13, color="Water_mass", axes=c(1, 3)) +
  scale_colour_manual(values = mass_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_13 <- PCoA_water_bc_0m_13 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_13)
ggsave(file.path(results_dir_m,"PCoA_water_bc_13.pdf"), plot = PCoA_water_bc_0m_13, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

#Generates PCoA plots all axes 2-3
ord_bc_0m_23 <- ordinate(water, "PCoA", "bray")
PCoA_water_bc_0m_23 <- plot_ordination(water, ord_bc_0m_13, color="Water_mass", axes=c(2, 3)) +
  scale_colour_manual(values = mass_pallet) + 
  geom_point(size=2.5, alpha=1)+ theme_bw() + 
  theme(text = element_text(size = 15))
PCoA_water_bc_0m_23 <- PCoA_water_bc_0m_23 + theme(legend.position = "none", aspect.ratio=1)
print(PCoA_water_bc_0m_23)
ggsave(file.path(results_dir_m,"PCoA_water_bc_23.pdf"), plot = PCoA_water_bc_0m_23, 
       device = 'pdf', width = 1500, height = 1000, units = 'px')

################################################################################

# Select core samples.

water_core <- prune_taxa(taxa_sums(water_core) > 0, water_core)

# Extract the OTU abundance matrix
water_core_otu_mat <- as(otu_table(water_core), "matrix")

# Convert to data frame
water_core_otu_df <- as.data.frame(water_core_otu_mat)

# Save to file (tab-separated, use CSV for comma-separated)
write.table(water_core_otu_df, file = file.path(results_dir_m, "Core_groups_otu_table.tsv"), sep = "\t", quote = FALSE, col.names = NA)


###LefSe all vs all core groups.
lef_out <- run_lefse(water_core, group = "Microbiome_core", norm = "CPM", taxa_rank = 'none',
                     kw_cutoff = 0.01, lda_cutoff = 3.5)
print(plot_ef_bar(lef_out))
lef_out
lef_out_OTUs_list <- lef_out@marker_table$feature
lef_out_taxonomy_tt <- water_core@tax_table[lef_out_OTUs_list, c('Order', 'Family', 'Genus')]
lef_taxonomy_df <- data.frame(as(tax_table(lef_out_taxonomy_tt), "matrix"))
lef_taxonomy_df$Full_taxonomy <- paste(lef_taxonomy_df$Order, lef_taxonomy_df$Family, lef_taxonomy_df$Genus, rownames(lef_taxonomy_df), sep='/')
lef_taxonomy_df$Full_taxonomy
lef_out@marker_table$feature <- lef_taxonomy_df$Full_taxonomy
lefse_table <- lef_out@marker_table
write.table(lefse_table, sep="\t", file.path(results_dir_m, "LefSe_taxa_core_vs_core_lda_3_5.tsv"))
plot_ef_bar(lef_out)
pdf(file.path(results_dir_m, "LefSe_taxa_core_vs_core_lda_3_5.pdf"), width = 10, height = 8)
plot(plot_ef_bar(lef_out))
dev.off()

# Split phylosec by core groups. Retain all non-zero OTUs.
MW_water_core <- subset_samples(water_core, Microbiome_core %in% c("MW core"))
MW_water_core <- prune_taxa(taxa_sums(MW_water_core) > 0, MW_water_core)
MW_water_core_otu_list <- taxa_names(MW_water_core)
MW_water_core_otu_list

KC_water_core <- subset_samples(water_core, Microbiome_core %in% c("KC core"))
KC_water_core <- prune_taxa(taxa_sums(KC_water_core) > 0, KC_water_core)
KC_water_core_otu_list <- taxa_names(KC_water_core)
KC_water_core_otu_list

<<<<<<< HEAD
CDW_water_core <- subset_samples(water_core, Microbiome_core %in% c("CDW core"))
CDW_water_core <- prune_taxa(taxa_sums(CDW_water_core) > 0, CDW_water_core)
CDW_water_core_otu_list <- taxa_names(CDW_water_core)
CDW_water_core_otu_list

CDW_BSBW_water_core <- subset_samples(water_core, Microbiome_core %in% c("CDW+BSBW core"))
CDW_BSBW_water_core <- prune_taxa(taxa_sums(CDW_BSBW_water_core) > 0, CDW_BSBW_water_core)
CDW_BSBW_water_core_otu_list <- taxa_names(CDW_BSBW_water_core)
CDW_BSBW_water_core_otu_list

BSBW_water_core <- subset_samples(water_core, Microbiome_core %in% c("BSBW core"))
BSBW_water_core <- prune_taxa(taxa_sums(BSBW_water_core) > 0, BSBW_water_core)
BSBW_water_core_otu_list <- taxa_names(BSBW_water_core)
BSBW_water_core_otu_list
=======
HAW_1_water_core <- subset_samples(water_core, Microbiome_cluster_core %in% c("CDW core"))
HAW_1_water_core <- prune_taxa(taxa_sums(HAW_1_water_core) > 0, HAW_1_water_core)
HAW_1_water_core_otu_list <- taxa_names(HAW_1_water_core)
HAW_1_water_core_otu_list

HAW_2_water_core <- subset_samples(water_core, Microbiome_cluster_core %in% c("CDW+BSBW core"))
HAW_2_water_core <- prune_taxa(taxa_sums(HAW_2_water_core) > 0, HAW_2_water_core)
HAW_2_water_core_otu_list <- taxa_names(HAW_2_water_core)
HAW_2_water_core_otu_list


# Plot a Venn diagram showing OTU sets intersections between core groups of water samples.
venn.plot <- venn.diagram(
  x = list(MW_core = MW_water_core_otu_list, KC_core = KC_water_core_otu_list, 
           CDW_core = CDW_water_core_otu_list, CDW_BSBW_core = CDW_BSBW_water_core_otu_list, 
           BSBW_core = BSBW_water_core_otu_list),
  filename = NULL,
  fill = c("pink2",  "slateblue", "darkkhaki","tomato3", "paleturquoise3"), # Your preferred colors
  fontfamily = "sans",
  cat.fontfamily = "sans",
  alpha = 0.5,
  cex = 1.5,
  category.cex = 1.5,
  resolution=300
)

pdf(file.path(results_dir_m, "Venn_taxa_core_vs_core.pdf"), width = 4.5, height = 4.5)
grid.draw(venn.plot)
dev.off()

# Retain core OTUs present in all samples of a set.
MW_water_core_otu_core <- filter_taxa(MW_water_core, function(x) all(x > 0), prune = TRUE)
MW_water_core_otu_core_list <- taxa_names(MW_water_core_otu_core)
MW_water_core_otu_core_list

KC_water_core_otu_core <- filter_taxa(KC_water_core, function(x) all(x > 0), prune = TRUE)
KC_water_core_otu_core_list <- taxa_names(KC_water_core_otu_core)
KC_water_core_otu_core_list

CDW_water_core_otu_core <- filter_taxa(CDW_water_core, function(x) all(x > 0), prune = TRUE)
CDW_water_core_otu_core_list <- taxa_names(CDW_water_core_otu_core)
CDW_water_core_otu_core_list

CDW_BSBW_water_core_otu_core <- filter_taxa(CDW_BSBW_water_core, function(x) all(x > 0), prune = TRUE)
CDW_BSBW_water_core_otu_core_list <- taxa_names(CDW_BSBW_water_core_otu_core)
CDW_BSBW_water_core_otu_core_list

BSBW_water_core_otu_core <- filter_taxa(BSBW_water_core, function(x) all(x > 0), prune = TRUE)
BSBW_water_core_otu_core_list <- taxa_names(BSBW_water_core_otu_core)
BSBW_water_core_otu_core_list

# Plot a Venn diagram showing core OTU sets intersections between core groups of water samples.
venn.plot <- venn.diagram(
  x = list(MW_core = MW_water_core_otu_core_list, KC_core = KC_water_core_otu_core_list, 
           CDW_core = CDW_water_core_otu_core_list, 
           CDW_BSBW_core = CDW_BSBW_water_core_otu_core_list, BSBW_core = BSBW_water_core_otu_core_list),
  filename = NULL,
  fill = c("pink2",  "slateblue", "darkkhaki","tomato3", "paleturquoise3"), # Your preferred colors
  fontfamily = "sans",
  cat.fontfamily = "sans",
  alpha = 0.5,
  cex = 1.5,
  category.cex = 1.5,
  resolution=300
)

pdf(file.path(results_dir_m, "Venn_taxa_core_vs_core_core_otus.pdf"), width = 4.5, height = 4.5)
grid.draw(venn.plot)
dev.off()


###SpeSpeNet (If the network is overlayed with kmeans cluster and positive correlations are used)
#https://utrecht-university.shinyapps.io/SpeSpeNet_v1/

#Get initial files for site

#0m core
otu_mat <- as(otu_table(water_core_0m), "matrix")
if(taxa_are_rows(water_core_0m)) {
  otu_mat <- t(otu_mat)
}
otu_mat <- t(otu_mat)
write.table(otu_mat, file = file.path(results_dir_m, "SpeSpeNet_otu_table_0m.csv"), sep = "\t", quote = FALSE, row.names = TRUE)

tax_mat <- as(tax_table(water_core_0m), "matrix")
write.table(tax_mat, file = file.path(results_dir_m, "SpeSpeNet_tax_table_0m.csv"), sep = "\t", quote = FALSE, row.names = TRUE)

metadata <- as(sample_data(water_core_0m), "data.frame")
write.table(metadata, file = file.path(results_dir_m, "SpeSpeNet_sample_metadata_0m.csv"), sep = "\t", quote = FALSE, row.names = TRUE)

#50m core
otu_mat <- as(otu_table(water_core_50m), "matrix")
if(taxa_are_rows(water_core_50m)) {
  otu_mat <- t(otu_mat)
}
otu_mat <- t(otu_mat)
write.table(otu_mat, file = file.path(results_dir_m, "SpeSpeNet_otu_table_50m.csv"), sep = "\t", quote = FALSE, row.names = TRUE)

tax_mat <- as(tax_table(water_core_50m), "matrix")
write.table(tax_mat, file = file.path(results_dir_m, "SpeSpeNet_tax_table_50m.csv"), sep = "\t", quote = FALSE, row.names = TRUE)

metadata <- as(sample_data(water_core_50m), "data.frame")
write.table(metadata, file = file.path(results_dir_m, "SpeSpeNet_sample_metadata_50m.csv"), sep = "\t", quote = FALSE, row.names = TRUE)

#all core
otu_mat <- as(otu_table(water_core), "matrix")
if(taxa_are_rows(water_core)) {
  otu_mat <- t(otu_mat)
}
otu_mat <- t(otu_mat)
write.table(otu_mat, file = file.path(results_dir_m, "SpeSpeNet_otu_table_all.csv"), sep = "\t", quote = FALSE, row.names = TRUE)

tax_mat <- as(tax_table(water_core), "matrix")
write.table(tax_mat, file = file.path(results_dir_m, "SpeSpeNet_tax_table_all.csv"), sep = "\t", quote = FALSE, row.names = TRUE)

metadata <- as(sample_data(water_core), "data.frame")
write.table(metadata, file = file.path(results_dir_m, "SpeSpeNet_sample_metadata_all.csv"), sep = "\t", quote = FALSE, row.names = TRUE)

# Create picture with SpeSpeNet results, define plotting themes:

themes <- list("Dark" = dark_theme_grey(),
               "Classic" = theme_bw(),
               "Minimal" = theme_minimal())

lineCol <- list("Dark" = "white",
                "Classic" = "black",
                "Minimal" = "black")

use_theme <- theme(axis.ticks = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   legend.title = element_text(size = 28),
                   legend.text = element_text(size = 30),
                   legend.position = c(0, 0),
                   legend.justification = c("left", "bottom"),
                   legend.background = element_rect(fill='transparent'),
                   legend.margin = margin(6, 6, 6, 6))

#0m core
tidy.net <- readRDS(file = "C:/Your/path/SpeSpeNet/0m_results/Tidygraph_object_2025-11-05_0m.rds")

node.data <- tidy.net%>%activate(nodes)%>%data.frame()
clus.data <- node.data$cluster

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)

set.seed(1)

plot.clust <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = as.factor(cluster)),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  guides(fill = guide_legend("Cluster", override.aes = list(size = 15))) +
  scale_fill_manual(values = as.vector(polychrome(length(unique(clus.data))+4)[-c(1,2,4,5)])) +
  themes[["Classic"]] +
  use_theme

pdf("C:/Your/path/SpeSpeNet/0m_results/SpeSpeNet_0m.pdf", 
    width = 10, height = 4)
plot(plot.clust)
dev.off()


#50m core
tidy.net <- readRDS(file = "C:/Your/path/SpeSpeNet/50m_results/Tidygraph_object_2025-11-05_50m.rds")

node.data <- tidy.net%>%activate(nodes)%>%data.frame()
clus.data <- node.data$cluster

rem_nodes <- which(degree(tidy.net)==0)
tidy.net <- delete_vertices(tidy.net,rem_nodes)

set.seed(1)

plot.clust <- ggraph(tidy.net, layout = 'fr') +
  geom_edge_arc0(
    width = 0.1,
    alpha = 1,
    strength = 0,
    color = lineCol[["Classic"]]) +
  geom_point(aes(x,y,size = nodeSize,fill = as.factor(cluster)),pch=21,color = "black") +
  scale_size_continuous(range = c(1, 12), guide = 'none') +
  guides(fill = guide_legend("Cluster", override.aes = list(size = 15))) +
  scale_fill_manual(values = as.vector(polychrome(length(unique(clus.data))+4)[-c(1,2,4,5)])) +
  themes[["Classic"]] +
  use_theme

pdf("C:/Your/path/SpeSpeNet/50m_results/SpeSpeNet_50m.pdf", 
    width = 10, height = 4)
plot(plot.clust)
dev.off()


###dbRDA
# Get matrix
comm_matrix <- as(otu_table(water), "matrix")
if (taxa_are_rows(otu_table(water))) { comm_matrix <- t(comm_matrix) }
env_data <- as(sample_data(water), "data.frame")
env_data$Depth_m <- as.numeric(as.character(env_data$Depth_m))
env_data$Temperature_C <- as.numeric(as.character(env_data$Temperature_C))
env_data$Salinity_PSU <- as.numeric(as.character(env_data$Salinity_PSU))
env_data$Water_mass <- as.factor(env_data$Water_mass)

# Standardization
comm_transf <- decostand(comm_matrix, method = "hellinger")
# Distance
diss <- vegdist(comm_transf, method = "bray")

# dbRDA Temperature
mod_temp <- dbrda(diss ~ Temperature_C, data = env_data)
anova(mod_temp, permutations = 999)
RsquareAdj(mod_temp)

# dbRDA Salinity
mod_sal <- dbrda(diss ~ Salinity_PSU, data = env_data)
anova(mod_sal, permutations = 999)
RsquareAdj(mod_sal)

# dbRDA Depth
mod_depth <- dbrda(diss ~ Depth_m, data = env_data)
anova(mod_depth, permutations = 999)
RsquareAdj(mod_depth)

# dbRDA physical parameters together
mod_phys <- dbrda(diss ~ Temperature_C + Salinity_PSU + Depth_m, data = env_data)
anova(mod_phys, permutations = 999)
RsquareAdj(mod_phys)

# dbRDA Water mass
mod_water <- dbrda(diss ~ Water_mass, data = env_data)
anova(mod_water, permutations = 999)
RsquareAdj(mod_water)

# Partial dbRDA — contribution of Water_mass beyond physical parameters
mod_partial_water <- dbrda(diss ~ Water_mass + Condition(Temperature_C + Salinity_PSU + Depth_m), data = env_data)
anova(mod_partial_water, permutations = 999)
RsquareAdj(mod_partial_water)

# Partial dbRDA — contribution of physical parameters beyond Water_mass
mod_partial_phys <- dbrda(diss ~ Temperature_C + Salinity_PSU + Depth_m + Condition(Water_mass), data = env_data)
anova(mod_partial_phys, permutations = 999)
RsquareAdj(mod_partial_phys)

###PERMANOVA
env_data$Microbiome_cluster <- as.factor(env_data$Microbiome_cluster)
dist_matrix <- vegdist(comm_matrix, method = "bray")
adonis2(dist_matrix ~ Microbiome_cluster, data = env_data, permutations = 999)