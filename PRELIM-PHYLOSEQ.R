# Phyloseq analysis for CERF 2025

# Setup----
# load libraries
library(phyloseq); library(ggplot2); library(ggplot2); library(data.table)
library(dplyr); 

# Graph formatting
pretty.theme <- function(){
  theme_bw() +
    theme(axis.text.x=element_text(size = 18, color="black"),
          axis.text.y=element_text(size = 18, color="black"),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          text=element_text(size=18),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA),
          axis.line = element_line(colour = "black"))
}

# Import----
# changed "-" to "_"
# change file names; need to reorder so metadata match otu table
mat = read.table("PHYLOSEQ-INPUT/COUNT-MATRIX-6.txt", header = TRUE, sep = "\t", row.names = 1)
tax = read.table("PHYLOSEQ-INPUT/TAXA-2.txt", header = TRUE, sep = "\t", row.names = 1)
meta = read.table("PHYLOSEQ-INPUT/METADATA-3.txt", header = TRUE, sep = "\t", row.names = 1)

dim(tax)
dim(mat)
dim(meta)

mat = as.matrix(mat)
tax = as.matrix(tax)

OTU <- otu_table(as.matrix(mat), taxa_are_rows = TRUE)
TAX = tax_table(tax)
META = sample_data(meta)

phy = phyloseq(OTU,TAX,META)
phy

# Filter----
phy %>%
  subset_taxa(Family != "Mitochondria",
              Order != "Chloroplast") -> phy.f
phy
phy.f 

# Filter out samples with super low reads
# Filter out samples with low sequence count
phy.f.final <- subset_samples(phy.f, sample_sums(phy.f) >=1000)
phy.f # 70 samples
phy.f.final # 70 samples

# relative abundance
per.f.final = transform_sample_counts(phy.f.final, function (x) x/sum(x)*100)

# ordination----
BC_distance <- phyloseq::distance(per.f.final, "bray") 
bcOrd <- ordinate(per.f.final, "NMDS", BC_distance)
plot_scree(bcOrd)

# Are there other environmental variables you want to look for?
p1 <- plot_ordination(per.f.final, bcOrd, type = "samples") +
  geom_point(aes(fill = Treatment), shape = 21, color = "black", stroke = 0.5, size = 5, alpha = 0.9) +
  scale_fill_viridis_d(option = "C") +
  pretty.theme() +
  labs(fill = "Treatment")
p1

p2 <- plot_ordination(per.f.final, bcOrd, type = "samples") +
  geom_point(aes(color = Treatment, shape = Location), size = 5, alpha = 0.9) +
  scale_color_viridis_d(option = "D") +
  pretty.theme() +
  labs(color = "Treatment", shape = "Location")
p2

\# Alpha diversity----
# But need to rarefy
alpha_plot <- plot_richness(per.f.final, x = "Treatment", measures = "Shannon") +
  geom_boxplot(aes(fill = Treatment)) +
  pretty.theme() +
  scale_fill_viridis_d(option = "C") +
  labs(y="Shannon Diversity", x = "Treatment") +
  theme(strip.background = element_blank(),
        strip.text = element_blank())
alpha_plot

# Barplot----
# Top 20 most abundant by genus
# Aggregate at the Genus level
phy_family <- tax_glom(per.f.final, taxrank = "Family")
taxa_sums <- taxa_sums(phy_family)
top20_taxa <- names(sort(taxa_sums, decreasing = TRUE))[1:20]
phy_top20 <- prune_taxa(top20_taxa, phy_family)

df_top20 <- psmelt(phy_top20)

pal_20 <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#FB9A99", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02",
  "#7570B3", "#c22131", "#66A61E", "#E6AB02", "#A6761D",
  "#d98302", "#A9A9F5", "#15425e", "#B2DF8A", "#1F78B4"
)

# Check the colors
bar_plot <- ggplot(df_top20, aes(x = Treatment, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal_20) +
  labs(x = "Depth (cm)", y = "Relative Abundance", fill = "Family") +
  pretty.theme()
# theme(legend.position = "none")
bar_plot
