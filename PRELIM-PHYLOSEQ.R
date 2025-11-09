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

p1 <- plot_ordination(per.f.final, bcOrd, type = "samples") +
  geom_point(aes(fill = Treatment), shape = 21, color = "black", stroke = 0.5, size = 5, alpha = 0.9) +
  scale_fill_viridis_d(option = "C") +
  pretty.theme() +
  labs(fill = "Treatment")
p1

p2 <- plot_ordination(per.f.final, bcOrd, type = "samples") +
  geom_point(aes(color = Treatment, shape = Timepoint), size = 5, alpha = 0.9) +
  scale_color_viridis_d(option = "C") +
  pretty.theme() +
  labs(fill = "Treatment")
p2
