args     = commandArgs(trailingOnly=TRUE)

infile  <- args[1]
outdir  <- args[2]

if (is.na(infile)) { # local testing
  infile <- "core_gene_alignment.fasta"
  outdir       <- "."
}

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)
library(viridisLite)

# Install biostrings if needed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("Biostrings")
}

library(Biostrings)

text_base_size   <- 7
ggplot_text_size <- text_base_size / ggplot2::.pt # Trick to make geom_text etc. same size and theme
theme_set(theme_cowplot(font_size = text_base_size, rel_small = 1, rel_tiny = 1, rel_large = 1))


msa <- readDNAStringSet(infile)
n_sequences <- length(msa)

x   <- consensusMatrix(msa)

colmax <- apply(X = x, MARGIN = 2, FUN = max)
i      <- x["-",] > 0          # columns with indels
r1     <- n_sequences - colmax # Number of sequences with minor alleles
r1[i]  <- -x[16,i]             # indels positions are set to -number of seqs with indels
rm(msa,colmax,i)

pd2 <- tibble(alleles=r1) %>%
  mutate(position=row_number()) %>% 
  {.}

rm(r1,x)

# 
# Now we need to calculate the density of SNPs surrounding the current snp.
#

neighbors_snps <- 10
neighbors_indels <- 100

pd <- pd2 %>%
  filter(alleles>0) %>% # Remove indels (-1)  and fixed positions (0)
  mutate(lead_10 = lead(position, n=neighbors_snps),
         lag_10  = lag(position,  n=neighbors_snps),
         snp_density = lead_10-lag_10) %>%
  mutate(snp_density = ifelse(is.na(snp_density), 2*(lead_10-position), snp_density)) %>%
  mutate(snp_density = ifelse(is.na(snp_density), 2*(position-lag_10), snp_density)) %>%
  select(-lead_10, -lag_10) %>%
  mutate(noise=runif(n=n(),min = -1, max=1)) %>%
  mutate(y=noise/snp_density) %>% 
  {.}

plot1 <- ggplot(pd, aes(x=position, y=y, color=alleles)) +
  geom_point(alpha=1, size=0.2, show.legend = TRUE) +
  ggtitle(paste(nrow(pd), "SNP positions in alignment of", n_sequences,"sequences.")) +
  xlab("Alignment position") + 
  scale_color_viridis_c("Minor count", direction = -1, option = "cividis") + 
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  NULL

pd <- pd2 %>%
  filter(alleles<0) %>% # Only look at indels
  mutate(alleles = abs(alleles)) %>%
  mutate(lead_10 = lead(position, n=neighbors_indels),
         lag_10  = lag(position, n=neighbors_indels),
         snp_density = lead_10-lag_10) %>%
  mutate(snp_density = ifelse(is.na(snp_density), 2*(lead_10-position), snp_density)) %>%
  mutate(snp_density = ifelse(is.na(snp_density), 2*(position-lag_10), snp_density)) %>%
  select(-lead_10, -lag_10) %>%
  mutate(noise=runif(n=n(),min = -1, max=1)) %>%
  mutate(y=noise/snp_density) %>% 
  {.}


plot2 <- ggplot(pd, aes(x=position, y=y, color=alleles)) +
  geom_point(alpha=1, size=0.2, show.legend = TRUE) +
  ggtitle(paste(nrow(pd), "INDEL positions in alignment of", n_sequences,"sequences.")) +
  xlab("Alignment position") + 
  scale_color_viridis_c("INDEL count", direction = -1, option = "cividis") + 
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  NULL

final <- plot_grid(plot1,plot2, ncol = 1)

ggsave(filename = paste(outdir, "/figure.alignment.alleles.png", sep=""), plot = final, width = 120, height=80, units = "mm", dpi = 300)
