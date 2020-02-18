### de novo screening of 9-mer ORFs to identify best predicted immunogenic peptides 
### Author: Chloe H. Lee

library(viridis)
library(ggridges)
library( RColorBrewer)

# import 9-mer peptides
peptides.dframe <- read_csv("data/nCoV_9mer_peptides.csv")

# predict immunogenicity by iPred
peptides_ipred <- peptides.dframe %>% mutate_if(is.factor, as.character) %>% 
  select(id, antigen.epitope = peptide) %>%
  predict_imm() 

peptides_ipred %>%
  ggplot(aes(x = imm.prob, fill = id)) +
  geom_density(alpha = 0.5) +
  scale_fill_brewer("source protein", palette = "Set1") +
  scale_x_continuous("P(Immunogenic)") +
  scale_y_continuous("Density") +
  theme_bw() +
  theme(aspect = 1, 
        legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# the number of immunogenic peptides by ORFs 
cols <- c(brewer.pal(6, "Dark2"), brewer.pal(6,"Set1"))
ggplot(peptides_ipred, aes(x = `imm.prob`, y = `id`, fill = id)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, gradient_lwd = 1.) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0)) + scale_fill_manual(values = cols) +
  scale_x_continuous("P(Immunogenic)") +
  scale_y_discrete("Density") +
  theme_bw() + theme(legend.position = "none")
  # arrange with corresponding open reading frames 


# import NetMHCpan 4.0 results on different HLA alleles 
library(readxl)
NetMHCpan_out_nCol_all_A0101 <- read_excel("data/NetMHCpan_out_nCol_all_A0101.xlsx")
NetMHCpan_out_nCol_all_A0201 <- read_excel("data/NetMHCpan_out_nCol_all_A0201.xlsx")
NetMHCpan_out_nCol_all_B0702 <- read_excel("data/NetMHCpan_out_nCol_all_B0702.xlsx")
NetMHCpan_out_nCol_all_B4001 <- read_excel("data/NetMHCpan_out_nCol_all_B4001.xlsx")
NetMHCpan_out_nCol_all_C0702 <- read_excel("data/NetMHCpan_out_nCol_all_C0702.xlsx")
 
netMHC_df <- select(NetMHCpan_out_nCol_all_A0101, Peptide, ID, A0101.Rank = Rank, A0101.NB = NB, HLA) %>% 
  inner_join(select(NetMHCpan_out_nCol_all_A0201,A0201.Rank = Rank,Peptide, A0201.NB = NB, HLA), by = "Peptide") %>%
  inner_join(select(NetMHCpan_out_nCol_all_B0702,B0702.Rank = Rank,Peptide, B0702.NB = NB, HLA), by = "Peptide") %>%
  inner_join(select(NetMHCpan_out_nCol_all_B4001,B4001.Rank = Rank, Peptide, B4001.NB = NB, HLA), by = "Peptide") %>%
  inner_join(select(NetMHCpan_out_nCol_all_C0702,C0702.Rank = Rank, Peptide, C0702.NB = NB, HLA), by = "Peptide") %>%
  mutate(bind = A0101.NB + A0201.NB + B0702.NB + B4001.NB + C0702.NB) %>%
  arrange(desc(bind))


# merge iPred with NetMHCpan score
peptides_ipred_netMHC <- peptides_ipred %>% 
  inner_join(select(netMHC_df, -ID), by = c("antigen.epitope" = "Peptide")) %>%
  arrange(desc(bind))

peptides_ipred_netMHC %>% filter(imm.prob >= 0.5) %>%
  group_by(bind) %>%
  summarize(n_epitopes = n())

plot(peptides_ipred_netMHC$bind, peptides_ipred_netMHC$imm.prob)
table(peptides_ipred_netMHC$bind)

library(ggpubr)
ggscatter(peptides_ipred_netMHC, x = "bind", y = "imm.prob", 
          size = 1, color = "white",add.params = list(color = "red"),
          xlab = "# HLA types predicted to bind", ylab = "iPred immunogenicity score") + 
  geom_jitter(width = 0.2, colour = "darkred")









