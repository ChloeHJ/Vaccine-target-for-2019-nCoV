### Sequence similarity of 2019-nCoV peptides with characterized immunogenic peptides 
### Author: Chloe H. Lee

# load library
library(seqinr)
library(plyr)
library(dnar)
library(tidyverse)

dir <- getwd()

# read immunogenic peptide list from IEDB
ImmuData_iedb <- read.csv("IEDB_epitope_Tcellpositive_human.csv") %>% filter(Object.Type == "Linear peptide") 

# read nCoV sequence fasta file 
setwd("data/") 
paths <- dir(pattern = "\\.fasta$")
names(paths) <- basename(paths)
all_seq <- lapply(paths, function(i) read.fasta(i, as.string = TRUE, seqtype = "AA")); str_length(all_seq)
orf.list <- list(
  envelop_protein = all_seq$nCov_envelop_protein.fasta$BBW89519.1[1],
  membrane_glycoprotein = all_seq$nCov_membrane_glycoprotein.fasta$BBW89520.1[1],
  nucleocapsid_phosphoprotein = all_seq$nCov_nucleocapsid_phosphoprotein.fasta$BBW89524.1[1],
  ORF10_protein = all_seq$nCov_ORF10_protein.fasta$BBW89525.1[1],
  orf1ab_polyprotein = all_seq$nCov_orf1ab_polyprotein.fasta$BBW89516.1[1],
  ORF3a_protein = all_seq$nCov_ORF3a_protein.fasta$BBW89518.1[1],
  ORF6_protein = all_seq$nCov_ORF6_protein.fasta$BBW89521.1[1],
  ORF7a_protein = all_seq$nCov_ORF7a_protein.fasta$BBW89522.1[1],
  ORF8_protein = all_seq$nCov_ORF8_protein.fasta$BBW89523.1[1],
  surface_glycoprotein = all_seq$nCov_surface_glycoprotein.fasta$BBW89517.1[1]
)
orf.df <- as.data.frame(t(as.data.frame(orf.list))); colnames(orf.df) <- "Peptide"
setwd(dir)


####### identify exact matches with immunogenic IEDB peptides ##########
all_iedb_grep.list <- list()
for(row in 1:nrow(ImmuData_iedb)){
  iedb <- as.character(ImmuData_iedb$Description)[row]
  all_iedb_grep.list[[row]] <- all_seq %>% str_extract_all(iedb); 
  names(all_iedb_grep.list[[row]]) <- gsub(".fasta", "", names(all_seq))
}
names(all_iedb_grep.list) <- ImmuData_iedb$Description

library(reshape2)
all_iedb_grep.df <- melt(all_iedb_grep.list); colnames(all_iedb_grep.df) <- c("iedb_peptides", "nCol_ORFs", "nCol_pattern")
all_iedb_grep.df <- all_iedb_grep.df %>% 
  inner_join(ImmuData_iedb, by = c("iedb_peptides" = "Description"))  %>%  
  filter(nchar(nCol_pattern) >5)

########## local alignment between nCoV ORF and IEDB peptides #################
library(Biostrings)
iedb_scores <- matrix(NA, ncol = length(orf.list), nrow = nrow(ImmuData_iedb), 
                 dimnames = list(as.character(ImmuData_iedb$Description), as.character(names(orf.list))))
iedb_pattern <- matrix(NA, ncol = length(orf.list), nrow = nrow(ImmuData_iedb), 
                  dimnames = list(as.character(ImmuData_iedb$Description), as.character(names(orf.list))))
iedb_subject <- matrix(NA, ncol = length(orf.list), nrow = nrow(ImmuData_iedb), 
                  dimnames = list(as.character(ImmuData_iedb$Description), as.character(names(orf.list))))

for(nCol in 1:nrow(orf.df)){
  for (iedb in 1:nrow(ImmuData_iedb)){
    aa1 <- AAString(orf.df$Peptide[nCol])
    aa2 <- AAString(ImmuData_iedb$Description[iedb])
    align <- pairwiseAlignment(aa1, aa2, type = "local", substitutionMatrix = "BLOSUM62", gapOpening = 5, gapExtension = 5)
    iedb_scores[iedb, nCol] <- align@score
    iedb_pattern[iedb, nCol] <- as.character(align@pattern)
    iedb_subject[iedb, nCol] <- as.character(align@subject)
  }
}

iedb_scores_table <- melt(iedb_scores); colnames(iedb_scores_table) <- c("iedb_peptides", "nCol_orf", "value")
iedb_pattern_table <- melt(iedb_pattern); colnames(iedb_pattern_table) <- c("iedb_peptides", "nCol_orf", "pattern")
iedb_subject_table <- melt(iedb_subject); colnames(iedb_subject_table) <- c("iedb_peptides", "nCol_orf", "subject")


all_iedb_nCol_table <- inner_join(iedb_scores_table, iedb_pattern_table, by = c("iedb_peptides", "nCol_orf")) %>% 
  inner_join(iedb_subject_table, by = c("iedb_peptides", "nCol_orf")) %>% 
  mutate_if(is.factor, as.character)
all_iedb_nCol_table <- all_iedb_nCol_table %>% mutate(normalizedValue = value/nchar(as.character(iedb_peptides))) 


# subset those with non-exact matches 
nonmatch_all_iedb_nCol_table <- all_iedb_nCol_table[is.na(ifelse(
  all_iedb_nCol_table$iedb_peptides==all_iedb_nCol_table$subject & 
    all_iedb_nCol_table$subjec == all_iedb_nCol_table$pattern,0,NA)),]
summary(nonmatch_all_iedb_nCol_table$normalizedValue)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.1197  0.9091  1.1500  1.2359  1.4667  5.6667 

match_all_iedb_nCol_table <- all_iedb_nCol_table[which(ifelse(
  all_iedb_nCol_table$iedb_peptides==all_iedb_nCol_table$subject & 
    all_iedb_nCol_table$subjec == all_iedb_nCol_table$pattern,0,NA) ==0),]
summary(match_all_iedb_nCol_table$normalizedValue)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 4.000   4.611   4.944   4.975   5.444   5.733 

# distribution of normalized alignment score 
library(ggforce)
ggplot() + 
  geom_histogram(aes(x = normalizedValue), dplyr::mutate(nonmatch_all_iedb_nCol_table, z = FALSE), bins = 50, fill="darkred") + 
  geom_histogram(aes(x = normalizedValue), dplyr::mutate(nonmatch_all_iedb_nCol_table, z = TRUE), bins = 300, fill="darkred") + 
  facet_zoom(xlim = c(4, 6), ylim = c(0, 8), zoom.data = z,
             horizontal = FALSE) + 
  theme(zoom.y = element_blank(), validate = FALSE) + theme_classic() + 
  labs(title = "Distribution of normalized alignment score for non-exact matching nCol peptides") +
  ylab("Normalized alignment score")  + xlab("Counts") 

ggplot(match_all_iedb_nCol_table, aes(x=normalizedValue)) + 
  geom_histogram( fill="darkred", bins = 60) +
  scale_x_discrete(drop=FALSE) + xlim(c(0, 6)) +
  theme(zoom.y = element_blank(), validate = FALSE) + theme_classic() + 
  labs(title = "Distribution of normalized alignment score for exact matching nCol peptides") +
  ylab("Normalized alignment score")  + xlab("Counts") 


# merging data information - match
match_all_iedb_nCol_table <- match_all_iedb_nCol_table[-which(nchar(match_all_iedb_nCol_table$iedb_peptides) <=3 ),]
match_all_iedb_nCol_epitope_table <- match_all_iedb_nCol_table %>% #28
  inner_join(ImmuData_iedb, by = c("iedb_peptides" = "Description")) 

# merging data information - nonmatch 
nonmatch_all_iedb_nCol_table <- nonmatch_all_iedb_nCol_table[-which(nchar(nonmatch_all_iedb_nCol_table$iedb_peptides) <=3 ),]
nonmatch_all_iedb_nCol_epitope_table <- nonmatch_all_iedb_nCol_table %>%    #361802
  inner_join(ImmuData_iedb, by = c("iedb_peptides" = "Description"))


