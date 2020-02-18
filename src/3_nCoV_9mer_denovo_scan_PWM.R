### Predict immunogenicity of 9-mer peptides from 2019-nCoV by positional weight matrix 
### Author: Chloe H. Lee

# load library
library(seqinr)
library(plyr)
library(dnar)

# import 9-mer peptides
peptides.dframe <- read_csv("data/nCoV_9mer_peptides.csv")
pep <- as.character(peptides.dframe$peptide)

# compute PWM score 
load("data/pwm.RData")
compute_score <- function(seq, pwm, trim_pos = c(1, 3, 4, 5, 6, 7, 8)){
  if(!is.null(trim_pos)){
    seqMat<- do.call(rbind,strsplit(seq,''))[, trim_pos] #seq <- apply(seqMat,1,paste,collapse="")
    pwm <- pwm[, trim_pos]
  }
  
  ids<-1:length(rownames(pwm) )
  names(ids)<-rownames(pwm)
  probs<-matrix(indexMatrix(ids[as.vector(seqMat)],rep(1:ncol(seqMat),each=nrow(seqMat)),pwm),nrow=nrow(seqMat))
  scores <-apply(probs,1,sum)
  score_df <- data.frame(Peptide = seq, Score = scores, Normalized_score = scores/7)
  
  return(score_df)
} 

binding_score <- compute_score(pep, binding_pwm)
activation_score <- compute_score(pep, activation_pwm)
killing_score <- compute_score(pep, killing_pwm)

# statistics of binding, activation and killing scores 
summary(binding_score$Normalized_score)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2713  0.4648  0.5362  0.5369  0.6043  0.9289 
summary(activation_score$Normalized_score)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1589  0.4134  0.4830  0.4875  0.5575  0.8996 
summary(killing_score$Normalized_score)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1573  0.4489  0.5120  0.5113  0.5726  0.8722 

library(plyr)
score_df <-  rename(binding_score[, c(1, 3)], replace = c("Normalized_score" = "Binding_score"))  %>% 
  join(rename(activation_score[, c(1, 3)], c("Normalized_score" = "Activation_score")), by = "Peptide") %>%
  join(rename(killing_score[, c(1, 3)], c("Normalized_score" = "Killing_score")), by = "Peptide") %>% 
  distinct(Peptide, .keep_all = TRUE) %>%
  mutate(AvgScore = rowMeans(.[,-1])) %>% arrange(desc(AvgScore))

score_df <- score_df %>% 
  mutate(multScore = Binding_score * Activation_score * Killing_score) %>% 
  mutate(geoMean = multScore^(1/3)) 

length(which(score_df$geoMean >= 0.5))
length(which(score_df$geoMean >= 0.6))
length(which(score_df$geoMean >= 0.7))
length(which(score_df$geoMean >= 0.8))


# pwm score statistics
library(ggplot2)
library(ggpubr)
library("gridExtra")

s1 <- ggplot(score_df, aes(x=Binding_score)) + 
  geom_histogram(binwidth=0.01, position="identity", fill="darkorange")+ 
  labs(title="Binding_score",x="Probability", y = "Counts")+
  theme_classic() + theme_pubclean()

s2 <- ggplot(score_df, aes(x=Activation_score)) + 
  geom_histogram(binwidth=0.01, position="identity", fill="darkorange")+ 
  labs(title="Activation_score",x="Probability", y = "Counts")+
  theme_classic() + theme_pubclean()

s3 <- ggplot(score_df, aes(x=Killing_score)) + 
  geom_histogram(binwidth=0.01, position="identity", fill="darkorange")+ 
  labs(title="Killing_score",x="Probability", y = "Counts")+
  theme_classic() + theme_pubclean()

s4 <- ggplot(score_df, aes(x=multScore)) + 
  geom_histogram(binwidth=0.01, position="identity", fill="darkgreen")+ 
  labs(title="Binding * Activation * Killing Score",x="Probability", y = "Counts")+
  theme_classic() + theme_pubclean()

s5 <- ggplot(score_df, aes(x=geoMean)) + 
  geom_histogram(binwidth=0.01, position="identity", fill="darkgreen")+ 
  labs(title="Geometric mean",x="Probability", y = "Counts")+
  theme_classic() + theme_pubclean()

grid.arrange(s1, s2, s3, s4, s5, ncol = 3, nrow = 2)



