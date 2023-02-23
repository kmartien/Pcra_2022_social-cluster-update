rm(list=ls())
library(tidyverse)
library(strataG)

setwd("/Users/Shared/KKMDocuments/Documents/Karen/Structure/Pseudorca/Data/Database files")

msat.dat <- read.csv("AS163.genotypes.csv")
hap.dat <- read.csv("AS163.haplotypes.csv")
seqs <- read.fasta("Pcra unique haps.fasta")
id.key <- read.csv("Pcra.id.key.csv")
clust <- read.csv("social.cluster.assignments.csv")
names(id.key)[2] <- "AnimalID"
names(clust)[1] <- "PhotoID"

dat.df <- left_join(hap.dat, id.key[,c(2,5)], by = "AnimalID") %>% 
  left_join(clust, by = "PhotoID") %>% 
  left_join(msat.dat[,c(1,6:37)], by = "AnimalID") %>%
  data.frame()

msat.g <- df2gtypes(dat.df,ploidy = 2, strata.col = 8, loc.col = 9)
msat.diff.estimates <- pairwiseTest(msat.g)
write.csv(pairwiseSummary(msat.diff.estimates)[,c(1,4,5,20,22,7)], file = "nuc.diff.estimates.csv")

hap.dat <- read.csv("/Users/Shared/KKMDocuments/Documents/Karen/Structure/Pseudorca/Manuscript drafts/social cluster update/soc.clust.all.hap.data.csv")
hap.g <- df2gtypes(hap.dat, ploidy = 1, strata.col = 2, loc.col = 4, sequences = seqs)
hap.diff.estimates <- pairwiseTest(hap.g)
write.csv(pairwiseSummary(hap.diff.estimates)[,c(1,4,5,10,7)], file = "mtdna.diff.estimates.csv")
t(table(dat.df[,c(6,8)]))

