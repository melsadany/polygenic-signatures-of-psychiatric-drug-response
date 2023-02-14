source("pipeline/01_23andme-prep.R")
source("pipeline/02_tx-impute.R")
source("pipeline/03_single-drug-recommendation.R")
source("pipeline/04_drug-comb-recommendation.R")

raw.file <- readline(prompt = "Enter the path for your 23andme data: ")
# raw.file <- "../data/phased_genome_statistical_20220923172100_Muhammad_Elsadany_v5_Full_20220925180737.txt"

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

prep_23andme_file(file.23.path = raw.file, output.path = "../data")
impute_tx_from_23(genotypes.file.path = "../data/genotypes-scored.tsv")
recommend_drug()
plot_drug_recomm(drugs.to.plot = c("sertraline", "methylphenidate","trazodone", "venlafaxine", "melatonin"))
combine_safe_drugs(med.of.interest = "sertraline")
