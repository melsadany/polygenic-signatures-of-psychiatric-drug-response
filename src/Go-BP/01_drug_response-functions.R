################################################################################
#                       pipeline for drug response prediction                  #
################################################################################
# rm(list = ls())
gc()
# .libPaths("/Users/msmuhammad/workbench/miniconda3/envs/tximpute2/lib/R/library")
# source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(tidyverse, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/tximpute/lib/R/library")
library(doMC, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/tximpute/lib/R/library")
library(readr, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/tximpute/lib/R/library")
library(data.table, lib.loc = "/Dedicated/jmichaelson-wdata/msmuhammad/workbench/miniconda3/envs/tximpute/lib/R/library")
################################################################################
################################################################################
# functions you need
################################################################################
############################ tx imputation section #############################
################################################################################
tissue <- "Brain_Frontal_Cortex_BA9"
impute.tx <- function(genotypes, weights, threads) {
  # the genotypes matrix is expected to have participants as rownames and genotypes as colnames
  # the weights matrix is expected to have 3 columns: variant, gene, weight
  ge <- intersect(colnames(genotypes), weights$variant)
  if (length(ge)==0) {
    print("no genotypes found to have weight")
    return(NULL)
  }
  filt.genotypes <- data.frame(lapply(genotypes, function(x) as.numeric(x)))[,-1]
  rownames(filt.genotypes) <- genotypes$IID
  colnames(filt.genotypes) <- colnames(genotypes)[-1]
  filt.genotypes <- filt.genotypes[,ge]
  gc()
  
  filt.weights <- weights %>%
    filter(variant %in% ge)
  gc()
  registerDoMC(cores = threads)

  imputed <- foreach(j=1:length(unique(filt.weights$gene)), .combine = cbind) %dopar% {
    # j=1
    gen <- unique(filt.weights$gene)[j]
    gene.weights <- filt.weights %>% filter(gene%in%gen) %>%
      pivot_wider(names_from = "variant", values_from = "weight") %>%
      column_to_rownames("gene")
    imp <- as.matrix(filt.genotypes[,colnames(gene.weights)]) %*% t(gene.weights)
    return(imp)
  }
  return(imputed)
}
################################################################################
############################ drug response section #############################
################################################################################
# drug response function using signature matrix
jaccard.index <- function(subject, drug) {
  # I'm assuming genes are in the first column of each of the input vectors/dataframes
  # I'm also assuming that the subject ID and the drug name are the column names of the second column
  df <- right_join(subject, drug)
  colnames(df) <- c("gene", "subject", "drug")
  df <- df %>%
    mutate(subject = ifelse(is.na(subject), 0, subject)) %>%
    mutate(UU = ifelse(subject>0 & drug>0, 1, 0)) %>%
    mutate(DD = ifelse(subject<0 & drug<0, 1, 0)) %>%
    mutate(UD = ifelse(subject>0 & drug<0, 1, 0)) %>%
    mutate(DU = ifelse(subject<0 & drug>0, 1, 0))
  UU <- df %>% filter(subject > 0 | drug > 0)
  DD <- df %>% filter(subject < 0 | drug < 0)
  UD <- df %>% filter(subject > 0 | drug < 0)
  DU <- df %>% filter(subject < 0 | drug > 0)
  jaccard <- ((sum(df$UU)/nrow(UU))+(sum(df$DD)/nrow(DD))-(sum(df$UD)/nrow(UD))-(sum(df$DU)/nrow(DU)))/2
  df2 <- data.frame(drug = jaccard)
  colnames(df2)[1] <- colnames(drug)[2]
  return(df2)
}
predict.response <- function(samples.tx, drug.sig, approach = 1, threads=6, set) {
  # assume a samples.tx is dataframe with genes as rownames
  # assume a drug.sig is a datframe with genes as rownames
  if (approach == 1) {
    genes.data <- data.frame(gene = rownames(drug.sig), drug = 1, participant_tx=0) %>%
      mutate(participant_tx =ifelse(gene %in% rownames(samples.tx), 1,0))
    if(set == "targets") {
      genes.data <- genes.data %>%
        filter(gene %in% c("SLC6A3","SLC6A2","HTR1A","CES1A1a"))
      drug.sig <- drug.sig %>% as.data.frame() %>%
        rownames_to_column("gene") %>%
        filter(gene %in% c("SLC6A3","SLC6A2","HTR1A","CES1A1a")) %>%
        column_to_rownames("gene")
    }
    samples.missing <- genes.data %>% filter(participant_tx==0) %>% dplyr::select(gene)
    samples.tx.full <- rbind(samples.tx, matrix(nrow = length(samples.missing[,1]), 
                                                ncol = ncol(samples.tx), 
                                                dimnames = list(samples.missing[,1], colnames(samples.tx)), 0))
    samples.tx.full <- samples.tx.full[rownames(drug.sig),]
    predicted.response <- cor(as.matrix(samples.tx.full), as.matrix(drug.sig))
    return(predicted.response)
  } else if (approach == 2) {
    int.genes <- intersect(rownames(samples.tx), rownames(drug.sig))
    drug.sig.filt <- drug.sig[int.genes,]
    samples.tx.filt <- samples.tx %>% as.data.frame() %>%
      rownames_to_column("gene") %>%
      filter(gene %in% int.genes) %>%
      column_to_rownames("gene")
    registerDoMC(cores=threads)
    predicted.response <- foreach( i=1:ncol(samples.tx.filt), .combine = rbind) %dopar% {
      # i=5
      participant <- colnames(samples.tx.filt)[i]
      response <- jaccard.index(subject = samples.tx.filt%>%as.data.frame()%>%
                                  rownames_to_column("gene")%>%dplyr::select(gene, i+1),
                                drug = drug.sig%>%as.data.frame()%>%rownames_to_column("gene")%>%
                                  dplyr::select(gene,2))
      ret <- data.frame(IID = participant, drug_corr = response)
      # print(ret)
      return(ret)
    }
    return(predicted.response)
  }
}
################################################################################
################################################################################
################################################################################
################################################################################
Go.BP <- function(from, genes_set, correct, cores, scale, method, genotypes_path) {
  ################################################################################
  # read genotypes for tissue/celltype
  ### choose
  if (from == "tissue") {
    tissue <- "Brain_Frontal_Cortex_BA9"
    genotypes <- fread(file = genotypes_path, header = T, nThread = cores)
    genotypes <- genotypes[-1,-1]
    gc()
    print(paste0("Done with: ", "reading genotypes file for tissue"))
    # get tissue weights
    tissue.weights <- read.table(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/projects/tx-imputation/UTMOST-GTEx-model-weights/tmp/rsid-for-", tissue), row.names = 1)
    ready.weights <- tissue.weights %>% dplyr::select(variant = ID_02_UTMOST, gene, weight) %>%
      distinct(variant, gene, .keep_all = T) %>%
      filter(variant %in% colnames(genotypes))
    gc()
    print(paste0("Done with: ", "reading weights file for tissue"))
  } else if (from == "celltype") {
    genotypes <- fread(file = genotypes_path, header = T, nThread = cores)
    genotypes <- genotypes[-1,-1]
    gc()
    print(paste0("Done with: ", "reading genotypes file for celltype"))
    # get celltype weights
    celltype <- "Excitatory"
    celltype.weights <- read_rds(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/", celltype,"-weights-fdr-sig.rds"))
    ready.weights <- celltype.weights %>%
      filter(FDR<0.05) %>%
      dplyr::select(variant=ID_37, gene, weight=beta) %>%
      distinct(variant, gene, .keep_all = T) %>%
      filter(variant %in% colnames(genotypes))
    rm(celltype.weights)
    gc()
    print(paste0("Done with: ", "reading weights file for celltype"))
  }
  ################################################################################
  # drug genes or all
  if (genes_set == "targets") {
    # impute tx for all genes
    imputed.tx <- impute.tx(genotypes = genotypes, 
                            weights = ready.weights%>%filter(gene %in% c("SLC6A3","SLC6A2","HTR1A","CES1A1a")), 
                            threads = cores)
    gc()
    print(paste0("Done with: ", "imputing tx"))
    imputed.tx <- imputed.tx %>% as.data.frame() %>% dplyr::select(any_of(c("SLC6A3","SLC6A2","HTR1A","CES1A1a")))
    print(paste0("Done with: ", "selecting targets"))
    if (is.null(imputed.tx)) {
      print("exit")
      return(NULL)
    }
  }else if (genes_set == "all") {
    # impute tx for all genes
    imputed.tx <- impute.tx(genotypes = genotypes, 
                            weights = ready.weights, threads = cores) %>% as.data.frame()
    gc()
    print(paste0("Done with: ", "imputing tx"))
    # imputed.tx <- imputed.tx %>% as.data.frame() 
    print(paste0("Done with: ", "keeping all genes"))
    if (is.null(imputed.tx)) {
      print("exit")
      return(NULL)
    }
  }
  ################################################################################
  # raw impute or PCs corrected
  if (correct == F) {
    tx.corrected <- imputed.tx %>% as.data.frame() 
    print(paste0("Done with: ", "not correcting tx for genetic PCs"))
  } else if (correct == T) {
    geno.pcs <- read_tsv("/Dedicated/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/PCA/all/PCs.tsv") %>%
      filter(IID %in% rownames(imputed.tx))
    rownames(geno.pcs) = geno.pcs$IID
    # correct for genetic pcs
    registerDoMC(cores)
    tx.corrected <- foreach(i=1:ncol(imputed.tx), .combine = cbind) %dopar% {
      # i=1
      gene <- colnames(imputed.tx)[i]
      df <- data.frame(geno.pcs[,1:6],exp=imputed.tx[rownames(geno.pcs),i])
      rownames(df) <- rownames(geno.pcs)
      df$corr <- residuals(glm(as.numeric(exp) ~ pc_01 + pc_02 + pc_03 + pc_04 + pc_05, data = df))
      colnames(df)[8] <- gene
      ret <- df%>%dplyr::select(8)
      return(ret)
    }
    print(paste0("Done with: ", "correcting imputed tx for genetic PCs"))
  }
  ################################################################################
  # corr or jaccard AND raw person or scaled
  cmap.drug.sig <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/cmap.of.int.rds")
  mph.sig <- cmap.drug.sig["methylphenidate",]%>%as.data.frame()%>%rename(methylphenidate=1)
  print(paste0("Done with: ", "reading mph signature"))
  if (scale == T ) {
    if (method == 1) {
      m <- predict.response(samples.tx = if (ncol(tx.corrected)>1) {t(scale(tx.corrected))} else {t(tx.corrected)}, 
                            drug.sig = scale(mph.sig), approach = method, threads = cores, set = genes_set) %>%
        as.data.frame() %>%
        rownames_to_column("IID") %>%
        rename(m = 2)
    }else {
      m <- predict.response(samples.tx = if (ncol(tx.corrected)>1) {t(scale(tx.corrected))} else {t(tx.corrected)}, 
                            drug.sig = scale(mph.sig), approach = method, threads = cores, set = genes_set) %>%
        as.data.frame() %>%
        rename(m = 2)
    }
    print(paste0("Done with: ", "predicting drug response with scaling option and approach ", method))
  } else if(scale == F){
    if (method == 1) {
      m <- predict.response(samples.tx = t(tx.corrected), drug.sig = scale(mph.sig), approach = method, threads = cores, set = genes_set) %>%
        as.data.frame() %>%
        rownames_to_column("IID") %>%
        rename(m = 2)
    }else {
      m <- predict.response(samples.tx = t(tx.corrected), drug.sig = scale(mph.sig), approach = method, threads = cores, set = genes_set) %>%
        as.data.frame() %>%
        rename(m = 2)
    }
    print(paste0("Done with: ", "predicting drug response without scaling and approach ", method))
  }
  return(m)
}
################################################################################
example <- 'Go.BP(from = "tissue", genes_set = "targets", correct = T, cores = 18, scale = T, method = 1, 
                  genotypes_path = "data/derivatives/spark-mph-samples-genotypes-LC-merged-Brain_Frontal_Cortex_BA9.xmat.gz")'
