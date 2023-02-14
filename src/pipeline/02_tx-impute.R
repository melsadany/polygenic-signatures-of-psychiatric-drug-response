source("https://raw.githubusercontent.com/melsadany/workbench/master/msmuhammad-source.R", local = T)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)
library(foreach)
library(doMC)
library(data.table)

impute_tx_from_23 <- function(genotypes.file.path = "../data/genotypes-scored.tsv", 
                              output.path = "../data/derivatives/imputed-tx", 
                              model.weights.dir = "../data/modeldata/UTMOST-GTEx-model-weights",
                              tissue="Brain_Anterior_cingulate_cortex_BA24",
                              threads = 4){
  start.time <- Sys.time()
  system(paste0("mkdir -p ", output.path))
  tissue.extracted.genotypes <- paste0(output.path, "genotypes-for-tissue-", tissue, ".xmat.gz")
  weights.dir <- model.weights.dir
  # pull tissue weights
  tissue.weights <- read.csv(file = paste0(weights.dir, "/", tissue, ".db.csv"), header = TRUE, row.names = 1)
  tissue.weights <- unique(tissue.weights)
  ##### 
  # rsid annotation section
  tmp.dir <- paste0(weights.dir, "/tmp")
  system(paste0("mkdir -p ", tmp.dir))
  
  # step 1
  tmp.rsid.filename <- paste0(tmp.dir, "/rsid-for-", tissue)
  tissue.rsid.info.filt <- read.table(tmp.rsid.filename, row.names = 1)
  
  genotypes.file <- read_tsv(genotypes.file.path) %>%
    column_to_rownames("ID_37")
  dim(genotypes.file)
  gc()
  
  
  #####
  # start of transcript imputation section 
  tissue.rsid.info.filt.2 <- tissue.rsid.info.filt[is.element(tissue.rsid.info.filt$ID_02_UTMOST, rownames(genotypes.file)),]
  
  # filter genotypes matrix
  genotypes.file.2 <- genotypes.file %>%
    rownames_to_column("ID_37") %>%
    filter(ID_37 %in% tissue.rsid.info.filt.2$ID_02_UTMOST) %>%
    column_to_rownames("ID_37")
  gc()
  
  # create a matrix for all tissue weights 
  # that matrix has genotypes found in cohort as rows and all genes affected by them as columns
  # if gene is not affected, it has a zero value
  # if it is, replace the zero with its weight from the tissue weights df
  
  tissue.rsid.info.filt.3 <- as.data.frame(matrix(0, nrow = length(unique(tissue.rsid.info.filt.2$ID_02_UTMOST)), ncol = length(unique(tissue.rsid.info.filt.2$gene))))
  colnames(tissue.rsid.info.filt.3) <- unique(tissue.rsid.info.filt.2$gene)
  rownames(tissue.rsid.info.filt.3) <- unique(tissue.rsid.info.filt.2$ID_02_UTMOST)
  
  c <- threads
  registerDoMC(c)
  
  for (i in 1:nrow(tissue.rsid.info.filt.2)) {
    tissue.rsid.info.filt.3[tissue.rsid.info.filt.2$ID_02_UTMOST[i], tissue.rsid.info.filt.2$gene[i]] <- tissue.rsid.info.filt.2$weight[i]
    # print(i)
  }
  print(paste("done with creating a matrix for tissue weights of",
              nrow(tissue.rsid.info.filt.3), "rows and", ncol(tissue.rsid.info.filt.3),
              "col", sep = " "))
  tissue.weights.all.matrix <- as.matrix(tissue.rsid.info.filt.3)
  rownames(tissue.weights.all.matrix) <- rownames(tissue.rsid.info.filt.3)
  colnames(tissue.weights.all.matrix) <- colnames(tissue.rsid.info.filt.3)
  
  
  tmp.genotypes.file.2 <- t(genotypes.file.2) %>% as.data.frame()
  
  # reorder the genotypes matrix to match ordering of tissue weights
  genotypes.file.r <- tmp.genotypes.file.2 %>%
    select(rownames(tissue.weights.all.matrix)) 
  # multiply matrices
  if (ncol(genotypes.file.2) > 0) {
    genotypes.file.r.matrix <- as.matrix(genotypes.file.r)
    colnames(genotypes.file.r.matrix) <- colnames(genotypes.file.r)
    if (any(is.element(rownames(tissue.weights.all.matrix), colnames(genotypes.file.r.matrix)))) {
      library(dplyr)
      print(paste0("found ", ncol(genotypes.file.r.matrix), " genotypes influencing expression of ", ncol(tissue.weights.all.matrix), " genes"))
      # this filtered.2 doesn't have column names (fixed)
      genotypes.file.matrix.filtered.2 <- matrix(as.numeric(genotypes.file.r.matrix), ncol = ncol(genotypes.file.r.matrix))
      colnames(genotypes.file.matrix.filtered.2) <- colnames(genotypes.file.r.matrix)
      print("I'm about to multiply matrices now")
      imputed.tx <- genotypes.file.matrix.filtered.2 %*% tissue.weights.all.matrix
      colnames(imputed.tx) <- colnames(tissue.weights.all.matrix)
      print("done with multiplication step, phew")
    }else {
      print("no genotypes found")
      imputed.tx <- matrix(0, ncol = ncol(tissue.weights.all.matrix), nrow = nrow(tmp.genotypes.file.2))
      rownames(imputed.tx) <- genotypes.file$IID
      colnames(imputed.tx) <- colnames(tissue.weights.all.matrix)
    }
  }else {
    print("no genotypes found")
    imputed.tx <- matrix(0, ncol = ncol(tissue.weights.all.matrix), nrow = nrow(tmp.genotypes.file.2))
    rownames(imputed.tx) <- genotypes.file$IID
    colnames(imputed.tx) <- colnames(tissue.weights.all.matrix)
  }
  
  imputed.tissue.tx.fname <- paste0(output.path, "/imputed-tx-of-", tissue)
  write_rds(imputed.tx, file = paste0(imputed.tissue.tx.fname, ".RDS"))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste0("Done imputing transcription in: ", time.taken))
  
}
