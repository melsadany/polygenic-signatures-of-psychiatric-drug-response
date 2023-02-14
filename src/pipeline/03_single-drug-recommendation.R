source("https://raw.githubusercontent.com/melsadany/workbench/master/msmuhammad-source.R", local = T)


recommend_drug <- function(imputed.tx.dir.path="../data/derivatives/imputed-tx",
                           tissue="Brain_Anterior_cingulate_cortex_BA24",
                           output.path="../data/derivatives/drug-recomm",
                           drugs.sig="../data/modeldata/LINCS/cp_mean_coeff_mat_tsv.rds",
                           drugs.metadata="../data/modeldata/LINCS/small_mol_metadata.rds",
                           meds.of.interest=c(),
                           KG.drug.ref.path="../data/modeldata/1KG/1KG-drugs-corr-cmap.rds",
                           threads=4){
  system(paste0("mkdir -p ", output.path))
  tissue.tx <- read_rds(paste0(imputed.tx.dir.path, "/imputed-tx-of-", tissue, ".RDS"))
  smol_meta <- read_rds(drugs.metadata) %>%
    mutate(pert_name=tolower(pert_name))
  mdrug <- read_rds(drugs.sig)
  rownames(mdrug) <- tolower(rownames(mdrug))
  
  if (length(meds.of.interest)>0) {
    meds <- tolower(meds.of.interest)
    drugs.metadata <- smol_meta %>%
      filter(pert_name %in% meds.of.interest) %>%
      distinct(pert_name, .keep_all = T)
    cmap.of.int <- mdrug[which(rownames(mdrug) %in% drugs.metadata$pert_name),]
    rm(mdrug)
    rm(smol_meta)
    gc()
  }else{
    meds <- rownames(mdrug)
    drugs.metadata <- smol_meta %>%
      dplyr::filter(pert_name %in% meds) %>%
      distinct(pert_name, .keep_all = T)
    cmap.of.int <- mdrug[which(rownames(mdrug) %in% drugs.metadata$pert_name),]
    rm(mdrug)
    rm(smol_meta)
    gc()
  }
  gene.int <- intersect(colnames(cmap.of.int), colnames(tissue.tx))
  tissue.tx <- tissue.tx[,gene.int]
  gc()
  cmap.of.int <- cmap.of.int[,gene.int]
  all(colnames(cmap.of.int) == colnames(tissue.tx))
  
  # reference correlations to 1KG correlations
  kg.drug.corr <- read_rds(KG.drug.ref.path)
  colnames(kg.drug.corr) <- tolower(colnames(kg.drug.corr))
  
  X1 <- tissue.tx %>% as.matrix()
  Y1 <- cmap.of.int %>% scale(scale = T, center = T) %>% t() %>% as.matrix()
  drug.corr <- cor(X1, Y1, method = "spearman")
  drug.corr <- -drug.corr
  
  
  library(doMC)
  registerDoMC(cores = threads)
  subject <- "you"
  df <- rbind(kg.drug.corr, drug.corr)
  df2 <- scale(df, scale = T, center = T)
  rownames(df2) <- c(rownames(kg.drug.corr), subject)
  ret <- as.data.frame(t(df2[subject,]))
  rownames(ret) <- subject
  drug.corr.ref <- ret %>% 
    t() %>%
    as.data.frame() %>%
    rownames_to_column("med")

  write_rds(drug.corr.ref, paste0(output.path, "/drug-corr-KG-ref_", tissue, ".rds"))
  write_tsv(drug.corr.ref, paste0(output.path, "/drug-corr-KG-ref_", tissue, ".tsv"))
}


plot_drug_recomm <- function(tissue.drug.recomm.path="../data/derivatives/drug-recomm/drug-corr-KG-ref_Brain_Anterior_cingulate_cortex_BA24.rds",
                             drugs.to.plot,
                             output.figs.path="../figs",
                             tissue="Brain_Anterior_cingulate_cortex_BA24"){
  drug.corr <- read_rds(tissue.drug.recomm.path) %>% column_to_rownames("med")
  
  drug.corr.scaled <- scale(drug.corr, scale = T, center = T)
  personal.drug.corr <- drug.corr.scaled[,"you"] %>% t() %>% as.data.frame()
  colnames(personal.drug.corr) <- rownames(drug.corr.scaled)
  library(fmsb)
  # you're adding 2 rows as your limits. Make sure you maintain the same order
  data <- rbind(rep(round(max(drug.corr.scaled), digits = 1),ncol(personal.drug.corr)), 
                rep(round(min(drug.corr.scaled), digits = 1),ncol(personal.drug.corr)), 
                personal.drug.corr)
  
  if (length(drugs.to.plot)>0) {
    meds <- drugs.to.plot
    data2 <- data[,is.element(colnames(data), meds)]
    svg(paste0(output.figs.path, "/drug-corr_EUR-1KG-ref_", tissue, ".svg"))
    radarchart(data2, axistype=1,
               #custom polygon
               pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.5), plwd=1.5, 
               #custom the grid
               cglcol="grey", cglty=1, axislabcol="grey", 
               caxislabels=seq(round(min(drug.corr.scaled)),round(max(drug.corr.scaled)),2),
               cglwd=0.8,
               #custom labels
               vlcex=0.8
    ) +
      mtext(paste0("Expected drug response\ntissue: ", tissue, "\nreferenced to 1KG"), cex = 0.8, font = 2)
    dev.off()
    
  }else{
    meds <- c("bupropion", "venlafaxin", "trazodone", "sertraline", "isocarboxazid", "fluoxetine", "escitalopram",
              "risperidone", "loxapine",
              "methylphenidate", "clonidine", 
              "propofol", "ketamine")
    data2 <- data[,is.element(colnames(data), meds)]
    svg(paste0(output.figs.path, "/drug-corr_EUR-1KG-ref_", tissue, ".svg"))
    radarchart(data2, axistype=1,
               #custom polygon
               pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.5), plwd=1.5, 
               #custom the grid
               cglcol="grey", cglty=1, axislabcol="grey", 
               caxislabels=seq(round(min(drug.corr.scaled)),round(max(drug.corr.scaled)),2),
               cglwd=0.8,
               #custom labels
               vlcex=0.8
    ) +
      mtext(paste0("Expected drug response\ntissue: ", tissue, "\nreferenced to 1KG"), cex = 0.8, font = 2)
    dev.off()
    
  }
}


