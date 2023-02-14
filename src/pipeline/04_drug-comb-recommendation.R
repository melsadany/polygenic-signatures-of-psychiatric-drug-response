source("https://raw.githubusercontent.com/melsadany/workbench/master/msmuhammad-source.R", local = T)


combine_safe_drugs <- function(drug.corr.path="../data/derivatives/drug-recomm/drug-corr-KG-ref_Brain_Anterior_cingulate_cortex_BA24.tsv",
                               med.of.interest=c(),
                               ddinteractions.path="../data/modeldata/DDInter/ddi-mtx.tsv",
                               output.path="../data/derivatives/drug-recomm"){
  drug.corr <- read_tsv(drug.corr.path) %>%
    mutate(med = tolower(med))
  ddi <- read_tsv(ddinteractions.path) %>%
    column_to_rownames("meds") 
  meds.int <-intersect(colnames(ddi), unique(drug.corr$med)) 
  ddi.filt <- ddi[meds.int,meds.int] 
  ###
  # filter meds.recom
  meds.corr.ref <- drug.corr %>% 
    column_to_rownames("med") %>% 
    t() %>%
    as.data.frame()
  meds.corr.ref.filt <- meds.corr.ref[,meds.int]
  all(colnames(ddi.filt) == colnames(meds.corr.ref.filt))
  ###
  
  if (length(med.of.interest)>0) {
    med <- med.of.interest
  }else{
    med <- colnames(meds.corr.ref.filt)[max.col(meds.corr.ref.filt)]
  }
  ###
  # safe combs that don't interact with best positive responsive drug by participant
  safe.combs <- matrix(nrow = length(meds.int), ncol = nrow(meds.corr.ref.filt), 0)
  rownames(safe.combs) <- meds.int
  colnames(safe.combs) <- rownames(meds.corr.ref.filt)
  
  for (i in 1:ncol(safe.combs)) {
    # p <- colnames(safe.combs)[i]
    # get a vector of drugs interacting with best recommended drug
    # then change the values to 1 for drugs that don't interact with it
    rep.v <- data.frame(inter = ddi.filt[,med%>%as.character()]) %>%
      mutate(inter = ifelse(inter == 0, 1, NA))
    safe.combs[,i] <- rep.v$inter %>% as.numeric()
  }
  # you now have a list of drug that are good to prescribed together
  ###
  # you need to filter the list above to only contain ones you expected to be working by patient
  safe.recom.combs <- safe.combs
  safe.recom.combs.cont <- safe.combs
  for (i in 1:ncol(safe.combs)) {
    p <- colnames(safe.combs)[i]
    # get a vector of drug recommendation for the participant
    # then change the values to 1 for drugs that have positive response
    rec.v <- data.frame(meds.corr.ref.filt[p,] %>%
                          t() %>% as.data.frame() %>%
                          dplyr::rename("recomm" = p) %>%
                          rownames_to_column("med") %>%
                          mutate(tx_rec = ifelse(recomm > 0, 1, NA)))
    
    new.v <- rec.v %>%
      mutate(ddi_v = safe.combs[,i]) %>%
      mutate(final_rec = ifelse(ddi_v == 1 & tx_rec == 1, 1, NA)) %>%
      mutate(final_rec_val = ifelse(final_rec == 1, recomm, NA))
    safe.recom.combs[,i] <- new.v$final_rec %>% as.numeric()
    safe.recom.combs.cont[,i] <- new.v$final_rec_val %>% as.numeric()
  }
  
  write_tsv(safe.recom.combs.cont%>%as.data.frame()%>%rownames_to_column("meds"), paste0(output.path, "/safe-recom-meds-combs.tsv"))
  
}