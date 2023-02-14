source("https://raw.githubusercontent.com/melsadany/workbench/master/msmuhammad-source.R", local = T)
library(snpStats)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)


prep_23andme_file <- function(file.23.path, 
                              output.path="../data"){
  geno.data <- read_tsv(file.23.path, skip = 19)
  geno.filtered <- geno.data %>%
    filter(!grepl("i", `# rsid`)) %>%
    mutate(rsid = `# rsid`) %>%
    select(-`# rsid`)
  gpos <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, geno.filtered$rsid, ifnotfound = "drop")
  seqlevelsStyle(gpos) <- "UCSC"
  z <- inferRefAndAltAlleles(gpos, BSgenome.Hsapiens.UCSC.hg19)
  mcols(gpos) <- cbind(mcols(gpos), z)
  rsid.info <- as.data.frame(gpos)
  rsid.info$alt_alleles <- substr(rsid.info$alt_alleles, 0, 1)
  rsid.info <- rsid.info %>%
    mutate(ID_37 = paste0(seqnames, ":", pos, ":", ref_allele, ":", alt_alleles))
  rsid.info <- rsid.info %>%
    mutate(CHR = seqnames) %>%
    mutate(POS = pos) %>%
    mutate(REF = ref_allele) %>%
    mutate(ALT = alt_alleles) %>%
    mutate(rsid = RefSNP_id) %>%
    column_to_rownames("RefSNP_id") %>%
    select(c(rsid, CHR, POS, REF, ALT, ID_37))
  genotypes.filt <- merge(rsid.info, geno.filtered, by = "rsid") %>%
    dplyr::rename(A01_23 = allele2) %>%
    dplyr::rename(A02_23 = allele1) %>%
    dplyr::rename(CHR_23 = chromosome) %>%
    dplyr::rename(POS_23 = position) 
  
  genotypes.filt.2 <- genotypes.filt %>%
    mutate(ALT = toupper(ALT)) %>%
    mutate(ref_a01_match = ifelse(REF==A01_23, T, F)) %>% 
    mutate(ref_a02_match = ifelse(REF==A02_23, T, F)) %>% 
    mutate(alt_a01_match = ifelse(ALT==A01_23, T, F)) %>%
    mutate(alt_a02_match = ifelse(ALT==A02_23, T, F)) %>%
    mutate(score = ifelse(ref_a01_match == T & ref_a02_match == T, 0, NA)) %>%
    mutate(score = ifelse(alt_a01_match == T & alt_a02_match == T, 2, score)) %>%
    mutate(score = ifelse(alt_a01_match == F & alt_a02_match == F & ref_a01_match == F & ref_a02_match == F, 2, score)) %>%
    mutate(score = ifelse(ref_a01_match != ref_a02_match, 1, score)) %>%
    filter(!duplicated(ID_37))
  genotype.file <- genotypes.filt.2 %>%
    select(c(ID_37,score)) %>%
    dplyr::rename(sample = score)
  
  write_tsv(genotype.file, paste0(output.path, "/genotypes-scored.tsv"))
}
