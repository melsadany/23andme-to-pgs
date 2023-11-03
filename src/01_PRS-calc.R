################################################################################
#                                PRS calculation tut                           #
################################################################################
rm(list = ls()); gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(bigsnpr)
################################################################################
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/scratch/PGS"
setwd(project.dir)
args <- commandArgs(trailingOnly = T)
################################################################################
source("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/customized-functions/calculate_ldpred_pgs.R")
library(bigsnpr);set.seed(123);library(data.table);library(magrittr);library(foreach);library(doMC)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = "tmp-data", fname = "map_hm3_plus.rds"))
###
# loop over gwas sumstats
gwas.sumstats.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/data/gwas-sumstats/ALL"
sumstats.meta <- data.frame(phenotype = list.files(gwas.sumstats.dir, full.names = F))
registerDoMC(cores = 10)
i <- as.numeric(args)
# foreach(i=1:nrow(sumstats.meta)) %dopar% {
  ###
  # 3. Load and transform the summary statistic file
  sumstats <- bigreadr::fread2(paste0(gwas.sumstats.dir, "/", sumstats.meta$phenotype[i])) 
  # Filter out hapmap SNPs
  sumstats.f <- sumstats[sumstats$rsid%in% info$rsid,] %>%
    filter(a1 != a0)
  rm(sumstats);gc()
  ###
  plink.rds.path <- "/Dedicated/jmichaelson-wdata/msmuhammad/scratch/PGS/data/derivatives/my_23-merged-w-1KG.hm3.rds"
  print(paste0("starting for: ", sumstats.meta$phenotype[i]))
  pgs <- calculate_ldpred_pgs(sumstats = sumstats.f,
                              plink_rds = plink.rds.path,
                              sd_y = 0,
                              pheno_name = sub(".txt", "", sumstats.meta$phenotype[i]),
                              n_core = 10,
                              build = 'hg19',
                              output_path = "/Dedicated/jmichaelson-wdata/msmuhammad/scratch/PGS/data/derivatives/pgs")
  gc()
  print(paste0("Done for: ", sub(".txt", "", sumstats.meta$phenotype[i])))
  # system("mkdir -p /Dedicated/jmichaelson-wdata/msmuhammad/scratch/PGS/data/derivatives/pgs")
  # write_rds(pgs, paste0("/Dedicated/jmichaelson-wdata/msmuhammad/scratch/PGS/data/derivatives/pgs/",
  #                       sub(".txt", "", sumstats.meta$phenotype[i]), ".rds"))
  # return(pgs)
# }
# write_rds(pgs.you, "/Dedicated/jmichaelson-wdata/msmuhammad/scratch/PGS/data/derivatives/your-raw-pgs.rds")
################################################################################
