##############################################################
#                        PRS calculation tut                 #
##############################################################


info <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  dir = "tmp-data", fname = "map_hm3_ldpred2.rds"))
str(info)


# Load packages bigsnpr and bigstatsr
library(bigsnpr)

# Read from bed/bim/fam, it generates .bk and .rds files.
snp_readBed("tmp-data/public-data3.bed")

# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("tmp-data/public-data3.rds")
# See how the file looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection
NCORES <- nb_cores()
# Read external summary statistics
sumstats <- bigreadr::fread2("tmp-data/public-data3-sumstats.txt")
str(sumstats)

# We split the genotype data using part of the data to choose hyper-parameters and 
# another part of the data to evaluate statistical properties of polygenic risk score 
# such as AUC. Here we consider that there are 350 individuals to be used as 
# validation set to tune hyper-parameters for LDpred2-grid. The other 153 individuals 
# are used as test set to evaluate the final models.

set.seed(1)
ind.val <- sample(nrow(G), 350)
ind.test <- setdiff(rows_along(G), ind.val)


# Matching variants between genotype data and summary statistics
sumstats$n_eff <- sumstats$N
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)  # use rsid instead of pos


# Computing LDpred2 scores genome-wide
# First, you need to compute correlations between variants. 
# POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)
# To avoid downloading "large" files, this has been precomputed
POS2 <- obj.bigSNP$map$genetic.dist

tmp <- tempfile(tmpdir = "tmp-data")

for (chr in 1:22) {
  
  # print(chr)
  
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'G'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                   infos.pos = POS2[ind.chr2]
                   #, ncores = NCORES
                   )
  
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

file.size(corr$sbk) / 1024^3  # file size in GB

# LDpred2-inf: infinitesimal model




















