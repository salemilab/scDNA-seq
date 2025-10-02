
## Initialize stable working environment and store time of initiation
rm(list=ls())


ptm <- proc.time()
# List of packages for session
.packages <-  c("parallel", "optparse", "phylotools", "dplyr") # May need to incorporate code for familyR (https://rdrr.io/github/emillykkejensen/familyR/src/R/get_children.R) i fno longer supported.
#.github_packages <- c("mrc-ide/skygrowth", "GuangchuangYu/treeio")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install CRAN packages (if not already installed)
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# .inst_github <- .packages %in% installed.packages()
# ## Install GitHub packages(if not already installed)
# if(length(.github_packages[!.inst_github]) > 0) try(remotes::install_github(.github_packages[!.inst_github]))
# if(length(.github_packages[!.inst_github]) > 0) try(devtools::install_github(.github_packages[!.inst_github]))

# Load packages into session 
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) BiocManager::install(.packages[!.inst], type = "source", checkBuilt = TRUE)
#.inst_github <- .github_packages %in% installed.packages()
## Install GitHub packages(if not already installed)
#if(length(.github_packages[!.inst_github]) > 0) try(remotes::install_github(.github_packages[!.inst_github]))
#if(length(.github_packages[!.inst_github]) > 0) try(devtools::install_github(.github_packages[!.inst_github]))
## Will need to remove install section if using on cluster ###################################

# Load packages into session 
lapply(.packages, require, character.only=TRUE)
#lapply(gsub(".+\\/(.+)", "\\1", .github_packages), require, character.only=TRUE)

numCores = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
if(is.na(numCores)) {
  numCores=detectCores()
}

option_list = list(
  make_option(c("-m", "--macaque"), type="character", default=NA, help="macaque ID"),
  make_option(c("-d", "--dpi"), type="numeric", default=NA, help="sample collection time point in days post-infection"),
  make_option(c("-t", "--tissue"), type="character", default=NA, help="sampled tissue")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

fa_file=list.files(pattern="\\d.fa")

fa=read.FASTA(fa_file)

fa_headers=names(fa)

fa_barcodes=data.frame(barcode=gsub("Consensus_([A-Z]+-1).+", "\\1", fa_headers),
                       isolate=paste0(opt$dpi, opt$macaque, opt$tissue, "D_", 
                                      sprintf("%04d", seq.int(length(fa_headers)))))

write.csv(fa_barcodes, paste(opt$dpi, opt$macaque, opt$tissue, "header_barcode.csv", sep="_"),
          quote=F, row.names=F)

names(fa) = fa_barcodes$isolate

write.FASTA(fa, paste(opt$dpi, opt$macaque, opt$tissue, "renamed.fa", sep="_"))




