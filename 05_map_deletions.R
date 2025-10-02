
rm(list=ls())

.packages=c("parallel", "phylotools", "dplyr", "tidyverse", "ape", 
            "cluster", "dendextend", "treeio", "ggtree", 
            "scales", "cowplot", "fpc", "pbmcapply", "data.table",
            "png", "patchwork")

lapply(.packages, require, character.only=T, quietly = TRUE)

`%notin%` <- Negate(`%in%`)

numCores = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
if(is.na(numCores)) {
  numCores=detectCores()
}

if("optparse" %in% installed.packages()) {
  library(optparse)
} else {
  install.packages("optparse")
}

.option_list = list( 
  make_option(c("-c", "--clust"), type="character",  default="N",
              help="Option to choose clustering method (default) or simply top 10 most representative sequences"),
  make_option(c("-t", "--tol.level"), type="numeric", default=10,
              help="number of clusters to choose from"),
  make_option(c("-p", "--panel"), type="character", default="CO243",
              help="primer panel"),
  make_option(c("-f", "--full"), type="character", default="N",
              help="specifies whether working with single sample or full alignment")
  
)

.opt_parser = OptionParser(option_list=.option_list)
opt = parse_args(.opt_parser)


print(opt)

## Note: will need to change values of cluster sizes if change the tol.level!

## Test dataset ##############################################
# seqs <- list(
#   seq1="AAAAAAAAAA----------AAAAAAAAAA----------AAAAAAAAAA",
#   seq2="AAAAAAAAA-----------AAAAAAAAAA-----------AAAAAAAAA",
#   seq3="----------AAAAAAAAAA----------AAAAAAAAAA----------",
#   seq4="---------AAAAAAAAAAA----------AAAAAAAAAAA---------",
#   seq5="AAAA---AAA----------AAA---AAA-----------AAA---AAAA",
#   seq6="AAAA---AAAA--------AAAA---AAAA----------AAA---AAAA",
#   seq7="AAAAAAAAAA----------AAA---AAA-----------AAA---AAAA",
#   seq8="AAAAAAAAAA----------AAA---AAA-----------AAAAAAAAAA"
# )
############################################################### 

#setwd("/Users/macbook/Dropbox (UFL)/SIV/scDNA-seq/DNA+protein/Figures_LED/R_make_figs/TrimmedR2s/SIV_alignments")

fa_file = list.files(path="aln", pattern=".aln$", full.names=T)
seqs = read.dna(fa_file, format="fasta", as.matrix=FALSE, as.character=TRUE)

if (opt$full=="Y"){
#  seqs=seqs[!grepl("R_", names(seqs))]  
  metadata_file = list.files(pattern="metadata")
    meta = read.delim(metadata_file, sep='\t', header=F, col.names=c("id", "panel", "seq")) 
}

seqs=mclapply(seqs, toupper, mc.cores=numCores)

new_seqs=mclapply(seqs, function(x) {
  x=x[1:9501]
  x=str_c(x, collapse="")
  xl=nchar(str_match(x, "^[-N]{200,}"))
  tryCatch(if(!is.na(xl)) {
    x=gsub(paste0("^([-N]{",xl,",})"), paste(rep("-", xl), collapse=""), x)
  } else{ 
    x=x }, error=function(e) NULL)
  xl2=nchar(str_match(x, "N{200,}-*$"))
  tryCatch(if(!is.na(xl2)) {
    x=gsub(paste0("(N{",xl2,",}-*)$"), paste(rep("-", xl2), collapse=""), x)
  } else {
    x=x }, error=function(e) NULL)
  # Below runs out of memory!
  # xl_internal=lapply(regmatches(x, gregexpr("N{400,}", x)), nchar)
  # xl_internal[[1]]=sort(xl_internal[[1]], decreasing=TRUE)
  # tryCatch(if(!is.na(xl_internal)) {
  #   i=1
  #   while(i<=length(xl_internal[[1]])) {
  #     len=xl_internal[[1]][i]
  #     x=gsub(paste0("(N{",len,",})"), paste(rep("-", len), collapse=""), x)
  #     i=i+1
  #   } #End while loop
  # } else {
  #   x=x }, error=function(e) NULL)
  
  xl_internal=gregexpr("N{400,}", x)
  tryCatch(if(!is.na(xl_internal)) {
    coords=data.frame(start=xl_internal[[1]],
                      size=lapply(xl_internal, attributes)[[1]]$match.length) %>%
      mutate(end=start+size-1)
    x=str_split(x, "")[[1]]
    i=1
    while(i<=nrow(coords)) {
      start=coords$start[i]
      end=coords$end[i]
      size=coords$size[i]
      x[start:end] = rep("-", size)
      i=i+1
    } # End while loop
    x=str_c(x, collapse="")
  } else {
     x=x }, error=function(e) NULL)
  return(x)
}, mc.cores=numCores)

 
#coded_seqs <- lapply(seqs, strsplit, NULL) ## Currently looks like "ACGTACGT", but needs to look like "A" "C" "G" "T"....
#coded_seqs <- lapply(seqs, unlist) # Even though read.dna gives us a list, it's a list of lists, or nested list, and we need it to be a list of vectors.
#ref <- coded_seqs[names(coded_seqs)=="KU892415.2"] #Pull out reference sequence
#ref <- ref[[1]] # We want a vector instead of a list
#coded_seqs <- coded_seqs[names(coded_seqs)!="KU892415.2"]# Now remove reference from sequences.

ref=new_seqs[[1]]
filtered_seqs=new_seqs[-1]

filtered_seqs = lapply(filtered_seqs, strsplit, NULL) ## Currently looks like "ACGTACGT", but needs to look like "A" "C" "G" "T"....
filtered_seqs = lapply(filtered_seqs, unlist) # Even though read.dna gives us a list, it's a list of lists, or nested list, and we need it to be a list of vectors.
filtered_seqs = mclapply(filtered_seqs, function(x) {
  xl = length(x[x!="-" & x!="N"])
   if(xl > 400) { # At least two amplicons present (min size is now 200)
    x=x
  } else {
    x=NULL
  }
  
  return(x)
}, mc.cores=numCores)
filtered_seqs = Filter(Negate(is.null), filtered_seqs)

basename=gsub("_filtered.fa.aln", "", paste0(fa_file))

fa_headers=names(filtered_seqs)

# barcodes=read.csv(list.files(pattern="header_barcode.csv")) %>%
#   dplyr::filter(isolate %in% fa_headers) %>%
#   dplyr::mutate(barcode=gsub("-1", "", barcode))
# 
# write.csv(barcodes, "barcodes_final.csv", row.names=F, quote=F)

hypermuts=read.table("hypermut_seqs.txt", sep='\t', header=T) %>%
  dplyr::filter(P.val<=0.05)

remove_spaces <- function(df) {
  df[] <- lapply(df, function(x) gsub(" ", "", x))
  return(df)
}

hypermuts=remove_spaces(hypermuts)
hypermuts = as.character(hypermuts$isolate)

print("Masking hypermutated sites and saving final alignment...")
filtered_seqs = pbmclapply(1:length(filtered_seqs), function(seq) {
  x=filtered_seqs[seq][[1]]
  if(names(filtered_seqs[seq]) %in% hypermuts) {
      for (i in 2:length(ref)) { # For each base in the ref sequence (starting with the second one)...
        if(isTRUE(grepl("GG|GA", paste0(ref[i-1], ref[i])) &
        grepl("AG|AA", paste0(x[i-1], x[i]) ) ) ) { # If the corresponding nucleotide pair in sequence x is hypermutated....
            x[i-1]="N"
            x[i]="N"
        } # End if-else statement
      } # End for loop
  } # End if-else statement
#  names(x) = names(filtered_seqs[seq])
  return(x)
}, mc.cores=numCores)
names(filtered_seqs) = fa_headers

write.dna(filtered_seqs, paste0(basename, "_final.fa"), format="fasta")
write.dna(filtered_seqs[names(filtered_seqs) %notin% hypermuts], 
        paste0(basename, "_final_hypermutsremoved.fa"), format="fasta")







################################################################################
# Determine coordinates of bases not covered by primers
print("Reading in amplicon panel data...")
if(opt$panel=="CO243") {
  panel_summary=tryCatch(read.table("/orange/salemi/share/data/Tapestri_Runs/panels/hPSC_CO243_2/rheMac10_CO243.2.amplicons",
                           header=F, sep='\t', col.names = c("genome", "start", "end", "amplicon")), error=function(e) {
                             NULL })
  if (is.null(panel_summary)) {
    panel_summary = read.table("../CO243.amplicons.bed", header=F, sep='\t', col.names = c("genome", "start", "end", "amplicon"))
  }
  
} else{
  if(opt$panel=="UF.HIV") {
  panel_summary=tryCatch(read.table("/blue/salemi/share/tapestri/panels/SIV_UF.HIV1/UF.amplicons.bed",
                                    header=F, sep='\t', col.names = c("genome", "start", "end", "amplicon")), error=function(e) {
                                      NULL })
  if (is.null(panel_summary)) {
    panel_summary = read.table("../UF.amplicons.bed", header=F, sep='\t', col.names = c("genome", "start", "end", "amplicon"))
  }
  } else {print("Incorrect Tapestri panel specified. Using default of CO243.")
  }
}

# Amps 31 and 35 don't seem to be working anywhere, so remove...
panel_summary=filter(panel_summary, amplicon!= "UF_AMP31" & amplicon!="UF_AMP35" & genome=="SIV_consensus2")
panel_summary$missing_start=NA
panel_summary$missing_end=NA
panel_summary$missing_start[1]=1
panel_summary$missing_end[1]=panel_summary$start[1]-1

for (i in 2:nrow(panel_summary)) {
  panel_summary$missing_start[i]=panel_summary$end[i-1]+1
  panel_summary$missing_end[i]=panel_summary$start[i]-1
} 

panel_summary=panel_summary[1:(nrow(panel_summary)-2),]
################################################################################                        
print("Reading in gene annotation data...")

gene_anno = tryCatch(
  read.table("/blue/salemi/share/tapestri/annotated_ref.txt", 
             header = T, sep='\t'), error = function(e) {
               NULL
             })
if (is.null(gene_anno)) {
  gene_anno = read.table("../annotated_ref.txt", header = T, sep='\t')
}

print("Searching first for point mutations...")
point_muts=mclapply(filtered_seqs, function(x) {
  start_coords = gene_anno$SIV[gene_anno$start=="start"]
  stop_coords = gene_anno$SIV[gene_anno$end=="end"]
  stop_codons = list(c("TAA"), c("TAG"), c("TGA"))
  
  result = mclapply(3:length(x), function(site) {
    codon=paste0(c(x[site-2], x[site-1], x[site]), collapse="")
    if( grepl("[-N]", codon) ) {
      return("pass") 
    } else {
      if (all(c(site-2, site-1, site) %in% start_coords) &
          codon != "ATG") {
        return("point")
      } else {
        if (all(c(site-2, site-1, site) %in% stop_coords) &
            codon %notin% stop_codons) {
          return("point")
        } else {
          return("pass")
        }
      }
    }
  }, mc.cores=numCores)
  
  if ("point" %in% result) {
    result = "point"
  } else {
    result = "pass"
  }
  
  #Premature stop codon in nef
  if ( paste0(x[9479:9481], collapse="") == "TAA" ) {
    result = "point"
  }
  
  return(result)
}, mc.cores=numCores)

# 5 = Not covered
# 0 = 2-LTR
# 1 = Base
# 2 = N
# 3 = hyper-A
# 4 = '-'
# 6 = UF.HIV primer

print("Categorizing 2-LTR and hypermutants...")
# First code non-included regions as '5
coded_seqs = mclapply(filtered_seqs, function(x) {
  for (i in 1:nrow(panel_summary)) {
    start=panel_summary$missing_start[i]
    end=panel_summary$missing_end[i]
    x[start:end] = gsub(".", "5", x[start:end])
  }
  return(x)
}, mc.cores=numCores)

# Next code 2-LTR and individual bases
coded_seqs <- pbmclapply(coded_seqs, function(x) { 
  x <- x[1:9501] 

  
  ## RRE is located at coordinates 8277-8392 in our reference (alignment in RRE folder) and used by Bruner to probe for hypermutants.
  ## Below is separate code to analyze percentages of these mutations in the entire RRE,
  ## Or, alternatively, just the two mutations used by Bruner to distinguish hypermutants.
  ## Either way, neither method will work because we have all ambiguous sites in this region.

  # First approach:
  
  # RRE_start=8277
  # RRE_end=8392
  
  
  # Extract RRE sequence region from both reference and sequence x
  #ref_RRE=ref[RRE_start:RRE_end]
  #x_RRE=x[RRE_start:RRE_end]
  
  # Create empty score sheets
  # score_ref=NA
  #  score_x=NA
  # for (i in 2:length(ref_RRE)) { # For each base in the RRE sequence (starting with the second one)...
  #   score_ref[i]=as.integer(grepl("gg|ga", paste0(ref_RRE[i-1], ref_RRE[i]))) #If the GA or AG is present in the reference sequence score it as 1...
  #   if(isTRUE(grepl("gg|ga", paste0(ref_RRE[i-1], ref_RRE[i])) &
  #      grepl("ag|aa", paste0(x_RRE[i-1], x_RRE[i]) ) ) ) { # If the corresponding nucleotide pair in sequence x is hypermutated....
  #     score_x[i]=1 # Score it as 1...
  #   } else {
  #     score_x[i]=0 # If not, score it as 0.
  #   }
  # }  
  # 
  #  if(sum(score_x, na.rm=T) > 0.50*sum(score_ref, na.rm=T)) { ## If the sum of the scores for sequence x is >90% of the sum of scores for ref...
  #    x[RRE_start:RRE_end]=rep("2", length(x[RRE_start:RRE_end])) ## Replace entire RRE string with "2"
  #  } else {
  #    x[RRE_start:RRE_end]=x[RRE_start:RRE_end] # If not, then keep string as is
  #  }
  
  
  ## Second approach using probe stretch instead (see figure in paper):
  
  # RRE_start=8322
  # RRE_end=8357
  # 
  # ref_RRE=ref[RRE_start:RRE_end]
  # x_RRE=x[RRE_start:RRE_end]
  # 
  # 
  # if(isTRUE(grepl("ag|aa", paste0(x_RRE[4], x_RRE[5]) , ignore.case=T) &
  #           grepl("ag|aa", paste0(x_RRE[33], x_RRE[34]), ignore.case=T )) ) { # If the corresponding nucleotide pairs in sequence x is hypermutated....
  #   x[RRE_start:RRE_end]=rep("3", length(x[RRE_start:RRE_end])) ## Replace entire RRE string with "2"
  # } else {
  #   x[RRE_start:RRE_end]=x[RRE_start:RRE_end] # If not, then keep string as is
  # }

  ## AMP1 (2-LTR) as separate category
  
  AMP1_start=7
  AMP1_end=250


  x[AMP1_start:AMP1_end] <- gsub("[A-Za-z]", "0", x[AMP1_start:AMP1_end]) # Replace all bases in 2-LTR junction with  "0"
  
  
  # match <- paste(rep("N", 200), collapse="")
  #for (i in  200:length(x)) {
  #  if (isTRUE(grepl(match, paste(x[(i-199):i], collapse=""), ignore.case=T))) {
  #    x[(i-199):i] <- rep("2", 200)
  #  } else {
  #    x[(i-199):i] <-x[(i-199):i]
  #  }
  #}
  
  
  x <- gsub("N", "2", x) # Now replace all Ns 
  x <- gsub("[A-Za-z]", "1", x) # Now replace all bases 
  x <- gsub("-", "4", x) # And replace all gaps 
  
  
  return(x)
}, mc.cores=numCores) 

# Next code hypermutants 
if(length(hypermuts)!=0) {
  nm = names(coded_seqs)
  coded_seqs = mclapply(seq_along(coded_seqs), function(x) {
    
    if(names(coded_seqs)[x] %in% hypermuts) {
      coded_seqs[[x]] = gsub("1", "3", coded_seqs[[x]]) 
    } else {
      coded_seqs[[x]] = coded_seqs[[x]]
    }
    return(coded_seqs[[x]])
  }, mc.cores=numCores)
  names(coded_seqs) = nm
}

print("Searching for sequences belonging to original [bad] primers...")
# Next code sequences resulting from bad primers
if(opt$full=="Y"){
  nm = names(coded_seqs)
  coded_seqs = mclapply(seq_along(coded_seqs), function(x) {
    nm_meta=names(coded_seqs)[x]
    panel = meta$panel[nm_meta==meta$id]
    if(panel=="UF.HIV") {
      coded_seqs[[x]] = gsub("1", "6", coded_seqs[[x]]) 
    } else {
      coded_seqs[[x]] = coded_seqs[[x]]
    }
    return(coded_seqs[[x]])
  }, mc.cores=numCores)
  names(coded_seqs) = nm
} # End if-else statement


if(opt$clust=="Y") {
  print("Clustering data...")
  seqs_df <- as.data.frame(do.call(rbind, coded_seqs), stringsAsFactors = T) %>%
    mutate_all(~factor(., ordered = F, levels=c("4", "3", "2", "1", "0"))) # Translate categorical data into a data frame for dissimilarity matrix calculation
  
  dist_mat <- daisy(as.data.frame(seqs_df)) # Create distance matrix of categorical data
  
  
  clust_rawdata <- hclust(as.dist(dist_mat,
                                  diag = TRUE,
                                  upper = FALSE)) # Create hierarchical cluster tree
  
  
  
  
  pdf(file="dendrogram.pdf")
  plot(clust_rawdata, labels=FALSE)
  #abline(h = opt$tol.level, col = 2, lty = 2)
  dev.off()
  
  max_clusters = 100
  max_list <- seq(2, max_clusters, 1)
  clust_rawdata_sil <- mclapply(seq_along(max_list), function(x) {
    return(data.frame(num.clust=max_list[x],
                      si.width=summary(silhouette(cutree(clust_rawdata, k=max_list[x]), as.dist(dist_mat)))$avg.width))
    
  }, mc.cores=numCores)
  clust_rawdata_sil = do.call(rbind, clust_rawdata_sil)
  
  print("Plotting sihouette for various cluster numbers to find optimal threshold")
  pdf(file="optimal_cluster_num.pdf")
  ggplot(data = clust_rawdata_sil,
         aes(x=num.clust, y=si.width)) +
    geom_point()+
    geom_line()+
    labs(x = "Number of clusters", y = "Average silhouette width")
  dev.off()
  
  max_clust=clust_rawdata_sil$num.clust[clust_rawdata_sil$si.width == max(clust_rawdata_sil$si.width)]
  
  print(paste0("The number of clusters with greatest average si.width is ", max_clust))
  
  if(isTRUE(max_clust >= opt$tol.level)) {
    opt$tol.level=max_clust}
  
  
  print("Printing silhouette statistics")
  ## Urge you to look up silhouette for clustering
  
  clust_rawdata_sil <- silhouette(cutree(clust_rawdata, k=opt$tol.level) ,
                                  as.dist(dist_mat))
  
  rownames(clust_rawdata_sil) <- names(coded_seqs)
  
  #write.csv(clust_rawdata_sil, "rawdata.csv", quote=F)
  
  
  pdf(file=paste0("silhouette_", opt$tol.level, ".pdf"))
  plot(clust_rawdata_sil)
  dev.off()
  
  
  ## If you want to test multiple clustering methods, you can fill this table to keep track of performance
  if(isTRUE(exists("summary_widths", envir=globalenv()))) {
    summary_widths=summary_widths
  } else {
    summary_widths <- data.frame(method=NA,
                                 hclust_method=NA,
                                 SI_avg_width=NA,
                                 tol.level=NA)
  }
  
  summary_widths <- rbind(summary_widths,
                          data.frame(method="clust_rawdata",
                                     hclust_method="average",
                                     SI_avg_width=summary(clust_rawdata_sil)$avg.width,
                                     tol.level=opt$tol.level))
  
  ## After filling table, be sure to choose the best method for downstream analysis
  best_method <- filter(summary_widths,
                        SI_avg_width==max(summary_widths$SI_avg_width, na.rm=T))$method
  best_method_sil <- paste0(best_method, "_sil")
  
  
  print("Manipulating data for plotting")
  ## Create dataframes with info on cluster assignment for individual sequences and
  ## Sequence representatives for each cluster
  cluster_assign <- as.data.frame(get(best_method_sil)[,c(1,3)] )%>%
    rownames_to_column(., var="seq")
  
  write.csv(cluster_assign, paste0("cluster_assignment_", opt$tol.level, ".csv"), quote=F, row.names=F)
  
  cluster_reps <- cluster_assign %>%
    group_by(cluster) %>%
    filter(sil_width == max(sil_width)) %>%
    filter(row_number()==1)
  
  cluster_size <- cluster_assign %>%
    group_by(cluster) %>%
    summarize(size=length(seq))
  
  hist(cluster_size$size)
  
  
  dend.tree <- get(best_method) %>%
    #  as.dendrogram() %>%
    as.phylo()
  
  new.dend.tree <- drop.tip(dend.tree,
                            names(coded_seqs)[names(coded_seqs) %notin% cluster_reps$seq])
  
  
  cluster_size_tbl <- merge(cluster_size, cluster_reps) %>%
    rename(label=seq) %>%
    mutate(cluster_intervals=cut(size, 5, dig.lab=2)) %>%
    mutate(cluster_size=cut(size, 5, labels=F))
  
  new.dend.tree <- as.treedata(as_tibble(new.dend.tree)) %>%
    left_join(., cluster_size_tbl)
  
  
  new.dend.tree@extraInfo <- new.dend.tree@extraInfo %>%
    mutate(cluster_size = replace_na(cluster_size, as.integer(6))) %>%
    mutate(cluster_size = factor(cluster_size, levels=c("6", "1", "2","3","4", "5")))
  
  sizes = setNames( c(0.25, cut(1:5, 5, labels=F)/2), levels(new.dend.tree@extraInfo$cluster_size))
  
  print("Generating cluster tree")
  t1 <- ggtree(new.dend.tree, aes(size=as.factor(cluster_size))) +
    scale_size_manual(values=sizes, labels = c("Backbone",
                                               levels(cluster_size_tbl$cluster_intervals)), name="Cluster size") +
                                              # "(1,81]", "(81,160]", "(160,240]", "(240,320]", "(320,400]"), name="Cluster size") +
    theme(legend.position = c(0.2,0.9),
          legend.background = element_rect(fill=alpha('white', 0)),
          text=element_text(size=18),
          legend.text = element_text(size = 12),
          plot.margin = margin(1,1,1,1.5, "cm"),
          legend.title = element_text(size=16))
  
  ## Pull out sequences for mapping and order according to tree
  is_tip <- new.dend.tree@phylo$edge[,2] <= length(new.dend.tree@phylo$tip.label)
  ordered_tips <- new.dend.tree@phylo$edge[is_tip, 2]
  ordered_tips <- new.dend.tree@phylo$tip.label[ordered_tips]
  
  rep_seqs <- do.call(bind_rows, coded_seqs[names(coded_seqs) %in% cluster_reps$seq])
  genomic_idx <- match(names(rep_seqs), ordered_tips)
  rep_seqs  <- rep_seqs[ , genomic_idx] 
  
  d3 <- rep_seqs %>% rownames_to_column('Site') %>% pivot_longer(cols = -Site) %>%
    mutate(value=factor(value, levels=c("0", "1", "2", "3", "4", "5")),
           Site=as.integer(Site))
  print("Potting tree with scatter plot representing sequence data")
  cols=c("0"="magenta", "1"="darkblue", "2"="grey", "3"="coral", "4"="khaki")
  
  p1 <-  filter(d3, Site <= 250) %>%
    ggplot(aes(x=name, y=Site, color=value)) +
    geom_point(shape=19, size=1) +
    theme_minimal() +
    coord_flip() +
    ylab("2-LTR") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="none",
          text=element_text(size=18)) +
    scale_color_manual(values=cols, limits=levels(d3$value))
  
  
  p2 <-  filter(d3, Site>250, Site<=9498) %>%
    ggplot(aes(x=name, y=Site-125, color=value)) +
    geom_point(shape=19, size=1) +
    scale_y_continuous(breaks= pretty_breaks()) +
    theme_minimal() +
    coord_flip() +
    ylab("Site") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=18)) +
    scale_color_manual(labels=c("2-LTR", "Observed", "Ambiguous", "Hypermutated", "Missing"),
                       name="",
                       values=cols,
                       limits=levels(d3$value)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  
  #multiplot(t1, p1, nrow=1)
  
  final <- plot_grid(t1, p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 1/20, 0.75))
  
  
  ggsave(plot=final, file=paste0("clustered_defective_profiles_", opt$tol.level, ".pdf"), height=11, width=8.5, units="in")
  
  #                  name=c(rep("A", 10), rep("B", 10)),
  #                  value=as.factor(c(1,1,1,0,0,0,0,2,2,1,1,0,0,0,0,0,0,2,2,1))
  # )
  # 
  # d4 %>%
  #   ggplot(aes(x=name, y=Site, color=value)) + 
  #   geom_point(shape=19, size=10) +
  #   scale_y_continuous(breaks= pretty_breaks()) +
  #   theme_minimal() +
  #   coord_flip() +
  #   theme(axis.text.y = element_blank(),
  #         axis.title.y = element_blank(),
  #         panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank())
  
  
  ##############################################################################
}


############################################################################################################################################

amb_thresh=0.98

print("Categorizing defects...")
defects <- mclapply(seq_along(coded_seqs), function(seq) {
  x=coded_seqs[[seq]]
  
#  five_prime=c(x[255:512],x[624:886])
  five_prime=x[255:532]
#  five_prime=c(x[256:532, x[625:886]])
#  three_prime=x[9031:9501] # because missing AMP36!
 # three_prime=c(x[8799:9025],x[9232:9501])
  three_prime=x[9232:9501]
  internal=x[917:8789]

  five_prime_amb=length(five_prime[five_prime=="4"])/length(five_prime[five_prime!="5"])
  three_prime_amb=length(three_prime[three_prime=="4"])/length(three_prime[three_prime!="5"])
#  internal_amb=length(internal[internal=="4"])/length(internal[internal!="5"])
  internal_col=str_c(internal, collapse="")
  xl_internal=lapply(regmatches(internal_col, gregexpr("[425]{1200,}", internal_col)), nchar) 
  internal_amb=length(xl_internal[[1]])
  
  psi=x[1035:1186]
  psi_amb=length(psi[psi=="4"])/length(psi[psi!="5"])


  # hypermut = ifelse("3" %in% x &
  #                     five_prime_amb<0.99 & 
  #                     three_prime_amb<0.99, 1, 0)
  
  # hypermut_del = ifelse("3" %in% x & 
  #                         hypermut==0, 1, 0)
  
  hypermut = ifelse("3" %in% x, 1, 0)
  
  point = ifelse((hypermut==0 &
                    #hypermut_del==0 &
                    five_prime_amb < amb_thresh & 
                    three_prime_amb < amb_thresh &
                    internal_amb == 0) &
                   (point_muts[[seq]]=="point" |
                      psi_amb>=0.99), 1, 0)
  
  intact = ifelse(hypermut==0 &
                    #hypermut_del==0 &
                    five_prime_amb < amb_thresh & 
                    three_prime_amb < amb_thresh &
                    internal_amb == 0 &
                    point==0, 1, 0)
  
  unknown = ifelse(hypermut==0 &
                     point==0 &
                     intact==0 &
                     "6" %in% x, 1, 0)
  
  internal_del = ifelse(hypermut==0 &
                          unknown==0 &
                          five_prime_amb < amb_thresh & 
                          three_prime_amb < amb_thresh &
                          internal_amb >= 1, 1, 0) 
  
  five_prime_del = ifelse(hypermut==0 &
                            #hypermut_del==0 &
                            unknown==0 &
                            internal_del==0 &
                            five_prime_amb >= amb_thresh & 
                            three_prime_amb < amb_thresh, 1, 0)
  
  three_prime_del = ifelse(hypermut==0 &
                             #hypermut_del==0 &
                             unknown==0 &
                             internal_del==0 &
                             three_prime_amb >= amb_thresh & 
                             five_prime_amb < amb_thresh, 1, 0)
  
  both_del = ifelse(hypermut==0 &
                      #hypermut_del==0 &
                      unknown==0 &
                      internal_del==0 &
                      five_prime_amb >= amb_thresh & 
                      three_prime_amb >= amb_thresh, 1, 0)

  i=251
  f=0
  while(i<9000) {
    if(isTRUE((x[i]=="4" | x[i]=="5") & x[i]!="6")) {
      f=f+1
      i=i+1
    } else{
      i=9001
    }
  }

  rev_x = rev(x)
  i=1
  r=0
  while(i<9000) {
    if(isTRUE((rev_x[i]=="4" | rev_x[i]=="5") & x[i]!="6")) {
      r=r+1
      i=i+1
    } else{
      i=9001
    }
  }

  result=data.frame(five_prime_del=five_prime_del,
                    three_prime_del=three_prime_del,
                    both_del=both_del,
                    internal_del=internal_del,
                    unknown=unknown,
                    hypermut=hypermut,
                    #hypermut_del=hypermut_del,
                    point=point,
                    intact=intact,
                    five_prime_del_length=f,
                    three_prime_del_length=r)

}, mc.cores=numCores)
defects <- rbindlist(defects, idcol="isolate")
defects$isolate=names(coded_seqs)

write.csv(defects, "sequence_profiles.csv", quote=F, row.names=F)
## Arrange deletions by deletion type and then arrange sequences by same order
# then plot using the same code above

if(opt$full=="Y") {
d_order = c("intact", "point", 
            "five_prime_del","three_prime_del", "internal_del",  "both_del",
            "hypermut", "unknown")
} else {
  d_order = c("intact", "point", 
              "five_prime_del","three_prime_del", "internal_del",  "both_del",
              "hypermut")
}

defects2 = defects %>%
  gather(del_type, value, five_prime_del:intact) %>%
  dplyr::filter(value==1) %>%
  dplyr::mutate(del_type = factor(del_type, levels=d_order))

name_levels = group_by(defects2, del_type) %>%
  dplyr::arrange(desc(three_prime_del_length), desc(five_prime_del_length), .by_group=T) %>%
  dplyr::ungroup() %>%
  dplyr::select(isolate)

  
defects2 =group_by(defects2, del_type)  %>% dplyr::group_split()

## Now just need to pick representative sequences (using distance matrix for each category?)

print("Picking representative sequences...")
rep_seqs_dist = do.call(rbind,
                        lapply(defects2, function(x) {
    init = seqs[names(seqs) %in% x$isolate]
    dist = dist.dna(as.DNAbin(init), as.matrix=2)
    rownames(dist) = x$isolate
    top = head(
      data.frame(isolate=rownames(dist), rs=rowSums(dist, na.rm=T), ms=rowMeans(dist, na.rm=T)) %>%
      dplyr::arrange(desc(rs)), n=10)
    return(top)
  }))

rep_seqs = coded_seqs[names(coded_seqs) %in% rep_seqs_dist$isolate]
#write.FASTA(as.DNAbin(rep_nc), file=paste0(x$del_type[1], "_seqs.fasta"))


name_levels = filter(name_levels, isolate %in% names(rep_seqs))
name_levels = as.vector(name_levels$isolate)


rep_seqs_df  = as.data.frame(t(do.call(rbind, rep_seqs)))

if(opt$full=="Y") {
    levs=c("0", "1", "2", "3", "4", "5", "6")
} else {
    levs=c("0", "1", "2", "3", "4", "5")
}

d4 <- rep_seqs_df %>% rownames_to_column('Site') %>% pivot_longer(cols = -Site) %>%
  mutate(value=factor(value, levels=levs),
         Site=as.integer(Site),
         name=factor(name, levels=name_levels)) 

rep_seqs_dist$isolate = factor(rep_seqs_dist$isolate, levels=name_levels)
  


print("Plotting tree with scatter plot representing sequence data")
if(opt$full=="Y"){
    cols=c("0"="darkred", "1"="darkblue", "2"="darkgrey", "3"="coral", "4"="khaki", "5"="lightgrey", "6"="black")
    labs = c("2-LTR", "Observed", "Ambiguous", "Hypermutated*", "Missing", "Not covered", "Unknown")
    } else {
        cols=c("0"="darkred", "1"="darkblue", "2"="darkgrey", "3"="coral", "4"="khaki", "5"="lightgrey")
        labs = c("2-LTR", "Observed", "Ambiguous", "Hypermutated*", "Missing", "Not covered")
    }

p3 <-  filter(d4, Site <= 250, Site>=6) %>%
  ggplot(aes(x=name, y=Site, color=value)) +
  geom_point(shape=19, size=1, show.legend=TRUE) +
  theme_minimal() +
  coord_flip() +
  ylab("2-LTR") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
 #       legend.position="none",
        text=element_text(size=14)) +
  scale_color_manual(labels=c("2-LTR", "Observed", "Ambiguous", "Hypermutated*", "Missing", "Not covered"),
                     name="",
                     values=cols,
                     limits=c("0"),
                     drop=FALSE) +
  guides(colour = guide_legend(override.aes = list(size=10)))

legend1=suppressWarnings(get_legend(p3)) #2-LTR legend

p3 <-  filter(d4, Site <= 250, Site>=6) %>%
  ggplot(aes(x=name, y=Site, color=value)) +
  geom_point(shape=19, size=1) +
  theme_minimal() +
  coord_flip() +
  ylab("2-LTR") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        text=element_text(size=14)) +
  scale_color_manual(labels=c("2-LTR", "Observed", "Ambiguous", "Hypermutated*", "Missing", "Not covered"),
                     name="",
                     values=cols,
                     limits=levels(d4$value),
                     drop=FALSE)


p4 <-  filter(d4, Site>250, Site<=9501) %>%
  ggplot(aes(x=name, y=Site-125, color=value)) +
  geom_point(shape=19, size=1, show.legend=T) +
#  scale_y_continuous(breaks= pretty_breaks()) +
  theme_minimal() +
  coord_flip() +
  ylab("Site") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=14),
        legend.direction="horizontal") +
  scale_color_manual(labels=c("Observed", "Ambiguous", "Hypermutated*", "Missing", "Not covered"),
                     name="",
                     values=cols,
                     limits=levels(d4$value),
                     drop=FALSE) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  ylim(c(8100,8500))

legend2=suppressWarnings(get_legend(p4))

p4 <-  filter(d4, Site>250, Site<=9501) %>%
  ggplot(aes(x=name, y=Site-125, color=value)) +
  geom_point(shape=19, size=1) +
  scale_y_continuous(breaks= pretty_breaks()) +
  theme_minimal() +
  coord_flip() +
  ylab("Site") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=16),
        legend.position="none",
        legend.direction="horizontal") +
  scale_color_manual(labels=c("2-LTR", "Observed", "Ambiguous", "Hypermutated*", "Missing", "Not covered"),
                     name="",
                     values=cols,
                     limits=levels(d4$value),
                     drop=FALSE) +
  guides(colour = guide_legend(override.aes = list(size=10)))

p5 = dplyr::select(rep_seqs_dist, isolate, ms) %>% 
  ggplot(aes(x=isolate, y=1, fill=ms)) + 
  ylim(0.5, 1.5) +
  geom_tile() + 
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=12),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.title=element_text(hjust=0.5)) +
  scale_fill_gradient2("Mean pairwise distance\n (subs/site)", 
                       low=muted("darkblue"), mid="lightgrey", high=muted("darkred"), 
                       midpoint=max(rep_seqs_dist$ms)/2,
                       na.value="white",
                       limits=c(0,0.08))
 #                      breaks=c(0.01, 0.025, 0.04))

legend3 = suppressWarnings(get_legend(p5)) # Legend for pariwise distance

p5 = dplyr::select(rep_seqs_dist, isolate, ms) %>% 
  ggplot(aes(x=isolate, y=1, fill=ms)) + 
  ylim(0.5, 1.5) +
  geom_tile() + 
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=12),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                legend.position = "none") +
         scale_fill_gradient2("Mean pairwise distance\n (subs/site)", 
                             low=muted("darkblue"), mid="lightgrey", high=muted("darkred"), 
                             midpoint=max(rep_seqs_dist$ms)/2,
                             limits=c(0,0.08),
                             na.value="white", breaks=c(0.01, 0.025, 0.04, 0.06))
 

fp=plot_grid(
  plot_grid(p3, p4, p5, legend3, align = "h", nrow = 1, 
            rel_widths = c(0.75, 9.5, 1,3)),
  plot_grid(legend1, legend2, NULL, align="h", nrow=1, rel_widths=c(1.5, 8.5, 4)), 
  rel_heights=c(10,2), ncol=1)

layout <- "
ABBBBBBBBC
ABBBBBBBBC
ABBBBBBBBC
ABBBBBBBBC
ABBBBBBBBC
ABBBBBBBBC
ABBBBBBBBC
########DD
E#FFFFF#DD
########DD
"
fp = p3 + p4 + p5 + legend3 + legend1 +legend2 + 
  plot_layout(design=layout)
fp
ggsave(plot=fp, "defective_profiles.png", width=11, height=6, units="in")


animals=unique(gsub("\\d+(.{4}).+", "\\1", defects$isolate))
if(length(animals)>1) {
pie_data <-  mutate(defects, animal=gsub("\\d+(.{4}).+", "\\1", isolate)) %>%
  gather(del_type, value, five_prime_del:intact)  %>%
  dplyr::group_by(animal) %>%
  dplyr::mutate(total=n())

pie_data=group_by(pie_data, animal, del_type) %>%
  dplyr::summarize(prop = sum(value)/total) %>%
  dplyr::arrange(desc(prop)) %>%
  distinct()

pie_data$del_type = factor(pie_data$del_type, 
                           levels=d_order) 
} else {
  pie_data <-  mutate(defects, animal=gsub("\\d+(.{4}).+", "\\1", isolate)) %>%
    gather(del_type, value, five_prime_del:intact)  %>%
    #  dplyr::group_by(del_type) %>%
    dplyr::group_by(animal, del_type) %>%
    dplyr::summarize(total=sum(value)) 
  
  pie_data<- pie_data %>%
    dplyr::arrange(desc(total)) %>%
    dplyr::mutate(prop = total / sum(pie_data$total) *100) %>%
    mutate(pct = prop / sum(prop))
  
  #  mutate(ypos = cumsum(prop) - (0.5*prop) )
  
  pie_data$del_type = factor(pie_data$del_type, 
                             levels=d_order)
}


#pie_data$ypos = c(60, 13, 8, 3, 98)
# Basic piechart
p6=ggplot(pie_data, aes(x="", y=prop, fill=del_type)) +
  geom_bar(stat="identity", width=1, color="black") +
#  coord_polar("y") +
  theme_void() + 
  #  theme(legend.position="none") +
  # geom_text(aes(y = ypos, label = paste(round(prop,1), "%")),  color = "white", size=6) +
  scale_fill_manual(name="",
                    values=c("darkblue", "violetred4", "seagreen", "magenta","powderblue", "purple4", "coral"),
                    labels = c("Intact",  "Psi deletion or point mutation", "5' Deletion", 
                               "3' Deletion", "Large internal deletions", 
                               "3' + 5' Deletion", "Hypermutated*")) +
  theme(
#  axis.ticks=element_blank(),
#        axis.text.y=element_blank(),
#        axis.text.x=element_blank(),
        axis.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position="right") +
#  geom_text(data=subset(pie_data, pct> 0),
#            aes(x = 1.75, label=scales::percent(pct, accuracy = .1)), 
#            position = position_stack(vjust = 0.5)) +
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) +
  facet_grid(~animal)

ggsave(plot=p6, "bar_chart.png", width=5.5, height=4, units="in")

# p6_bar=ggplot(pie_data, aes(x=animal, y=pct, fill=del_type)) +
#  geom_bar(stat="identity", width=1, color="white", linewidth=0.25) +
#  theme_minimal() + 
#  scale_fill_manual(name="",
#                    values=c("darkblue", "violetred4", "seagreen", "magenta","powderblue", "purple4", "coral"),
#                    labels = c("Intact",  "Psi deletion or point mutation", "5' Deletion", 
#                               "3' Deletion", "Large internal deletion", 
#                               "3' + 5' Deletion", "Hypermutated*")) +
#  xlab("Animal")+
#  ylab("% Genomes") +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        legend.text=element_text(size=16),
#        legend.position="right") +
#  guides(fill=guide_legend(ncol=1,byrow=TRUE))

# ggsave(plot=p6_bar, "pie_per_animal.png", width=6, height=4, units="in")

dist_data<- dplyr::select(defects, -hypermut, -intact, -point) %>%
  gather(del_type, value, five_prime_del:both_del) %>%
  dplyr::filter(value==1) %>%
  gather(dist_type, length, five_prime_del_length, three_prime_del_length)
write.csv(dist_data, "deletion_lengths.csv", row.names=F, quote=F)

p7=ggplot(dist_data) +
  geom_density(aes(x=length, fill=dist_type)) +
  facet_wrap(del_type~dist_type, ncol=2, scales="free_y",
             labeller = labeller(del_type = 
                                   c("both_del" = "3' + 5' Deletion",
                                     "five_prime_del" = "5' Deletion",
                                     "three_prime_del" = "3' Deletion"),
                                 dist_type=
                                   c("five_prime_del_length" = "5' Deletion Length Distribution",
                                     "three_prime_del_length" = "3' Deletion Length Distribution"))) +
  labs(x="Deletion length", y="Kernel Density") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        text=element_text(size=18),
        axis.text.y = element_blank()) +
  scale_fill_manual(values=c("seagreen", "magenta"))
# scale_y_continuous(labels = scientific)

tryCatch(ggsave(plot=p7, "distribution_plots.png", width=7, height=4, units="in"),
         error=function(e) NULL)


#require(ggpubr)
#pdf("deletion_distributions.pdf", height=11, width=8.5)
#ggarrange(p1,p2, ncol=1)
#dev.off()

# c1 = plot_grid(p3, p4, align="hv", axis = "l", nrow = 1, rel_widths = c(0.05, 1))
# #gsave(plot=c1, "defective_profiles_a.pdf", height=4, width=11, units="in")
# c2 = plot_grid(p5, p6, align="hv", axis = "l",  nrow = 1, rel_widths = c(0.8, 1))
# #gsave(plot=c2, "defective_profiles_bc.pdf", height=4, width=11, units="in")
# 
# final=plot_grid(c1, c2, nrow=2, align="v", axis="r", rel_heights = c(1,1.5), labels=c("a", "b", "c"))
# img <- readPNG(system.file("img", "Rlogo.png", package="png"))
# 
# ggsave(plot=final, "defective_profiles_abc.pdf", height=6.5, width=11, units="in")


#   write.csv(group_by(dist_data, del_type, dist_type) %>%
#   summarize(mean_length = mean(length),
#             med_length=median(length)), "distance_distributions.csv", quote=F, row.names=F)
# 
# # Find peaks in deletion data ############################################################
# 
# findPeaks = function (x, thresh = 0) {
#   pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0) + 2
#   if (!missing(thresh)) {
#     pks[x[pks - 1] - x[pks] > thresh]
#   }
#   else pks
# }
# 
# ## Plot peaks for 5+3 deletions only  
# p_select=ggplot(dist_data[dist_data$del_type=="both_del",]) +
#   geom_density(aes(x=length, fill=dist_type)) +
#   labs(x="Deletion length", y="Kernel Density") +
#   theme_minimal() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "none",
#         text=element_text(size=18),
#         axis.text.y = element_blank()) +
#   scale_fill_manual(values=c("steelblue4", "steelblue3"))
# 
# plot_build = ggplot_build(p_select)
# 
# #plot_build$data[[1]]$x[findPeaks(plot_build$data[[1]]$density)]
# 
# ## Now export peaks for all deletion types:
# 
# p_select2=ggplot(dist_data) +
#   geom_density(aes(x=length, fill=dist_type)) +
#   labs(x="Deletion length", y="Kernel Density") +
#   theme_minimal() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "none",
#         text=element_text(size=18),
#         axis.text.y = element_blank()) +
#   scale_fill_manual(values=c("steelblue4", "steelblue3")) +
#   facet_wrap(~del_type)
# 
# plot_build2 = ggplot_build(p_select2)
# 
# df = plot_build2$data[[1]] %>%
#   dplyr::select(PANEL, fill, x, density) %>%
#   dplyr::rename(del_length=x) %>%
#   dplyr::mutate(fill = ifelse(fill=="steelblue4", "five_prime_del_length", "three_prime_del_length"),
#          PANEL = ifelse(PANEL==1, "both_del", ifelse(PANEL==2, "5_prime_del", "3_prime_del"))) %>%
#   dplyr::group_by(PANEL, fill) %>%
#   dplyr::group_split(.keep=T)
# 
# peaks = do.call(rbind, lapply(df,   function(x) {
#   x = arrange(x, PANEL, fill, del_length)
#   peak_length = findPeaks(x$density)
#   if (length(peak_length)==0) {
#     peak_length=0
#   }
#   return(data.frame(group=x$PANEL[1], del_type=x$fill[1], peak_length = peak_length))
# })) %>%
#   arrange(group, del_type)
# 
# write.csv(peaks, "peak_deletion_lengths.csv", row.names=F, quote=F)

