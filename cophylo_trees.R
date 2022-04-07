require(dplyr)
require(ape)
require(treeio)
require(ggtree)
require(phytools)
require(parallel)

numCores=detectCores()

setwd("../lmaps_treefiles")

#gp120 coordinates for custom ref: 6832 - 8172

animals=c("HG02", "JA41", "JP70") # Read from list in future

cophyloTree <- function(animal) {
  
  SIV_tree <- list.files(path=animal, pattern="all.fa.treefile")
  
  SIV_tree <- multi2di(read.tree(paste0(animal, "/", SIV_tree)))
  SIV_tree <- root(SIV_tree, outgroup="KU892415.1", resolve.root=T)
  SIV_tree$node.label <- as.numeric(if_else(grepl("\\d+\\/\\d+", SIV_tree$node.label), gsub("(\\d+)\\/\\d+", "\\1", SIV_tree$node.label), ""))
  assign(paste0(animal, "_SIV_tree"), SIV_tree, envir = globalenv())
  
  
  env_tree <- list.files(path=animal, pattern="gp120.fa.treefile")
  
  env_tree <- multi2di(read.tree(paste0(animal, "/", env_tree)))
  env_tree <- root(env_tree, outgroup="KU892415.1", resolve.root=T)
  env_tree$node.label <- as.numeric(if_else(grepl("\\d+\\/\\d+", env_tree$node.label), gsub("(\\d+)\\/\\d+", "\\1", env_tree$node.label), ""))
  assign(paste0(animal, "_env_tree"), env_tree, envir = globalenv())
  
  
  association <- cbind(SIV_tree$tip.label, SIV_tree$tip.label)[-1,]
  
  
  #co_trees <- cophylo(SIV_tree, TRIM5a_tree, assoc = association) 
  co_trees <- cophylo(SIV_tree, env_tree, assoc = association) 
  co_trees$trees[[1]]$node.label <- as.numeric(co_trees$trees[[1]]$node.label)
  co_trees$trees[[2]]$node.label <- as.numeric(co_trees$trees[[2]]$node.label)
  
  p1<-character(length(co_trees$trees[[1]]$node.label))
  #p2<-character(length(TRIM5a_tree$node.label))
  
  p2<-character(length(co_trees$trees[[2]]$node.label))
  
  #set up ranges that correspond to node support levels
  p1[co_trees$trees[[1]]$node.label >= 90] <- "black"
  p1[co_trees$trees[[1]]$node.label < 90] <- "grey"
  p1[is.na(co_trees$trees[[1]]$node.label)] <- "white"

  
  p2[co_trees$trees[[2]]$node.label >= 90] <- "black"
  p2[co_trees$trees[[2]]$node.label < 90] <- "grey"
  p2[is.na(co_trees$trees[[2]]$node.label)] <- "white"
  
  
  
  #Color according to genome coveragetiplabels.cophylo(pie=to.matrix(x[obj$trees[[1]]$tip.label],c("a","b")),
  SIV_meta <- data.frame(ID=co_trees$trees[[1]]$tip.label) %>%
    mutate(tissue=if_else(grepl("K_", ID), "Lymph Node", 
                          ifelse(grepl("KU892415.1", ID), "Outgroup", "PBMC")),
           time=if_else(grepl("119", ID), "cART", 
                        ifelse(grepl("KU892415.1", ID), "Outgroup", "Post-cART")),
           asr=paste(time, tissue, sep=" "))
  
  assign(paste0(animal, "_SIV_meta"), SIV_meta, envir = globalenv())
  
  
  
  env_meta <- data.frame(ID=co_trees$trees[[2]]$tip.label) %>%
    mutate(tissue=if_else(grepl("K_", ID), "Lymph Node", 
                          ifelse(grepl("KU892415.1", ID), "Outgroup", "PBMC")),
           time=if_else(grepl("119", ID), "cART", 
                        ifelse(grepl("KU892415.1", ID), "Outgroup", "Post-cART")),
           asr=paste(time, tissue, sep=" "))
  
  assign(paste0(animal, "_env_meta"), env_meta, envir = globalenv())
  
  
  
  
  p3<-setNames(character(length(co_trees$trees[[1]]$tip.label)), co_trees$trees[[1]]$tip.label)
  
  p3[SIV_meta$time == "cART"] <- "aquamarine"
  p3[SIV_meta$time != "cART"] <- "red"
  p3[SIV_meta$time == "Outgroup"] <- "black"
  
  p4<-setNames(character(length(co_trees$trees[[2]]$tip.label)), co_trees$trees[[2]]$tip.label)
  
  p4[env_meta$time == "cART"] <- "aquamarine"
  p4[env_meta$time != "cART"] <- "red"
  p3[env_meta$time == "Outgroup"] <- "black"
  
  
  ## Shapes
  
  p5<-setNames(numeric(length(co_trees$trees[[1]]$tip.label)), co_trees$trees[[1]]$tip.label)
  
  p5[SIV_meta$tissue == "Lymph Node"] <- as.numeric(21)
  p5[SIV_meta$tissue != "Lymph Node"] <- as.numeric(22)
  p3[SIV_meta$tissue == "Outgroup"] <- as.numeric(21)
  
  
  p6<-setNames(numeric(length(co_trees$trees[[2]]$tip.label)), co_trees$trees[[2]]$tip.label)
  
  p6[env_meta$tissue == "Lymph Node"] <- as.numeric(21)
  p6[env_meta$tissue != "Lymph Node"] <- as.numeric(22)
  p3[env_meta$tissue == "Outgroup"] <- as.numeric(21)
  
  #fn=function(x) abs(x)
  pdf(file=paste0(animal, "_cophylo_env.pdf"), 8.5, 11)
  plot(co_trees,link.type="curved",link.lwd=4,
       link.lty="solid",link.col=make.transparent("blue",0.25))
  tiplabels.cophylo(bg=p3, pch=p5, cex=1.5, which="left")
  tiplabels.cophylo(bg=p4,pch=p6, cex=1.5, which="right")
  nodelabels.cophylo(pch=21,bg=p1,cex=1, which="left")
  nodelabels.cophylo(pch=21,bg=p2,cex=1, which="right")
  add.scale.bar(x=0.03, y=0, cex = 0.7, font = 2, col = "black")
  dev.off() 
} # End function

lapply(animals, cophyloTree)