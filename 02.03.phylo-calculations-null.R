####recommended machine: 03#####

####packages ----
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("picante"))
suppressPackageStartupMessages(library("phytools"))

sessionInfo()

args <- commandArgs(trailingOnly = T)
input_dir <- args[1]
output <- args[2]
SR.start <- as.numeric(args[3])
SR.end <- as.numeric(args[4])

####load data ----
suppressWarnings(load(file.path(input_dir, "DT_sPlot3.0.RData")))
suppressWarnings(load(file.path(input_dir, "01.tree", "PhyloTree_sPlot3.0.RData")))
suppressWarnings(load(file.path(input_dir, "pruned.trees.traitsp.RData")))
#select sPlot species in PhyloTree

sPlot.species <- DT2 %>% 
  distinct(Species) 

PD3 <- function(mylist) {
  spp1 <- mylist %>%
    pivot_wider(names_from = Species, 
                values_from = Abundance) %>% 
    as.data.frame()
  
  pruned.tree <- prune.sample(spp1, sPlot.phylo.tree$scenario.3)
  
  ####PD calculation----
  
  PDnull <- data.frame(RQE.phyl.null = NA,
                       MPD.null = NA)
  
  if(length(colnames(spp1))>1) {
    #PDobs <- pd(spp1, mylist$pruned.tree[[1]], include.root = FALSE)
    PDnull$MPD.null <- mpd(spp1, cophenetic.phylo(pruned.tree))
    #PDobs$MNTD <- mntd(spp1, cophenetic.phylo(mylist$pruned.tree[[1]]))
    if(is.ultrametric(pruned.tree)==TRUE){
      invisible(capture.output(PDnull$RQE.phyl.null <- raoD(spp1, pruned.tree)$Dkk))
    } else {
      pruned.tree <- prune.sample(spp1, sPlot.phylo.tree$scenario.3) %>% 
                      force.ultrametric(., method = "nnls")
      invisible(capture.output(PDnull$RQE.phyl.null <- raoD(spp1, pruned.tree)$Dkk))
    }
  }
  return(PDnull)
}

RQE.out <- list()
MPD.out <- list()

spec <- c(SR.start:SR.end)

for(i in 1:length(spec)) {
PD.out <- replicate(499,
  sPlot.in.Tree %>% 
  sample_n(spec[i]) %>% 
  as.data.frame() %>%  
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Abundance = 1), simplify = FALSE) %>% 
  map_dfr(., PD3)

RQE.out[[i]] <- PD.out$RQE.phyl.null %>% as.data.frame()
colnames(RQE.out[[i]]) <- paste0(spec[i])
MPD.out[[i]] <- PD.out$MPD.null %>% as.data.frame() 
colnames(MPD.out[[i]]) <- paste0(spec[i])
}

RQE.phyl.null <- do.call(cbind, RQE.out)
MPD.null <- do.call(cbind, MPD.out)

list.out <- list(MPD.null, RQE.phyl.null)

saveRDS(list.out, file = output)
