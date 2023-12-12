####recommended machine: 03#####

####packages ----
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("picante"))

sessionInfo()

args <- commandArgs(trailingOnly = T)
input_dir <- args[1]
output <- args[2]
start <- as.numeric(args[3])
end <- as.numeric(args[4])

####load data ----
suppressWarnings(load(file.path(input_dir, "sPlotOpen.RData")))
#suppressWarnings(load(file.path(input_dir, "header_sPlot3.0.RData")))
suppressWarnings(load(file.path(input_dir, "01.tree", "PhyloTree_sPlot3.0.RData")))
suppressWarnings(load(file.path(input_dir, "pruned.trees.traitsp.sPlotOpen.RData")))
#select sPlot species in PhyloTree

sPlot.species <- DT2.oa %>% 
  distinct(Species) 

PD3 <- function(mylist) {
  spp1 <- mylist$mydata[[1]] %>%
    filter(Species %in% mylist$pruned.tree[[1]]$tip.label) %>% 
    group_by(Species) %>% 
    dplyr::summarise(Abundance = (sum(Abundance)>0) * 1) %>% 
    pivot_wider(names_from = Species, 
                values_from = Abundance) %>% 
    as.data.frame()
  
  #length(mylist$pruned.tree[[1]])
  
  #  ####PD calculation----
  
  if(length(colnames(spp1))>1) {
                    PDobs <- pd(spp1, mylist$pruned.tree[[1]], include.root = FALSE)
                    PDobs$MPD <- mpd(spp1, cophenetic.phylo(mylist$pruned.tree[[1]]))
                    PDobs$MNTD <- mntd(spp1, cophenetic.phylo(mylist$pruned.tree[[1]]))
          if(is.ultrametric(mylist$pruned.tree[[1]])==TRUE){
                    invisible(capture.output(PDobs$RaoD.phyl <- raoD(spp1, mylist$pruned.tree[[1]])$Dkk))
          } else {
                    PDobs$RaoD.phyl <- NA
                  }
    } else if(length(colnames(spp1))==1){
                    PDobs <- pd(spp1, mylist$pruned.tree[[1]], include.root = TRUE)
                    PDobs$RaoD.phyl <- NA
                    PDobs$MPD <- NA
                    PDobs$MNTD <- NA
    } else {
                    PDobs <- data.frame(PD = NA, SR =NA)
                    PDobs$RaoD.phyl <- NA
                    PDobs$MPD <- NA
                    PDobs$MNTD <- NA
    }
  return(PDobs)
}

PlotIDs <- DT2.oa %>% filter(Species %in% sPlot.in.Tree$Species) %>% 
                   distinct(PlotObservationID) %>% pull(PlotObservationID)

if(end < length(PlotIDs)) {
  PlotIDs <- PlotIDs[start:end]
} else if(end >= length(PlotIDs)) {PlotIDs <- PlotIDs[start:length(PlotIDs)]}


PD.out <- DT2.oa %>% 
          filter(PlotObservationID %in% PlotIDs) %>% 
          dplyr::select(PlotObservationID, Species, Original_abundance) %>% 
          filter(Species %in% sPlot.in.Tree$Species) %>% 
          mutate(Species = gsub(" ", "_", Species)) %>%
          group_by(PlotObservationID) %>% 
          dplyr::summarise(mydata=list(data.frame(Species=Species,
                                                  Abundance=Original_abundance))) %>% 
          left_join(header.oa %>%
                       dplyr::select(PlotObservationID, GIVD_ID),
                     by="PlotObservationID") %>%
           left_join(PD.out0 %>%
                       dplyr::select(GIVD_ID, pruned.tree),
                     by="GIVD_ID") %>%
          split(.$PlotObservationID) %>% 
          map_dfr(., PD3) %>% 
          mutate(PlotObservationID = PlotIDs)

saveRDS(PD.out, file = output)
