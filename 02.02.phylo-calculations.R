#### recommended machine: 03#####

#### packages ----
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("picante"))

sessionInfo()

args <- commandArgs(trailingOnly = T)
input_dir <- args[1]
output <- args[2]
start <- as.numeric(args[3])
end <- as.numeric(args[4])

#### load data ----
suppressWarnings(load(file.path(input_dir, "DT_sPlot3.0.RData")))
suppressWarnings(load(file.path(input_dir, "header_sPlot3.0.RData")))
suppressWarnings(load(file.path(input_dir, "01.tree", "PhyloTree_sPlot3.0.RData")))
suppressWarnings(load(file.path(input_dir, "pruned.trees.traitsp.RData")))
# select sPlot species in PhyloTree

sPlot.species <- DT2 %>%
  distinct(Species)

# sPlot.in.Tree <- sPlot.species %>%
#   filter(Species %in%
#            gsub("_", " ", unique(sPlot.phylo.tree$scenario.3$tip.label)))
#
#            # (tree.species %>%
#            #               filter(status != "fail to bind") %>%
#            #               pull(species)))

PD3 <- function(mylist) {
  spp1 <- mylist$mydata[[1]] %>%
    filter(Species %in% mylist$pruned.tree[[1]]$tip.label) %>%
    group_by(Species) %>%
    dplyr::summarise(Abundance = (sum(Abundance) > 0) * 1) %>%
    pivot_wider(
      names_from = Species,
      values_from = Abundance
    ) %>%
    as.data.frame()

  # length(mylist$pruned.tree[[1]])

  #  ####PD calculation----

  if (length(colnames(spp1)) > 1) {
    PDobs <- pd(spp1, mylist$pruned.tree[[1]], include.root = FALSE)
    PDobs$MPD <- mpd(spp1, cophenetic.phylo(mylist$pruned.tree[[1]]))
    PDobs$MNTD <- mntd(spp1, cophenetic.phylo(mylist$pruned.tree[[1]]))
    if (is.ultrametric(mylist$pruned.tree[[1]]) == TRUE) {
      invisible(capture.output(PDobs$RaoD.phyl <- raoD(spp1, mylist$pruned.tree[[1]])$Dkk))
    } else {
      PDobs$RaoD.phyl <- NA
    }
  } else if (length(colnames(spp1)) == 1) {
    PDobs <- pd(spp1, mylist$pruned.tree[[1]], include.root = TRUE)
    PDobs$RaoD.phyl <- NA
    PDobs$MPD <- NA
    PDobs$MNTD <- NA
  } else {
    PDobs <- data.frame(PD = NA, SR = NA)
    PDobs$RaoD.phyl <- NA
    PDobs$MPD <- NA
    PDobs$MNTD <- NA
  }
  return(PDobs)
}

PlotIDs <- DT2 %>%
  filter(Species %in% sPlot.in.Tree$Species) %>%
  filter(between(PlotObservationID, start, end)) %>%
  distinct(PlotObservationID) %>%
  pull(PlotObservationID)

PD.out <- DT2 %>%
  filter(between(PlotObservationID, start, end)) %>%
  dplyr::select(PlotObservationID, Species, Abundance) %>%
  filter(Species %in% sPlot.in.Tree$Species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  group_by(PlotObservationID) %>%
  dplyr::summarise(mydata = list(data.frame(
    Species = Species,
    Abundance = Abundance
  ))) %>%
  left_join(
    header %>%
      dplyr::select(PlotObservationID, `GIVD ID`),
    by = "PlotObservationID"
  ) %>%
  left_join(
    PD.out0 %>%
      dplyr::select(`GIVD ID`, pruned.tree),
    by = "GIVD ID"
  ) %>%
  split(.$PlotObservationID) %>%
  map_dfr(., PD3) %>%
  mutate(PlotObservationID = PlotIDs)

saveRDS(PD.out, file = output)
