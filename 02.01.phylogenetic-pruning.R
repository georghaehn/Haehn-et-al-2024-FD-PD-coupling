#### recommended machine: 02 ----

#### packages ----
library("tidyverse")
library("tidyr")
library("picante")
library("future")
library("furrr")

#### load sPlot data ----
load("~/share/groups/sPlot/releases//sPlot3.0/DT_sPlot3.0.RData")
load("~/share/groups/sPlot/releases//sPlot3.0/header_sPlot3.0.RData")
load("02.data/PhyloTree_sPlot3.0.Rdata")
load("~/share/groups/sPlot/releases//sPlot3.0/Traits_CWMs_sPlot3.RData")

# select sPlot species present in the Phylogeny and with known traits

sPlot.species <- DT2 %>%
  distinct(Species)

species.trait <- sPlot.traits %>%
  select(Species, SLA_mean, PlantHeight_mean, SpecificRootLength_mean) %>%
  filter(complete.cases(.)) %>%
  pull(Species)

sPlot.in.Tree <- sPlot.species %>%
  filter(Species %in% (tree.species %>%
    filter(status != "fail to bind") %>%
    pull(species))) %>%
  filter(Species %in% species.trait)

pruning <- function(myspecies, tree.to.prune) {
  spp1 <- matrix(
    data = 1,
    nrow = 1,
    ncol = length(unique(myspecies)),
    dimnames = list(
      "mydataset",
      unique(myspecies)
    )
  )

  pruned.tree <- prune.sample(spp1, tree.to.prune)
  return(pruned.tree)
}

### Build the pruned tree once for each dataset ----

PD.out0 <- DT2 %>%
  left_join(
    header %>%
      dplyr::select(PlotObservationID, `GIVD ID`),
    by = "PlotObservationID"
  ) %>%
  filter(Species %in% sPlot.in.Tree$Species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  group_by(`GIVD ID`) %>%
  dplyr::summarize(pruned.tree = list(pruning(
    myspecies = Species,
    tree.to.prune = sPlot.phylo.tree$scenario.3
  )))

save(PD.out0, sPlot.in.Tree, file = "02.data/pruned.trees.traitsp.RData")

rm(list = ls())

#### redo for sPlotOpen----
#### load data ----
load("02.data/sPlotOpen.RData")
load("02.data/PhyloTree_sPlot3.0.Rdata")
load("02.data/Traits_CWMs_sPlot3.RData")

# select species in phylogeny and with known traits

sPlot.species <- DT2.oa %>%
  distinct(Species)

species.trait <- sPlot.traits %>%
  select(Species, SLA_mean, PlantHeight_mean, SpecificRootLength_mean) %>%
  filter(complete.cases(.)) %>%
  pull(Species)

sPlot.in.Tree <- sPlot.species %>%
  filter(Species %in% (tree.species %>%
    filter(status != "fail to bind") %>%
    pull(species))) %>%
  filter(Species %in% species.trait)

pruning <- function(myspecies, tree.to.prune) {
  spp1 <- matrix(
    data = 1,
    nrow = 1,
    ncol = length(unique(myspecies)),
    dimnames = list(
      "mydataset",
      unique(myspecies)
    )
  )

  pruned.tree <- prune.sample(spp1, tree.to.prune)
  return(pruned.tree)
}

### Build the pruned tree once for each dataset ----

PD.out0 <- DT2.oa %>%
  left_join(
    header.oa %>%
      dplyr::select(PlotObservationID, GIVD_ID),
    by = "PlotObservationID"
  ) %>%
  filter(Species %in% sPlot.in.Tree$Species) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  group_by(GIVD_ID) %>%
  dplyr::summarize(pruned.tree = list(pruning(
    myspecies = Species,
    tree.to.prune = sPlot.phylo.tree$scenario.3
  )))

save(PD.out0, sPlot.in.Tree, file = "02.data/pruned.trees.traitsp.sPlotOpen.RData")
