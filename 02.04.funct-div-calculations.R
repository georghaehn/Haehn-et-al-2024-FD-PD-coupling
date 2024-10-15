#### recommended machine: 03#####

#### packages ----
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("FD"))
suppressPackageStartupMessages(library("funrar"))

sessionInfo()

args <- commandArgs(trailingOnly = T)
input_dir <- args[1]
output <- args[2]
start <- as.numeric(args[3])
end <- as.numeric(args[4])

#### load data ----
suppressWarnings(load(file.path(input_dir, "DT_sPlot3.0.RData")))
suppressWarnings(load(file.path(input_dir, "Traits_CWMs_sPlot3.RData")))
suppressWarnings(load(file.path(input_dir, "pruned.trees.traitsp.RData")))

rm(CWM)

traits <- sPlot.traits %>%
  select(c(
    "Species",
    "SLA_mean",
    "PlantHeight_mean",
    "SpecificRootLength_mean"
  ))

raoQ <- function(abund, trait, Hill = TRUE, scale = FALSE, method = "default") {
  abund <- as.matrix(abund)

  anames <- colnames(abund)[order(colnames(abund))]

  abund <- abund[, anames, drop = FALSE]

  trait <- as.matrix(trait)

  trait <- trait[anames, anames]

  if (ncol(abund) != ncol(trait)) {
    stop("Not all species in the abundance matrix appear in the trait matrix!")
  }

  abund <- abund / rowSums(abund)

  if (method == "default") {
    Q <- apply(abund, 1, function(x) crossprod(x, trait %*% x))
  }

  if (method == "divc") {
    Q <- apply(abund, 1, function(x) x %*% trait^2 %*% (x / 2 / sum(x)^2))
  }

  if (Hill == TRUE) Q <- 1 / (1 - Q)

  if (scale == TRUE) Q <- Q / max(Q)

  names(Q) <- rownames(abund)

  return(Q)
}

FD3 <- function(mylist) {
  spp1 <- mylist$mydata[[1]] %>%
    select(Species, Abundance) %>%
    arrange(Species) %>%
    group_by(Species) %>%
    # arrange(desc("Species")) %>%
    dplyr::summarise(Abundance = (sum(Abundance, na.rm = T) > 0) * 1) %>%
    filter(Abundance > 0) %>%
    pivot_wider(
      names_from = Species,
      values_from = Abundance
    ) %>%
    as.data.frame()

  FD.out <- data.frame(FDis.SLA = NA)

  trait.SLA <- mylist$mydata[[1]] %>%
    select(Species, SLA_mean) %>%
    distinct(Species, .keep_all = T) %>%
    arrange(Species) %>%
    filter(complete.cases(.))

  if (nrow(trait.SLA) > 1 &&
    ncol(spp1) > 1) {
    spp1.SLA <- spp1[, colnames(spp1) %in%
      trait.SLA$Species]

    trait.SLA <- trait.SLA %>%
      filter(Species %in% colnames(spp1.SLA)) %>%
      column_to_rownames("Species") %>%
      compute_dist_matrix()

    if (identical(rownames(trait.SLA), colnames(spp1.SLA)) &&
      sum(trait.SLA) > 0) {
      SLA <- dbFD(trait.SLA, spp1.SLA,
        w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F
      )

      FD.out$FDis.SLA <- SLA$FDis
      FD.out$RQE.SLA <- raoQ(spp1.SLA, trait.SLA)
    }
  }

  trait.HEIGHT <- mylist$mydata[[1]] %>%
    select(Species, PlantHeight_mean) %>%
    distinct(Species, .keep_all = T) %>%
    arrange(Species) %>%
    filter(complete.cases(.))


  if (nrow(trait.HEIGHT) > 1 &&
    ncol(spp1) > 1) {
    spp1.HEIGHT <- spp1[, colnames(spp1) %in%
      trait.HEIGHT$Species]

    trait.HEIGHT <- trait.HEIGHT %>%
      filter(Species %in% colnames(spp1.HEIGHT)) %>%
      column_to_rownames("Species") %>%
      compute_dist_matrix()

    if (identical(rownames(trait.HEIGHT), colnames(spp1.HEIGHT)) &&
      sum(trait.HEIGHT) > 0) {
      HEIGHT <- dbFD(trait.HEIGHT, spp1.HEIGHT,
        w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F
      )

      FD.out$FDis.HEIGHT <- HEIGHT$FDis
      FD.out$RQE.HEIGHT <- raoQ(spp1.HEIGHT, trait.HEIGHT)
    }
  }

  trait.ROOT <- mylist$mydata[[1]] %>%
    select(Species, SpecificRootLength_mean) %>%
    distinct(Species, .keep_all = T) %>%
    arrange(Species) %>%
    filter(complete.cases(.))

  if (nrow(trait.ROOT) > 1 &&
    ncol(spp1) > 1) {
    spp1.ROOT <- spp1[, colnames(spp1) %in%
      trait.ROOT$Species]

    trait.ROOT <- trait.ROOT %>%
      filter(Species %in% colnames(spp1.ROOT)) %>%
      column_to_rownames("Species") %>%
      compute_dist_matrix()

    if (identical(rownames(trait.ROOT), colnames(spp1.ROOT)) &&
      sum(trait.ROOT) > 0) {
      ROOT <- dbFD(trait.ROOT, spp1.ROOT,
        w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F
      )

      FD.out$FDis.ROOT <- ROOT$FDis
      FD.out$RQE.ROOT <- raoQ(spp1.ROOT, trait.ROOT)
    }
  }

  trait.MULTI <- mylist$mydata[[1]] %>%
    select(-Abundance) %>%
    distinct(Species, .keep_all = T) %>%
    arrange(Species) %>%
    filter(complete.cases(.))


  if (nrow(trait.MULTI) > 1 &&
    ncol(spp1) > 1) {
    spp1.MULTI <- spp1[, colnames(spp1) %in%
      trait.MULTI$Species]

    trait.MULTI <- trait.MULTI %>%
      filter(Species %in% colnames(spp1.MULTI)) %>%
      column_to_rownames("Species") %>%
      compute_dist_matrix()

    if (identical(rownames(trait.MULTI), colnames(spp1.MULTI)) &&
      sum(trait.MULTI) > 0) {
      MULTI <- dbFD(trait.MULTI, spp1.MULTI,
        w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F
      )

      FD.out$FDis.MULTI <- MULTI$FDis
      FD.out$RQE.MULTI <- raoQ(spp1.MULTI, trait.MULTI)
    }
  }

  return(FD.out)
}

PlotIDs <- DT2 %>%
  filter(between(PlotObservationID, start, end)) %>%
  filter(Species %in% sPlot.in.Tree$Species) %>%
  group_by(PlotObservationID) %>%
  filter(n() > 1) %>%
  distinct(PlotObservationID) %>%
  pull(PlotObservationID)

FD.out <- DT2 %>%
  filter(between(PlotObservationID, start, end)) %>%
  filter(Species %in% sPlot.in.Tree$Species) %>%
  dplyr::select(PlotObservationID, Species, Abundance) %>%
  left_join(traits, by = "Species") %>%
  group_by(PlotObservationID) %>%
  filter(n() > 1) %>%
  dplyr::summarise(mydata = list(data.frame(
    Species = Species,
    Abundance = Abundance,
    SLA_mean = SLA_mean,
    PlantHeight_mean = PlantHeight_mean,
    SpecificRootLength_mean = SpecificRootLength_mean
  ))) %>%
  split(.$PlotObservationID) %>%
  map_dfr(., FD3) %>%
  mutate(PlotObservationID = PlotIDs)


saveRDS(FD.out, file = output)
