#### recommended machine: 03#####

#### packages ----
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("picante"))
suppressPackageStartupMessages(library("phytools"))

sessionInfo()

args <- commandArgs(trailingOnly = T)
input_dir <- args[1]
output <- args[2]
SR.start <- as.numeric(args[3])
SR.end <- as.numeric(args[4])

#### load data ----
# suppressWarnings(load(file.path(input_dir, "DT_sPlot3.0.RData")))
suppressWarnings(load(file.path(input_dir, "01.tree", "PhyloTree_sPlot3.0.RData")))
suppressWarnings(load(file.path(input_dir, "pruned.trees.traitsp.RData")))
Species.Biome <- readRDS(file.path(input_dir, "species-list-biome.RDS"))

# get biomes

biom <- Species.Biome %>%
  pull(Biome) %>%
  unique()

# function to calculate PD
PD3 <- function(mylist) {
  spp1 <- mylist %>%
    pivot_wider(
      names_from = Species,
      values_from = Abundance
    ) %>%
    as.data.frame()

  pruned.tree <- prune.sample(spp1, sPlot.phylo.tree$scenario.3)

  #### PD calculation----

  PDnull <- data.frame(
    RQE.phyl.null = NA,
    MPD.null = NA
  )

  if (length(colnames(spp1)) > 1) {
    # PDobs <- pd(spp1, mylist$pruned.tree[[1]], include.root = FALSE)
    PDnull$MPD.null <- mpd(spp1, cophenetic.phylo(pruned.tree))
    # PDobs$MNTD <- mntd(spp1, cophenetic.phylo(mylist$pruned.tree[[1]]))
    if (is.ultrametric(pruned.tree) == TRUE) {
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

for (i in 1:length(spec)) {
  RQE.bio <- list()
  MPD.bio <- list()

  for (j in 1:length(biom)) {
    PD.out <- replicate(499,
      Species.Biome %>%
        filter(
          Biome == biom[j],
          Species %in% sPlot.in.Tree$Species
        ) %>%
        sample_n(spec[i], weight = Species_count) %>%
        as.data.frame() %>%
        mutate(Species = gsub(" ", "_", Species)) %>%
        mutate(Abundance = 1) %>%
        dplyr::select(Species, Abundance),
      simplify = FALSE
    ) %>%
      map_dfr(., PD3)

    RQE.bio[[j]] <- PD.out$RQE.phyl.null %>%
      as.data.frame() %>%
      `colnames<-`("RQE.phyl.null") %>%
      mutate(
        Biome = biom[j],
        SR = spec[i]
      )

    MPD.bio[[j]] <- PD.out$MPD.null %>%
      as.data.frame() %>%
      `colnames<-`("MPD.null") %>%
      mutate(
        Biome = biom[j],
        SR = spec[i]
      )
  }


  RQE.out[[i]] <- do.call(rbind, RQE.bio)

  MPD.out[[i]] <- do.call(rbind, MPD.bio)
}

RQE.phyl.null <- do.call(rbind, RQE.out)
MPD.null <- do.call(rbind, MPD.out)

list.out <- list(MPD.null, RQE.phyl.null)

saveRDS(list.out, file = output)
