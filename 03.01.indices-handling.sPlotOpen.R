#### recommended machine: 01####

#### packages ----
library("tidyverse")
library("tidyr")
library("data.table")
library("raster")
library("ncdf4")
library("sf")

#### Phylogenetic Data----

# load phylo data
PD.names <- list.files("../03.data/eve-calculations/PD-calculation-sPlotOpen")

PD.list <- list()

for (i in 1:length(PD.names)) {
  PD.list[[i]] <- readRDS(paste0("../03.data/eve-calculations/PD-calculation-sPlotOpen/", PD.names[i]))
}

df.PD <- bind_rows(PD.list)

saveRDS(df.PD, file = "../03.data/01a.PD-values-sPlotOpen.Rds")

#### Functional Data----

# load data
FD.names <- list.files("../03.data/eve-calculations/FD-calculation-sPlotOpen")

FD.list <- list()

for (i in 1:length(FD.names)) {
  FD.list[[i]] <- readRDS(paste0("../03.data/eve-calculations/FD-calculation-sPlotOpen/", FD.names[i]))
}

df.FD <- bind_rows(FD.list)

saveRDS(df.FD, file = "../03.data/01a.FD-values-sPlotOpen.Rds")
