#### recommended machine: 01####


#### packages ----
library("tidyverse")
library("tidyr")
library("data.table")
library("raster")
library("ncdf4")
library("sf")
library("rgdal")
library("RStoolbox")
library("downloader")

#### Phylogenetic Data----

# load phylo data
PD.names <- list.files("02.data/eve-calculations/PD-calculation")

PD.list <- list()

for (i in 1:length(PD.names)) {
  PD.list[[i]] <- readRDS(paste0("02.data/eve-calculations/PD-calculation/", PD.names[i]))
}

df.PD <- bind_rows(PD.list)

# load phylo sqrt data
PD.sqrt.names <- list.files("02.data/eve-calculations/phylo-calculations-sqrt")

PD.sqrt.list <- list()

for (i in 1:length(PD.sqrt.names)) {
  PD.sqrt.list[[i]] <- readRDS(paste0("02.data/eve-calculations/phylo-calculations-sqrt/", PD.sqrt.names[i]))
}

df.PD.sqrt <- bind_rows(PD.sqrt.list)

# load phylo null data

PD.null.names <- list.files("02.data/eve-calculations/PD-null-calculation", pattern = "*.Rds")

PD.null.list <- list()

for (i in 1:length(PD.null.names)) {
  PD.null.list[[i]] <- readRDS(paste0("02.data/eve-calculations/PD-null-calculation/", PD.null.names[i]))

  names(PD.null.list[[i]]) <- c("MPD-null", "RQE-phyl-null", "RQE-phyl-sqrt-null", "MPD-sqrt-null")
}

list.MPD.null <- lapply(PD.null.list, `[[`, "MPD-null")
list.RQE.phyl.null <- lapply(PD.null.list, `[[`, "RQE-phyl-null")
list.MPD.sqrt.null <- lapply(PD.null.list, `[[`, "MPD-sqrt-null")
list.RQE.phyl.sqrt.null <- lapply(PD.null.list, `[[`, "RQE-phyl-sqrt-null")

df.MPD.null <- bind_cols(list.MPD.null)
df.RQE.phyl.null <- bind_cols(list.RQE.phyl.null)
df.MPD.sqrt.null <- bind_cols(list.MPD.sqrt.null)
df.RQE.phyl.sqrt.null <- bind_cols(list.RQE.phyl.sqrt.null)

PD.null <- apply(df.MPD.null, 2, mean) %>%
  as.data.frame() %>%
  rename(., mean.MPD = `.`) %>%
  mutate(median.MPD = apply(df.MPD.null, 2, median) %>% unname()) %>%
  mutate(sd.MPD = apply(df.MPD.null, 2, sd) %>% unname()) %>%
  mutate(mean.MPD.sqrt = apply(df.MPD.sqrt.null, 2, mean) %>% unname()) %>%
  mutate(median.MPD.sqrt = apply(df.MPD.sqrt.null, 2, median) %>% unname()) %>%
  mutate(sd.MPD.sqrt = apply(df.MPD.sqrt.null, 2, sd) %>% unname()) %>%
  mutate(mean.RQE.phyl = apply(df.RQE.phyl.null, 2, mean) %>% unname()) %>%
  mutate(median.RQE.phyl = apply(df.RQE.phyl.null, 2, median) %>% unname()) %>%
  mutate(sd.RQE.phyl = apply(df.RQE.phyl.null, 2, sd) %>% unname()) %>%
  mutate(mean.RQE.phyl.sqrt = apply(df.RQE.phyl.sqrt.null, 2, mean) %>% unname()) %>%
  mutate(median.RQE.phyl.sqrt = apply(df.RQE.phyl.sqrt.null, 2, median) %>% unname()) %>%
  mutate(sd.RQE.phyl.sqrt = apply(df.RQE.phyl.sqrt.null, 2, sd) %>% unname()) %>%
  mutate(SR.null = rownames(.))

lm(mean.RQE.phyl ~ median.RQE.phyl, data = PD.null) %>% summary() # distribution of expected value is symmetric
lm(mean.MPD ~ median.MPD, data = PD.null) %>% summary() # not sure here

PD.indices <- df.PD
PD.null.matrix <- list(df.MPD.null, df.RQE.phyl.null, df.MPD.sqrt.null, df.RQE.phyl.sqrt.null)
names(PD.null.matrix) <- c("MPD-null", "RQE-phyl-null", "MPD-sqrt-null", "RQE-phyl-sqrt-null")

PD.sqrt.indices <- df.PD.sqrt %>% dplyr::select(RQEP.sqrt = RaoD.phyl, MPD.sqrt = MPD, PlotObservationID)

save(PD.indices, PD.sqrt.indices, PD.null, PD.null.matrix, file = "02.data/03.PD-values.RData")

#### Functional Data----

# load data
FD.names <- list.files("02.data/eve-calculations/FD-calculation")

FD.list <- list()

for (i in 1:length(FD.names)) {
  FD.list[[i]] <- readRDS(paste0("02.data/eve-calculations/FD-calculation/", FD.names[i]))
}

df.FD <- bind_rows(FD.list)

# load FD null data

FD.null.names <- list.files("02.data/eve-calculations/FD-null-calculation", pattern = "*.Rds")

FD.null.list <- list()

for (i in 1:length(FD.null.names)) {
  FD.null.list[[i]] <- readRDS(paste0("02.data/eve-calculations/FD-null-calculation/", FD.null.names[i]))
}

list.FDis.SLA.null <- lapply(FD.null.list, `[[`, "FDis.SLA.null")
list.RQE.SLA.null <- lapply(FD.null.list, `[[`, "RQE.SLA.null")

df.FDis.SLA.null <- bind_cols(list.FDis.SLA.null)
df.RQE.SLA.null <- bind_cols(list.RQE.SLA.null)

list.FDis.HEIGHT.null <- lapply(FD.null.list, `[[`, "FDis.HEIGHT.null")
list.RQE.HEIGHT.null <- lapply(FD.null.list, `[[`, "RQE.HEIGHT.null")

df.FDis.HEIGHT.null <- bind_cols(list.FDis.HEIGHT.null)
df.RQE.HEIGHT.null <- bind_cols(list.RQE.HEIGHT.null)

list.FDis.ROOT.null <- lapply(FD.null.list, `[[`, "FDis.ROOT.null")
list.RQE.ROOT.null <- lapply(FD.null.list, `[[`, "RQE.ROOT.null")

df.FDis.ROOT.null <- bind_cols(list.FDis.ROOT.null)
df.RQE.ROOT.null <- bind_cols(list.RQE.ROOT.null)

list.FDis.MULTI.null <- lapply(FD.null.list, `[[`, "FDis.MULTI.null")
list.RQE.MULTI.null <- lapply(FD.null.list, `[[`, "RQE.MULTI.null")

df.FDis.MULTI.null <- bind_cols(list.FDis.MULTI.null)
df.RQE.MULTI.null <- bind_cols(list.RQE.MULTI.null)

FD.null <- apply(df.FDis.SLA.null, 2, mean) %>%
  as.data.frame() %>%
  rename(., mean.FDis.SLA = `.`) %>%
  # mutate(quant05.FDis.SLA = apply(df.FDis.SLA.null, 2, quantile, probs = .05) %>% unname()) %>%
  # mutate(quant95.FDis.SLA = apply(df.FDis.SLA.null, 2, quantile, probs = .95) %>% unname()) %>%
  mutate(sd.FDis.SLA = apply(df.FDis.SLA.null, 2, sd) %>% unname()) %>%
  mutate(median.FDis.SLA = apply(df.FDis.SLA.null, 2, median) %>% unname()) %>%
  mutate(SR.null = rownames(.)) %>%
  mutate(mean.RQE.SLA = apply(df.RQE.SLA.null, 2, mean) %>% unname()) %>%
  # mutate(quant05.RQE.SLA = apply(df.RQE.SLA.null, 2, quantile, probs = .05) %>% unname()) %>%
  # mutate(quant95.RQE.SLA = apply(df.RQE.SLA.null, 2, quantile, probs = .95) %>% unname()) %>%
  mutate(sd.RQE.SLA = apply(df.RQE.SLA.null, 2, sd) %>% unname()) %>%
  mutate(median.RQE.SLA = apply(df.RQE.SLA.null, 2, median) %>% unname()) %>%
  mutate(mean.FDis.HEIGHT = apply(df.FDis.HEIGHT.null, 2, mean) %>% unname()) %>%
  # mutate(quant05.FDis.HEIGHT = apply(df.FDis.HEIGHT.null, 2, quantile, probs = .05) %>% unname()) %>%
  # mutate(quant95.FDis.HEIGHT = apply(df.FDis.HEIGHT.null, 2, quantile, probs = .95) %>% unname()) %>%
  mutate(sd.FDis.HEIGHT = apply(df.FDis.HEIGHT.null, 2, sd) %>% unname()) %>%
  mutate(median.FDis.HEIGHT = apply(df.FDis.HEIGHT.null, 2, median) %>% unname()) %>%
  mutate(mean.RQE.HEIGHT = apply(df.RQE.HEIGHT.null, 2, mean) %>% unname()) %>%
  # mutate(quant05.RQE.HEIGHT = apply(df.RQE.HEIGHT.null, 2, quantile, probs = .05) %>% unname()) %>%
  # mutate(quant95.RQE.HEIGHT = apply(df.RQE.HEIGHT.null, 2, quantile, probs = .95) %>% unname()) %>%
  mutate(sd.RQE.HEIGHT = apply(df.RQE.HEIGHT.null, 2, sd) %>% unname()) %>%
  mutate(median.RQE.HEIGHT = apply(df.RQE.HEIGHT.null, 2, median) %>% unname()) %>%
  mutate(mean.FDis.ROOT = apply(df.FDis.ROOT.null, 2, mean) %>% unname()) %>%
  # mutate(quant05.FDis.ROOT = apply(df.FDis.ROOT.null, 2, quantile, probs = .05) %>% unname()) %>%
  # mutate(quant95.FDis.ROOT = apply(df.FDis.ROOT.null, 2, quantile, probs = .95) %>% unname()) %>%
  mutate(sd.FDis.ROOT = apply(df.FDis.ROOT.null, 2, sd) %>% unname()) %>%
  mutate(median.FDis.ROOT = apply(df.FDis.ROOT.null, 2, median) %>% unname()) %>%
  mutate(mean.RQE.ROOT = apply(df.RQE.ROOT.null, 2, mean) %>% unname()) %>%
  # mutate(quant05.RQE.ROOT = apply(df.RQE.ROOT.null, 2, quantile, probs = .05) %>% unname()) %>%
  # mutate(quant95.RQE.ROOT = apply(df.RQE.ROOT.null, 2, quantile, probs = .95) %>% unname()) %>%
  mutate(sd.RQE.ROOT = apply(df.RQE.ROOT.null, 2, sd) %>% unname()) %>%
  mutate(median.RQE.ROOT = apply(df.RQE.ROOT.null, 2, median) %>% unname()) %>%
  mutate(mean.FDis.MULTI = apply(df.FDis.MULTI.null, 2, mean) %>% unname()) %>%
  # mutate(quant05.FDis.MULTI = apply(df.FDis.MULTI.null, 2, quantile, probs = .05) %>% unname()) %>%
  # mutate(quant95.FDis.MULTI = apply(df.FDis.MULTI.null, 2, quantile, probs = .95) %>% unname()) %>%
  mutate(sd.FDis.MULTI = apply(df.FDis.MULTI.null, 2, sd) %>% unname()) %>%
  mutate(median.FDis.MULTI = apply(df.FDis.MULTI.null, 2, median) %>% unname()) %>%
  mutate(mean.RQE.MULTI = apply(df.RQE.MULTI.null, 2, mean) %>% unname()) %>%
  # mutate(quant05.RQE.MULTI = apply(df.RQE.MULTI.null, 2, quantile, probs = .05) %>% unname()) %>%
  # mutate(quant95.RQE.MULTI = apply(df.RQE.MULTI.null, 2, quantile, probs = .95) %>% unname()) %>%
  mutate(sd.RQE.MULTI = apply(df.RQE.MULTI.null, 2, sd) %>% unname()) %>%
  mutate(median.RQE.MULTI = apply(df.RQE.MULTI.null, 2, median) %>% unname())

lm(mean.RQE.MULTI ~ median.RQE.MULTI, data = FD.null) %>% summary() # distribution of expected value is symmetric
lm(mean.FDis.MULTI ~ median.FDis.MULTI, data = FD.null) %>% summary() # distribution of expected value is symmetric

FD.indices <- df.FD
FD.null.matrix <- list(
  df.FDis.SLA.null,
  df.RQE.SLA.null,
  df.FDis.HEIGHT.null,
  df.RQE.HEIGHT.null,
  df.FDis.ROOT.null,
  df.RQE.ROOT.null,
  df.FDis.MULTI.null,
  df.RQE.MULTI.null
)

names(FD.null.matrix) <- c(
  "FDis.SLA.null",
  "RQE.SLA.null",
  "FDis.HEIGHT.null",
  "RQE.HEIGHT.null",
  "FDis.ROOT.null",
  "RQE.ROOT.null",
  "FDis.MULTI.null",
  "RQE.MULTI.null"
)

save(FD.indices, FD.null, FD.null.matrix, file = "02.data/03.FD-values.RData")
