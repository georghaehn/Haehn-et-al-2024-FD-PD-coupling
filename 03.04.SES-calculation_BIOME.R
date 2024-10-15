#### recommended machine: 01####


#### packages ----
library("tidyverse")
library("tidyr")
library("data.table")
library("raster")
library("ncdf4")
library("sf")
library("rgdal")


# load phylo null BIOME data

PD.null.names <- list.files("02.data/eve-calculations/PD-null-biome", pattern = "*.Rds")

PD.null.list <- list()

for (i in 1:length(PD.null.names)) {
  PD.null.list[[i]] <- readRDS(paste0("02.data/eve-calculations/PD-null-biome/", PD.null.names[i]))
}

PD.ll <- lapply(PD.null.list, FUN = function(x) {
  x[[1]]$var <- colnames(x[[1]])[1]
  x[[2]]$var <- colnames(x[[2]])[1]
  colnames(x[[1]])[1] <- "value"
  colnames(x[[2]])[1] <- "value"
  rbind(x[[1]], x[[2]])
})

PD.null <- do.call(rbind, PD.ll)

PD.null.biome <- PD.null %>%
  group_by(SR, Biome, var) %>%
  summarise(
    mean.bio = mean(value),
    sd.bio = sd(value)
  ) %>%
  as.data.frame() %>%
  pivot_wider(values_from = c(mean.bio, sd.bio), names_from = var)

# load phylo null BIOME weighted data

PD.null.w.names <- list.files("02.data/eve-calculations/PD-null-biome-weighted", pattern = "*.Rds")

PD.null.w.list <- list()

for (i in 1:length(PD.null.w.names)) {
  PD.null.w.list[[i]] <- readRDS(paste0("02.data/eve-calculations/PD-null-biome-weighted/", PD.null.w.names[i]))
}

PD.w.ll <- lapply(PD.null.w.list, FUN = function(x) {
  x[[1]]$var <- colnames(x[[1]])[1]
  x[[2]]$var <- colnames(x[[2]])[1]
  colnames(x[[1]])[1] <- "value"
  colnames(x[[2]])[1] <- "value"
  rbind(x[[1]], x[[2]])
})

PD.w.null <- do.call(rbind, PD.w.ll)

PD.null.w.biome <- PD.w.null %>%
  mutate(
    var = gsub("MPD.null", "MPD.null.weigh", var),
    var = gsub("RQE.phyl.null", "RQE.phyl.null.weigh", var)
  ) %>%
  group_by(SR, Biome, var) %>%
  summarise(
    mean.bio = mean(value),
    sd.bio = sd(value)
  ) %>%
  as.data.frame() %>%
  pivot_wider(values_from = c(mean.bio, sd.bio), names_from = var)

# load funct div null BIOME data

FD.null.names <- list.files("02.data/eve-calculations/FD-null-biome", pattern = "*.Rds")

FD.null.list <- list()

for (i in 1:length(FD.null.names)) {
  FD.null.list[[i]] <- readRDS(paste0("02.data/eve-calculations/FD-null-biome/", FD.null.names[i]))
}

FD.ll <- lapply(FD.null.list, FUN = function(x) {
  x[[1]]$var <- colnames(x[[1]])[1]
  x[[2]]$var <- colnames(x[[2]])[1]
  colnames(x[[1]])[1] <- "value"
  colnames(x[[2]])[1] <- "value"
  rbind(x[[1]], x[[2]])
})

FD.null <- do.call(rbind, FD.ll)

FD.null.biome <- FD.null %>%
  group_by(SR, Biome, var) %>%
  summarise(
    mean.bio = mean(value),
    sd.bio = sd(value)
  ) %>%
  as.data.frame() %>%
  pivot_wider(values_from = c(mean.bio, sd.bio), names_from = var)

# load funct div null BIOME weighted data

FD.null.w.names <- list.files("02.data/eve-calculations/FD-null-biome-weighted", pattern = "*.Rds")

FD.null.w.list <- list()

for (i in 1:length(FD.null.w.names)) {
  FD.null.w.list[[i]] <- readRDS(paste0("02.data/eve-calculations/FD-null-biome-weighted/", FD.null.w.names[i]))
}

FD.w.ll <- lapply(FD.null.w.list, FUN = function(x) {
  x[[1]]$var <- colnames(x[[1]])[1]
  x[[2]]$var <- colnames(x[[2]])[1]
  colnames(x[[1]])[1] <- "value"
  colnames(x[[2]])[1] <- "value"
  rbind(x[[1]], x[[2]])
})

FD.w.null <- do.call(rbind, FD.w.ll)

FD.null.w.biome <- FD.w.null %>%
  mutate(
    var = gsub("FDis.multi.null", "FDis.multi.null.weigh", var),
    var = gsub("RQE.multi.null", "RQE.multi.null.weigh", var)
  ) %>%
  group_by(SR, Biome, var) %>%
  summarise(
    mean.bio = mean(value),
    sd.bio = sd(value)
  ) %>%
  as.data.frame() %>%
  pivot_wider(values_from = c(mean.bio, sd.bio), names_from = var)
# read compiled data

sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-data.RDS")

data <- PD.null.biome %>%
  mutate(m = paste0(Biome, "_", SR)) %>%
  dplyr::select(-c(SR, Biome)) %>%
  left_join(FD.null.biome %>%
    mutate(m = paste0(Biome, "_", SR)) %>%
    dplyr::select(-c(SR, Biome))) %>%
  left_join(PD.null.w.biome %>%
    mutate(m = paste0(Biome, "_", SR)) %>%
    dplyr::select(-c(SR, Biome))) %>%
  left_join(FD.null.w.biome %>%
    mutate(m = paste0(Biome, "_", SR)) %>%
    dplyr::select(-c(SR, Biome))) %>%
  left_join(sPlot.data %>%
    mutate(m = paste0(sBiome, "_", SR))) %>%
  mutate(
    SES.RQEP.B = (RaoD.phyl - mean.bio_RQE.phyl.null) / sd.bio_RQE.phyl.null,
    SES.RQEF.B = (RQE.MULTI - mean.bio_RQE.multi.null) / sd.bio_RQE.multi.null,
    SES.RQEP.BW = (RaoD.phyl - mean.bio_RQE.phyl.null.weigh) / sd.bio_RQE.phyl.null.weigh,
    SES.RQEF.BW = (RQE.MULTI - mean.bio_RQE.multi.null.weigh) / sd.bio_RQE.multi.null.weigh
  )

sum(is.na(sPlot.data$SES.RQEP))
sum(is.na(sPlot.data$SES.RQEF))
sum(is.na(sPlot.data$SES.RQEP.B))
sum(is.na(sPlot.data$SES.RQEP.BW))
sum(is.na(sPlot.data$SES.RQEF.BW))

saveRDS(data, "02.data/03.sPlot.PD.FD.CD-PCA-BW-data.RDS")

### add phylogenetic kingdoms

data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-BW-data.RDS")

# read data with assigned phylokingdoms

pk <- readRDS("02.data/phylo-kingdom-data.RDS")

data <- left_join(data, pk %>%
  dplyr::select(PlotObservationID, Phyl.clust = Cluster),
by = "PlotObservationID"
)

table(data[, c("Country", "Phyl.clust")])

# load phylo null

PD.null.names <- list.files("02.data/eve-calculations/PD-null-phyl-clust", pattern = "*.Rds")

PD.null.list <- list()

for (i in 1:length(PD.null.names)) {
  PD.null.list[[i]] <- readRDS(paste0("02.data/eve-calculations/PD-null-phyl-clust/", PD.null.names[i]))
}

PD.ll <- lapply(PD.null.list, FUN = function(x) {
  x[[1]]$var <- colnames(x[[1]])[1]
  x[[2]]$var <- colnames(x[[2]])[1]
  colnames(x[[1]])[1] <- "value"
  colnames(x[[2]])[1] <- "value"
  rbind(x[[1]], x[[2]])
})

PD.null <- do.call(rbind, PD.ll)

PD.null.pk <- PD.null %>%
  group_by(SR, Phyl.clust, var) %>%
  summarise(
    mean.pk = mean(value),
    sd.pk = sd(value)
  ) %>%
  as.data.frame() %>%
  pivot_wider(values_from = c(mean.pk, sd.pk), names_from = var)

# load funct div null

FD.null.names <- list.files("02.data/eve-calculations/FD-null-phyl-clust", pattern = "*.Rds")

FD.null.list <- list()

for (i in 1:length(FD.null.names)) {
  FD.null.list[[i]] <- readRDS(paste0("02.data/eve-calculations/FD-null-phyl-clust/", FD.null.names[i]))
}

FD.ll <- lapply(FD.null.list, FUN = function(x) {
  x[[1]]$var <- colnames(x[[1]])[1]
  x[[2]]$var <- colnames(x[[2]])[1]
  colnames(x[[1]])[1] <- "value"
  colnames(x[[2]])[1] <- "value"
  rbind(x[[1]], x[[2]])
})

FD.null <- do.call(rbind, FD.ll)

FD.null.pk <- FD.null %>%
  group_by(SR, Phyl.clust, var) %>%
  summarise(
    mean.pk = mean(value),
    sd.pk = sd(value)
  ) %>%
  as.data.frame() %>%
  pivot_wider(values_from = c(mean.pk, sd.pk), names_from = var)

data1 <- PD.null.pk %>%
  mutate(m = paste0(Phyl.clust, "_", SR)) %>%
  dplyr::select(-c(SR, Phyl.clust)) %>%
  left_join(FD.null.pk %>%
    mutate(m = paste0(Phyl.clust, "_", SR)) %>%
    dplyr::select(-c(SR, Phyl.clust))) %>%
  left_join(data %>%
    mutate(m = paste0(Phyl.clust, "_", SR))) %>%
  mutate(
    SES.RQEP.PK = (RaoD.phyl - mean.pk_RQE.phyl.null) / sd.pk_RQE.phyl.null,
    SES.RQEF.PK = (RQE.MULTI - mean.pk_RQE.multi.null) / sd.pk_RQE.multi.null
  )

sum(is.na(data1$SES.RQEP))
sum(is.na(data1$SES.RQEF))
sum(is.na(data1$SES.RQEP.PK))
sum(is.na(data1$SES.RQEF.PK))

saveRDS(data1, "02.data/03.sPlot.PD.FD.CD-PCA-BW-PK-data.RDS")

### available dataset

data1 <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-BW-PK-data.RDS")

data.public <- data1 %>%
  dplyr::select(
    Longitude, Latitude,
    PD, MPD, MNTD, RaoD.phyl, RQEP.sqrt,
    MPD.sqrt, SES.RQEP, SES.RQEP.sqrt, SES.MPD, SES.MPD.sqrt,
    SES.RQEP.B, SES.RQEP.BW, SES.RQEP.PK,
    RQE.MULTI, FDis.MULTI, SES.RQEF, SES.FDis,
    SES.RQEF.B, SES.RQEF.BW, SES.RQEF.PK,
    mean.pk_MPD.null:sd.bio_RQE.multi.null.weigh,
    mean.MPD:sd.RQE.phyl.sqrt,
    mean.FDis.SLA:sd.RQE.MULTI
  ) %>%
  mutate(
    Longitude = round(jitter(round(Longitude, digits = 2)), digits = 2),
    Latitude = round(jitter(round(Latitude, digits = 2)), digits = 2)
  ) %>%
  filter(!is.na(SES.RQEF.PK)) %>%
  dplyr::select(-m)

write.csv(data.public, "02.data/sPlot_FD-PD_data.csv")

data1.open <- readRDS("02.data/01c.sPlotOpen.PD.FD.CD-data.RDS")

data.open <- data1.open %>%
  dplyr::select(
    Longitude, Latitude,
    PD, MPD, MNTD, RaoD.phyl, SES.RQEP, SES.MPD,
    FDis.SLA:RQE.MULTI,
    SES.RQEF, SES.FDis,
    mean.MPD:sd.RQE.phyl,
    mean.FDis.SLA:sd.RQE.MULTI
  )

write.csv(data.open, "02.data/sPlotOpen_FD-PD_data.csv")
