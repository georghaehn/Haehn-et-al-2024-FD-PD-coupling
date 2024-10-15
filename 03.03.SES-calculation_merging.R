#### recommended machine: 01####

#### packages ----
library("tidyverse")
library("tidyr")
library("data.table")

#### load data ----

load("02.data/header_sPlot3.0.RData")
load("02.data/DT_sPlot3.0.RData")
load("02.data/pruned.trees.traitsp.RData")

# add species richness for each plot, SR is based on species present in the phylogeny and TRY

Plots <- DT2 %>%
  left_join(., DT2 %>% group_by(PlotObservationID) %>%
    summarize_at(.vars = "Abundance", .funs = sum, na.rm = T) %>%
    rename("sum.Abundance" = "Abundance"),
  by = "PlotObservationID"
  ) %>%
  mutate(Abu = Abundance / sum.Abundance) %>%
  filter(Species %in% sPlot.in.Tree$Species) %>%
  group_by(PlotObservationID) %>%
  summarise_at(.vars = "Abu", .funs = sum, na.rm = T) %>%
  filter(Abu >= .5) %>%
  pull(PlotObservationID)

header <- header %>%
  filter(PlotObservationID %in% Plots) %>%
  left_join(
    DT2 %>%
      filter(Species %in% sPlot.in.Tree$Species) %>%
      group_by(PlotObservationID) %>%
      summarise(species.richness = n_distinct(Species)),
    by = "PlotObservationID"
  )

rm(DT2, PD.out0)

# the same for the sPlotOpen data

load("02.data/sPlotOpen.Rdata")

Plots.open <- DT2.oa %>%
  left_join(., DT2.oa %>% group_by(PlotObservationID) %>%
    summarize_at(.vars = "Original_abundance", .funs = sum, na.rm = T) %>%
    rename("sum.Abundance" = "Original_abundance"),
  by = "PlotObservationID"
  ) %>%
  mutate(Abu = Original_abundance / sum.Abundance) %>%
  filter(Species %in% sPlot.in.Tree$Species) %>%
  group_by(PlotObservationID) %>%
  summarise_at(.vars = "Abu", .funs = sum, na.rm = T) %>%
  filter(Abu >= .5) %>%
  pull(PlotObservationID)

header.oa <- header.oa %>%
  filter(PlotObservationID %in% Plots.open) %>%
  left_join(
    DT2.oa %>%
      filter(Species %in% sPlot.in.Tree$Species) %>%
      group_by(PlotObservationID) %>%
      summarise(species.richness = n_distinct(Species)),
    by = "PlotObservationID"
  )

rm(DT2.oa, CWM_CWV.oa, metadata.oa, reference.oa)

load("02.data/01.PD-values.RData")
PD.open <- readRDS("02.data/01a.PD-values-sPlotOpen.Rds")
load("02.data/01.FD-values.RData")
FD.open <- readRDS("02.data/01a.FD-values-sPlotOpen.Rds")
clim.var <- readRDS("02.data/01b.climate-data.RDS")
clim.open <- readRDS("02.data/01b.climate-data-sPlotOpen.RDS")

rm(FD.null.matrix, PD.null.matrix)

# merge data
df1 <- header %>%
  left_join(., clim.var %>%
    dplyr::select(-c("Longitude", "Latitude")),
  by = "PlotObservationID"
  ) %>%
  left_join(., PD.indices,
    by = "PlotObservationID"
  ) %>%
  left_join(., PD.null %>%
    mutate(SR.null = as.numeric(SR.null)),
  by = c("species.richness" = "SR.null")
  ) %>%
  left_join(., FD.indices,
    by = "PlotObservationID"
  ) %>%
  left_join(., FD.null %>%
    mutate(SR.null = as.numeric(SR.null)),
  by = c("species.richness" = "SR.null")
  )
rm(header)

df.open <- header.oa %>%
  left_join(., clim.open %>%
    dplyr::select(-c("Longitude", "Latitude")),
  by = "PlotObservationID"
  ) %>%
  left_join(., PD.open,
    by = "PlotObservationID"
  ) %>%
  left_join(., PD.null %>%
    mutate(SR.null = as.numeric(SR.null)),
  by = c("species.richness" = "SR.null")
  ) %>%
  left_join(., FD.open,
    by = "PlotObservationID"
  ) %>%
  left_join(., FD.null %>%
    mutate(SR.null = as.numeric(SR.null)),
  by = c("species.richness" = "SR.null")
  )

rm(header.oa)

#### indices ----
df1$SES.RQEP <- (df1$RaoD.phyl - df1$mean.RQE.phyl) / df1$sd.RQE.phyl
df1$SES.RQEF <- (df1$RQE.MULTI - df1$mean.RQE.MULTI) / df1$sd.RQE.MULTI
df1$SES.MPD <- (df1$MPD - df1$mean.MPD) / df1$sd.MPD
df1$SES.FDis <- (df1$FDis.MULTI - df1$mean.FDis.MULTI) / df1$sd.FDis.MULTI

df.open$SES.RQEP <- (df.open$RaoD.phyl - df.open$mean.RQE.phyl) / df.open$sd.RQE.phyl
df.open$SES.RQEF <- (df.open$RQE.MULTI - df.open$mean.RQE.MULTI) / df.open$sd.RQE.MULTI
df.open$SES.MPD <- (df.open$MPD - df.open$mean.MPD) / df.open$sd.MPD
df.open$SES.FDis <- (df.open$FDis.MULTI - df.open$mean.FDis.MULTI) / df.open$sd.FDis.MULTI

df1.1 <- df1 %>%
  filter(
    !is.na(df1$SES.RQEP),
    !is.na(df1$SES.RQEF),
    !is.na(df1$SES.MPD),
    !is.na(df1$SES.FDis),
    !is.na(Latitude),
    !is.na(Longitude),
    species.richness > 1
  ) %>%
  filter(SES.RQEF < 5e+10) %>% # exclude outlier
  as.data.frame()

df.open.1 <- df.open %>%
  filter(
    !is.na(df.open$SES.RQEP),
    !is.na(df.open$SES.RQEF),
    !is.na(df.open$SES.MPD),
    !is.na(df.open$SES.FDis),
    !is.na(Latitude),
    !is.na(Longitude),
    species.richness > 1
  ) %>%
  filter(SES.RQEF < 5e+10) %>%
  as.data.frame()

saveRDS(df1.1, "02.data/01c.sPlot.PD.FD.CD-data.RDS")
saveRDS(df.open.1, "02.data/01c.sPlotOpen.PD.FD.CD-data.RDS")
