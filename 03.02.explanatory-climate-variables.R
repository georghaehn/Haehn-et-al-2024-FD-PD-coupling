####recommended machine: 01####

####packages ----
library("tidyverse")
library("tidyr")
library("dplyr")
library("data.table")
library("raster")
library("ncdf4")
library("sf")
library("rgdal")
library("RStoolbox")
library("downloader")

#### CHELSA data ----

CHELSA.stack <- stack(paste("02.data/CHELSA_clim-data/CHELSA_bio10_", 
                            c(1:19), ".tif", sep=""))

load("02.data/header_sPlot3.0.RData")

sPlot.spa <- SpatialPointsDataFrame(coords = header %>% 
                                      filter(!is.na(Longitude)) %>% 
                                      filter(!is.na(Latitude)) %>% 
                                      dplyr::select(Longitude, Latitude), 
                                    data = data.frame(PlotObservationID = header %>% 
                                                        filter(!is.na(Longitude)) %>% 
                                                        filter(!is.na(Latitude)) %>% 
                                                        pull(PlotObservationID)),
                                    proj4string = crs("+proj=longlat +datum=WGS84 +no_defs"))


chelsa.world <- raster::extract(CHELSA.stack, sPlot.spa)

chelsa.world <- as.data.frame(chelsa.world) %>% 
  rename_at(.vars=vars(starts_with("CHELSA")),
            .funs=~stringr::str_sub(., -8,-1))

colnames(chelsa.world)[1:9] <- c("bio10_1", "bio10_2", "bio10_3", "bio10_4", "bio10_5", "bio10_6",
                                 "bio10_7", "bio10_8", "bio10_9")

#load climate data
piThresh <- readRDS("02.data/clim/StableClim_V1.0.1/StableClim_piControl_thresholds.RDS")
pastReg <- readRDS("02.data/clim/StableClim_V1.0.1/StableClim_past_RegionalTemperatureRegressions.RDS")
#futReg <- readRDS("02.data/clim/StableClim_V1.0.1/StableClim_HistoricalRCP_RegionalTemperatureRegressions.RDS")
# Past regression rasters
pastRastTrend <- brick("02.data/clim/StableClim_V1.0.1/ncdf/regressions/StableClim_Regression_past_ts.nc",
                       varname = "ts_trend")
pastRastVar <- brick("02.data/clim/StableClim_V1.0.1/ncdf/regressions/StableClim_Regression_past_ts.nc",
                     varname = "ts_variability")
pastRastSNR <- brick("02.data/clim/StableClim_V1.0.1/ncdf/regressions/StableClim_Regression_past_ts.nc",
                     varname = "ts_snr")

# Geopackage containing the regions
regionalShp <- read_sf("02.data/clim/StableClim_V1.0.1/gpkg/StableClim_VectorData.gpkg", layer = "WallaceRegions")

periods <- piThresh$`All Periods [signed]`
thresh <- periods[RegionType == "Global" & Region == "Land/Sea", ][["90%"]] ## threshold in °C/Year
thresh

pastReg <- pastReg[["Global.Land/Sea"]]
#futReg <- futReg[["rcp26.Global.Land"]]

pastReg_RCC <- pastReg[pastReg$Slope <= thresh, ]
#futReg_RCC <- futReg[futReg$Slope >= thresh, ]

# make Start names the same format as the layer names
layerIDX <- gsub("-", "X\\.", pastReg_RCC$Start)
layerIDX <- which(names(pastRastTrend) %in% layerIDX)
# Subset the raster data, and calculate median values through time
rcc_past_trend <- calc(readAll(pastRastTrend[[layerIDX]]), median)
rcc_past_var <- calc(readAll(pastRastVar[[layerIDX]]), median)
rcc_past_snr <- calc(readAll(pastRastSNR[[layerIDX]]), median)

world_raster <- rasterize(as(regionalShp, "Spatial"),
                          y = raster(res = res(rcc_past_var)),
                          getCover = TRUE)

world_trend <- mask(rcc_past_trend*100, world_raster)
world_var <- mask(rcc_past_var, world_raster)
world_snr <- mask(rcc_past_snr, world_raster)

plot(world_var)

saveRDS(world_var, file = "02.data/03.climate-variability.RDS")

stable.clim <- readRDS("02.data/03.climate-variability.RDS")

stable.clim.world <- raster::extract(stable.clim, sPlot.spa)

sPlot.out <- sPlot.spa@data %>% 
  as.tbl() %>% 
  bind_cols(as.data.frame(sPlot.spa@coords)) %>% 
  bind_cols(chelsa.world) %>% 
  bind_cols(stable.clim.world)

colnames(sPlot.out)

colnames(sPlot.out)[colnames(sPlot.out)=="...23"] <- "stable.clim"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_1"] <- "mean.annual.air.temp"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_2"] <- "mean.diurnal.air.temp.range"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_3"] <- "isothermality"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_4"] <- "temp.season"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_5"] <- "mean.daily.max.air.temp.warm.m"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_6"] <- "mean.daily.min.air.temp.cold.m"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_7"] <- "annual.range.air.temp"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_8"] <- "mean.daily.air.temp.wet.qu"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_9"] <- "mean.daily.air.temp.dry.qu"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_10"] <- "mean.daily.air.temp.warm.qu"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_11"] <- "mean.daily.air.temp.cold.qu"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_12"] <- "annual.prec"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_13"] <- "prec.wet.m"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_14"] <- "prec.dry.m"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_15"] <- "prec.season"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_16"] <- "mean.monthly.prec.wet.qu"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_17"] <- "mean.monthly.prec.dry.qu"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_18"] <- "mean.monthly.prec.warm.qu"
colnames(sPlot.out)[colnames(sPlot.out)=="bio10_19"] <- "mean.monthly.prec.cold.qu"

colnames(sPlot.out)

saveRDS(sPlot.out, "02.data/03.climate-data.RDS")

rm(list = ls())

#### calc standardized effect size ----
####load data ----

load("02.data/header_sPlot3.0.RData")
load("02.data/DT_sPlot3.0.RData")
load("02.data/pruned.trees.traitsp.RData")

#add species richness for each plot, SR is based on species present in the phylogeny and TRY

Plots <- DT2 %>% 
  left_join(., DT2 %>% group_by(PlotObservationID) %>% 
              summarize_at(.vars = "Abundance", .funs = sum, na.rm = T) %>% 
              rename("sum.Abundance" = "Abundance"),
            by = "PlotObservationID") %>% 
  mutate(Abu = Abundance/sum.Abundance) %>% 
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
    by = "PlotObservationID")

rm(DT2, PD.out0)

load("02.data/03.PD-values.RData")
load("02.data/03.FD-values.RData")
clim.var <- readRDS("02.data/03.climate-data.RDS")


rm(FD.null.matrix, PD.null.matrix)

#merge data
df1 <- header %>% 
  left_join(., clim.var %>% 
              dplyr::select(-c("Longitude", "Latitude")),
            by = "PlotObservationID") %>% 
  left_join(., PD.indices,
            by = "PlotObservationID") %>% 
  left_join(., PD.sqrt.indices,
            by = "PlotObservationID") %>% 
  left_join(., PD.null %>% 
              mutate(SR.null = as.numeric(SR.null)),
            by = c("species.richness" = "SR.null")) %>% 
  left_join(., FD.indices,
            by = "PlotObservationID") %>% 
  left_join(., FD.null %>% 
              mutate(SR.null = as.numeric(SR.null)),
            by = c("species.richness" = "SR.null")) 
rm(header)

#### indices ----
df1$SES.RQEP <- (df1$RaoD.phyl - df1$mean.RQE.phyl)/df1$sd.RQE.phyl
df1$SES.RQEP.sqrt <- (df1$RQEP.sqrt - df1$mean.RQE.phyl.sqrt)/df1$sd.RQE.phyl.sqrt
df1$SES.RQEF <- (df1$RQE.MULTI - df1$mean.RQE.MULTI)/df1$sd.RQE.MULTI
df1$SES.MPD <- (df1$MPD - df1$mean.MPD)/df1$sd.MPD
df1$SES.MPD.sqrt <- (df1$MPD.sqrt - df1$mean.MPD.sqrt)/df1$sd.MPD.sqrt
df1$SES.FDis <- (df1$FDis.MULTI - df1$mean.FDis.MULTI)/df1$sd.FDis.MULTI

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
  filter(SES.RQEF < 5e+10) %>%  #exclude outlier
  as.data.frame() %>%
  rename("plot.size" = "Relevé area (m²)")


library("dggridR")

dggs <- dgconstruct(area=30*30, metric=T, resround='down')

df1.1$cell <- dgGEO_to_SEQNUM(dggs, df1.1$Longitude, df1.1$Latitude)$seqnum

saveRDS(df1.1, "02.data/03.sPlot.PD.FD.CD-data.RDS")

rm(list = ls())

#### PCA ----

library("dplyr")
library("sp")
library("tidyr")
library("dismo")
library("gbm")
library("mgcv")
library("ggplot2")
library("ggfortify")
library("corrplot")
library("rgl")
library("ggrepel")

sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-data.RDS")

clim.PCA <- sPlot.data %>% 
  as.data.frame() %>% 
  dplyr::select(c(mean.annual.air.temp:mean.monthly.prec.cold.qu)) %>% 
  rename("Mean annual air temperature" = "mean.annual.air.temp") %>%
  rename("Mean diurnal air temperature range" = "mean.diurnal.air.temp.range") %>%
  rename("Isothermality" = "isothermality") %>%
  rename("Temperature seasonality" = "temp.season") %>%
  rename("Mean daily maximum air temperature \n warmest month" = "mean.daily.max.air.temp.warm.m") %>% 
  rename("Mean daily minimum air temperature \n coldest month" = "mean.daily.min.air.temp.cold.m") %>% 
  rename("Annual range air temperature" = "annual.range.air.temp") %>% 
  rename("Mean daily air temperature \n wettest quarter" = "mean.daily.air.temp.wet.qu") %>% 
  rename("Mean daily air temperature \n driest quarter" = "mean.daily.air.temp.dry.qu") %>% 
  rename("Mean daily air temperature \n warmest quarter" = "mean.daily.air.temp.warm.qu") %>% 
  rename("Mean daily air temperature \n coldest quarter" = "mean.daily.air.temp.cold.qu") %>% 
  rename("Annual precipitation" = "annual.prec") %>% 
  rename("Precipitation wettest month" = "prec.wet.m") %>% 
  rename("Precipitation driest month" = "prec.dry.m") %>% 
  rename("Precipitation seasonality" = "prec.season") %>% 
  rename("Mean monthly precipitation \n wettest quarter" = "mean.monthly.prec.wet.qu") %>% 
  rename("Mean monthly precipitation \n driest quarter" = "mean.monthly.prec.dry.qu") %>% 
  rename("Mean monthly precipitation \n warmest quarter" = "mean.monthly.prec.warm.qu") %>% 
  rename("Mean monthly precipitation \n coldest quarter" = "mean.monthly.prec.cold.qu")
#SES.RQEP, SES.RQEF, SES.MPD, SES.FDis))

PCA.clim <- prcomp(clim.PCA, center = T, scale = T)

expl.var <- round((PCA.clim$sdev^2/sum(PCA.clim$sdev^2)*100),1)
### ENVFIT - first 5 PC axes explain >90% of variation
set.seed(15)
sa.chelsa <- sample(1:nrow(clim.PCA), 10000)
pca.envfit.chelsa <- vegan::envfit(PCA.clim$x[sa.chelsa,], clim.PCA %>%
                                     slice(sa.chelsa),
                                   choices = c(1:5))
print(pca.envfit.chelsa)

vec.sp.df <- round(as.data.frame(pca.envfit.chelsa$vectors$arrows*sqrt(pca.envfit.chelsa$vectors$r)) *10 ,
                   digits = 2)
vec.sp.df$clim <- names(pca.envfit.chelsa$vectors$r)
write.csv(vec.sp.df, "__Submission/Figures/PCA-tab.csv")

vec.sp.df <- as.data.frame(pca.envfit.chelsa$vectors$arrows*sqrt(pca.envfit.chelsa$vectors$r)) *10
vec.sp.df$Var <- paste0("Clim", 1:nrow(vec.sp.df))
vec.sp.df$Names <- rownames(pca.envfit.chelsa$vectors$arrows)

(p1 <- ggplot(data= PCA.clim$x %>% 
                as.data.frame() %>% 
                as.tbl()) + #  %>% 
    #sample_frac(0.05)) + 
    geom_hex(aes(x=PC1, y=PC2), col=gray(0.5), lwd=0.1) +
    geom_segment(data=vec.sp.df,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.1, "cm")), alpha=0.9,
                 colour=gray(0.1)) + 
    geom_text_repel(data=vec.sp.df, aes(x=PC1, y=PC2, label=Var), size=10, segment.alpha = 1/3
                    #position = position_dodge(width = 1) 
    )+
    coord_fixed() +
    scale_x_continuous(limits = c(-10, 32), name = paste("PC1 (", expl.var[1], "%)", sep="")) + 
    scale_y_continuous(name = paste("PC2 (", expl.var[2], "%)", sep="")) + 
    scale_fill_distiller(palette = "Spectral",
                         trans = "log",
                         limits = c(1, 3000),
                         breaks = c(10, 100, 1000)) + 
    theme_bw() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size=36),
          legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          axis.title.x = element_text(size = 40),
          axis.text.x = element_text(size = 36),
          axis.title.y = element_text(size = 40),
          axis.text.y = element_text(size = 36)) )


ggsave("03.results/03.PCA-clim_PC1-PC2.png", p1, height=15, width=20, units="in", dpi=300)

paste0(vec.sp.df$Var, ":", vec.sp.df$Names)


sPlot.out <- sPlot.data %>%
  bind_cols(as.data.frame(PCA.clim$x) %>%
              dplyr::select(PC1:PC5) %>%
              mutate_at(.vars = vars(PC2), #flip PCA axes pointing to negative direction
                        .funs = list(~ -.)))


saveRDS(sPlot.out, "02.data/03.sPlot.PD.FD.CD-PCA-data.RDS")

