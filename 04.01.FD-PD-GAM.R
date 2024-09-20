####recommended machine: 02####

library("tidyverse")
library("sp")
library("tidyr")
library("dismo")
library("gbm")
library("mgcv")
library("marginaleffects")
library("viridis")
library("doSNOW")

sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-BW-data.RDS")

##funct ~ phyl sPlot open ----

sPlot.oa <- readRDS("02.data/03.01.sPlotOpen-header-FD-PD.RDS")

coords <- sPlot.oa  %>% dplyr::select(Longitude, Latitude)
indices <- sPlot.oa %>%  
  dplyr::select(
    PDQ.PA,
    FDQ.ALL.PA,
    FDQ.MULTI.PA,
    FDQ.SLA.PA,
    FDQ.HEIGHT.PA,
    FDQ.ROOT.PA,
    FDQ.SM.PA,
    FDQ.LDMC.PA,
    FDQ.N.PA,
    FDQ.P.PA,
    FDQ.CHRO.PA)

data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                   data = indices,
                                   proj4string = CRS("+proj=eck4"))

# spatial.data.wgs has to be a spatial dataframe
# orig observed data

threads <- 16

cl <- makeCluster(threads)
registerDoSNOW(cl)

result <- foreach(i=2:ncol(data.wgs@data)) %dopar% {
  
  form <- paste(colnames(data.wgs@data)[i], "~ PDQ.PA + s(Longitude, Latitude, bs = \"sos\")")
  
  mgcv::gam(as.formula(form), family = "gaussian", method = "REML", data = data.wgs)
  
}
stopCluster(cl)

for(i in 1:length(result)) {
  print(summary(result[[i]]))
}

saveRDS(result, "02.data/sPlotOpen-PA-GAMS.Rds")

### Abundance data
# 
# coords <- sPlot.oa  %>% dplyr::select(Longitude, Latitude)
# indices <- sPlot.oa %>%  
#   dplyr::select(
#     PDQ.ABU,
#     FDQ.ALL.ABU,
#     FDQ.SUB.ABU,
#     FDQ.MULTI.ABU,
#     FDQ.SLA.ABU,
#     FDQ.HEIGHT.ABU,
#     FDQ.ROOT.ABU,
#     FDQ.SM.ABU,
#     FDQ.LDMC.ABU,
#     FDQ.GF.ABU,
#     FDQ.N.ABU,
#     FDQ.P.ABU,
#     FDQ.CHRO.ABU)
# 
# data.wgs <- SpatialPointsDataFrame(coords = coords, 
#                                    data = indices,
#                                    proj4string = CRS("+proj=eck4"))
# 
# # spatial.data.wgs has to be a sABUtial dataframe
# # orig observed data
# 
# threads <- 16
# 
# cl <- makeCluster(threads)
# registerDoSNOW(cl)
# 
# result <- foreach(i=2:ncol(data.wgs@data)) %dopar% {
#   
#   form <- paste(colnames(data.wgs@data)[i], "~ PDQ.ABU + s(Longitude, Latitude, bs = \"sos\")")
#   
#   mgcv::gam(as.formula(form), family = "gaussian", method = "REML", data = data.wgs)
#   
# }
# stopCluster(cl)
# 
# for(i in 1:length(result)) {
#   print(summary(result[[i]]))
# }
# 
# saveRDS(result, "02.data/sPlotOpen-ABU-GAMS.Rds")

## RQEF ~ RQEP sPlot3 ----
# 
# #for - nonfro -> loose 284698 plots with no classification
# data.for <- sPlot.data %>% filter(is.forest == T)
# coords.for <- data.for %>% dplyr::select(Longitude, Latitude)
# indices.for <- data.for %>%  
#   dplyr::select(
#     RQEP = RaoD.phyl, RQEF = RQE.MULTI, 
#     RQEP.sqrt,
#     SES.RQEP, SES.RQEF,
#     SES.RQEP.sqrt)
# 
# data.for.wgs <- SpatialPointsDataFrame(coords = coords.for, 
#                                    data = indices.for,
#                                    proj4string = crs("+proj=eck4"))
# 
# data.nonfor <- sPlot.data %>% filter(is.forest == F)
# coords.nonfor <- data.nonfor %>% dplyr::select(Longitude, Latitude)
# indices.nonfor <- data.nonfor %>%  
#   dplyr::select(
#     RQEP = RaoD.phyl, RQEF = RQE.MULTI, 
#     RQEP.sqrt,
#     SES.RQEP, SES.RQEF,
#     SES.RQEP.sqrt)
# 
# data.nonfor.wgs <- SpatialPointsDataFrame(coords = coords.nonfor, 
#                                        data = indices.nonfor,
#                                        proj4string = crs("+proj=eck4"))
# 
# 
coords <- sPlot.data %>% 
  filter(!is.na(Longitude)) %>% 
  dplyr::select(Longitude, Latitude)
indices <- sPlot.data %>% 
  filter(!is.na(Longitude)) %>%
  dplyr::select(
    RaoD.phyl,
    RQE.MULTI,
    RQEP.sqrt,
    SES.RQEP,
    SES.RQEF,
    SES.RQEP.sqrt)

data.wgs <- SpatialPointsDataFrame(coords = coords,
                                   data = indices,
                                   proj4string = CRS("+proj=eck4"))

# spatial.data.wgs has to be a spatial dataframe
# rm(indices, coords, sPlot.data, data.for, data.nonfor, indices.for, indices.nonfor, coords.for, coords.nonfor)
# gc()

# orig observed data
model1 <- mgcv::gam(RQE.MULTI ~ RaoD.phyl + s(Longitude, Latitude, 
                                            bs = "sos"), 
                    family = "gaussian", 
                    method = "REML", 
                    data = data.wgs)

saveRDS(model1, "02.data/04.GAM-RQEF-RQEP.Rds")

sum1 <- summary(model1)

rm(model1)

# orig observed functional diversity, phyl div based on sqrt(distance)
model2 <- mgcv::gam(RQE.MULTI ~ RQEP.sqrt + s(Longitude, Latitude, 
                                            bs = "sos"), 
                    family = "gaussian", 
                    method = "REML", 
                    data = data.wgs)

saveRDS(model2, "02.data/04.GAM-RQEF-RQEP.sqrt.Rds")

sum2 <- summary(model2)
        
rm(model2)

model2.1 <- mgcv::gam(RQE.MULTI ~ 1 + s(Longitude, Latitude, 
                                                bs = "sos"), 
                      family = "gaussian", 
                      method = "REML", 
                      data = data.wgs)

saveRDS(model2.1, "02.data/04.GAM-RQEF-RQEP.Smooth.Rds")

summary(model2.1)

rm(model2.1)

# standardized effect size of FD + PD
model3 <- mgcv::gam(SES.RQEF ~ SES.RQEP + s(Longitude, Latitude, 
                                            bs = "sos"), 
                    family = "gaussian", 
                    method = "REML", 
                    data = data.wgs)

saveRDS(model3, "02.data/04.GAM-SES-RQEF-RQEP.Rds")

sum3 <- summary(model3)

rm(model3)

model3.1 <- mgcv::gam(SES.RQEF ~ 1 + s(Longitude, Latitude, 
                                            bs = "sos"), 
                    family = "gaussian", 
                    method = "REML", 
                    data = data.wgs)

saveRDS(model3.1, "02.data/04.GAM-SES-RQEF-smooth.Rds")

sum3 <- summary(model3.1)

rm(model3)

# standardized effect size of FD and PD based on sqrt(distance)
model4 <- mgcv::gam(SES.RQEF ~ SES.RQEP.sqrt + s(Longitude, Latitude, 
                                            bs = "sos"), 
                    family = "gaussian", 
                    method = "REML", 
                    data = data.wgs)

saveRDS(model4, "02.data/04.GAM-SES-RQEF-RQEP.sqrt.Rds")

sum4 <- summary(model4)

rm(model4)

sum1
sum2
sum3
sum4

# SES based on biomes and weighted

coords <- sPlot.data %>% 
  filter(!is.na(Longitude)) %>% 
  dplyr::select(Longitude, Latitude)
indices <- sPlot.data %>% 
  filter(!is.na(Longitude)) %>% 
  dplyr::select(
    SES.RQEF.B,
    SES.RQEF.BW,
    SES.RQEP.B,
    SES.RQEP.BW)

data.wgs <- SpatialPointsDataFrame(coords = coords,
                                   data = indices,
                                   proj4string = CRS("+proj=eck4"))

modb <- mgcv::gam(SES.RQEF.B ~ SES.RQEP.B + s(Longitude, Latitude, 
                                      bs = "sos"), 
              family = "gaussian", 
              method = "REML", 
              data = data.wgs)

saveRDS(modb, "02.data/04.GAM-SES-RQEF-RQEP.bio.Rds")

modb1 <- mgcv::gam(SES.RQEF.B ~ 1 + s(Longitude, Latitude, 
                                              bs = "sos"), 
                  family = "gaussian", 
                  method = "REML", 
                  data = data.wgs)

saveRDS(modb1, "02.data/04.GAM-SES-RQEF-RQEP.bio.smooth.Rds")

modbw <- mgcv::gam(SES.RQEF.BW ~ SES.RQEP.BW + s(Longitude, Latitude, 
                                      bs = "sos"), 
              family = "gaussian", 
              method = "REML", 
              data = data.wgs)

saveRDS(modbw, "02.data/04.GAM-SES-RQEF-RQEP.bio.weigh.Rds")

modbw.1 <- mgcv::gam(SES.RQEF.BW ~ 1 + s(Longitude, Latitude, 
                                         bs = "sos"), 
                     family = "gaussian", 
                     method = "REML", 
                     data = data.wgs)

saveRDS(modbw.1, "02.data/04.GAM-SES-RQEF.bio.weigh-smooth.Rds")


## gams for phyl cluster

data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-BW-PK-data.RDS")
colnames(data)

coords <- data  %>% 
  filter(!is.na(SES.RQEF.PK)) %>% 
  dplyr::select(Longitude, Latitude)

indices <- data %>%  
  dplyr::select(
    SES.RQEP.PK,
    SES.RQEF.PK) %>% 
  filter(!is.na(SES.RQEF.PK))

data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                   data = indices,
                                   proj4string = CRS("+proj=eck4"))

modpk <- mgcv::gam(SES.RQEF.PK ~ SES.RQEP.PK + s(Longitude, Latitude, 
                                                 bs = "sos"), 
                   family = "gaussian", 
                   method = "REML", 
                   data = data.wgs)

summary(modpk)

saveRDS(modpk, "02.data/04.GAM-SES-RQEF-RQEP.PK.Rds")

modpk.1 <- mgcv::gam(SES.RQEF.PK ~ 1 + s(Longitude, Latitude, 
                                         bs = "sos"), 
                     family = "gaussian", 
                     method = "REML", 
                     data = data.wgs)

summary(modpk.1)

saveRDS(modpk.1, "02.data/04.GAM-SES-RQEF.PK-smooth.Rds")


# make the predictions for each model to get the residuals
mod.FD.PD <- readRDS("02.data/04.GAM-RQEF-RQEP.Rds")

summary(mod.FD.PD)

pred.FD.PD <- predictions(mod.FD.PD)

saveRDS(pred.FD.PD, "02.data/04.GAM-RQEF-RQEP_predictions.Rds")

mod.FD.PD.sqrt <- readRDS("02.data/04.GAM-RQEF-RQEP.sqrt.Rds")

summary(mod.FD.PD.sqrt)

pred.FD.PD.sqrt <- predictions(mod.FD.PD.sqrt)

saveRDS(pred.FD.PD.sqrt, "02.data/04.GAM-RQEF-RQEP.sqrt_predictions.Rds")

mod.FD.PD.SES <- readRDS("02.data/04.GAM-SES-RQEF-RQEP.Rds")

summary(mod.FD.PD.SES)

pred.FD.PD.SES <- predictions(mod.FD.PD.SES)

saveRDS(pred.FD.PD.SES, "02.data/04.GAM-SES-RQEF-RQEP_predictions.Rds")

mod.FD.PD.sqrt.SES <- readRDS("02.data/04.GAM-SES-RQEF-RQEP.sqrt.Rds")

summary(mod.FD.PD.sqrt.SES)

pred.FD.PD.sqrt.SES <- predictions(mod.FD.PD.sqrt.SES)

saveRDS(pred.FD.PD.sqrt.SES, "02.data/04.GAM-SES-RQEF-RQEP.sqrt_predictions.Rds")

#####

sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-BW-data.RDS") %>% 
  as.data.frame() %>% 
  filter(!is.na(Longitude),
         !is.na(SES.RQEF.BW))

# check relationship of SES.PD and devisions

dev <- readRDS("02.data/perc-fern-sPlot.RDS")

spp <- sPlot.data %>% 
  left_join(dev,
            by = "PlotObservationID")

lmh <- lm(SES.RQEP.BW ~ Other, spp)

summary(lmh)

lmm <- lm(SES.RQEP.BW ~ Other, spp %>% 
            filter(Other != 100))
summary(lmm)

lmmf <- lm(SES.RQEP.BW ~ Fern, spp %>% 
             filter(Other != 100))
summary(lmmf)

lmmg <- lm(SES.RQEP.BW ~ Gymnosperm, spp %>% 
             filter(Other != 100))
summary(lmmg)

lmml <- lm(SES.RQEP.BW ~ Lycopods, spp %>% 
             filter(Other != 100))
summary(lmml)

nrow(spp %>% filter(Other == 100)) / nrow(spp) * 100

lmhb <- lm(SES.RQEP.B ~ Other, spp)

summary(lmhb)

lmhg <- lm(SES.RQEP ~ Other, spp)

summary(lmhg)

# remove everything and load the predcitions in the environment
rm(list = ls())

pred.FD.PD <- readRDS("02.data/04.GAM-RQEF-RQEP_predictions.Rds")

pred.FD.PD.sqrt <- readRDS("02.data/04.GAM-RQEF-RQEP.sqrt_predictions.Rds")

#pred.FD.PD.sqrt.for <- readRDS("02.data/04.GAM-SES-RQEF-RQEP.for_predictions.Rds")

#pred.FD.PD.sqrt.nonfor <- readRDS("02.data/04.GAM-SES-RQEF-RQEP.nonfor_predictions.Rds")

pred.FD.PD.SES <- readRDS("02.data/04.GAM-SES-RQEF-RQEP_predictions.Rds")

pred.FD.PD.sqrt.SES <- readRDS("02.data/04.GAM-SES-RQEF-RQEP.sqrt_predictions.Rds")

# load sPlot data
sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-BW-PK-data.RDS") %>% 
  as.data.frame() %>% 
  filter(!is.na(Longitude),
         !is.na(SES.RQEF.BW))

# predict FD from PD based on the sequence between minimum and maximum of PD, 0.1 steps

#reload mods
mod.FD.PD <- readRDS("02.data/04.GAM-RQEF-RQEP.Rds")
mod.FD.PD.sqrt <- readRDS("02.data/04.GAM-RQEF-RQEP.sqrt.Rds")
mod.FD.PD.smooth <- readRDS("02.data/04.GAM-RQEF-RQEP.Smooth.Rds")
mod.FD.PD.SES <- readRDS("02.data/04.GAM-SES-RQEF-RQEP.Rds")
mod.FD.PD.sqrt.SES <- readRDS("02.data/04.GAM-SES-RQEF-RQEP.sqrt.Rds")
mod.MPD <- readRDS("02.data/04.MPD-GAM.Rds")
mod.SES.B <- readRDS("02.data/04.GAM-SES-RQEF-RQEP.bio.Rds")
mod.SES.BW <- readRDS("02.data/04.GAM-SES-RQEF-RQEP.bio.weigh.Rds")
mod.SES.PK <- readRDS("02.data/04.GAM-SES-RQEF-RQEP.PK.Rds")


FD.PD.preds <- predictions(mod.FD.PD, newdata = datagrid(
  RQEP = seq(
    min(sPlot.data$RaoD.phyl), 
    max(sPlot.data$RaoD.phyl), 
    1)
))

FD.PD.sqrt.preds <- predictions(mod.FD.PD.sqrt, newdata = datagrid(
  RQEP.sqrt = seq(
    min(sPlot.data$RQEP.sqrt), 
    max(sPlot.data$RQEP.sqrt), 
    .1)
))

FD.PD.SES.preds <- predictions(mod.FD.PD.SES, newdata = datagrid(
  SES.RQEP = seq(
    min(sPlot.data$SES.RQEP), 
    max(sPlot.data$SES.RQEP), 
    .1)
))

FD.PD.sqrt.SES.preds <- predictions(mod.FD.PD.sqrt.SES, newdata = datagrid(
  SES.RQEP.sqrt = seq(
    min(sPlot.data$SES.RQEP.sqrt), 
    max(sPlot.data$SES.RQEP.sqrt), 
    .1)
))

MPD.SES.preds <- predictions(mod.MPD, newdata = datagrid(
  SES.MPD = seq(
    min(sPlot.data$SES.MPD), 
    max(sPlot.data$SES.MPD), 
    .1)
))

B.SES.preds <- predictions(mod.SES.B, newdata = datagrid(
  SES.RQEP.B = seq(
    min(sPlot.data$SES.RQEP.B), 
    max(sPlot.data$SES.RQEP.B), 
    .1)
))

BW.SES.preds <- predictions(mod.SES.BW, newdata = datagrid(
  SES.RQEP.BW = seq(
    min(sPlot.data$SES.RQEP.BW), 
    max(sPlot.data$SES.RQEP.BW), 
    .1)
))

PK.SES.preds <- predictions(mod.SES.PK, newdata = datagrid(
  SES.RQEP.PK = seq(
    min(sPlot.data$SES.RQEP.PK), 
    max(sPlot.data$SES.RQEP.PK), 
    .1)
))

# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# 
# d.RQE <- data.frame(x = pred.FD.PD$RQEP,
#                     y = pred.FD.PD$RQEF-pred.FD.PD$predicted)
# p.RQE <- ggplot() +
#   geom_point(data = d.RQE, aes(x = x, y = y,
#                                color = atan(range01(d.RQE$y)/range01(d.RQE$x))),
#              alpha = atan(range01(d.RQE$y) + range01(d.RQE$x))) +
#   geom_line(data = FD.PD.preds, aes(x = RaoD.phyl, y = estimate-RQE.MULTI), linewidth = 2) +
#   labs(y = expression(paste("Residuals FD"[Q])), x = expression(paste("PD"[Q]))) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 40),
#     axis.text.x = element_text(size = 36),
#     axis.title.y = element_text(size = 40),
#     axis.text.y = element_text(size = 36),
#     legend.position = "none"
#   ) +
#   #scale_fill_distiller(palette = "Spectral") +
#   scale_color_viridis_c() +
#   #scale_y_continuous(breaks = c(1, 1.5, 2), limits = c(0.9, 2.2)) +
#   geom_text(
#     data = data.frame(
#       xpos = Inf,
#       ypos = Inf,
#       annotateText = "R-sq. (adj) = 0.185\np-value < 0.001",
#       hjustvar = 1.1,
#       vjustvar = 1.1
#     ), 
#     aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText), nudge_x =  -1.1, size = 10)
# 
# 
# d.RQE.sqrt <- data.frame(x = pred.FD.PD.sqrt$RQEP.sqrt,
#                          y = pred.FD.PD.sqrt$RQEF-pred.FD.PD.sqrt$predicted)
# p.RQE.sqrt <- ggplot() +
#   geom_point(data = d.RQE.sqrt, aes(x = x, y = y,
#                                     color = atan(range01(d.RQE.sqrt$y)/range01(d.RQE.sqrt$x))),
#              alpha = atan(range01(d.RQE.sqrt$y) + range01(d.RQE.sqrt$x))) +
#   #geom_ribbon(data = FD.PD.sqrt.preds, aes(x = RQEP.sqrt, ymin = conf.low, ymax = conf.high), fill = "grey" , alpha = 0.8) +
#   geom_line(data = FD.PD.sqrt.preds, aes(x = RQEP.sqrt, y = estimate-RQEF), linewidth = 2) +
#   labs(y = "", x = expression(paste("PD"[Q]," based on square root distances"))) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 40),
#     axis.text.x = element_text(size = 36),
#     axis.title.y = element_text(size = 40),
#     axis.text.y = element_blank(),
#     legend.position = "none"
#   ) +
#   scale_color_viridis_c() +
#   #scale_y_continuous(breaks = c(1, 1.5, 2), limits = c(0.9, 2.2)) +
#   geom_text(
#     data = data.frame(
#       xpos = Inf,
#       ypos = Inf,
#       annotateText = "R-sq. (adj) = 0.362\np-value < 0.001",
#       hjustvar = 1.1,
#       vjustvar = 1.1
#     ), 
#     aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText), nudge_x =  -1.1, size = 10)


res <- residuals.gam(mod.FD.PD.smooth)

SES.RQE.S <- readRDS("02.data/04.GAM-SES-RQEF-smooth.Rds")
summary(SES.RQE.S)

res.S <- residuals.gam(SES.RQE.S)


(p.FPD <- ggplot() +
    geom_hex(aes(y = sPlot.data$RQE.MULTI, 
                 x = sPlot.data$RaoD.phyl),
             bins = 40) +
    #geom_point(aes(y = pred.orig.RQEP$SES.RQEP-pred.orig.RQEP$predicted, x = pred.orig.RQEP$PC1), color = "grey80", alpha = 0.5) +
    geom_line(data = FD.PD.preds, 
              aes(x = RQEP, 
                  y = estimate), 
              linewidth = 2, lty = 1) +
    labs(y = expression(paste("FD"[Q])), 
         x = expression(paste("PD"[Q])) ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40),
      legend.title=element_text(size=36), 
      legend.text=element_text(size=36),
      legend.position = "right",
      #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
      legend.key.height = unit(1.1, "cm"), 
      legend.key.width = unit(1.1, "cm"),
      plot.title = element_text(hjust = .5, size = 70)
    ) +
    # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
    #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
    #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 100000),
                        breaks = c(10, 100, 1000, 10000)) 
)

(p.FPD.sqrt <- ggplot() +
    geom_hex(aes(y = sPlot.data$RQE.MULTI, 
                 x = sPlot.data$RQEP.sqrt),
             bins = 40) +
    #geom_point(aes(y = pred.orig.RQEP$SES.RQEP-pred.orig.RQEP$predicted, x = pred.orig.RQEP$PC1), color = "grey80", alpha = 0.5) +
    geom_line(data = FD.PD.sqrt.preds, 
              aes(x = RQEP.sqrt, 
                  y = estimate), 
              linewidth = 2, lty = 1) +
    labs(y = "", 
         x = expression(paste("PD"[Q]," square root distances")) ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40),
      legend.title=element_text(size=36), 
      legend.text=element_text(size=36),
      legend.position = "right",
      #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
      legend.key.height = unit(1.1, "cm"), 
      legend.key.width = unit(1.1, "cm"),
      plot.title = element_text(hjust = .5, size = 70)
    ) +
    # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
    #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
    #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 100000),
                        breaks = c(10, 100, 1000, 10000)) 
)

(p.SES.sqrt <- ggplot() +
    geom_hex(aes(y = sPlot.data$SES.RQEF, 
                 x = sPlot.data$SES.RQEP.sqrt),
             bins = 40) +
    #geom_point(aes(y = pred.orig.RQEP$SES.RQEP-pred.orig.RQEP$predicted, x = pred.orig.RQEP$PC1), color = "grey80", alpha = 0.5) +
    geom_line(data = FD.PD.sqrt.SES.preds, 
              aes(x = SES.RQEP.sqrt, 
                  y = estimate), 
              linewidth = 2, lty = 1) +
    labs(y = bquote(atop("SES.FD"[Q],
                         "Global species pool")), 
         x = bquote(atop("SES.PD"[Q]~"square root distances",
                         "Global species pool"))  ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40),
      legend.title=element_text(size=36), 
      legend.text=element_text(size=36),
      legend.position = "right",
      #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
      legend.key.height = unit(1.1, "cm"), 
      legend.key.width = unit(1.1, "cm"),
      plot.title = element_text(hjust = .5, size = 70)
    ) +
    # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
    #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
    #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 100000),
                        breaks = c(10, 100, 1000, 10000)) 
)


(p.SES0 <- ggplot() +
  geom_hex(aes(y = sPlot.data$SES.RQEF, 
               x = sPlot.data$SES.RQEP),
           bins = 40) +
  geom_ribbon(data = FD.PD.SES.preds,
              aes(ymin = conf.low, ymax = conf.high,
                  x = SES.RQEP), fill = "grey70", alpha = 0.5) +
  #geom_point(aes(y = pred.orig.RQEP$SES.RQEP-pred.orig.RQEP$predicted, x = pred.orig.RQEP$PC1), color = "grey80", alpha = 0.5) +
  geom_line(data = FD.PD.SES.preds, 
            aes(x = SES.RQEP, 
                y = estimate), 
            linewidth = 2, lty = 1) +
  labs(y = bquote(atop("SES.FD"[Q],
                       "Global species pool")), 
       x = bquote(atop("SES.PD"[Q],
                       "Global species pool")) ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 44),
    axis.text.x = element_text(size = 40),
    axis.title.y = element_text(size = 44),
    axis.text.y = element_text(size = 40),
    legend.title=element_text(size=36), 
    legend.text=element_text(size=36),
    legend.position = "right",
    #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
    legend.key.height = unit(1.1, "cm"), 
    legend.key.width = unit(1.1, "cm"),
    plot.title = element_text(hjust = .5, size = 70)
  ) +
  # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
  #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
  #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
  scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                      limits = c(1, 100000),
                      breaks = c(10, 100, 1000, 10000)) 
)

(p.MPD <- ggplot() +
    geom_hex(aes(y = sPlot.data$SES.FDis, 
                 x = sPlot.data$SES.MPD),
             bins = 40) +
    geom_ribbon(data = MPD.SES.preds,
                aes(ymin = conf.low, ymax = conf.high,
                    x = SES.MPD), fill = "grey70", alpha = 0.5) +
    #geom_point(aes(y = pred.orig.RQEP$SES.RQEP-pred.orig.RQEP$predicted, x = pred.orig.RQEP$PC1), color = "grey80", alpha = 0.5) +
    geom_line(data = MPD.SES.preds, 
              aes(x = SES.MPD, 
                  y = estimate), 
              linewidth = 2, lty = 1) +
    labs(y = bquote(atop("SES.FDis",
                         "Global species pool")), 
         x = bquote(atop("SES.MPD",
                         "Global species pool")), ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40),
      legend.title=element_text(size=36), 
      legend.text=element_text(size=36),
      legend.position = "right",
      #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
      legend.key.height = unit(1.1, "cm"), 
      legend.key.width = unit(1.1, "cm"),
      plot.title = element_text(hjust = .5, size = 70)
    ) +
    # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
    #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
    #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 100000),
                        breaks = c(10, 100, 1000, 10000)) 
)

(p.SES.B <- ggplot() +
    geom_hex(aes(y = sPlot.data$SES.RQEF.B, 
                 x = sPlot.data$SES.RQEP.B),
             bins = 40) +
    geom_ribbon(data = B.SES.preds,
                aes(ymin = conf.low, ymax = conf.high,
                    x = SES.RQEP.B), fill = "grey70", alpha = 0.5) +
    #geom_point(aes(y = pred.orig.RQEP$SES.RQEP-pred.orig.RQEP$predicted, x = pred.orig.RQEP$PC1), color = "grey80", alpha = 0.5) +
    geom_line(data = B.SES.preds, 
              aes(x = SES.RQEP.B, 
                  y = estimate), 
              linewidth = 2, lty = 1) +
    labs(y = bquote(atop("SES.FD"[Q],
                         "Biome species pool")), 
         x = bquote(atop("SES.PD"[Q],
                         "Biome species pool")) ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40),
      legend.title=element_text(size=36), 
      legend.text=element_text(size=36),
      legend.position = "none",
      #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
      legend.key.height = unit(1.1, "cm"), 
      legend.key.width = unit(1.1, "cm"),
      plot.title = element_text(hjust = .5, size = 70)
    ) +
    # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
    #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
    #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 110000),
                        breaks = c(10, 100, 1000, 10000)) 
)

(p.SES.BW <- ggplot() +
    geom_hex(aes(y = sPlot.data$SES.RQEF.BW, 
                 x = sPlot.data$SES.RQEP.BW),
             bins = 40) +
    geom_ribbon(data = BW.SES.preds,
                aes(ymin = conf.low, ymax = conf.high,
                    x = SES.RQEP.BW), fill = "grey70", alpha = 0.5) +
    #geom_point(aes(y = pred.orig.RQEP$SES.RQEP-pred.orig.RQEP$predicted, x = pred.orig.RQEP$PC1), color = "grey80", alpha = 0.5) +
    geom_line(data = BW.SES.preds, 
              aes(x = SES.RQEP.BW, 
                  y = estimate), 
              linewidth = 2, lty = 1) +
    labs(y = bquote(atop("SES.FD"[Q],
                         "Biome weighted species pool")), 
         x = bquote(atop("SES.PD"[Q],
                         "Biome weighted species pool")) ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40),
      legend.title=element_text(size=36), 
      legend.text=element_text(size=36),
      legend.position = "none",
      #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
      legend.key.height = unit(1.1, "cm"), 
      legend.key.width = unit(1.1, "cm"),
      plot.title = element_text(hjust = .5, size = 70)
    ) +
    # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
    #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
    #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 110000),
                        breaks = c(10, 100, 1000, 10000)) 
)

(p.SES.PK <- ggplot() +
    geom_hex(aes(y = sPlot.data$SES.RQEF.PK, 
                 x = sPlot.data$SES.RQEP.PK),
             bins = 40) +
    geom_ribbon(data = PK.SES.preds,
                aes(ymin = conf.low, ymax = conf.high,
                    x = SES.RQEP.PK), fill = "grey70", alpha = 0.5) +
    #geom_point(aes(y = pred.orig.RQEP$SES.RQEP-pred.orig.RQEP$predicted, x = pred.orig.RQEP$PC1), color = "grey80", alpha = 0.5) +
    geom_line(data = BW.SES.preds, 
              aes(x = SES.RQEP.BW, 
                  y = estimate), 
              linewidth = 2, lty = 1) +
    labs(y = bquote(atop("SES.FD"[Q],
                         "Phytogeo species pool")), 
         x = bquote(atop("SES.PD"[Q],
                         "Phytogeo species pool")) ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40),
      legend.title=element_text(size=36), 
      legend.text=element_text(size=36),
      legend.position = "none",
      #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
      legend.key.height = unit(1.1, "cm"), 
      legend.key.width = unit(1.1, "cm"),
      plot.title = element_text(hjust = .5, size = 70)
    ) +
    # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
    #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
    #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 110000),
                        breaks = c(10, 100, 1000, 10000)) 
)

style <- "AABB \n CCDD \n EEFF \n GGHI"

p.out <- p.FPD + p.FPD.sqrt + p.SES0 + p.SES.sqrt + p.MPD + p.SES.B +
  p.SES.PK + guide_area() + plot_layout(guides = "collect", design = style) +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 40))

ggsave("__Submission/Figures/04.GAM.FD-PD.png", p.out, 
       height=35, width=25, units="in", dpi=600, bg = "white", limitsize = FALSE)

ggsave("__Submission/Figures/04.GAM.FD-PD.pdf", p.out, 
       height=35, width=25, units="in", dpi=600, bg = "white", limitsize = FALSE)


(p.SES <- ggplot() +
  geom_hex(aes(y = res.S, 
               x = pred.FD.PD.SES$SES.RQEP),
           bins = 40) +
  #geom_point(aes(y = pred.orig.RQEP$SES.RQEP-pred.orig.RQEP$predicted, x = pred.orig.RQEP$PC1), color = "grey80", alpha = 0.5) +
  geom_line(data = FD.PD.SES.preds, 
            aes(x = SES.RQEP, y = estimate), 
            linewidth = 2, lty = 1) +
  labs(y = expression(paste("Residuals SES.FD"[Q])), 
       x = expression(paste("SES.PD"[Q])) ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 44),
    axis.text.x = element_text(size = 40),
    axis.title.y = element_text(size = 44),
    axis.text.y = element_text(size = 40),
    legend.title=element_text(size=36), 
    legend.text=element_text(size=36),
    legend.position = "right",
    #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
    legend.key.height = unit(1.1, "cm"), 
    legend.key.width = unit(1.1, "cm"),
    plot.title = element_text(hjust = .5, size = 70)
  ) +
  # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
  #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
  #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
  scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                      limits = c(1, 100000),
                      breaks = c(10, 100, 1000, 10000)) 
)


library("patchwork")

p.out <- (p.FPD + p.FPD.sqrt) /
  (p.SES0 + p.SES.sqrt) /
  (p.MPD + p.SES.B) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 40))

ggsave("__Submission/Figures/04.GAM.FD-PD.png", p.out, 
       height=30, width=25, units="in", dpi=600, bg = "white", limitsize = FALSE)




res.S <- residuals.gam(mod.SES.BW.s)

data <- sPlot.data %>%
  filter(!is.na(SES.RQEP.BW),
         !is.na(SES.RQEF.BW))

d.RQE.SES <- data.frame(x = sPlot.data$SES.RQEP,
                        y = sPlot.data$SES.RQEF#-pred.FD.PD.SES$predicted
                        )

n <- 1000000
set.seed(999)
x <- sort(runif(n,min(sPlot.data$SES.RQEF.BW),max(sPlot.data$SES.RQEF.BW)))

err <- rnorm(n, mean=mean(sPlot.data$SES.RQEF.BW), sd=sd(sPlot.data$SES.RQEF.BW))/1.96
hist(err)
quantile(err, probs = c(0.025, 0.975))
y <- x+err

model1 <- lm(y~x)
summary(model1)
conf <- predict(model1,interval="confidence")
pred <- predict(model1,interval="prediction")

data1 <- data.frame(x,y, pred)

ggplot(data = data1) +
  geom_line(aes(x = x, y = upr)) 
  #geom_line(aes(x = x, y = lwr))

library("ggpubr")
library("png")

ggplot(data1, aes(x=x, y=y))+
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70")+
 # geom_point() +
  geom_smooth(method="lm", formula="y~x", fill="lightblue")


bg <- expand.grid(x = c(1:1000),
                  y = c(1:1000)) %>% 
  ggplot(aes(x,y,fill=atan(y/x), alpha=atan(y+x)))+
  geom_raster(alpha = 0.4) +
  scale_fill_viridis() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        legend.position="none") +
  coord_fixed() 

ggsave("999.bg.png", bg, 
       height=15, width=15, units="in", dpi=300, bg = "white")

img <- readPNG("999.bg-crop.png")

inter.low <- mean(data1$lwr-data1$x)
inter.high <- mean(data1$upr-data1$x)

yy <- data %>% 
  filter(SES.RQEF.BW > inter.low + SES.RQEP.BW &
           SES.RQEF.BW < inter.high + SES.RQEP.BW)

hist(yy$SES.RQEP.BW)
quantile(yy$SES.RQEP.BW, probs = c(0.025, 0.975))


p.RQE.SES <- ggplot() +
  # geom_point(data = d.RQE.SES, aes(x = x, y = y,
  #               color = atan(range01(d.RQE.SES$y)/range01(d.RQE.SES$x))),
  #               alpha = atan(range01(d.RQE.SES$y) + range01(d.RQE.SES$x))
  #                ) +
  background_image(img) +
  stat_binhex(data = d.RQE.SES, aes(x = x, y = y, 
                                    fill = ..count..)) + 
  scale_fill_gradient(name = "count", trans = "log",
                      low = alpha("grey20", 0), high = "grey20",
                         #limits = c(1, 100000),
                         breaks = c(10, 100, 1000, 10000)) +
  geom_ribbon(data = data1, aes(ymin = lwr, ymax = upr, x = x), fill = "grey70", alpha = 0.3) +
  geom_smooth(data = data1, aes(x = x, y = y), method="lm", formula="y~x", color = "black") + 
  geom_line(data = BW.SES.preds, aes(x = SES.RQEP.BW, y = estimate), color = "blue", linewidth = 2) +
  # geom_line(aes(x = c(quantile(yy$SES.RQEP.BW, probs = c(0.025, 0.975))[1],
  #                     quantile(yy$SES.RQEP.BW, probs = c(0.025, 0.975))[1]), 
  #               y = c(-Inf,
  #                     inter.high + quantile(yy$SES.RQEP.BW, probs = c(0.025, 0.975))[1])),
  #           linetype = 2,
  #           linewidth = 1.5) + 
  # geom_line(aes(x = c(quantile(yy$SES.RQEP.BW, probs = c(0.025, 0.975))[2],
  #                     quantile(yy$SES.RQEP.BW, probs = c(0.025, 0.975))[2]), 
  #               y = c(-Inf,
  #                     inter.high + quantile(yy$SES.RQEP.BW, probs = c(0.025, 0.975))[2])),
  #           linetype = 2,
#           linewidth = 1.5) + 
  labs(y = expression(paste("SES.FD"[Q])), x = expression(paste("SES.PD"[Q]))) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(size = 36),
    axis.title.y = element_text(size = 40),
    axis.text.y = element_text(size = 36),
    legend.position=c(.9,.65),
    legend.key.height = unit(1.1, "cm"), 
    legend.key.width = unit(0.5, "cm"),
    legend.title=element_text(size=36), 
    legend.text=element_text(size=36),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) + 
  scale_color_viridis_c()
  # geom_text(
  #   data = data.frame(
  #     xpos = Inf,
  #     ypos = Inf,
  #     annotateText = 
  #     "R-sq. (adj) = 0.058\np-value < 0.001",
  #     hjustvar = 1.1,
  #     vjustvar = 1.1
  #   ), 
  #   aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText), nudge_x =  -1.1, size = 10)


c <- sPlot.data %>% 
  filter(SES.RQEF.BW > inter.low + SES.RQEP.BW &
           SES.RQEF.BW < inter.high + SES.RQEP.BW) %>% 
  pull(PlotObservationID)
  #nrow() / nrow(sPlot.data)*100

d.FD <- sPlot.data %>% 
  filter(SES.RQEF.BW > inter.high + SES.RQEP.BW) %>% 
  pull(PlotObservationID)
  #nrow() / nrow(sPlot.data)*100

d.PD <- sPlot.data %>% 
  filter(SES.RQEF.BW < inter.low + SES.RQEP.BW) %>% 
  pull(PlotObservationID)
  #nrow() / nrow(sPlot.data)*100

data2 <- sPlot.data %>% 
  mutate(Status = ifelse(PlotObservationID %in% c, "Coupling",
                         ifelse(PlotObservationID %in% d.FD, "Decoupling higher FD",
                                ifelse(PlotObservationID %in% d.PD, "Decoupling higher PD", NA))))

saveRDS(data2, "02.data/04.sPlot.PD.FD.CD-PCA-BW-Status-data.RDS")

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

log.th <- sPlot.data %>% 
    dplyr::select(SES.RQEF.BW, SES.RQEP.BW) %>% 
    mutate(x = range01(SES.RQEP.BW),
           y = range01(SES.RQEF.BW) ) %>% 
    filter(SES.RQEF.BW > inter.low + SES.RQEP.BW &
           SES.RQEF.BW < inter.high + SES.RQEP.BW) %>% 
    mutate(x = ifelse(x == 0, 0.001, x),
           y = ifelse(y == 0, 0.001, y)) %>%
    mutate(rel.BW = y/x) %>%
    mutate(rel.log.BW = log(rel.BW))

min(log.th$rel.log.BW)
max(log.th$rel.log.BW)

library("ggdendro")
library("viridis")
library("patchwork")

(top.left <- mtcars %>% 
    as.data.frame() %>% 
    mutate(car = rownames(mtcars)) %>% 
    filter(car == "Merc 280C" | car == "Mazda RX4" |
             car == "Merc 240D" | car == "Merc 230") %>% 
    dplyr::select(-car) %>% 
    dist() %>% 
    hclust() %>% 
    as.dendrogram() %>% 
    dendro_data(., type = "rectangle") %>%
    segment() %>% 
    ggplot() +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1.5) + 
    geom_segment(aes(x = 2.5, y = 150, xend = 2.5, yend = 65), linewidth = 1.5) +
    #geom_point(aes(x = 2, y = -4), shape = 17, colour = "red", size = 10) +
    #geom_point(aes(x = 1, y = -4), shape = 19, colour = "darkgreen", size = 10) +
    labs(title = "Decoupling with FD > PD",
         subtitle = paste0(round(d.FD, 2), "% Observations")) +
    coord_flip() + 
    scale_y_reverse(expand = c(0, 0),
                    limits = c(150, -150/4)) +
    scale_x_continuous(limits = c(0.9,4.1)) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 32),
          plot.subtitle = element_text(hjust = 0.5, size = 30),
          #panel.border = element_rect(colour = "darkgrey", fill = NA, linewidth = 2),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)

(mid <- mtcars %>% 
    as.data.frame() %>% 
    mutate(car = rownames(mtcars)) %>% 
    filter(car == "Merc 280C" | car == "Mazda RX4" |
             car == "Merc 240D" | car == "Merc 230") %>% 
    dplyr::select(-car) %>% 
    dist() %>% 
    hclust() %>% 
    as.dendrogram() %>% 
    dendro_data(., type = "rectangle") %>%
    segment() %>% 
    ggplot() +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1.5) + 
    geom_segment(aes(x = 2.5, y = 70, xend = 2.5, yend = 65), linewidth = 1.5) +
    #geom_point(aes(x = 2, y = -4), shape = 17, colour = "red", size = 10) +
    #geom_point(aes(x = 1, y = -4), shape = 19, colour = "darkgreen", size = 10) +
    labs(title = expression(paste("Coupling with FD" %~~% "PD")),
         subtitle = paste0(round(c, 2), "% Observations")) +
    coord_flip() + 
    scale_y_reverse(expand = c(0, 0),
                    limits = c(70, -70/4)) +
    scale_x_continuous(limits = c(0.9,4.1)) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 32),
          plot.subtitle = element_text(hjust = 0.5, size = 30),
          #panel.border = element_rect(colour = "darkgrey", fill = NA, linewidth = 2),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)

(bottom.right <- mtcars %>% 
    as.data.frame() %>% 
    mutate(car = rownames(mtcars)) %>% 
    filter(car == "Merc 280C" | car == "Mazda RX4" |
             car == "Merc 240D" | car == "Merc 230") %>% 
    dplyr::select(-car) %>% 
    dist() %>% 
    hclust() %>% 
    as.dendrogram() %>% 
    dendro_data(., type = "rectangle") %>%
    segment() %>% 
    ggplot() +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1.5) + 
    geom_segment(aes(x = 2.5, y = 70, xend = 2.5, yend = 65), linewidth = 1.5) +
    #geom_point(aes(x = 2, y = -4), shape = 17, colour = "red", size = 10) +
    #geom_point(aes(x = 1, y = -4), shape = 19, colour = "darkgreen", size = 10) +
    labs(title = "Decoupling with FD < PD",
         subtitle = paste0(round(d.PD, 2), "% Observations")
    ) +
    coord_flip() + 
    scale_y_reverse(expand = c(0, 0),
                    limits = c(70, -70/4)) +
    scale_x_continuous(limits = c(0.9,4.1)) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 32),
          plot.subtitle = element_text(hjust = 0.5, size = 30),
          #panel.border = element_rect(colour = "darkgrey", fill = NA, linewidth = 2),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)

out <- p.RQE.SES + 
  inset_element(top.left, left = 0.05, bottom = 0.75, right = 0.29, top = 0.99) + #top left
  inset_element(mid,  left = 0.71, bottom = 0.75, right = 0.95, top = 0.99) + #top right
  inset_element(bottom.right, left = 0.71, bottom = 0.01, right = 0.95, top = 0.25)  #bottom right 

ggsave("__Submission/Figures/04.GAM.FD-PD.BG.png", out, 
       height=15, width=20, units="in", dpi=600, bg = "white")

save.image("03.results/_temp.04.FD-PD-GAM-SES.RData")


library("patchwork")
library("dggridR")
library("rnaturalearth")


dggs <- dgconstruct(area=50*50, metric=T, resround='down')

pred.FD.PD.SES$cell <- dgGEO_to_SEQNUM(dggs, 
                                       pred.FD.PD.SES$Longitude, 
                                       pred.FD.PD.SES$Latitude)$seqnum
pred.FD.PD.SES$SES <- res.S

d.phyl <- pred.FD.PD.SES %>% 
  group_by(cell) %>% 
  summarise(smooth = mean(SES, na.rm = T))

grid.SES <- dgcellstogrid(dggs, d.phyl$cell) %>%
  st_as_sf() %>% 
  mutate(cell = d.phyl$cell) %>% 
  mutate(smooth = d.phyl$smooth) %>% 
  st_transform("+proj=eck4") %>% 
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))

world <- ne_countries(returnclass = "sf")%>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()

(p.phyl.geo <- grid.SES %>% 
    .[-which.max(st_area(grid.SES)),] %>% 
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    # geom_sf(data = bb, col = "grey20", fill = NA) +
    coord_sf(crs = "+proj=eck4") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=36),
          legend.text=element_text(size=36),
          #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          legend.position = "right", 
          plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    geom_sf(aes(color = smooth, fill = smooth), lwd = 0) +
    scale_fill_viridis_c(option = "plasma",
                         breaks = c(-4,0,4)) +
    scale_color_viridis_c(option = "plasma",
                          breaks = c(-4,0,4))
    # scale_fill_gradient2(low = "#91bfdb",
    #                      mid = "#ffffbf",
    #                      high = "#fc8d59",
    #                      # limits = c(-5,5),
    #                      breaks = c(-4,0,4)
    # )+
    # # midpoint = 0) +
    # scale_color_gradient2(low = "#91bfdb",
    #                       mid = "#ffffbf",
    #                       high = "#fc8d59",
    #                       breaks = c(-4,0,4)
    # )
  # guides(color = guide_colourbar(title.position="top", title.hjust = 0.5),
  #        size = guide_legend(title.position="top", title.hjust = 0.5))
)

library("patchwork")

layout <- "AABBC"

p.out.SES <- p.SES + p.phyl.geo + guide_area() +
  plot_layout(design = layout, heights = c(1,1), guides = "collect") +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 40))

ggsave("__Submission/Figures/04.S-GAM_FD-PD.RES.png", p.out.SES, 
       height=10, width=25, units="in", dpi=600, bg = "white")

dggs <- dgconstruct(area=50*50, metric=T, resround='down')

pred.FD.PD$cell <- dgGEO_to_SEQNUM(dggs, 
                                       pred.FD.PD$Longitude, 
                                       pred.FD.PD$Latitude)$seqnum
pred.FD.PD$S <- res

d.phyl <- pred.FD.PD %>% 
  group_by(cell) %>% 
  summarise(smooth = mean(S, na.rm = T))

grid <- dgcellstogrid(dggs, d.phyl$cell) %>%
  st_as_sf() %>% 
  mutate(cell = d.phyl$cell) %>% 
  mutate(smooth = d.phyl$smooth) %>% 
  st_transform("+proj=eck4") %>% 
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))

world <- ne_countries(returnclass = "sf")%>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()

(p.phyl.geo <- grid %>% 
    .[-which.max(st_area(grid)),] %>% 
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    # geom_sf(data = bb, col = "grey20", fill = NA) +
    coord_sf(crs = "+proj=eck4") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=36),
          legend.text=element_text(size=36),
          #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          legend.position = "right", 
          plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    geom_sf(aes(color = smooth, fill = smooth), lwd = 0) +
    scale_fill_viridis_c(option = "plasma",
                         limits = c(-0.5,0.5),
                         breaks = c(-0.4,0,0.4)) +
    scale_color_viridis_c(option = "plasma",
                          limits = c(-0.5,0.5),
                          breaks = c(-0.4,0,0.4)) 
    # scale_fill_gradient2(low = "#91bfdb",
    #                      mid = "#ffffbf",
    #                      high = "#fc8d59",
    #                      limits = c(-0.5,0.5),
    #                      breaks = c(-0.4,0,0.4)
    # )+
    # # midpoint = 0) +
    # scale_color_gradient2(low = "#91bfdb",
    #                       mid = "#ffffbf",
    #                       high = "#fc8d59",
    #                       limits = c(-0.5,0.5),
    #                       breaks = c(-0.4,0,0.4)
    # )
  # guides(color = guide_colourbar(title.position="top", title.hjust = 0.5),
  #        size = guide_legend(title.position="top", title.hjust = 0.5))
)

# pred.FD.PD.sqrt$cell <- dgGEO_to_SEQNUM(dggs, pred.FD.PD.sqrt$Longitude, pred.FD.PD.sqrt$Latitude)$seqnum
# 
# d.phyl.sqrt <- pred.FD.PD.sqrt %>% 
#   group_by(cell) %>% 
#   summarise(smooth = mean(predicted, na.rm = T))
# 
# grid.sqrt <- dgcellstogrid(dggs, d.phyl.sqrt$cell) %>%
#   st_as_sf() %>% 
#   mutate(cell = d.phyl.sqrt$cell) %>% 
#   mutate(smooth = d.phyl.sqrt$smooth) %>% 
#   st_transform("+proj=eck4") %>% 
#   st_wrap_dateline(options = c("WRAPDATELINE=YES"))
# 
# (p.phyl.sqrt.geo <- grid.sqrt %>% 
#     .[-which.max(st_area(grid.sqrt)),] %>% 
#     ggplot() +
#     geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
#     #geom_sf(data = bb, col = "grey20", fill = NA) +
#     coord_sf(crs = "+proj=eck4") +
#     theme_minimal() +
#     theme(axis.text = element_blank(),
#           axis.title = element_blank(),
#           legend.title=element_text(size=36),
#           legend.text=element_text(size=36),
#           legend.background = element_rect(linewidth=0.1, linetype="solid", colour = 1),
#           legend.key.height = unit(1.1, "cm"),
#           legend.key.width = unit(1.1, "cm"),
#           legend.position = "right") +
#     geom_sf(aes(color = smooth), lwd = 1) +
#     scale_color_gradient2(low = "#91bfdb",
#                          mid = "#ffffbf",
#                          high = "#fc8d59",
#                          limits = c(1,2),
#                          breaks = c(1.1,1.5,1.9),
#                          midpoint = 1.5)
#     # scale_fill_gradient2(low = "#91bfdb",
#     #                       mid = "#ffffbf",
#     #                       high = "#fc8d59",
#     #                       limits = c(1,2),
#     #                       breaks = c(1.1,1.5,1.9),
#     #                       midpoint = 1.5)
#     # guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5),
#     #        size = guide_legend(title.position="top", title.hjust = 0.5))
# )
layout <- "AABB
           CCDD"
# 
# p.out <- p.FPD +
#                  
#   p.FPD.sqrt + 
#                 
#   p.SES.sqrt + 
#                 
#   p.phyl.geo + 
#                 
#   #guide_area() +
#   plot_layout(design = layout, heights = c(0.8, 0.8), guides = "collect") +
#   plot_annotation(tag_levels = "A")  & 
#   theme(plot.tag = element_text(size = 40))
# 
# 
# ggsave("__Submission/Figures/S04.FD-PD-GAM_sqrt.png", p.out,
#        height=25, width=30, units="in", dpi=600, bg = "white")

rm(list = ls())
### smooths
sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-BW-PK-data.RDS")
sPlot.data <- sPlot.data %>% as.data.frame() %>% filter(!is.na(Longitude))

mod.FD.S <- readRDS("02.data/04.GAM-RQEF-RQEP.Smooth.Rds")
summary(mod.FD.S)

mod.SES.glob <- readRDS("02.data/04.GAM-SES-RQEF-smooth.Rds")
summary(mod.SES.glob)

mod.SES.bio <- readRDS("02.data/04.GAM-SES-RQEF-RQEP.bio.smooth.Rds")
summary(mod.SES.bio)

mod.SES.biowe <- readRDS("02.data/04.GAM-SES-RQEF.bio.weigh-smooth.Rds")
summary(mod.SES.biowe)

mod.SES.pk <- readRDS("02.data/04.GAM-SES-RQEF.PK-smooth.Rds")
summary(mod.SES.pk)

sPlot.data$res.FD <- residuals.gam(mod.FD.S, type = "response")
sPlot.data$res.SES.glob <- residuals.gam(mod.SES.glob, type = "response")

sPlot.data <- sPlot.data %>% 
  filter(!is.na(SES.RQEF.BW),
         !is.na(SES.RQEP.BW))

sPlot.data$res.SES.bio <- residuals.gam(mod.SES.bio, type = "response")
sPlot.data$res.SES.biowe <- residuals.gam(mod.SES.biowe, type = "response")
sPlot.data$res.SES.pk <- residuals.gam(mod.SES.pk, type = "response")

d.RQEF <- sPlot.data %>% 
  group_by(cell) %>% 
  summarise(s.FD = mean(res.FD, na.rm = T),
            s.Glob = mean(res.SES.glob, na.rm = T),
            s.Bio = mean(res.SES.bio, na.rm = T),
            s.Biowe = mean(res.SES.biowe, na.rm = T),
            s.PK = mean(res.SES.pk, na.rm = T)
            )

dggs <- dgconstruct(area=50*50, metric=T, resround='down')

grid.RQEF <- dgcellstogrid(dggs, d.RQEF$cell) %>%
  st_as_sf() %>% 
  mutate(cell = d.RQEF$cell) %>% 
  mutate(s.FD = d.RQEF$s.FD,
         s.Glob = d.RQEF$s.Glob,
         s.Bio = d.RQEF$s.Bio,
         s.Biowe = d.RQEF$s.Biowe,
         s.PK = d.RQEF$s.PK) %>%
  st_transform("+proj=eck4") %>% 
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))

world <- ne_countries(returnclass = "sf")%>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()


(FD.geo <- grid.RQEF %>% 
    mutate(area = as.numeric(st_area(.))) %>%
    filter(area <= quantile(area, c(0.999))) %>%  
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    # geom_sf(data = bb, col = "grey20", fill = NA) +
    coord_sf(crs = "+proj=eck4") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=36),
          legend.text=element_text(size=36),
          #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          legend.position = "bottom", 
          plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    geom_sf(aes(color = s.FD, fill = s.FD), lwd = 0) +
    scale_fill_viridis_c(option = "plasma",
                         breaks = c(-0.3,0,0.3),
                         labels = c("0.3", "0", "0.3"),
                         name = expression(paste("FD"[Q]))) +
    scale_color_viridis_c(option = "plasma",
                          breaks = c(-0.3,0,0.3),
                          labels = c("0.3", "0", "0.3"),
                          name = expression(paste("FD"[Q]))) +
    guides(
      fill = "none",
      color = guide_colourbar(
        #title = expression(paste("Higher\nSES.PD"[Q], symbol("\256"), "Higher\nSES.FD"[Q])),
        title.position="top", title.hjust = 0.5
      )) )

(glob.geo <- grid.RQEF %>% 
    mutate(area = as.numeric(st_area(.))) %>%
    filter(area <= quantile(area, c(0.999))) %>%  
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    # geom_sf(data = bb, col = "grey20", fill = NA) +
    coord_sf(crs = "+proj=eck4") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=36),
          legend.text=element_text(size=36),
          #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          legend.position = "bottom", 
          plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    geom_sf(aes(color = s.Glob, fill = s.Glob), lwd = 0) +
    scale_fill_viridis_c(option = "plasma",
                         breaks = c(-4,0,4),
                         name = expression(paste("SES.FD"[Q], "(Global)"))) +
    scale_color_viridis_c(option = "plasma",
                          breaks = c(-4,0,4),
                          name = expression(paste("SES.FD"[Q], "(Global)"))) +
    guides(
      fill = "none",
      color = guide_colourbar(
        #title = expression(paste("Higher\nSES.PD"[Q], symbol("\256"), "Higher\nSES.FD"[Q])),
        title.position="top", title.hjust = 0.5
      )) )

(bio.geo <- grid.RQEF %>% 
    mutate(area = as.numeric(st_area(.))) %>%
    filter(area <= quantile(area, c(0.999))) %>%  
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    # geom_sf(data = bb, col = "grey20", fill = NA) +
    coord_sf(crs = "+proj=eck4") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=36),
          legend.text=element_text(size=36),
          #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          legend.position = "bottom", 
          plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    geom_sf(aes(color = s.Bio, fill = s.Bio), lwd = 0) +
    scale_fill_viridis_c(option = "plasma",
                         breaks = c(-5, 0, 5),
                         name = expression(paste("SES.FD"[Q], "(Biome)"))) +
    scale_color_viridis_c(option = "plasma",
                          breaks = c(-5, 0, 5),
                          name = expression(paste("SES.FD"[Q], "(Biome)")))+
    guides(
      fill = "none",
      color = guide_colourbar(
        #title = expression(paste("Higher\nSES.PD"[Q], symbol("\256"), "Higher\nSES.FD"[Q])),
        title.position="top", title.hjust = 0.5
      )) )

(biowe.geo <- grid.RQEF %>% 
    mutate(area = as.numeric(st_area(.))) %>%
    filter(area <= quantile(area, c(0.999))) %>%  
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    # geom_sf(data = bb, col = "grey20", fill = NA) +
    coord_sf(crs = "+proj=eck4") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=36),
          legend.text=element_text(size=36),
          #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          legend.position = "bottom", 
          plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    geom_sf(aes(color = s.Biowe, fill = s.Biowe), lwd = 0) +
    scale_fill_viridis_c(option = "plasma",
                         breaks = c(-5,0,5),
                         name = expression(paste("SES.FD"[Q], "(Biome, weighted)"))) +
    scale_color_viridis_c(option = "plasma",
                          breaks = c(-5,0,5),
                          name = expression(paste("SES.FD"[Q], "(Biome, weighted)")))+
    guides(
      fill = "none",
      color = guide_colourbar(
        #title = expression(paste("Higher\nSES.PD"[Q], symbol("\256"), "Higher\nSES.FD"[Q])),
        title.position="top", title.hjust = 0.5
      )) )

(pk.geo <- grid.RQEF %>% 
    mutate(area = as.numeric(st_area(.))) %>%
    filter(area <= quantile(area, c(0.999))) %>%  
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    # geom_sf(data = bb, col = "grey20", fill = NA) +
    coord_sf(crs = "+proj=eck4") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=36),
          legend.text=element_text(size=36),
          #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          legend.position = "bottom", 
          plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    geom_sf(aes(color = s.PK, fill = s.PK), lwd = 0) +
    scale_fill_viridis_c(option = "plasma",
                         breaks = c(-5,0,5),
                         name = expression(paste("SES.FD"[Q], "(Phytogeo)"))) +
    scale_color_viridis_c(option = "plasma",
                          breaks = c(-5,0,5),
                          name = expression(paste("SES.FD"[Q], "(Phytogeo)")))+
    guides(
      fill = "none",
      color = guide_colourbar(
        #title = expression(paste("Higher\nSES.PD"[Q], symbol("\256"), "Higher\nSES.FD"[Q])),
        title.position="top", title.hjust = 0.5
      )) )

library("patchwork")

layout = "AABB
          CCDD
          EEFF"

p.s <- FD.geo + glob.geo +
  bio.geo + biowe.geo +
  pk.geo +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 40))

ggsave("__Submission/Figures/04.smooths.png", p.s, 
       height=12, width=15, units="in", dpi=600, bg = "white", limitsize = FALSE)

ggsave("__Submission/Figures/04.smooths.pdf", p.s, 
       height=12, width=15, units="in", dpi=600, bg = "white", limitsize = FALSE)

#sPlot Open

#Presence-Absence
models <- readRDS("02.data/sPlotOpen-PA-GAMS.Rds")
sPlot.oa <- readRDS("02.data/03.01.sPlotOpen-header-FD-PD.RDS")

RQE <- c("FDQ.ALL.PA",
         "FDQ.MULTI.PA",
         "FDQ.SLA.PA",
         "FDQ.HEIGHT.PA",
         "FDQ.ROOT.PA",
         "FDQ.SM.PA",
         "FDQ.LDMC.PA",
         "FDQ.N.PA",
         "FDQ.P.PA",
         "FDQ.CHRO.PA")

preds <- list()
  
for(i in 1:length(models)) {
  preds[[i]] <- marginaleffects::predictions(models[[i]], newdata = datagrid(
    PDQ.PA = seq(
      min(sPlot.oa$PDQ.PA, na.rm = TRUE), 
      max(sPlot.oa$PDQ.PA, na.rm = TRUE), 
      10)
  ))
}

saveRDS(preds, "02.data/04.oa.preds-PA-GAM.RDS")
preds <- readRDS("02.data/04.oa.preds-PA-GAM.RDS")

size <-  data.frame(Response = c("Chromosome n",
                                 "Plant height",
                                 "LDMC",
                                 "Leaf N",
                                 "Leaf P",
                                 "SRL",
                                 "SLA",
                                 "Seed mass",
                                 "All traits"),
                       Width = c(rep(3, times = 8), 3.5) )

df.preds <- bind_rows(preds) %>% 
  mutate(Response = rep(RQE, each = 31)) %>% 
  mutate(
    Response = gsub("FDQ.CHRO.PA", "Chromosome n", Response),
    Response = gsub("FDQ.HEIGHT.PA", "Plant height", Response),
    Response = gsub("FDQ.LDMC.PA", "LDMC", Response),
    Response = gsub("FDQ.N.PA", "Leaf N", Response),
    Response = gsub("FDQ.P.PA", "Leaf P", Response),
    Response = gsub("FDQ.ROOT.PA", "SRL", Response),
    Response = gsub("FDQ.SLA.PA", "SLA", Response),
    Response = gsub("FDQ.SM.PA", "Seed mass", Response),
    Response = gsub("FDQ.ALL.PA", "All traits", Response)
  ) %>% 
  filter(!Response %in% c("FDQ.MULTI.PA","FDQ.SUB.PA","Growthform")) %>% 
  left_join(size)


cols <- c("Chromosome n" = "mediumvioletred",
          "Growthform" = "seagreen1",
          "Plant height" = "forestgreen",
          "LDMC" = "deepskyblue1",
          "Leaf N" = "royalblue",
          "Leaf P" = "darkgoldenrod1",
          "SRL" = "sienna",
          "SLA" = "firebrick1",
          "Seed mass" = "mediumorchid1",
          "All traits" = "black"
          )

linewidth <- c("Chromosome n" = 2,
          "Plant height" = 2,
          "LDMC" = 2,
          "Leaf N" = 2,
          "Leaf P" = 2,
          "SRL" = 2,
          "SLA" = 2,
          "Seed mass" = 2,
          "All traits" = 4
)



(p.oa <- ggplot() +
  #geom_point(aes(y = sPlot.oa$RQE.ALL, x = sPlot.oa$RaoD.phyl), color = "grey80", alpha = 0.5) +
  geom_line(data = df.preds, aes(x = PDQ.PA, y = estimate, color = Response, size = Response)) +
  labs(y = expression("FD"[Q]), 
       x = expression("PD"[Q]),
       colour = expression(paste("FD"[Q], " based on")),
       title = "Relationship of FD and PD based on presence-absence data") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(size = 36),
    axis.title.y = element_text(size = 40),
    axis.text.y = element_text(size = 36),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 40),
    plot.title = element_text(size = 46, face = "bold", hjust = 0.5)
  ) +
  scale_color_manual(values = cols) +
  scale_size_manual(values = linewidth) +
  scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
)

ggsave("03.results/04.FD-PD-PA.oa-GAM.png", p.oa,
       height=15, width=25, units="in", dpi=300, bg = "white")



# #Abundance
# models<- readRDS("02.data/sPlotOpen-ABU-GAMS.Rds")
# 
# RQE <- c("FDQ.ALL.ABU",
#          "FDQ.SUB.ABU",
#          "FDQ.MULTI.ABU",
#          "FDQ.SLA.ABU",
#          "FDQ.HEIGHT.ABU",
#          "FDQ.ROOT.ABU",
#          "FDQ.SM.ABU",
#          "FDQ.LDMC.ABU",
#          "FDQ.GF.ABU",
#          "FDQ.N.ABU",
#          "FDQ.P.ABU",
#          "FDQ.CHRO.ABU")
# 
# preds.abu <- list()
# 
# for(i in 1:length(models)) {
#   preds.abu[[i]] <- marginaleffects::predictions(models[[i]], newdata = datagrid(
#     PDQ.ABU = seq(
#       min(sPlot.oa$PDQ.ABU, na.rm = TRUE), 
#       max(sPlot.oa$PDQ.ABU, na.rm = TRUE), 
#       10)
#   ))
# }
# 
# saveRDS(preds.abu, "02.data/04.oa.preds-ABU-GAM.RDS")
# preds.abu <- readRDS("02.data/04.oa.preds-ABU-GAM.RDS")
# 
# df.preds.abu <- bind_rows(preds.abu) %>% 
#   mutate(Response = rep(RQE, each = 59)) %>% 
#   mutate(
#     Response = gsub("FDQ.CHRO.ABU", "Chromosome n", Response),
#     Response = gsub("FDQ.GF.ABU", "Growthform", Response),
#     Response = gsub("FDQ.HEIGHT.ABU", "Plant height", Response),
#     Response = gsub("FDQ.LDMC.ABU", "LDMC", Response),
#     Response = gsub("FDQ.N.ABU", "Leaf N", Response),
#     Response = gsub("FDQ.P.ABU", "Leaf P", Response),
#     Response = gsub("FDQ.ROOT.ABU", "SRL", Response),
#     Response = gsub("FDQ.SLA.ABU", "SLA", Response),
#     Response = gsub("FDQ.SM.ABU", "Seed mass", Response),
#     Response = gsub("FDQ.ALL.ABU", "All traits", Response)
#   ) %>% 
#   filter(!Response %in% c("FDQ.MULTI.ABU","FDQ.SUB.ABU"))
# 
# 
# (p.abu <- ggplot() +
#     #geom_point(aes(y = sPlot.oa$RQE.ALL, x = sPlot.oa$RaoD.phyl), color = "grey80", alpha = 0.5) +
#     geom_line(data = df.preds.abu, aes(x = PDQ.ABU, y = estimate, color = Response), size = 2) +
#     labs(y = expression("FD"[Q]), 
#          x = expression("PD"[Q]),
#          colour = expression(paste("FD"[Q], " based on")),
#          title = "Relationship of FD and PD based on abundance data") +
#     theme_bw() +
#     theme(
#       axis.title.x = element_text(size = 40),
#       axis.text.x = element_text(size = 36),
#       axis.title.y = element_text(size = 40),
#       axis.text.y = element_text(size = 36),
#       legend.text = element_text(size = 36),
#       legend.title = element_text(size = 40),
#       plot.title = element_text(size = 46, face = "bold", hjust = 0.5)
#     ) +
#     scale_color_manual(values = cols) +
#     scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
# )
# 
# ggsave("03.results/04.FD-PD-ABU.oa-GAM.png", p.abu,
#        height=15, width=25, units="in", dpi=300, bg = "white")
#deprecated:
# p.hist.in <- ggplot() +
#   geom_point(aes(x = pred.FD.PD.sqrt$RQEP.sqrt, y = pred.FD.PD.sqrt$RQEF), color = "grey", alpha = 0.3) +
#   geom_ribbon(data = FD.PD.sqrt.preds, aes(x = RQEP.sqrt, ymin = conf.low, ymax = conf.high), fill = "grey" , alpha = 0.8) +
#   geom_line(data = FD.PD.sqrt.preds, aes(x = RQEP.sqrt, y = predicted)) +
#   labs(y = "Residuals Functional Entropy", x = "Phylogenetic Entropy based on \nsquare root distances") +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 40),
#     axis.text.x = element_text(size = 36),
#     axis.title.y = element_text(size = 40),
#     axis.text.y = element_text(size = 36)
#   ) +
#   geom_text(
#     data = data.frame(
#       xpos = Inf,
#       ypos = Inf,
#       annotateText = "R-sq. (adj) = 0.362\np-value < 0.001",
#       hjustvar = 1.1,
#       vjustvar = 1.1
#     ), 
#     aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText), nudge_x =  -1.1, size = 10)
# 
# 
# p.hist <- ggExtra::ggMarginal(p.hist.in, type = "histogram")
# 
# ggsave("03.results/04.FD-PD-GAM-hist.png", p.hist,
#        height=20, width=20, units="in", dpi=600, bg = "white")
# 
# out <- sPlot.data %>% 
#   filter(RQE.MULTI == 2, RaoD.phyl < 100)
# 
# out2 <- sPlot.data %>% 
#   filter(between(RQE.MULTI, 1.799999999999999, 1.800000000001))
# 
# out3 <- sPlot.data %>% 
#   filter(SES.RQEF < -5)
# 
# out4 <- sPlot.data %>% 
#   filter(between(SES.RQEF, 5.5, 7), SES.RQEP < 0)
# 
# #plot for and nonfor
# 
# p.hist.in.for <- ggplot() +
#   geom_point(aes(x = pred.FD.PD.sqrt.for$RQEP.sqrt, y = pred.FD.PD.sqrt.for$RQEF), color = "grey", alpha = 0.3) +
#   geom_ribbon(data = FD.PD.sqrt.for.preds, aes(x = RQEP.sqrt, ymin = conf.low, ymax = conf.high), fill = "grey" , alpha = 0.8) +
#   geom_line(data = FD.PD.sqrt.for.preds, aes(x = RQEP.sqrt, y = predicted)) +
#   labs(y = "Residuals Functional Entropy", x = "Phylogenetic Entropy based on \nsquare root distances"
#        #title = "FOREST"
#        ) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 40),
#     axis.text.x = element_text(size = 36),
#     axis.title.y = element_text(size = 40),
#     axis.text.y = element_text(size = 36)
#     #plot.title = element_text(size = 44, colour = "darkgreen", vjust = 0.2)
#   ) +
#   scale_y_continuous(breaks = c(1.25, 1.50, 1.75, 2.00), limits = c(1, 2.1)) +
#   scale_x_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 17)) +
#   geom_text(
#     data = data.frame(
#       xpos = Inf,
#       ypos = Inf,
#       annotateText = "R-sq. (adj) = 0.294\np-value < 0.001",
#       hjustvar = 1.1,
#       vjustvar = 1.1
#     ), 
#     aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText), nudge_x =  -1.1, size = 10)+
#   geom_text(
#     data = data.frame(
#       xpos = -Inf,
#       ypos = -Inf,
#       annotateText = "FOREST",
#       hjustvar = 0,
#       vjustvar = -0.1
#     ), 
#     aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText), size = 20, color = "darkgreen")
# 
# 
# p.hist.for <- ggExtra::ggMarginal(p.hist.in.for, type = "histogram")
# 
# p.hist.in.nonfor <- ggplot() +
#   geom_point(aes(x = pred.FD.PD.sqrt.nonfor$RQEP.sqrt, y = pred.FD.PD.sqrt.nonfor$RQEF), color = "grey", alpha = 0.3) +
#   geom_ribbon(data = FD.PD.sqrt.nonfor.preds, aes(x = RQEP.sqrt, ymin = conf.low, ymax = conf.high), fill = "grey" , alpha = 0.8) +
#   geom_line(data = FD.PD.sqrt.nonfor.preds, aes(x = RQEP.sqrt, y = predicted)) +
#   labs(y = "Residuals Functional Entropy", x = "Phylogenetic Entropy based on \nsquare root distances"
#        #title = "NON-FOREST"
#        ) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 40),
#     axis.text.x = element_text(size = 36),
#     axis.title.y = element_text(size = 40),
#     axis.text.y = element_text(size = 36)
#     #plot.title = element_text(size = 44, colour = "brown", vjust = -2)
#   ) +
#   scale_y_continuous(breaks = c(1.25, 1.50, 1.75, 2.00), limits = c(1, 2.1)) +
#   scale_x_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 17)) +
#   geom_text(
#     data = data.frame(
#       xpos = Inf,
#       ypos = Inf,
#       annotateText = "R-sq. (adj) = 0.414\np-value < 0.001",
#       hjustvar = 1.1,
#       vjustvar = 1.1
#     ), 
#     aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText), nudge_x =  -1.1, size = 10) +
#   geom_text(
#     data = data.frame(
#       xpos = -Inf,
#       ypos = -Inf,
#       annotateText = "NON-FOREST",
#       hjustvar = 0,
#       vjustvar = -0.1
#     ), 
#     aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText), size = 20, color = "brown")
# 
# 
# p.hist.nonfor <- ggExtra::ggMarginal(p.hist.in.nonfor, type = "histogram")
# 
# p.hist.out <- ggpubr::ggarrange(p.hist.for, p.hist.nonfor, ncol = 1, nrow = 2, align = "hv")
# 
# ggsave("03.results/04.FD-PD-GAM-hist-for.nonfor.png", p.hist.out,
#        height=30, width=20, units="in", dpi=600, bg = "white")
# 
