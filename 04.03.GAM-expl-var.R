####recommended machine: 01####

library("tidyverse")
library("sp")
library("tidyr")
library("dismo")
library("gbm")
library("mgcv")
library("marginaleffects")
library("viridis")
library("doSNOW")

sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-data.RDS")

exp <- c("PC1", 
         "PC2",
         "PC3",
         "PC4",
         "PC5",
         "plot.size", "Plants recorded", #data bias
         "sBiome", "is.forest", #Vegetation type
         "stable.clim" #Longterm climate stability (LGM)
)

sPlot.data$is.forest <- as.factor(sPlot.data$is.forest)
sPlot.data <- sPlot.data %>% as.data.frame() %>% filter(!is.na(plot.size),
                                                        !is.na(stable.clim),
                                                        !is.na(PC1),
                                                        !is.na(PC2),
                                                        !is.na(PC3),
                                                        !is.na(PC4),
                                                        !is.na(PC5),
                                                        !is.na(is.forest))

coords <- sPlot.data  %>% dplyr::select(Longitude, Latitude)
indices <- sPlot.data %>%  
  dplyr::select(
    SES.RQEP, exp)

data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                   data = indices,
                                   proj4string = crs("+proj=eck4"))

mod <- mgcv::gam(SES.RQEP ~ PC1 + is.forest + s(Longitude, Latitude, 
                                    bs = "sos"), 
                    family = "gaussian", 
                    method = "REML", 
                    data = data.wgs)

saveRDS(mod, "02.data/05.GAM_SES.RQEP-exp.RDS")

coords <- sPlot.data  %>% 
  dplyr::select(Longitude, Latitude)
indices <- sPlot.data %>% 
  dplyr::select(
    SES.RQEF, exp)

data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                   data = indices,
                                   proj4string = crs("+proj=eck4"))

mod <- mgcv::gam(SES.RQEF ~ stable.clim + PC2 + PC3 + PC5 + plot.size + is.forest + 
                   s(Longitude, Latitude, 
                     bs = "sos"), 
                 family = "gaussian", 
                 method = "REML", 
                 data = data.wgs)

saveRDS(mod, "02.data/05.GAM_SES.RQEF-exp.RDS")


coords <- sPlot.data  %>% 
  dplyr::select(Longitude, Latitude)
indices <- sPlot.data %>% 
  dplyr::select(
    SES.RQEF, exp)

data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                   data = indices,
                                   proj4string = CRS("+proj=eck4"))

forms <- c("SES.RQEF ~ PC2 + PC3 + PC5 + plot.size + is.forest+ s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEF ~ stable.clim + PC3 + PC5 + plot.size + is.forest+ s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEF ~ stable.clim + PC2 + PC5 + plot.size + is.forest+ s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEF ~ stable.clim + PC2 + PC3 + plot.size + is.forest+ s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEF ~ stable.clim + PC2 + PC3 + PC5 + is.forest+ s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEF ~ stable.clim + PC2 + PC3 + PC5 + plot.size + s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEF ~ 1 + s(Longitude, Latitude, bs = \"sos\")")

threads <- 7

cl <- parallel::makeCluster(threads)
registerDoSNOW(cl)

result <- foreach(i=1:length(forms)) %dopar% {
  
  mod <- mgcv::gam(as.formula(forms[i]), 
                   family = "gaussian", 
                   method = "REML", 
                   data = data.wgs)
  
  saveRDS(mod, paste0("02.data/05.GAM_SES.RQEF-exp_", i, ".RDS"))
}




sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-data.RDS")

exp <- c("PC1", 
         "PC2",
         "PC3",
         "PC4",
         "PC5",
         "plot.size", "Plants recorded", #data bias
         "sBiome", "is.forest", #Vegetation type
         "stable.clim" #Longterm climate stability (LGM)
)

sPlot.data$is.forest <- as.factor(sPlot.data$is.forest)
sPlot.data <- sPlot.data %>% as.data.frame() %>% filter(
                                                        !is.na(PC1),
                                                        !is.na(is.forest))

coords <- sPlot.data  %>% 
  dplyr::select(Longitude, Latitude)
indices <- sPlot.data %>% 
  dplyr::select(
    SES.RQEP, exp)

data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                   data = indices,
                                   proj4string = CRS("+proj=eck4"))

forms <- c("SES.RQEP ~ is.forest + s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEP ~ PC1 + s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEP ~ 1 + s(Longitude, Latitude, bs = \"sos\")")

threads <- 3

cl <- parallel::makeCluster(threads)
registerDoSNOW(cl)

result <- foreach(i=1:length(forms)) %dopar% {
  
  mod <- mgcv::gam(as.formula(forms[i]), 
                   family = "gaussian", 
                   method = "REML", 
                   data = data.wgs)
  
  saveRDS(mod, paste0("02.data/05.GAM_SES.RQEP-exp_", i, ".RDS"))
}


