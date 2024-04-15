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
                                   proj4string = CRS("+proj=eck4"))

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
                                   proj4string = CRS("+proj=eck4"))

mod <- mgcv::gam(SES.RQEF ~ stable.clim + PC2 + PC3 + PC5 + plot.size + is.forest + 
                   s(Longitude, Latitude, 
                     bs = "sos"), 
                 family = "gaussian", 
                 method = "REML", 
                 data = data.wgs)

saveRDS(mod, "02.data/05.GAM_SES.RQEF-exp.RDS")

mod <- mgcv::gam(SES.RQEF ~ stable.clim + PC2 + PC3 + PC5 + plot.size * is.forest + 
                   s(Longitude, Latitude, 
                     bs = "sos"), 
                 family = "gaussian", 
                 method = "REML", 
                 data = data.wgs)

saveRDS(mod, "02.data/05.GAM_SES.RQEF-exp-inter.RDS")

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



##### BW ----

sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-BW-data.RDS")

exp <- c("PC1", 
         "PC2",
         "PC3",
         "PC4",
         "PC5",
         "Plants recorded", #data bias
         "is.forest", #Vegetation type
         "stable.clim" #Longterm climate stability (LGM)
)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

sPlot.data$is.forest <- as.factor(sPlot.data$is.forest)
sPlot.data <- sPlot.data %>% 
  as.data.frame() %>% 
  filter(!is.na(stable.clim),
         !is.na(PC1),
         !is.na(PC2),
         !is.na(PC3),
         !is.na(PC4),
         !is.na(PC5),
         !is.na(is.forest),
         !is.na(SES.RQEF.BW))  %>% 
  mutate(SES.RQEP.BW = range01(SES.RQEP.BW),
         SES.RQEF.BW = range01(SES.RQEF.BW)) %>% 
  mutate(SES.RQEP.BW = ifelse(SES.RQEP.BW == 0, 0.001, SES.RQEP.BW),
         SES.RQEF.BW = ifelse(SES.RQEF.BW == 0, 0.001, SES.RQEF.BW)) %>%
  mutate(rel.BW = SES.RQEF.BW/SES.RQEP.BW) %>%
  mutate(rel.log.BW = log(rel.BW)) %>% 
  filter(rel.log.BW > -Inf,
         rel.log.BW < Inf)

coords <- sPlot.data  %>% dplyr::select(Longitude, Latitude)
indices <- sPlot.data %>%  
  dplyr::select(
    SES.RQEP.BW,
    SES.RQEF.BW, 
    rel.log.BW,
    exp)

data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                   data = indices,
                                   proj4string = CRS("+proj=eck4"))

forms <- c("SES.RQEP.BW ~ is.forest + PC1 + s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEP.BW ~ is.forest + s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEP.BW ~ PC1 + s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEP.BW ~ 1 + s(Longitude, Latitude, bs = \"sos\")",
           
           "SES.RQEF.BW ~ stable.clim + PC2 + PC5 + s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEF.BW ~ PC2 + PC5 + s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEF.BW ~ stable.clim + PC5 + s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEF.BW ~ stable.clim + PC2 + s(Longitude, Latitude, bs = \"sos\")",
           "SES.RQEF.BW ~ 1 + s(Longitude, Latitude, bs = \"sos\")",
           
           "rel.log.BW ~ is.forest + PC1 + s(Longitude, Latitude, bs = \"sos\")",
           "rel.log.BW ~ is.forest + s(Longitude, Latitude, bs = \"sos\")",
           "rel.log.BW ~ PC1 + s(Longitude, Latitude, bs = \"sos\")",
           "rel.log.BW ~ 1 + s(Longitude, Latitude, bs = \"sos\")"
           )

threads <- 4

cl <- parallel::makeCluster(threads)
registerDoSNOW(cl)

result <- foreach(i=1:length(forms)) %dopar% {
  
  mod <- mgcv::gam(as.formula(forms[i]), 
                   family = "gaussian", 
                   method = "REML", 
                   data = data.wgs)
  
  saveRDS(mod, paste0("02.data/05.GAM_", substr(forms[i],1, 11), "exp_", i, ".RDS"))
  
}
