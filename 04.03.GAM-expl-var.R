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
library("gratia")

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

rm(list = ls())

### by groups

sPlot.data <- readRDS("02.data/04.sPlot.PD.FD.CD-PCA-BW-Status-data.RDS")

sPlot.data$is.forest <- as.factor(sPlot.data$is.forest)
sPlot.data$Status <- factor(sPlot.data$Status, levels = c("Coupling", "Decoupling higher FD", "Decoupling higher PD"))

levels(sPlot.data$Status)
#"Coupling"             "Decoupling higher FD" "Decoupling higher PD"

sPlot.data <- sPlot.data %>% 
  as.data.frame() %>% 
  filter(!is.na(stable.clim),
         !is.na(PC1),
         !is.na(PC2),
         !is.na(PC3),
         !is.na(PC4),
         !is.na(PC5),
         !is.na(is.forest)) 

exp <- c("PC1", 
         "PC2",
         "PC3",
         "PC4",
         "PC5",
         "Plants recorded", #data bias
         "is.forest", #Vegetation type
         "stable.clim" #Longterm climate stability (LGM)
)

size.n = 10000
rep.n = 100

threads <- 10

cl <- parallel::makeCluster(threads)
registerDoSNOW(cl)

result <- foreach(i=1:rep.n) %dopar% {
  
  library("tidyverse")
  library("sp")
  library("mgcv")
  
  data.s <- sPlot.data %>% 
    group_by(Status) %>%
    slice_sample(n = size.n) %>%
    ungroup()
  
  coords <- data.s %>% dplyr::select(Longitude, Latitude)
  
  indices <- data.s %>%
    mutate(ST = as.numeric(Status)) %>% 
    dplyr::select(
      ST,
      Status,
      exp)
  
  data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                     data = indices,
                                     proj4string = CRS("+proj=eck4"))
  
  mod <- mgcv::gam(ST ~ is.forest + s(stable.clim) + s(PC1) + s(PC2) + s(PC5) +
                     s(Longitude, Latitude, 
                       bs = "sos"), 
                   family = ocat(R=3), 
                   method = "REML", 
                   data = data.wgs)
  
  out <- list(data.s, mod)
  
  saveRDS(out, paste0("02.data/04.GAM_status-rep-", i, ".RDS"))
  
  out
}

ll <- list()

for(i in 1:100) {
  
  ll[[i]] <- readRDS(paste0("02.data/04.GAM_status-rep-", i, ".RDS"))
  
}

ll1 <- list()


for(i in 1:length(ll)) {
  
  mod <- ll[[i]][[2]]
  df <- ll[[i]][[1]]
  
  sg <- summary(mod)
  
  dev.exp <- data.frame(
    Expl = "Dev. expl.",
    V1 = sg$dev.expl
  )  
  
  pv <- sg$s.table[,4] %>%
    t() %>% 
    as.data.frame() %>%
    `colnames<-`(c("D", "A", "B", "C", "latlong")) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("Expl") %>% 
    filter(Expl != "latlong")
  
  pfor <- sg$p.table[,4] %>% 
    t() %>% 
    as.data.frame() %>% 
    `colnames<-`(c("x", "E")) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("Expl") %>% 
    filter(Expl != "x")
  
  pout <- rbind(pv, pfor, dev.exp)
  
  pout$rep <- i
  
  thresh <- gratia::theta(mod) |>
    tibble::as_tibble() |>
    setNames(c("threshold"))
  
  thresh$rep <- i
  
  
  
  ### PC1 
  
  ds1 <- data_slice(mod, PC1 = evenly(df$PC1, n = 100))
  fv1 <- fitted_values(mod, data = ds1, scale = "link")
  fv1$rep <- i

  ### PC2
  
  ds2 <- data_slice(mod, PC2 = evenly(df$PC2, n = 100))
  fv2 <- fitted_values(mod, data = ds2, scale = "link") # <- link scale
  fv2$rep <- i

  ### PC5
  
  ds3 <- data_slice(mod, PC5 = evenly(df$PC5, n = 100))
  fv3 <- fitted_values(mod, data = ds3, scale = "link") # <- link scale
  fv3$rep <- i

  ### Stable climate
  
  ds4 <- data_slice(mod, stable.clim = evenly(df$stable.clim, n = 100))
  fv4 <- fitted_values(mod, data = ds4, scale = "link") # <- link scale
  fv4$rep <- i

  ### Forest
  
  ds5 <- data_slice(mod, is.forest = evenly(df$is.forest, n = 100))
  fv5 <- fitted_values(mod, data = ds5, scale = "link") # <- link scale
  fv5$rep <- i

  if(i == 1) {
    
    ll1[[1]] <- fv1
    ll1[[2]] <- fv2
    ll1[[3]] <- fv3
    ll1[[4]] <- fv4
    ll1[[5]] <- fv5
    ll1[[6]] <- thresh
    ll1[[7]] <- pout
    
  } else {
    
    ll1[[1]] <- rbind(ll1[[1]], fv1)
    ll1[[2]] <- rbind(ll1[[2]], fv2)
    ll1[[3]] <- rbind(ll1[[3]], fv3)
    ll1[[4]] <- rbind(ll1[[4]], fv4)
    ll1[[5]] <- rbind(ll1[[5]], fv5)
    ll1[[6]] <- rbind(ll1[[6]], thresh)
    ll1[[7]] <- rbind(ll1[[7]], pout)
    
  }
  
  names(ll1) <- c("PC1", "PC2", "PC5", "Stable.clim", "Forest", "Thresh", "pv")
}

saveRDS(ll1, "02.Data/04.GAM-status-ll1.RDS")


### changing the levels of the groups

sPlot.data <- readRDS("02.data/04.sPlot.PD.FD.CD-PCA-BW-Status-data.RDS")

sPlot.data$is.forest <- as.factor(sPlot.data$is.forest)
sPlot.data$Status <- factor(sPlot.data$Status, levels = c("Decoupling higher FD", "Coupling", "Decoupling higher PD"))

levels(sPlot.data$Status)
#"Coupling"             "Decoupling higher FD" "Decoupling higher PD"

sPlot.data <- sPlot.data %>% 
  as.data.frame() %>% 
  filter(!is.na(stable.clim),
         !is.na(PC1),
         !is.na(PC2),
         !is.na(PC3),
         !is.na(PC4),
         !is.na(PC5),
         !is.na(is.forest)) 

exp <- c("PC1", 
         "PC2",
         "PC3",
         "PC4",
         "PC5",
         "Plants recorded", #data bias
         "is.forest", #Vegetation type
         "stable.clim" #Longterm climate stability (LGM)
)

size.n = 10000
rep.n = 100

threads <- 10

cl <- parallel::makeCluster(threads)
registerDoSNOW(cl)

result <- foreach(i=1:rep.n) %dopar% {
  
  library("tidyverse")
  library("sp")
  library("mgcv")
  
  data.s <- sPlot.data %>% 
    group_by(Status) %>%
    slice_sample(n = size.n) %>%
    ungroup()
  
  coords <- data.s %>% dplyr::select(Longitude, Latitude)
  
  indices <- data.s %>%
    mutate(ST = as.numeric(Status)) %>% 
    dplyr::select(
      ST,
      Status,
      exp)
  
  data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                     data = indices,
                                     proj4string = CRS("+proj=eck4"))
  
  mod <- mgcv::gam(ST ~ is.forest + s(stable.clim) + s(PC1) + s(PC2) + s(PC5) +
                     s(Longitude, Latitude, 
                       bs = "sos"), 
                   family = ocat(R=3), 
                   method = "REML", 
                   data = data.wgs)
  
  out <- list(data.s, mod)
  
  saveRDS(out, paste0("02.data/04.GAM_status1-rep-", i, ".RDS"))
  
  out
}

ll <- list()

for(i in 1:100) {
  
  ll[[i]] <- readRDS(paste0("02.data/04.GAM_status1-rep-", i, ".RDS"))
  
}

summary(ll[[10]][[2]])

summary(ll[[1]][[1]]$Status)

ll1 <- list()


for(i in 1:length(ll)) {
  
  mod <- ll[[i]][[2]]
  df <- ll[[i]][[1]]
  
  sg <- summary(mod)
  
  dev.exp <- data.frame(
    Expl = "Dev. expl.",
    V1 = sg$dev.expl
  )  
  
  pv <- sg$s.table[,4] %>%
    t() %>% 
    as.data.frame() %>%
    `colnames<-`(c("D", "A", "B", "C", "latlong")) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("Expl") %>% 
    filter(Expl != "latlong")
  
  pfor <- sg$p.table[,4] %>% 
    t() %>% 
    as.data.frame() %>% 
    `colnames<-`(c("x", "E")) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("Expl") %>% 
    filter(Expl != "x")
  
  pout <- rbind(pv, pfor, dev.exp)
  
  pout$rep <- i
  
  thresh <- gratia::theta(mod) |>
    tibble::as_tibble() |>
    setNames(c("threshold"))
  
  thresh$rep <- i
  
  
  
  ### PC1 
  
  ds1 <- data_slice(mod, PC1 = evenly(df$PC1, n = 100))
  fv1 <- fitted_values(mod, data = ds1, scale = "link")
  fv1$rep <- i
  
  ### PC2
  
  ds2 <- data_slice(mod, PC2 = evenly(df$PC2, n = 100))
  fv2 <- fitted_values(mod, data = ds2, scale = "link") # <- link scale
  fv2$rep <- i
  
  ### PC5
  
  ds3 <- data_slice(mod, PC5 = evenly(df$PC5, n = 100))
  fv3 <- fitted_values(mod, data = ds3, scale = "link") # <- link scale
  fv3$rep <- i
  
  ### Stable climate
  
  ds4 <- data_slice(mod, stable.clim = evenly(df$stable.clim, n = 100))
  fv4 <- fitted_values(mod, data = ds4, scale = "link") # <- link scale
  fv4$rep <- i
  
  ### Forest
  
  ds5 <- data_slice(mod, is.forest = evenly(df$is.forest, n = 100))
  fv5 <- fitted_values(mod, data = ds5, scale = "link") # <- link scale
  fv5$rep <- i
  
  if(i == 1) {
    
    ll1[[1]] <- fv1
    ll1[[2]] <- fv2
    ll1[[3]] <- fv3
    ll1[[4]] <- fv4
    ll1[[5]] <- fv5
    ll1[[6]] <- thresh
    ll1[[7]] <- pout
    
  } else {
    
    ll1[[1]] <- rbind(ll1[[1]], fv1)
    ll1[[2]] <- rbind(ll1[[2]], fv2)
    ll1[[3]] <- rbind(ll1[[3]], fv3)
    ll1[[4]] <- rbind(ll1[[4]], fv4)
    ll1[[5]] <- rbind(ll1[[5]], fv5)
    ll1[[6]] <- rbind(ll1[[6]], thresh)
    ll1[[7]] <- rbind(ll1[[7]], pout)
    
  }
  
  names(ll1) <- c("PC1", "PC2", "PC5", "Stable.clim", "Forest", "Thresh", "pv")
}

saveRDS(ll1, "02.Data/04.GAM-status-ll2.RDS")


### changing the levels of the groups

sPlot.data <- readRDS("02.data/04.sPlot.PD.FD.CD-PCA-BW-Status-data.RDS")

sPlot.data$is.forest <- as.factor(sPlot.data$is.forest)
sPlot.data$Status <- factor(sPlot.data$Status, levels = c("Decoupling higher PD", "Coupling", "Decoupling higher FD"))

levels(sPlot.data$Status)
#"Coupling"             "Decoupling higher FD" "Decoupling higher PD"

sPlot.data <- sPlot.data %>% 
  as.data.frame() %>% 
  filter(!is.na(stable.clim),
         !is.na(PC1),
         !is.na(PC2),
         !is.na(PC3),
         !is.na(PC4),
         !is.na(PC5),
         !is.na(is.forest)) 

exp <- c("PC1", 
         "PC2",
         "PC3",
         "PC4",
         "PC5",
         "Plants recorded", #data bias
         "is.forest", #Vegetation type
         "stable.clim" #Longterm climate stability (LGM)
)

size.n = 10000
rep.n = 100

threads <- 10

cl <- parallel::makeCluster(threads)
registerDoSNOW(cl)

result <- foreach(i=1:rep.n) %dopar% {
  
  library("tidyverse")
  library("sp")
  library("mgcv")
  
  data.s <- sPlot.data %>% 
    group_by(Status) %>%
    slice_sample(n = size.n) %>%
    ungroup()
  
  coords <- data.s %>% dplyr::select(Longitude, Latitude)
  
  indices <- data.s %>%
    mutate(ST = as.numeric(Status)) %>% 
    dplyr::select(
      ST,
      Status,
      exp)
  
  data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                     data = indices,
                                     proj4string = CRS("+proj=eck4"))
  
  mod <- mgcv::gam(ST ~ is.forest + s(stable.clim) + s(PC1) + s(PC2) + s(PC5) +
                     s(Longitude, Latitude, 
                       bs = "sos"), 
                   family = ocat(R=3), 
                   method = "REML", 
                   data = data.wgs)
  
  out <- list(data.s, mod)
  
  saveRDS(out, paste0("02.data/04.GAM_status2-rep-", i, ".RDS"))
  
  out
}

library(gratia)

ll <- list()

for(i in 1:100) {
  
  ll[[i]] <- readRDS(paste0("02.data/04.GAM_status2-rep-", i, ".RDS"))
  
}

summary(ll[[10]][[2]])

summary(ll[[1]][[1]]$Status)

ll1 <- list()


for(i in 1:length(ll)) {
  
  mod <- ll[[i]][[2]]
  df <- ll[[i]][[1]]
  
  sg <- summary(mod)
  dev.exp <- data.frame(
    Expl = "Dev. expl.",
    V1 = sg$dev.expl
  )  
  pv <- sg$s.table[,4] %>%
    t() %>% 
    as.data.frame() %>%
    `colnames<-`(c("D", "A", "B", "C", "latlong")) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("Expl") %>% 
    filter(Expl != "latlong")
  
  pfor <- sg$p.table[,4] %>% 
    t() %>% 
    as.data.frame() %>% 
    `colnames<-`(c("x", "E")) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("Expl") %>% 
    filter(Expl != "x")
  
  pout <- rbind(pv, pfor, dev.exp)
  
  pout$rep <- i
  
  thresh <- gratia::theta(mod) |>
    tibble::as_tibble() |>
    setNames(c("threshold"))
  
  thresh$rep <- i
  
  
  
  ### PC1 
  
  ds1 <- data_slice(mod, PC1 = evenly(df$PC1, n = 100))
  fv1 <- fitted_values(mod, data = ds1, scale = "link")
  fv1$rep <- i
  
  ### PC2
  
  ds2 <- data_slice(mod, PC2 = evenly(df$PC2, n = 100))
  fv2 <- fitted_values(mod, data = ds2, scale = "link") # <- link scale
  fv2$rep <- i
  
  ### PC5
  
  ds3 <- data_slice(mod, PC5 = evenly(df$PC5, n = 100))
  fv3 <- fitted_values(mod, data = ds3, scale = "link") # <- link scale
  fv3$rep <- i
  
  ### Stable climate
  
  ds4 <- data_slice(mod, stable.clim = evenly(df$stable.clim, n = 100))
  fv4 <- fitted_values(mod, data = ds4, scale = "link") # <- link scale
  fv4$rep <- i
  
  ### Forest
  
  ds5 <- data_slice(mod, is.forest = evenly(df$is.forest, n = 100))
  fv5 <- fitted_values(mod, data = ds5, scale = "link") # <- link scale
  fv5$rep <- i
  
  if(i == 1) {
    
    ll1[[1]] <- fv1
    ll1[[2]] <- fv2
    ll1[[3]] <- fv3
    ll1[[4]] <- fv4
    ll1[[5]] <- fv5
    ll1[[6]] <- thresh
    ll1[[7]] <- pout
    
  } else {
    
    ll1[[1]] <- rbind(ll1[[1]], fv1)
    ll1[[2]] <- rbind(ll1[[2]], fv2)
    ll1[[3]] <- rbind(ll1[[3]], fv3)
    ll1[[4]] <- rbind(ll1[[4]], fv4)
    ll1[[5]] <- rbind(ll1[[5]], fv5)
    ll1[[6]] <- rbind(ll1[[6]], thresh)
    ll1[[7]] <- rbind(ll1[[7]], pout)
    
  }
  
  names(ll1) <- c("PC1", "PC2", "PC5", "Stable.clim", "Forest", "Thresh", "pv")
}

saveRDS(ll1, "02.Data/04.GAM-status-ll3.RDS")

## linear instead of smooth


sPlot.data <- readRDS("02.data/04.sPlot.PD.FD.CD-PCA-BW-Status-data.RDS")

sPlot.data$is.forest <- as.factor(sPlot.data$is.forest)
sPlot.data$Status <- factor(sPlot.data$Status, levels = c("Decoupling higher PD", "Coupling", "Decoupling higher FD"))

levels(sPlot.data$Status)
#"Decoupling higher PD" "Coupling"             "Decoupling higher FD"

sPlot.data <- sPlot.data %>% 
  as.data.frame() %>% 
  filter(!is.na(stable.clim),
         !is.na(PC1),
         !is.na(PC2),
         !is.na(PC3),
         !is.na(PC4),
         !is.na(PC5),
         !is.na(is.forest)) 

exp <- c("PC1", 
         "PC2",
         "PC3",
         "PC4",
         "PC5",
         "Plants recorded", #data bias
         "is.forest", #Vegetation type
         "stable.clim" #Longterm climate stability (LGM)
)

size.n = 10000
rep.n = 100

threads <- 10

cl <- parallel::makeCluster(threads)
registerDoSNOW(cl)

result <- foreach(i=1:rep.n) %dopar% {
  
  library("tidyverse")
  library("sp")
  library("mgcv")
  
  data.s <- sPlot.data %>% 
    group_by(Status) %>%
    slice_sample(n = size.n) %>%
    ungroup()
  
  coords <- data.s %>% dplyr::select(Longitude, Latitude)
  
  indices <- data.s %>%
    mutate(ST = as.numeric(Status)) %>% 
    dplyr::select(
      ST,
      Status,
      exp)
  
  data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                     data = indices,
                                     proj4string = CRS("+proj=eck4"))
  
  mod <- mgcv::gam(ST ~ is.forest + stable.clim + PC1 + PC2 + PC5 +
                     s(Longitude, Latitude, 
                       bs = "sos"), 
                   family = ocat(R=3), 
                   method = "REML", 
                   data = data.wgs)
  
  out <- list(data.s, mod)
  
  saveRDS(out, paste0("02.data/04.GAM_status3-rep-", i, ".RDS"))
  
  out
}

library(gratia)

ll <- list()

for(i in 1:100) {
  
  ll[[i]] <- readRDS(paste0("02.data/04.GAM_status3-rep-", i, ".RDS"))
  
}

summary(ll[[10]][[2]])

summary(ll[[1]][[1]]$Status)

ll1 <- list()


for(i in 1:length(ll)) {
  
  mod <- ll[[i]][[2]]
  df <- ll[[i]][[1]]
  
  sg <- summary(mod)
  dev.exp <- data.frame(
    Expl = "Dev. expl.",
    V1 = sg$dev.expl
  )  
  pv <- sg$p.table[,4] %>%
    t() %>% 
    as.data.frame() %>%
    `colnames<-`(c("x", "E", "D", "A", "B", "C")) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("Expl") %>% 
    filter(Expl != "x")
  
  # pfor <- sg$p.table[,4] %>% 
  #   t() %>% 
  #   as.data.frame() %>% 
  #   `colnames<-`(c()) %>% 
  #   t() %>% 
  #   as.data.frame() %>% 
  #   rownames_to_column("Expl") %>% 
  #   filter(Expl != "x")
  
  pout <- rbind(pv, dev.exp)
  
  pout$rep <- i
  
  thresh <- gratia::theta(mod) |>
    tibble::as_tibble() |>
    setNames(c("threshold"))
  
  thresh$rep <- i
  
  
  
  ### PC1 
  
  ds1 <- data_slice(mod, PC1 = evenly(df$PC1, n = 100))
  fv1 <- fitted_values(mod, data = ds1, scale = "link")
  fv1$rep <- i
  
  ### PC2
  
  ds2 <- data_slice(mod, PC2 = evenly(df$PC2, n = 100))
  fv2 <- fitted_values(mod, data = ds2, scale = "link") # <- link scale
  fv2$rep <- i
  
  ### PC5
  
  ds3 <- data_slice(mod, PC5 = evenly(df$PC5, n = 100))
  fv3 <- fitted_values(mod, data = ds3, scale = "link") # <- link scale
  fv3$rep <- i
  
  ### Stable climate
  
  ds4 <- data_slice(mod, stable.clim = evenly(df$stable.clim, n = 100))
  fv4 <- fitted_values(mod, data = ds4, scale = "link") # <- link scale
  fv4$rep <- i
  
  ### Forest
  
  ds5 <- data_slice(mod, is.forest = evenly(df$is.forest, n = 100))
  fv5 <- fitted_values(mod, data = ds5, scale = "link") # <- link scale
  fv5$rep <- i
  
  if(i == 1) {
    
    ll1[[1]] <- fv1
    ll1[[2]] <- fv2
    ll1[[3]] <- fv3
    ll1[[4]] <- fv4
    ll1[[5]] <- fv5
    ll1[[6]] <- thresh
    ll1[[7]] <- pout
    
  } else {
    
    ll1[[1]] <- rbind(ll1[[1]], fv1)
    ll1[[2]] <- rbind(ll1[[2]], fv2)
    ll1[[3]] <- rbind(ll1[[3]], fv3)
    ll1[[4]] <- rbind(ll1[[4]], fv4)
    ll1[[5]] <- rbind(ll1[[5]], fv5)
    ll1[[6]] <- rbind(ll1[[6]], thresh)
    ll1[[7]] <- rbind(ll1[[7]], pout)
    
  }
  
  names(ll1) <- c("PC1", "PC2", "PC5", "Stable.clim", "Forest", "Thresh", "pv")
}

saveRDS(ll1, "02.Data/04.GAM-status-ll4.RDS")


## for the whole dataset - unbalanced

sPlot.data <- readRDS("02.data/04.sPlot.PD.FD.CD-PCA-BW-Status-data.RDS")

sPlot.data$is.forest <- as.factor(sPlot.data$is.forest)
sPlot.data$Status <- factor(sPlot.data$Status, levels = c("Decoupling higher PD", "Coupling", "Decoupling higher FD"))

levels(sPlot.data$Status)
#"Decoupling higher PD" "Coupling"             "Decoupling higher FD"

sPlot.data <- sPlot.data %>% 
  as.data.frame() %>% 
  filter(!is.na(stable.clim),
         !is.na(PC1),
         !is.na(PC2),
         !is.na(PC3),
         !is.na(PC4),
         !is.na(PC5),
         !is.na(is.forest)) 

exp <- c("PC1", 
         "PC2",
         "PC3",
         "PC4",
         "PC5",
         "Plants recorded", #data bias
         "is.forest", #Vegetation type
         "stable.clim" #Longterm climate stability (LGM)
)


coords <- sPlot.data %>% dplyr::select(Longitude, Latitude)
  
indices <- sPlot.data %>%
  mutate(ST = as.numeric(Status)) %>% 
  dplyr::select(
    ST,
    Status,
    exp)
  
data.wgs <- SpatialPointsDataFrame(coords = coords, 
                                   data = indices,
                                   proj4string = CRS("+proj=eck4"))
  
mod <- mgcv::gam(ST ~ is.forest + stable.clim + PC1 + PC2 + PC5 +
                   s(Longitude, Latitude, 
                     bs = "sos"), 
                 family = ocat(R=3), 
                 method = "REML", 
                 data = data.wgs)
  
saveRDS(mod, "02.data/04.GAM_status-compl-data.RDS")

mod <- readRDS("02.data/04.GAM_status-compl-data.RDS")

summary(mod)  

sg <- summary(mod)
dev.exp <- data.frame(
  Expl = "Dev. expl.",
  V1 = sg$dev.expl
)  
pv <- sg$p.table[,4] %>%
  t() %>% 
  as.data.frame() %>%
  `colnames<-`(c("x", "E", "D", "A", "B", "C")) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Expl") %>% 
  filter(Expl != "x")

pout <- rbind(pv, dev.exp)

thresh <- gratia::theta(mod) |>
  tibble::as_tibble() |>
  setNames(c("threshold"))

df <- sPlot.data

### PC1 

ds1 <- data_slice(mod, PC1 = evenly(df$PC1, n = 1000))
fv1 <- fitted_values(mod, data = ds1, scale = "link")

### PC2

ds2 <- data_slice(mod, PC2 = evenly(df$PC2, n = 1000))
fv2 <- fitted_values(mod, data = ds2, scale = "link") # <- link scale

### PC5

ds3 <- data_slice(mod, PC5 = evenly(df$PC5, n = 1000))
fv3 <- fitted_values(mod, data = ds3, scale = "link") # <- link scale

### Stable climate

ds4 <- data_slice(mod, stable.clim = evenly(df$stable.clim, n = 1000))
fv4 <- fitted_values(mod, data = ds4, scale = "link") # <- link scale

### Forest

ds5 <- data_slice(mod, is.forest = evenly(df$is.forest, n = 1000))
fv5 <- fitted_values(mod, data = ds5, scale = "link") # <- link scale

ll1 <- list()

ll1[[1]] <- fv1
ll1[[2]] <- fv2
ll1[[3]] <- fv3
ll1[[4]] <- fv4
ll1[[5]] <- fv5
ll1[[6]] <- thresh
ll1[[7]] <- pout

names(ll1) <- c("PC1", "PC2", "PC5", "Stable.clim", "Forest", "Thresh", "pv")

saveRDS(ll1, "02.Data/04.GAM-status-ll5.RDS")
