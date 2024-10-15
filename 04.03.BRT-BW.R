#### recommended machine: 02#####

#### packages ----
library("tidyverse")
library("dismo")
library("gbm")
library("doSNOW")

sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-BW-data.RDS")

#### GBM ----
sPlot.data$is.forest <- as.factor(sPlot.data$is.forest)
sPlot.data <- sPlot.data %>%
  as.data.frame() %>%
  filter(
    !is.na(Longitude),
    !is.na(SES.RQEF.BW)
  )

exp <- c(
  "PC1",
  "PC2",
  "PC3",
  "PC4",
  "PC5",
  # "plot.size", #not relevant because of SES
  "Plants recorded", # data bias
  # "sBiome", #not relevant because of Biome species pools
  "is.forest", # Vegetation type
  "stable.clim" # Longterm climate stability (LGM)
)

threads <- 2

cl <- makeCluster(threads)
registerDoSNOW(cl)

result <- foreach(i = 1:2) %dopar% {
  if (i == 1) {
    library("dismo")
    library("gbm")

    set.seed(seed = 163)

    gbm.step(
      data = sPlot.data,
      gbm.x = exp,
      gbm.y = "SES.RQEF.BW",
      family = "gaussian",
      tree.complexity = 5,
      learning.rate = 0.01,
      bag.fraction = 0.5,
      max.trees = 20000,
      silent = F
    )
  } else {
    library("dismo")
    library("gbm")

    set.seed(seed = 163)

    gbm.step(
      data = sPlot.data,
      gbm.x = exp,
      gbm.y = "SES.RQEP.BW",
      family = "gaussian",
      tree.complexity = 5,
      learning.rate = 0.01,
      bag.fraction = 0.5,
      max.trees = 20000,
      silent = F
    )
  }
}

saveRDS(result, "02.data/04.BRT-results-SES.BW.RDS")

summary(result[[1]])

summary(result[[2]])


range01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

dat <- sPlot.data %>%
  filter(
    !is.na(SES.RQEP.BW),
    !is.na(SES.RQEF.BW)
  ) %>%
  mutate(
    SES.RQEP.BW = range01(SES.RQEP.BW),
    SES.RQEF.BW = range01(SES.RQEF.BW)
  ) %>%
  mutate(
    SES.RQEP.BW = ifelse(SES.RQEP.BW == 0, 0.001, SES.RQEP.BW),
    SES.RQEF.BW = ifelse(SES.RQEF.BW == 0, 0.001, SES.RQEF.BW)
  ) %>%
  mutate(rel.BW = SES.RQEF.BW / SES.RQEP.BW) %>%
  mutate(rel.log.BW = log(rel.BW))

#### GBM ----
dat$is.forest <- as.factor(dat$is.forest)
dat <- dat %>%
  filter(
    rel.log.BW > -Inf,
    rel.log.BW < Inf
  ) |>
  # mutate(rel.log = ifelse(rel.log == -Inf, -100, rel.log), #-Inf replaced by -100, otherwise the BRT will not run, same for Inf
  #        rel.log = ifelse(rel.log == Inf, 100, rel.log)) %>%
  as.data.frame()

set.seed(seed = 163)

exp <- c(
  "PC1",
  "PC2",
  "PC3",
  "PC4",
  "PC5",
  # "plot.size",
  "Plants recorded", # data bias
  # "sBiome",
  "is.forest", # Vegetation type
  "stable.clim" # Longterm climate stability (LGM)
)

model.GBM <- gbm.step(
  data = dat,
  gbm.x = exp,
  gbm.y = "rel.log.BW",
  family = "gaussian",
  tree.complexity = 5,
  learning.rate = 0.01,
  bag.fraction = 0.5,
  max.trees = 20000,
  silent = F
)



saveRDS(model.GBM, "02.data/04.BRT-results-SES.BW.rel.log.RDS")

summary(model.GBM)
