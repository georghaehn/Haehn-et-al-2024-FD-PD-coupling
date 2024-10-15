#### recommended machine: 03#####

#### packages ----
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("dismo"))
suppressPackageStartupMessages(library("gbm"))
sessionInfo()

args <- commandArgs(trailingOnly = T)
input_dir <- args[1]
output_dir <- args[2]
start <- args[3]

# sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-data.RDS")

sPlot.data <- readRDS(file.path(input_dir, "03.sPlot.PD.FD.CD-PCA-data.RDS"))


#### GBM ----
sPlot.data$sBiome <- as.factor(sPlot.data$sBiome)
sPlot.data$is.forest <- as.factor(sPlot.data$is.forest)
sPlot.data <- sPlot.data %>% as.data.frame()

set.seed(seed = 163)

div.indices <- c("SES.RQEP", "SES.RQEF", "RaoD.phyl", "RQE.MULTI", "rel.log")

exp <- c(
  "PC1",
  "PC2",
  "PC3",
  "PC4",
  "PC5",
  "plot.size", "Plants recorded", # data bias
  "sBiome", "is.forest", # Vegetation type
  "stable.clim" # Longterm climate stability (LGM)
)

# exp.variables <- c("mean.daily.air.temp.wet.qu",
#                    "mean.daily.air.temp.warm.qu",
#                    "annual.range.air.temp",
#                    "mean.monthly.prec.wet.qu",
#                    "mean.monthly.prec.warm.qu",
#                    "plot.size", "Plants recorded", #data bias
#                    "sBiome", "is.forest", #Vegetation type
#                    "stable.clim" #Longterm climate stability (LGM),
# )

i <- as.numeric(start)

if (i != 5) {
  model.GBM <- gbm.step(
    data = sPlot.data,
    gbm.x = exp,
    gbm.y = div.indices[i],
    family = "gaussian",
    tree.complexity = 5,
    learning.rate = 0.01,
    bag.fraction = 0.5,
    max.trees = 20000,
    silent = F
  )

  saveRDS(model.GBM, file.path(output_dir, paste0("04.BRT-results-", div.indices[i], ".RDS")))
} else {
  range01 <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }

  dat <- sPlot.data %>%
    mutate(
      SES.RQEP = range01(SES.RQEP),
      SES.RQEF = range01(SES.RQEF)
    ) %>%
    mutate(
      SES.RQEP = ifelse(SES.RQEP == 0, 0.001, SES.RQEP),
      SES.RQEF = ifelse(SES.RQEF == 0, 0.001, SES.RQEF)
    ) %>%
    mutate(rel = SES.RQEF / SES.RQEP) %>%
    mutate(rel.log = log(rel))


  dat$sBiome <- as.factor(dat$sBiome)
  dat$is.forest <- as.factor(dat$is.forest)
  dat <- dat %>%
    filter(
      rel.log > -Inf,
      rel.log < Inf
    ) %>%
    as.data.frame()

  set.seed(seed = 163)

  model.GBM <- gbm.step(
    data = dat,
    gbm.x = exp,
    gbm.y = div.indices[i],
    family = "gaussian",
    tree.complexity = 5,
    learning.rate = 0.01,
    bag.fraction = 0.5,
    max.trees = 20000,
    silent = F
  )

  saveRDS(model.GBM, file.path(output_dir, paste0("04.BRT-results-", div.indices[i], ".RDS")))
}
