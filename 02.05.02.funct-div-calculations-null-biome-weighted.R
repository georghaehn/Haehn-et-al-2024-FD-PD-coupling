####recommended machine: 03#####

####packages ----
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("FD"))
suppressPackageStartupMessages(library("funrar"))
sessionInfo()

args <- commandArgs(trailingOnly = T)
input_dir <- args[1]
output <- args[2]
SR.start <- as.numeric(args[3])
SR.end <- as.numeric(args[4])

####load data ----
suppressWarnings(load(file.path(input_dir, "Traits_CWMs_sPlot3.RData")))
suppressWarnings(load(file.path(input_dir, "pruned.trees.traitsp.RData")))

rm(CWM)

Species.Biome <- readRDS(file.path(input_dir, "species-list-biome.RDS"))

#get biomes

biom <- Species.Biome %>% 
  pull(Biome) %>% 
  unique()

traits <- sPlot.traits %>% 
  select(c("Species", 
           "SLA_mean", 
           "PlantHeight_mean",
           "SpecificRootLength_mean")) %>% 
  filter(complete.cases(.)) %>% 
  filter(Species %in% sPlot.in.Tree$Species)

raoQ <- function(abund, trait, Hill = TRUE, scale = FALSE, method = "default") {
  
  abund <- as.matrix(abund)
  
  anames <- colnames(abund)[order(colnames(abund))]
  
  abund <- abund[, anames, drop = FALSE]
  
  trait <- as.matrix(trait)
  
  trait <- trait[anames, anames]
  
  if(ncol(abund) != ncol(trait))
    
    stop("Not all species in the abundance matrix appear in the trait matrix!") 
  
  abund <- abund / rowSums(abund)
  
  if(method == "default") 
    
    Q <- apply(abund, 1, function(x) crossprod(x, trait %*% x))
  
  if(method == "divc")
    
    Q <- apply(abund, 1, function(x) x %*% trait^2 %*% (x/2/sum(x)^2))
  
  if(Hill == TRUE) Q <- 1/(1 - Q)
  
  if(scale == TRUE) Q <- Q / max(Q)
  
  names(Q) <- rownames(abund)
  
  return(Q)
  
}

FD3 <- function(mylist) {
  spp1 <- mylist %>%
    arrange(Species) %>% 
    pivot_wider(names_from = Species, 
                values_from = Abundance) %>% 
    as.data.frame()
  
  # trait.SLA <- traits %>% 
  #   select(Species, SLA_mean) %>% 
  #   filter(Species %in% colnames(spp1)) %>% 
  #   distinct(Species, .keep_all = T) %>% 
  #   arrange(Species) %>% 
  #   column_to_rownames("Species") %>% 
  #   filter(complete.cases(.)) %>% 
  #   compute_dist_matrix()
  # 
  # trait.HEIGHT <- traits %>% 
  #   select(Species, PlantHeight_mean) %>% 
  #   filter(Species %in% colnames(spp1)) %>% 
  #   distinct(Species, .keep_all = T) %>% 
  #   arrange(Species) %>% 
  #   column_to_rownames("Species") %>% 
  #   filter(complete.cases(.)) %>% 
  #   compute_dist_matrix()
  # 
  # trait.ROOT <- traits %>% 
  #   select(Species, SpecificRootLength_mean) %>% 
  #   filter(Species %in% colnames(spp1)) %>% 
  #   distinct(Species, .keep_all = T) %>% 
  #   arrange(Species) %>% 
  #   column_to_rownames("Species") %>% 
  #   filter(complete.cases(.)) %>% 
  #   compute_dist_matrix()
  
  trait.MULTI <- traits %>% 
    filter(Species %in% colnames(spp1)) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species) %>% 
    column_to_rownames("Species") %>% 
    filter(complete.cases(.)) %>% 
    compute_dist_matrix()
  
  
  if(length(rownames(trait.MULTI))>1 &&
     #length(rownames(trait.SLA))>1 &&
     #length(rownames(trait.ROOT))>1 &&
     #length(rownames(trait.HEIGHT))>1 &&
     length(colnames(spp1))>1
  ){
    FD.out <- data.frame(FDis.MULTI = NA)
    # SLA <- dbFD(trait.SLA, spp1, w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
    # FD.out$FDis.SLA <- SLA$FDis
    # FD.out$RQE.SLA <- raoQ(spp1, trait.SLA)
    # 
    # HEIGHT <- dbFD(trait.HEIGHT, spp1, w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
    # FD.out$FDis.HEIGHT <- HEIGHT$FDis
    # FD.out$RQE.HEIGHT <- raoQ(spp1, trait.HEIGHT)
    # 
    # ROOT <- dbFD(trait.ROOT, spp1, w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
    # FD.out$FDis.ROOT <- ROOT$FDis
    # FD.out$RQE.ROOT <- raoQ(spp1, trait.ROOT)
    
    MULTI <- dbFD(trait.MULTI, spp1, w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
    FD.out$FDis.MULTI <- MULTI$FDis
    FD.out$RQE.MULTI <- raoQ(spp1, trait.MULTI)
    
  } else {
    FD.out <- data.frame(FDis.SLA = NA)
  }
  
  return(FD.out)
}

# FDis.SLA <- list()
# RQE.SLA <- list()
# 
# FDis.HEIGHT <- list()
# RQE.HEIGHT <- list()
# 
# FDis.ROOT <- list()
# RQE.ROOT <- list()

FDis.bio <- list()
RQE.bio <- list()

FDis.out <- list()
RQE.out <- list()

spec <- c(SR.start:SR.end)
species <- traits$Species %>% as.data.frame() %>% rename(., Species = `.`)

for(i in 1:length(spec)) {
  
  for (j in 1:length(biom)) {
    
    FD.out <- replicate(499,
                        Species.Biome %>% 
                          filter(Biome == biom[j],
                                 Species %in% species$Species) %>% 
                          sample_n(spec[i], weight = Species_count) %>% 
                          as.data.frame() %>% 
                          #left_join(traits, by = "Species") %>% 
                          mutate(Abundance = 1) %>% 
                          dplyr::select(Species, Abundance), simplify = FALSE) %>% 
      map_dfr(., FD3)
    
    # FDis.SLA[[i]] <- FD.out$FDis.SLA %>% as.data.frame()
    # colnames(FDis.SLA[[i]]) <- paste0(spec[i])
    # RQE.SLA[[i]] <- FD.out$RQE.SLA %>% as.data.frame()
    # colnames(RQE.SLA[[i]]) <- paste0(spec[i])
    # 
    # FDis.HEIGHT[[i]] <- FD.out$FDis.HEIGHT %>% as.data.frame()
    # colnames(FDis.HEIGHT[[i]]) <- paste0(spec[i])
    # RQE.HEIGHT[[i]] <- FD.out$RQE.HEIGHT %>% as.data.frame()
    # colnames(RQE.HEIGHT[[i]]) <- paste0(spec[i])
    # 
    # FDis.ROOT[[i]] <- FD.out$FDis.ROOT %>% as.data.frame()
    # colnames(FDis.ROOT[[i]]) <- paste0(spec[i])
    # RQE.ROOT[[i]] <- FD.out$RQE.ROOT %>% as.data.frame()
    # colnames(RQE.ROOT[[i]]) <- paste0(spec[i])
    
    FDis.bio[[j]] <- FD.out$FDis.MULTI %>% 
      as.data.frame() %>%
      `colnames<-`("FDis.multi.null") %>% 
      mutate(Biome = biom[j],
             SR = spec[i])
    
    RQE.bio[[j]] <- FD.out$RQE.MULTI %>% 
      as.data.frame() %>% 
      `colnames<-`("RQE.multi.null") %>% 
      mutate(Biome = biom[j],
             SR = spec[i])
  }
  
  RQE.out[[i]] <- do.call(rbind, RQE.bio)
  
  FDis.out[[i]] <- do.call(rbind, FDis.bio)
}


# FDis.SLA.null <- do.call(cbind, FDis.SLA)
# RQE.SLA.null <- do.call(cbind, RQE.SLA)
# 
# FDis.HEIGHT.null <- do.call(cbind, FDis.HEIGHT)
# RQE.HEIGHT.null <- do.call(cbind, RQE.HEIGHT)
# 
# FDis.ROOT.null <- do.call(cbind, FDis.ROOT)
# RQE.ROOT.null <- do.call(cbind, RQE.ROOT)

RQE.MULTI.null <- do.call(rbind, RQE.out)
FDis.MULTI.null <- do.call(rbind, FDis.out)

list.out <- list(
  # FDis.SLA.null,
  # RQE.SLA.null,
  # FDis.HEIGHT.null,
  # RQE.HEIGHT.null, 
  # FDis.ROOT.null,
  # RQE.ROOT.null,
  FDis.MULTI.null,
  RQE.MULTI.null)

names(list.out) <- c(
  # "FDis.SLA.null",
  # "RQE.SLA.null",
  # "FDis.HEIGHT.null",
  # "RQE.HEIGHT.null", 
  # "FDis.ROOT.null",
  # "RQE.ROOT.null",
  "FDis.MULTI.null",
  "RQE.MULTI.null")

saveRDS(list.out, file = output)