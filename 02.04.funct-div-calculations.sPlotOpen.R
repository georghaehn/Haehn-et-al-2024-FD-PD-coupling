####recommended machine: 03#####

####packages ----
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("FD"))
suppressPackageStartupMessages(library("funrar"))

sessionInfo()

args <- commandArgs(trailingOnly = T)
input_dir <- args[1]
output <- args[2]
start <- as.numeric(args[3])
end <- as.numeric(args[4])

####load data ----
suppressWarnings(load(file.path(input_dir, "sPlotOpen.RData")))
suppressWarnings(load(file.path(input_dir, "Traits_CWMs_sPlot3.RData")))
suppressWarnings(load(file.path(input_dir, "pruned.trees.traitsp.sPlotOpen.RData")))

rm(CWM)

traits <- sPlot.traits %>% 
  dplyr::select(Species,
                GrowthForm, 
               LeafN_mean, 
               LeafP_mean, 
               SLA_mean, 
               PlantHeight_mean, 
               SeedMass_mean, 
               LDMC_mean,
               SpecificRootLength_mean, 
               Chromosome.n_mean) %>%
  filter(complete.cases(.))

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
  spp1 <- mylist$mydata[[1]] %>%
    dplyr::select(Species, Abundance) %>% 
    arrange(Species) %>%
    group_by(Species) %>% 
    #arrange(desc("Species")) %>% 
    dplyr::summarise(Abundance = (sum(Abundance, na.rm = T)>0) * 1) %>% 
    filter(Abundance > 0) %>% 
    pivot_wider(names_from = Species, 
                values_from = Abundance) %>% 
    as.data.frame()
  
  FD.out <- data.frame(PlotObservationID = mylist$PlotObservationID)
  
  trait.ALL <- mylist$mydata[[1]] %>% 
    dplyr::select(-Abundance) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
    
  
  if(nrow(trait.ALL)>1 &&
     ncol(spp1)>1) {
    
    spp1.ALL <-  spp1[,colnames(spp1) %in% 
                           trait.ALL$Species]
    
    trait.ALL <- trait.ALL %>% 
      filter(Species %in% colnames(spp1.ALL)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.ALL),colnames(spp1.ALL)) &&
       sum(trait.ALL)>0) {
      
    ALL <- dbFD(trait.ALL, spp1.ALL, 
                  w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
    
    FD.out$FDis.ALL <- ALL$FDis
    FD.out$RQE.ALL <- raoQ(spp1.ALL, trait.ALL)
    
    
    
    }
  }
  
  trait.SUB <- mylist$mydata[[1]] %>% 
    dplyr::select(Species,
                  SeedMass_mean, 
           LDMC_mean, 
           GrowthForm) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
  
  
  if(nrow(trait.SUB)>1 &&
     ncol(spp1)>1) {
    
    spp1.SUB <-  spp1[,colnames(spp1) %in% 
                        trait.SUB$Species]
    
    trait.SUB <- trait.SUB %>% 
      filter(Species %in% colnames(spp1.SUB)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.SUB),colnames(spp1.SUB)) &&
       sum(trait.SUB)>0) {
      
      SUB <- dbFD(trait.SUB, spp1.SUB, 
                  w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
      
      FD.out$FDis.SUB <- SUB$FDis
      FD.out$RQE.SUB <- raoQ(spp1.SUB, trait.SUB)
      
      
      
    }
  }

  
  trait.MULTI <- mylist$mydata[[1]] %>% 
    dplyr::select(Species,
                  SLA_mean,
                  SpecificRootLength_mean,
                  PlantHeight_mean) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
  
  
  if(nrow(trait.MULTI)>1 &&
     ncol(spp1)>1) {
    
    spp1.MULTI <-  spp1[,colnames(spp1) %in% 
                        trait.MULTI$Species]
    
    trait.MULTI <- trait.MULTI %>% 
      filter(Species %in% colnames(spp1.MULTI)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.MULTI),colnames(spp1.MULTI)) &&
       sum(trait.MULTI)>0) {
      
      MULTI <- dbFD(trait.MULTI, spp1.MULTI, 
                  w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
      
      FD.out$FDis.MULTI <- MULTI$FDis
      FD.out$RQE.MULTI <- raoQ(spp1.MULTI, trait.MULTI)
      
      
      
    }
  }
  
  trait.SLA <- mylist$mydata[[1]] %>% 
    dplyr::select(Species, SLA_mean) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
  
  if(nrow(trait.SLA)>1 &&
     ncol(spp1)>1) {
    
    spp1.SLA <-  spp1[,colnames(spp1) %in% 
                        trait.SLA$Species]
    
    trait.SLA <- trait.SLA %>% 
      filter(Species %in% colnames(spp1.SLA)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.SLA),colnames(spp1.SLA)) &&
       sum(trait.SLA)>0) {
      
      SLA <- dbFD(trait.SLA, spp1.SLA, 
                  w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
      
      FD.out$FDis.SLA <- SLA$FDis
      FD.out$RQE.SLA <- raoQ(spp1.SLA, trait.SLA)
    }
  }
  
  trait.HEIGHT <- mylist$mydata[[1]] %>% 
    dplyr::select(Species, PlantHeight_mean) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
  
  
  if(nrow(trait.HEIGHT)>1 &&
     ncol(spp1)>1) {
    
    spp1.HEIGHT <-  spp1[,colnames(spp1) %in% 
                           trait.HEIGHT$Species]
    
    trait.HEIGHT <- trait.HEIGHT %>% 
      filter(Species %in% colnames(spp1.HEIGHT)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.HEIGHT),colnames(spp1.HEIGHT)) &&
       sum(trait.HEIGHT)>0){
      
      HEIGHT <- dbFD(trait.HEIGHT, spp1.HEIGHT, 
                     w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
      
      FD.out$FDis.HEIGHT <- HEIGHT$FDis
      FD.out$RQE.HEIGHT <- raoQ(spp1.HEIGHT, trait.HEIGHT)
    }
  }
  
  trait.ROOT <- mylist$mydata[[1]] %>% 
    dplyr::select(Species, SpecificRootLength_mean) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
  
  if(nrow(trait.ROOT)>1 &&
     ncol(spp1)>1) {
    
    spp1.ROOT <-  spp1[,colnames(spp1) %in% 
                         trait.ROOT$Species]
    
    trait.ROOT <- trait.ROOT %>% 
      filter(Species %in% colnames(spp1.ROOT)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.ROOT),colnames(spp1.ROOT)) &&
       sum(trait.ROOT)>0) {
      
      ROOT <- dbFD(trait.ROOT, spp1.ROOT, 
                   w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
      
      FD.out$FDis.ROOT <- ROOT$FDis
      FD.out$RQE.ROOT <- raoQ(spp1.ROOT, trait.ROOT)
    }
  }
  
  trait.SM <- mylist$mydata[[1]] %>% 
    dplyr::select(Species, SeedMass_mean) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
  
  if(nrow(trait.SM)>1 &&
     ncol(spp1)>1) {
    
    spp1.SM <-  spp1[,colnames(spp1) %in% 
                         trait.SM$Species]
    
    trait.SM <- trait.SM %>% 
      filter(Species %in% colnames(spp1.SM)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.SM),colnames(spp1.SM)) &&
       sum(trait.SM)>0) {
      
      SM <- dbFD(trait.SM, spp1.SM, 
                   w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
      
      FD.out$FDis.SM <- SM$FDis
      FD.out$RQE.SM <- raoQ(spp1.SM, trait.SM)
    }
  }
  
  trait.LDMC <- mylist$mydata[[1]] %>% 
    dplyr::select(Species, LDMC_mean) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
  
  if(nrow(trait.LDMC)>1 &&
     ncol(spp1)>1) {
    
    spp1.LDMC <-  spp1[,colnames(spp1) %in% 
                         trait.LDMC$Species]
    
    trait.LDMC <- trait.LDMC %>% 
      filter(Species %in% colnames(spp1.LDMC)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.LDMC),colnames(spp1.LDMC)) &&
       sum(trait.LDMC)>0) {
      
      LDMC <- dbFD(trait.LDMC, spp1.LDMC, 
                   w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
      
      FD.out$FDis.LDMC <- LDMC$FDis
      FD.out$RQE.LDMC <- raoQ(spp1.LDMC, trait.LDMC)
    }
  }
  
  trait.GF <- mylist$mydata[[1]] %>% 
    dplyr::select(Species, GrowthForm) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
  
  if(nrow(trait.GF)>1 &&
     ncol(spp1)>1) {
    
    spp1.GF <-  spp1[,colnames(spp1) %in% 
                         trait.GF$Species]
    
    trait.GF <- trait.GF %>% 
      filter(Species %in% colnames(spp1.GF)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.GF),colnames(spp1.GF)) &&
       sum(trait.GF)>0) {
      
      GF <- dbFD(trait.GF, spp1.GF, 
                   w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
      
      FD.out$FDis.GF <- GF$FDis
      FD.out$RQE.GF <- raoQ(spp1.GF, trait.GF)
    }
  }
  
  trait.N <- mylist$mydata[[1]] %>% 
    dplyr::select(Species, LeafN_mean) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
  
  if(nrow(trait.N)>1 &&
     ncol(spp1)>1) {
    
    spp1.N <-  spp1[,colnames(spp1) %in% 
                       trait.N$Species]
    
    trait.N <- trait.N %>% 
      filter(Species %in% colnames(spp1.N)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.N),colnames(spp1.N)) &&
       sum(trait.N)>0) {
      
      N <- dbFD(trait.N, spp1.N, 
                 w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
      
      FD.out$FDis.N <- N$FDis
      FD.out$RQE.N <- raoQ(spp1.N, trait.N)
    }
  }
  
  trait.P <- mylist$mydata[[1]] %>% 
    dplyr::select(Species, LeafP_mean) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
  
  if(nrow(trait.P)>1 &&
     ncol(spp1)>1) {
    
    spp1.P <-  spp1[,colnames(spp1) %in% 
                         trait.P$Species]
    
    trait.P <- trait.P %>% 
      filter(Species %in% colnames(spp1.P)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.P),colnames(spp1.P)) &&
       sum(trait.P)>0) {
      
      P <- dbFD(trait.P, spp1.P, 
                   w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
      
      FD.out$FDis.P <- P$FDis
      FD.out$RQE.P <- raoQ(spp1.P, trait.P)
    }
  }
  
  trait.CHRO <- mylist$mydata[[1]] %>% 
    dplyr::select(Species, Chromosome.n_mean) %>% 
    distinct(Species, .keep_all = T) %>% 
    arrange(Species)
  
  if(nrow(trait.CHRO)>1 &&
     ncol(spp1)>1) {
    
    spp1.CHRO <-  spp1[,colnames(spp1) %in% 
                         trait.CHRO$Species]
    
    trait.CHRO <- trait.CHRO %>% 
      filter(Species %in% colnames(spp1.CHRO)) %>% 
      column_to_rownames("Species") %>%
      compute_dist_matrix()
    
    if(identical(rownames(trait.CHRO),colnames(spp1.CHRO)) &&
       sum(trait.CHRO)>0) {
      
      CHRO <- dbFD(trait.CHRO, spp1.CHRO, 
                   w.abun = F, calc.CWM = F, calc.FRic = F, calc.FDiv = F)
      
      FD.out$FDis.CHRO <- CHRO$FDis
      FD.out$RQE.CHRO <- raoQ(spp1.CHRO, trait.CHRO)
    }
  }
  
  return(FD.out)
}

PlotIDs <- DT2.oa %>% 
  filter(Species %in% sPlot.in.Tree$Species) %>% 
  group_by(PlotObservationID) %>% 
  filter(n()>1) %>% 
  distinct(PlotObservationID) %>% pull(PlotObservationID)

if(end < length(PlotIDs)) {
  PlotIDs <- PlotIDs[start:end]
  } else if(end >= length(PlotIDs)) {PlotIDs <- PlotIDs[start:length(PlotIDs)]}

FD.out <- DT2.oa %>% 
  filter(PlotObservationID %in% PlotIDs) %>% 
  filter(Species %in% sPlot.in.Tree$Species) %>% 
  filter(Species %in% traits$Species) %>%
  dplyr::select(PlotObservationID, Species, Original_abundance) %>% 
  left_join(traits, by = "Species") %>% 
  group_by(PlotObservationID) %>% 
  filter(n()>1) %>% 
  dplyr::summarise(mydata=list(data.frame(Species=Species,
                                          Abundance=Original_abundance, 
                                          GrowthForm = GrowthForm, 
                                          LeafN_mean = LeafN_mean, 
                                          LeafP_mean = LeafP_mean, 
                                          SLA_mean = SLA_mean, 
                                          PlantHeight_mean = PlantHeight_mean, 
                                          SeedMass_mean =  SeedMass_mean, 
                                          LDMC_mean = LDMC_mean,
                                          SpecificRootLength_mean = SpecificRootLength_mean, 
                                          Chromosome.n_mean = Chromosome.n_mean
                                          ))) %>% 
  split(.$PlotObservationID) %>% 
  map_dfr(., FD3)# %>% 
 # mutate(PlotObservationID = PlotIDs)

  
saveRDS(FD.out, file = output)


