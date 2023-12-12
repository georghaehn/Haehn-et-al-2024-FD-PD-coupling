####recommended machine: 02#####

####packages ----
library("tidyverse")
library("plyr")
library("tidyr")
library("dplyr")
library("V.PhyloMaker")
library("ape")
library("phytools")
library("tidytree")

####load sPlot data ----
load("~/data/sPlot/releases/sPlot3.0/Backbone3.0.RData")
load("~/data/sPlot/releases/sPlot3.0/DT_sPlot3.0.RData")

#select species, genus and family from backbone data
bb.species <- Backbone %>% 
             filter(!is.na(Name_correct), 
                    Rank_correct == "species") %>% 
             dplyr::select(species = Name_short, 
                           genus = Genus_correct, 
                           family = Family_correct, 
                           Rank_correct) %>% 
             distinct(species, .keep_all = T) %>% 
             filter(!(genus %in% lichen.genera), 
                    !(family %in% mushroom.families))

sPlot.species <- DT2 %>% 
                 distinct(Species) 

data <- left_join(sPlot.species, 
                  bb.species, 
                  by = c("Species" = "species")) #to add the genus and family information

rm(DT2, Backbone) #reduce needed RAM

####calculate tree with V.PhyloMaker package ----
tree3 =  phylo.maker(sp.list = data,
                     tree = GBOTB.extended,
                     nodes = nodes.info.1,
                     scenarios = "S3") #tip for a new genus binds to the 1/2 point of the family branch

####explore the tree ----

sum(tree3$species.list[[5]] == "bind") #39247
sum(tree3$species.list[[5]] == "fail to bind") #8860
sum(tree3$species.list[[5]] == "prune") #28805

bind <- tree3$species.list %>% as.data.frame() %>% 
        filter(status == "bind")

fail <- tree3$species.list %>% as.data.frame() %>% 
        filter(status == "fail to bind") %>% 
        mutate(species = gsub(" ", "_", species)) %>% 
        mutate(genus = ifelse(is.na(genus), species, genus)) %>% 
        mutate(genus = sub("\\_.*", "", genus),
               forcebind = NA)

####resolve the failed species ----
for(i in 1:nrow(fail)) {
      
      #get the indices in tip.label for every entry in the matching genus
      x <- grep(paste0(fail[i, "genus"], "_"), tree3$scenario.3$tip.label)
        
        #for x length 1 the most recent common ancestor (MRCA) is not defined, 
        #adding the new tip at 1/2 way along the terminal edge leading to the one species of the genus
        if(length(x) == 1){
                          tree3$scenario.3 <- bind.tip(tree3$scenario.3, paste0(fail[i, "species"]), where = x,
                          position = 0.5 * tree3$scenario.3$edge.length[which(tree3$scenario.3$edge[ ,2] == x)])
                          fail[i, "status"] <- "1/2 terminal level"
                                                      
                          } else if(length(x) > 1) {
                      
                                  #for x > 1 new species is bind to most recent common ancestor
                                  y <- findMRCA(tree3$scenario.3, tree3$scenario.3$tip.label[x])

                                  tree3$scenario.3 <- bind.tip(tree3$scenario.3, paste0(fail[i, "species"]), where=y)
                                  fail[i, "status"] <- "MRCA"

                                  } else fail[i, "forcebind"] <- "failed"
    }

#check which species remain unresolved
missing <- fail %>% filter(forcebind == "failed")

tree.species <- sp.tree %>% filter(status != "fail to bind") %>% 
                rbind(fail %>% select(-forcebind))

####save ----

sPlot.phylo.tree <- tree3
save(tree.species, sPlot.phylo.tree, file = "02.data/PhyloTree_sPlot3.0.Rdata")
