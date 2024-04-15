library(tidyverse)

#load TRY trait data
load("02.data/Traits_CWMs_sPlot3.RData")
load("02.data/PhyloTree_sPlot3.0.Rdata")


species <- sPlot.traits |> 
  filter(gsub(" ", "_", Species) %in% sPlot.phylo.tree$scenario.3$tip.label) |> 
  dplyr::select(Species,
                SLA_mean, 
                PlantHeight_mean, 
                SpecificRootLength_mean) |> 
  drop_na()

load("02.data/header_sPlot3.0.RData")
load("02.Data/DT_sPlot3.0.RData")

## filter DT by species present in traits and phylo

DT <- DT2 |> 
  filter(Species %in% species$Species)

## biomes
sum(is.na(header$sBiome))

DTT <- DT |> 
  group_by(PlotObservationID, Species) |> 
  summarise(Abundance_sum = sum(Abundance, na.rm = TRUE)) |> 
  left_join(
    header |> 
      dplyr::select(PlotObservationID, Biome = sBiome),
    by = "PlotObservationID"
  ) |> 
  group_by(Biome, Species) |> 
  summarise(Abundance_sum = sum(Abundance_sum, na.rm = TRUE),
            Species_count = n())

DTTT <- DTT |> 
  filter(!is.na(Biome))

saveRDS(DTTT, "02.Data/species-list-biome.RDS")

## number of of ferns, lycopods and/or gymnosperms per biome

ferns <- c(
           "Anemiaceae",
           "Aspleniaceae",
           "Athyriaceae",
           "Blechnaceae",
           "Cibotiaceae",
           "Culcitaceae",
           "Cyatheaceae",
           "Cystodiaceae",
           "Cystopteridaceae",
           "Davalliaceae",
           "Dennstaedtiaceae",
           "Desmophlebiaceae",
           "Dicksoniaceae",
           "Didymochlaenaceae",
           "Diplaziopsidaceae",
           "Dipteridaceae",
           "Dryopteridaceae",
           "Equisetaceae",
           "Gleicheniaceae",
           "Hemidictyaceae",
           "Hymenophyllaceae",
           "Hypodematiaceae",
           "Lindsaeaceae",
           "Lomariopsidaceae",
           "Lonchitidaceae",
           "Loxsomataceae",
           "Lygodiaceae",
           "Marattiaceae",
           "Marsileaceae",
           "Matoniaceae",
           "Metaxyaceae",
           "Nephrolepidaceae",
           "Oleandraceae",
           "Onocleaceae",
           "Ophioglossaceae",
           "Osmundaceae",
           "Plagiogyriaceae",
           "Polypodiaceae",
           "Psilotaceae",
           "Pteridaceae",
           "Rhachidosoraceae",
           "Saccolomataceae",
           "Salviniaceae",
           "Schizaeaceae",
           "Tectariaceae",
           "Thelypteridaceae",
           "Thyrsopteridaceae",
           "Woodsiaceae")

gymno <- c("Araucariaceae",
           "Boweniaceae",
           "Cephalotaxaceae",
           "Cupressaceae",
           "Cycadaceae",
           "Ephedraceae",
           "Ginkgoaceae",
           "Gnetaceae",
           "Pinaceae",
           "Podocarpaceae",
           "Taxaceae",
           "Taxodiaceae",
           "Welwitschiaceae",
           "Zamiaceae")

lyco <- c("Lycopodiaceae",
          "Isoetaceae",
          "Selaginellaceae")

st <- left_join(DTTT %>% as.data.frame(), tree.species,
                by = c("Species" = "species"))

st.out <- st %>% 
  mutate(devision = ifelse(family %in% ferns, "Fern",
                           ifelse(family %in% gymno, "Gymnosperm",
                                  ifelse(family %in% lyco, "Lycopods", "Other")))) %>% 
  group_by(Biome) %>% 
  mutate(n_biome = n()) %>% 
  ungroup() %>% 
  group_by(Biome, devision) %>% 
  summarise(n = n(),
            n_biome = mean(n_biome)) %>% 
  mutate(rel = n/n_biome*100) %>% 
  pivot_wider(id_cols = Biome, names_from = devision, values_from = rel)

DT22 <- DT %>% 
  left_join(tree.species, by = c("Species" = "species")) %>% 
  mutate(Devision =  ifelse(family %in% ferns, "Fern",
                            ifelse(family %in% gymno, "Gymnosperm",
                                   ifelse(family %in% lyco, "Lycopods", "Other")))) %>% 
  group_by(PlotObservationID, Devision) %>% 
  summarise(n = n())

DT222 <- DT22 %>% 
  left_join(DT22 %>% 
              group_by(PlotObservationID) %>% 
              summarise(nn = sum(n))) %>% 
  mutate(perc.n = n/nn*100) %>% 
  pivot_wider(id_cols = PlotObservationID, names_from = Devision, values_from = perc.n, values_fill = 0)

saveRDS(DT222, "02.data/perc-fern-sPlot.RDS")
