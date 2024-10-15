#### recommended machine: 01#####

suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("picante"))
suppressPackageStartupMessages(library("phytools"))
suppressPackageStartupMessages(library("doSNOW"))
suppressPackageStartupMessages(library("ggridges"))
suppressPackageStartupMessages(library("patchwork"))

#### load data ----
suppressWarnings(load("02.data/PhyloTree_sPlot3.0.RData"))
suppressWarnings(load("02.data/Traits_CWMs_sPlot3.RData"))
suppressWarnings(load("02.data/pruned.trees.traitsp.RData"))

species.trait <- sPlot.traits %>%
  dplyr::select(
    Species,
    GrowthForm,
    LeafN_mean,
    LeafP_mean,
    SLA_mean,
    PlantHeight_mean,
    SeedMass_mean,
    LDMC_mean,
    SpecificRootLength_mean,
    Chromosome.n_mean
  ) %>%
  filter(complete.cases(cur_data()))

species.both <- species.trait %>%
  filter(Species %in% (tree.species %>%
    filter(status != "fail to bind") %>%
    pull(species))) %>%
  mutate(Species = gsub(" ", "_", Species))

rm(CWM, sPlot.traits, try.combined.means)

# prune for each family
pruning <- function(myspecies, tree.to.prune) {
  spp1 <- matrix(
    data = 1,
    nrow = 1,
    ncol = length(unique(myspecies)),
    dimnames = list(
      "mydataset",
      unique(myspecies)
    )
  )

  pruned.tree <- prune.sample(spp1, tree.to.prune)
  return(pruned.tree)
}

PD.fam <- tree.species %>%
  mutate(species = gsub(" ", "_", species)) %>%
  filter(species %in% species.both$Species) %>%
  group_by(family) %>%
  dplyr::summarize(pruned.tree = list(pruning(
    myspecies = species,
    tree.to.prune = sPlot.phylo.tree$scenario.3
  )))

traits <- c( # "GrowthForm",
  "LeafN_mean",
  "LeafP_mean",
  "SLA_mean",
  "PlantHeight_mean",
  "SeedMass_mean",
  "LDMC_mean",
  "SpecificRootLength_mean",
  "Chromosome.n_mean"
)

sig <- list()

threads <- 16

cl <- makeCluster(threads)
registerDoSNOW(cl)

result <- foreach(j = 1:length(traits)) %dopar% {
  for (i in 1:nrow(PD.fam)) {
    mytree <- PD.fam$pruned.tree[[i]]

    trait <- species.both %>%
      dplyr::mutate(Species = gsub(" ", "_", Species)) %>%
      dplyr::filter(Species %in% mytree$tip.label) %>%
      tibble::column_to_rownames("Species") %>%
      dplyr::select(traits[j]) %>%
      as.matrix()

    tryCatch(
      {
        sig[[i]] <- data.frame(
          lambda = phytools::phylosig(mytree, trait, method = "lambda", test = T)$lambda,
          family = PD.fam[i, "family"],
          Response = traits[j]
        )
      },
      error = function(e) {
        sig[[i]] <- NA
      }
    )
  }

  sig
}
stopCluster(cl)

names(result) <- traits

out <- lapply(result, function(x) {
  bind_rows(x)
})

out1 <- bind_rows(out)

df <- out1 %>%
  mutate(
    Response = gsub("Chromosome.n_mean", "Chromosome n", Response),
    Response = gsub("PlantHeight_mean", "Plant height", Response),
    Response = gsub("LDMC_mean", "LDMC", Response),
    Response = gsub("LeafN_mean", "Leaf N", Response),
    Response = gsub("LeafP_mean", "Leaf P", Response),
    Response = gsub("SpecificRootLength_mean", "SRL", Response),
    Response = gsub("SLA_mean", "SLA", Response),
    Response = gsub("SeedMass_mean", "Seed mass", Response)
  )



cols <- c(
  "Chromosome n" = "mediumvioletred",
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


# (p.fam <- df %>%
#     ggplot(aes(x=lambda, group=Response, fill=Response, color = Response, y=after_stat(scaled))) +
#     geom_density_ (adjust=.7, alpha=.4, linewidth = 1.3) +
#     theme_bw() +
#     scale_fill_manual(values = cols) +
#     scale_color_manual(values = cols) +
#     geom_vline(xintercept = 1, color = "darkred", lty = 5, linewidth = 1.2) +
#     theme(
#       axis.title.x = element_blank(),
#       axis.text.x = element_text(size = 36),
#       axis.title.y = element_blank(),
#       axis.text.y = element_text(size = 36),
#       plot.title = element_text(size = 48, hjust = .5),
#       legend.position = "top",
#       legend.title = element_blank(),
#       legend.text = element_text(size = 36)
#     ) +
#     labs(y = "Density", x = expression(paste("Pagel`s ", lambda)), title = "Family Phylogeny") +
#     annotate(x = 1, y = +Inf ,label = "Brownian motion", vjust = 1.2, geom = "label", color = "darkred", size = 10, label.size = NA) +
#     scale_y_continuous(breaks = c(0, .5, 1), limits = c(0, 1.1), labels = c("0", .5, "1.0"), expand = c(0, 0)) +
#     scale_x_continuous(breaks = c(2, 4, 6), limits = c(0, 6), labels = c(2, 4, 6), expand = c(0, 0))
# )
df[nrow(df) + 1, "Response"] <- "All traits"
df[nrow(df) + 1, "Response"] <- "Z"
(p.l.fam <- ggplot(df, aes(x = lambda, y = Response, group = Response)) +
  geom_density_ridges(aes(fill = as.factor(Response)), alpha = 0.85, bandwidth = 0.1) +
  theme_bw() +
  scale_fill_manual(
    values = cols,
    guide = guide_legend(reverse = TRUE),
    drop = FALSE
  ) +
  geom_vline(xintercept = 1, color = "darkred", lty = 5, linewidth = 1.2) +
  theme(
    axis.title.x = element_text(size = 36, color = "black"),
    axis.text.x = element_text(size = 32),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 36, color = "black"),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 48, hjust = .5),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 36)
  ) +
  labs(y = "Density", x = "Pagel`s λ") +
  annotate(x = 1.45, y = +Inf, vjust = 1, label = "Brownian motion", geom = "label", color = "darkred", size = 10, label.size = NA) +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(0, 2.2), labels = c(0, 1, 2), expand = c(0, 0))
)

# Blombergs K
sigK <- list()

cl <- makeCluster(threads)
registerDoSNOW(cl)

result <- foreach(j = 1:length(traits)) %dopar% {
  for (i in 1:nrow(PD.fam)) {
    mytree <- PD.fam$pruned.tree[[i]]

    trait <- species.both %>%
      dplyr::mutate(Species = gsub(" ", "_", Species)) %>%
      dplyr::filter(Species %in% mytree$tip.label) %>%
      tibble::column_to_rownames("Species") %>%
      dplyr::select(traits[j]) %>%
      as.matrix()

    tryCatch(
      {
        sigK[[i]] <- data.frame(
          K = phytools::phylosig(mytree, trait, method = "K", test = T)$K,
          family = PD.fam[i, "family"],
          Response = traits[j]
        )
      },
      error = function(e) {
        sigK[[i]] <- NA
      }
    )
  }

  sigK
}
stopCluster(cl)

names(result) <- traits

outK <- lapply(result, function(x) {
  bind_rows(x)
})

outK1 <- bind_rows(outK) %>%
  filter(!is.na(family))

dfK <- outK1 %>%
  mutate(
    Response = gsub("Chromosome.n_mean", "Chromosome n", Response),
    Response = gsub("PlantHeight_mean", "Plant height", Response),
    Response = gsub("LDMC_mean", "LDMC", Response),
    Response = gsub("LeafN_mean", "Leaf N", Response),
    Response = gsub("LeafP_mean", "Leaf P", Response),
    Response = gsub("SpecificRootLength_mean", "SRL", Response),
    Response = gsub("SLA_mean", "SLA", Response),
    Response = gsub("SeedMass_mean", "Seed mass", Response)
  )

dfK[nrow(dfK) + 1, "Response"] <- "All traits"
dfK[nrow(dfK) + 1, "Response"] <- "Z"

(p.K.fam <- ggplot(dfK, aes(x = K, y = Response, group = Response)) +
  geom_density_ridges(aes(fill = as.factor(Response)), alpha = 0.85, bandwidth = 0.1) +
  theme_bw() +
  scale_fill_manual(
    values = cols,
    guide = guide_legend(reverse = TRUE),
    drop = FALSE
  ) +
  geom_vline(xintercept = 1, color = "darkred", lty = 5, linewidth = 1.2) +
  theme(
    axis.title.x = element_text(size = 36, color = "black"),
    axis.text.x = element_text(size = 32),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 48, hjust = .5),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 36)
  ) +
  labs(y = "Density", x = expression(paste("Blomberg`s ", italic("K")))) +
  annotate(x = 1.45, y = +Inf, vjust = 1, label = "Brownian motion", geom = "label", color = "darkred", size = 10, label.size = NA) +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(0, 2.2), labels = c(0, 1, 2), expand = c(0, 0))
)


# GAM oa

preds <- readRDS("02.data/04.oa.preds-PA-GAM.RDS")

size <- data.frame(
  Response = c(
    "Chromosome n",
    "Plant height",
    "LDMC",
    "Leaf N",
    "Leaf P",
    "SRL",
    "SLA",
    "Seed mass",
    "All traits"
  ),
  Width = c(rep(3, times = 8), 3.5)
)

RQE <- c(
  "FDQ.ALL.PA",
  "FDQ.MULTI.PA",
  "FDQ.SLA.PA",
  "FDQ.HEIGHT.PA",
  "FDQ.ROOT.PA",
  "FDQ.SM.PA",
  "FDQ.LDMC.PA",
  "FDQ.N.PA",
  "FDQ.P.PA",
  "FDQ.CHRO.PA"
)


df.preds <- bind_rows(preds) %>%
  mutate(Response = rep(RQE, each = 31)) %>%
  mutate(
    Response = gsub("FDQ.CHRO.PA", "Chromosome n", Response),
    Response = gsub("FDQ.GF.PA", "Growthform", Response),
    Response = gsub("FDQ.HEIGHT.PA", "Plant height", Response),
    Response = gsub("FDQ.LDMC.PA", "LDMC", Response),
    Response = gsub("FDQ.N.PA", "Leaf N", Response),
    Response = gsub("FDQ.P.PA", "Leaf P", Response),
    Response = gsub("FDQ.ROOT.PA", "SRL", Response),
    Response = gsub("FDQ.SLA.PA", "SLA", Response),
    Response = gsub("FDQ.SM.PA", "Seed mass", Response),
    Response = gsub("FDQ.ALL.PA", "All traits", Response)
  ) %>%
  filter(!Response %in% c("FDQ.MULTI.PA", "FDQ.SUB.PA", "Growthform")) %>%
  left_join(size)

linewidth <- c(
  "Chromosome n" = 2,
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
  geom_line(data = df.preds, aes(x = PDQ.PA, y = estimate, color = Response, size = Response)) +
  labs(
    y = expression("FD"[Q]),
    x = expression("PD"[Q]),
    colour = expression(paste("FD"[Q], " based on"))
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(size = 36),
    axis.title.y = element_text(size = 40),
    axis.text.y = element_text(size = 36),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 40),
    plot.title = element_text(size = 46, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  scale_color_manual(values = cols) +
  scale_size_manual(values = linewidth)
)


design <- "
1112
3344
"

(p.out <- p.oa + guide_area() + p.l.fam + theme(axis.text.y = element_blank()) +
  p.K.fam +
  plot_layout(design = design, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 40))
)

ggsave("__Submission/Figures/oa-Phylosig.png", p.out,
  height = 15, width = 25, units = "in", dpi = 600, bg = "white"
)
