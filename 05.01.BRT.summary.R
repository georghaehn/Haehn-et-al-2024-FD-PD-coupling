####recommended machine: 01#####

####packages ----
library("tidyverse")
library("tidyr")
library("data.table")
library("raster")
library("ncdf4")
library("sf")
library("gbm")
library("dismo")

BRTs <- list.files("02.data/eve-calculations/BRT-calculation")

BRT.RQEP <- readRDS(paste0("02.data/eve-calculations/BRT-calculation/", BRTs[2]))

summary(BRT.RQEP)

BRT.RQEF <- readRDS(paste0("02.data/eve-calculations/BRT-calculation/", BRTs[3]))

summary(BRT.RQEF)

BRT.SES.RQEF <- readRDS(paste0("02.data/eve-calculations/BRT-calculation/", BRTs[4]))

summary(BRT.SES.RQEF)

gbm.plot(BRT.SES.RQEF, n.plots = 9)

BRT.SES.RQEP <- readRDS(paste0("02.data/eve-calculations/BRT-calculation/", BRTs[5]))

summary(BRT.SES.RQEP)

gbm.plot(BRT.SES.RQEP, n.plots = 9)

BRT.log <- readRDS(paste0("02.data/eve-calculations/BRT-calculation/", BRTs[6]))

summary(BRT.log)

GBM.sum.df <- summary(BRT.SES.RQEF) %>% as.data.frame() %>%
  mutate(., Response = "SES Functional Entropy") %>%
  bind_rows(.,
            summary(BRT.SES.RQEP) %>% as.data.frame() %>%
              mutate(., Response = "SES Phylogenetic Entropy")) %>%
  bind_rows(.,
            summary(BRT.log) %>% as.data.frame() %>%
              mutate(., Response = "log"))
            
            
GBM.sum.df.clean <- GBM.sum.df %>%
  mutate(var = gsub("Plant_recorded", "Plants recorded", var)) %>% 
  mutate(var = gsub("Biome", "Biome", var)) %>% 
  mutate(var = gsub("is.forest", "Forest or Non-Forest", var)) %>% 
  mutate(var = gsub("sBiome", "Biome", var)) %>% 
  mutate(var = gsub("stable.clim", "Climate variability after LGM", var)) %>%
  mutate(var = gsub("PC1", "PC1 - Annual precipitation", var)) %>%
  mutate(var = gsub("PC2", "PC2 - Temp. coldest quarter/month", var)) %>%
  mutate(var = gsub("PC3", "PC3 - Annual range temperature", var)) %>%
  mutate(var = gsub("PC4", "PC4 - Isothermality", var)) %>%
  mutate(var = gsub("PC5", "PC5 - Precipitation seasonality", var)) %>%
  mutate(var = gsub("plot.size", "Plot size", var)) %>% 
  mutate(var = gsub("Releve_area", "Plot size", var)) %>% 
  mutate(var = gsub("`Plants recorded`", "Plants recorded", var)) 

(p3 <- GBM.sum.df.clean %>% 
  ggplot() +
  geom_bar(data = . %>% filter(Response == "log"),
           aes(x=reorder(var, desc(var)), y=-rel.inf, fill="log"), stat = "identity", position = position_dodge2()) +
  # geom_errorbar(data = . %>% filter(!is.na(phylogenetic)) %>% arrange(phyl.group),
  #               aes(x=var, ymin=phylogenetic-sd, ymax=phylogenetic+sd), 
  #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
  geom_bar(data = . %>% filter(Response != "log"),
           aes(x=reorder(var, desc(var)), y=rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
  # geom_errorbar(data = . %>% filter(!is.na(functional)) %>% arrange(funct.group),
  #               aes(x=var, y = -functional, ymin=-functional-sd, ymax=-functional+sd), 
  #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
  scale_fill_manual(labels = c(expression(paste("log(SES.FD"[Q]/"SES.PD"[Q],")")),
                               expression(SES.FD[Q]), 
                               expression(SES.PD[Q])),
                    values = c("darkcyan", "olivedrab3", "darkorchid4")) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(y = "Relative influence [%]") +
  geom_abline(slope = 0, intercept = 10, linetype = 2, color = "black", size = 1.2) +
  geom_abline(slope = 0, intercept = -10, linetype = 2, color = "black", size = 1.2) +
  coord_flip() +
  guides(fill=guide_legend(title="Response")) +
  scale_y_continuous(breaks=c(-50,-25,0,25,50),
                     labels=c("50","25", "0", "25", "50"),
                     limits=c(-50,50)) +
  theme(
    axis.title.x = element_text(size = 44),
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 44, color = "black"),
    legend.title=element_blank(), 
    legend.text=element_text(size=40),
    legend.position = "top",
  )
)

ggsave(filename = "__Submission/Figures/04.BRT-summary.rel.png", width = 30, height = 20, units = "in", dpi = 600, p3)

(p3.1 <- GBM.sum.df.clean %>% 
    filter(Response != "log") %>%
    ggplot() +
    geom_bar(data = . %>% filter(Response == "SES Functional Entropy"),
             aes(x=reorder(var, desc(var)), y=-rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
    # geom_errorbar(data = . %>% filter(!is.na(phylogenetic)) %>% arrange(phyl.group),
    #               aes(x=var, ymin=phylogenetic-sd, ymax=phylogenetic+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    geom_bar(data = . %>% filter(Response == "SES Phylogenetic Entropy"),
             aes(x=reorder(var, desc(var)), y=rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
    # geom_errorbar(data = . %>% filter(!is.na(functional)) %>% arrange(funct.group),
    #               aes(x=var, y = -functional, ymin=-functional-sd, ymax=-functional+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    scale_fill_manual(labels = c(expression(SES.FD[Q]), 
                                 expression(SES.PD[Q])),
                      values = c("olivedrab3", "darkorchid4")) +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    labs(y = "Relative influence [%]") +
    geom_abline(slope = 0, intercept = 10, linetype = 2, color = "black", size = 1.2) +
    geom_abline(slope = 0, intercept = -10, linetype = 2, color = "black", size = 1.2) +
    coord_flip() +
    guides(fill=guide_legend(title="Response")) +
    scale_y_continuous(breaks=c(-50,-25,0,25,50),
                       labels=c("50","25", "0", "25", "50"),
                       limits=c(-50,50)) +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.text.y = element_text(size = 44, color = "black"),
      legend.title=element_blank(), 
      legend.text=element_text(size=40),
      legend.position = "top",
    )
)

ggsave(filename = "__Submission/Figures/04.BRT-summary.png", width = 30, height = 20, 
       units = "in", dpi = 600, p3.1)

remotes::install_github("njtierney/treezy")

library("treezy")

partial_plot <- function(x, vars){
  if(length(vars) > 1) {
  df_box <- list("vector", length(vars))
  for (i in (1:length(vars))) {
    df_box[[i]] <- partial_dependence(x, vars[[i]])
  }
  
  df <- dplyr::bind_rows(df_box) %>%
    mutate(variable = ifelse(variable == "plot.size", "Plot size", 
                             ifelse(variable == "stable.clim", "Climate variability", variable)))
  } else(df <- partial_dependence(x, vars[[1]]))
  df_mean <- df %>% 
    dplyr::group_by(variable) %>% 
    dplyr::summarise(mean = mean(fitted_function))
  
  ggplot2::ggplot(data = df %>%
                    dplyr::group_by(variable) %>%
                    dplyr::filter(value < quantile(value, 0.975),
                                  value > quantile(value, 0.025)), 
                  ggplot2::aes(x = value, y = fitted_function)) + 
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~variable, ncol = 2, scales = "free_x") + 
    ggplot2::geom_hline(data = df_mean,
                        ggplot2::aes(yintercept = mean), 
                        colour = "red", linetype = "dashed",
                        alpha = 0.75)  +
    ggplot2::labs(x = "Variable Values",
                  y = "Model Predicted Values") +
    theme(
      axis.title.x = element_text(size = 40),
      axis.text.x = element_text(size = 40, color = "black"),
      axis.title.y = element_text(size = 40),
      axis.text.y = element_text(size = 36),
      legend.text = element_text(size = 40), 
      strip.text = element_text(size = 50, face = "bold")
    ) 
  
}

PD <- partial_plot(BRT.SES.RQEP,
                var = c("PC1"))

ggsave("__Submission/Figures/04.BRT-par-SES.PD.png", PD, 
       height=15, width=20, units="in", dpi=600, bg = "white")

# df <- partial_dependence(BRT.SES.RQEP, "is.forest") %>%
#   mutate(f = ifelse(value == 1, "Non-Forest", "Forest")) %>%
#   mutate(variable = "Vegetation")
# 
# df_mean <- df %>% 
#   dplyr::group_by(variable) %>% 
#   dplyr::summarise(mean = mean(fitted_function))
# 
# Fo <- ggplot(data = df, aes(x = f, y = fitted_function)) + 
#   geom_point() +
#   facet_wrap(~variable, ncol = 2, scales = "free_x") + 
#   geom_hline(data = df_mean,
#                aes(yintercept = mean), 
#                colour = "red", linetype = "dashed",
#                alpha = 0.75)  +
#   labs(x = "Variable Values",
#        y = "")


FD <- partial_plot(x = BRT.SES.RQEF,
                vars = c(
                  "stable.clim",
                  "PC2",
                  "PC3",
                  "PC5",
                  "plot.size")) + 
  theme(
    axis.title.x = element_text(size = 70),
    axis.text.x = element_text(size = 50, color = "black"),
    axis.title.y = element_text(size = 70),
    axis.text.y = element_text(size = 50), 
    strip.text = element_text(size = 50, face = "bold")
  )

ggsave("__Submission/Figures/04.BRT-par-SES.FD.png", FD, 
       height=45, width=30, units="in", dpi=600, bg = "white")

FDPD <- partial_plot(x = BRT.log,
                   vars = c(
                     "PC1")) + 
  theme(
    axis.title.x = element_text(size = 70),
    axis.text.x = element_text(size = 50, color = "black"),
    axis.title.y = element_text(size = 70),
    axis.text.y = element_text(size = 50), 
    strip.text = element_text(size = 50, face = "bold")
  )

ggsave("__Submission/Figures/04.BRT-par-SES.FDPD.log.png", FDPD, 
       height=15, width=20, units="in", dpi=600, bg = "white")
