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
library("patchwork")
library("scales")

# BRTs <- list.files("02.data/eve-calculations/BRT-calculation")
# 
# BRT.RQEP <- readRDS(paste0("02.data/eve-calculations/BRT-calculation/", BRTs[2]))
# 
# summary(BRT.RQEP)
# 
# BRT.RQEF <- readRDS(paste0("02.data/eve-calculations/BRT-calculation/", BRTs[3]))
# 
# summary(BRT.RQEF)
# 
# BRT.SES.RQEF <- readRDS(paste0("02.data/eve-calculations/BRT-calculation/", BRTs[4]))
# 
# summary(BRT.SES.RQEF)
# 
# gbm.plot(BRT.SES.RQEF, n.plots = 9)
# 
# BRT.SES.RQEP <- readRDS(paste0("02.data/eve-calculations/BRT-calculation/", BRTs[5]))
# 
# summary(BRT.SES.RQEP)
# 
# gbm.plot(BRT.SES.RQEP, n.plots = 9)
# 
# BRT.log <- readRDS(paste0("02.data/eve-calculations/BRT-calculation/", BRTs[6]))
# 
# summary(BRT.log)
# 
# GBM.sum.df <- summary(BRT.SES.RQEF) %>% as.data.frame() %>%
#   mutate(., Response = "SES Functional Entropy") %>%
#   bind_rows(.,
#             summary(BRT.SES.RQEP) %>% as.data.frame() %>%
#               mutate(., Response = "SES Phylogenetic Entropy")) %>%
#   bind_rows(.,
#             summary(BRT.log) %>% as.data.frame() %>%
#               mutate(., Response = "log"))
#             
#             
# GBM.sum.df.clean <- GBM.sum.df %>%
#   mutate(var = gsub("Plant_recorded", "Plants recorded", var)) %>% 
#   mutate(var = gsub("Biome", "Biome", var)) %>% 
#   mutate(var = gsub("is.forest", "Forest or Non-Forest", var)) %>% 
#   mutate(var = gsub("sBiome", "Biome", var)) %>% 
#   mutate(var = gsub("stable.clim", "Climate variability after LGM", var)) %>%
#   mutate(var = gsub("PC1", "PC1 - Annual precipitation", var)) %>%
#   mutate(var = gsub("PC2", "PC2 - Temp. coldest quarter/month", var)) %>%
#   mutate(var = gsub("PC3", "PC3 - Annual range temperature", var)) %>%
#   mutate(var = gsub("PC4", "PC4 - Isothermality", var)) %>%
#   mutate(var = gsub("PC5", "PC5 - Precipitation seasonality", var)) %>%
#   mutate(var = gsub("plot.size", "Plot size", var)) %>% 
#   mutate(var = gsub("Releve_area", "Plot size", var)) %>% 
#   mutate(var = gsub("`Plants recorded`", "Plants recorded", var)) 
# 
# (p3 <- GBM.sum.df.clean %>% 
#   ggplot() +
#   geom_bar(data = . %>% filter(Response == "log"),
#            aes(x=reorder(var, desc(var)), y=rel.inf, fill="log"), stat = "identity", position = position_dodge2()) +
#   # geom_errorbar(data = . %>% filter(!is.na(phylogenetic)) %>% arrange(phyl.group),
#   #               aes(x=var, ymin=phylogenetic-sd, ymax=phylogenetic+sd), 
#   #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
#   # geom_bar(data = . %>% filter(Response != "log"),
#   #          aes(x=reorder(var, desc(var)), y=rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
#   # geom_errorbar(data = . %>% filter(!is.na(functional)) %>% arrange(funct.group),
#   #               aes(x=var, y = -functional, ymin=-functional-sd, ymax=-functional+sd), 
#   #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
#   scale_fill_manual(labels = c(expression(paste("log(SES.FD"[Q]/"SES.PD"[Q],")")),
#                                expression(SES.FD[Q]), 
#                                expression(SES.PD[Q])),
#                     values = c("darkcyan", "olivedrab3", "darkorchid4")) +
#   theme_bw() +
#   theme(axis.title.y = element_blank()) +
#   labs(y = "Relative influence [%]") +
#   geom_abline(slope = 0, intercept = 10, linetype = 2, color = "black", size = 1.2) +
#   #geom_abline(slope = 0, intercept = -10, linetype = 2, color = "black", size = 1.2) +
#   coord_flip() +
#   guides(fill=guide_legend(title="Response")) +
#   scale_y_continuous(breaks=c(0,25,50),
#                      labels=c("0", "25", "50"),
#                      limits=c(0,50)) +
#   theme(
#     axis.title.x = element_text(size = 44),
#     axis.text.x = element_text(size = 40),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     legend.title=element_blank(), 
#     legend.text=element_text(size=40),
#     legend.position = "top",
#   )
# )
# 
# #ggsave(filename = "__Submission/Figures/04.BRT-summary.rel.png", width = 30, height = 20, units = "in", dpi = 600, p3)
# 
# (p3.1 <- GBM.sum.df.clean %>% 
#     filter(Response != "log") %>%
#     ggplot() +
#     geom_bar(data = . %>% filter(Response == "SES Functional Entropy"),
#              aes(x=reorder(var, desc(var)), y=-rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
#     # geom_errorbar(data = . %>% filter(!is.na(phylogenetic)) %>% arrange(phyl.group),
#     #               aes(x=var, ymin=phylogenetic-sd, ymax=phylogenetic+sd), 
#     #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
#     geom_bar(data = . %>% filter(Response == "SES Phylogenetic Entropy"),
#              aes(x=reorder(var, desc(var)), y=rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
#     # geom_errorbar(data = . %>% filter(!is.na(functional)) %>% arrange(funct.group),
#     #               aes(x=var, y = -functional, ymin=-functional-sd, ymax=-functional+sd), 
#     #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
#     scale_fill_manual(labels = c(expression(SES.FD[Q]), 
#                                  expression(SES.PD[Q])),
#                       values = c("olivedrab3", "darkorchid4")) +
#     theme_bw() +
#     theme(axis.title.y = element_blank()) +
#     labs(y = "Relative influence [%]") +
#     geom_abline(slope = 0, intercept = 10, linetype = 2, color = "black", size = 1.2) +
#     geom_abline(slope = 0, intercept = -10, linetype = 2, color = "black", size = 1.2) +
#     coord_flip() +
#     guides(fill=guide_legend(title="Response")) +
#     scale_y_continuous(breaks=c(-50,-25,0,25,50),
#                        labels=c("50","25", "0", "25", "50"),
#                        limits=c(-50,50)) +
#     theme(
#       axis.title.x = element_text(size = 44),
#       axis.text.x = element_text(size = 40),
#       axis.text.y = element_text(size = 44, color = "black"),
#       legend.title=element_blank(), 
#       legend.text=element_text(size=40),
#       legend.position = "top",
#     )
# )
# 
# 
# pp <- p3.1 + p3 + #plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = "A")  & 
#   theme(plot.tag = element_text(size = 40),
#         legend.position = "top")
# 
# 
# ggsave(filename = "__Submission/Figures/24-03_04.BRT-summary.png", width = 30, height = 15, 
#        units = "in", dpi = 600, pp)
# 
# remotes::install_github("njtierney/treezy")
# 
# library("treezy")
# 
# partial_plot <- function(x, vars){
#   if(length(vars) > 1) {
#   df_box <- list("vector", length(vars))
#   for (i in (1:length(vars))) {
#     df_box[[i]] <- partial_dependence(x, vars[[i]])
#   }
#   
#   df <- dplyr::bind_rows(df_box) %>%
#     mutate(variable = ifelse(variable == "plot.size", "Plot size", 
#                              ifelse(variable == "stable.clim", "Climate variability", variable)))
#   } else(df <- partial_dependence(x, vars[[1]]))
#   df_mean <- df %>% 
#     dplyr::group_by(variable) %>% 
#     dplyr::summarise(mean = mean(fitted_function))
#   
#   ggplot2::ggplot(data = df %>%
#                     dplyr::group_by(variable) %>%
#                     dplyr::filter(value < quantile(value, 0.975),
#                                   value > quantile(value, 0.025)), 
#                   ggplot2::aes(x = value, y = fitted_function)) + 
#     ggplot2::geom_line() +
#     ggplot2::facet_wrap(~variable, ncol = 2, scales = "free_x") + 
#     ggplot2::geom_hline(data = df_mean,
#                         ggplot2::aes(yintercept = mean), 
#                         colour = "red", linetype = "dashed",
#                         alpha = 0.75)  +
#     ggplot2::labs(x = "Variable Values",
#                   y = "Model Predicted Values") +
#     theme(
#       axis.title.x = element_text(size = 40),
#       axis.text.x = element_text(size = 40, color = "black"),
#       axis.title.y = element_text(size = 40),
#       axis.text.y = element_text(size = 36),
#       legend.text = element_text(size = 40), 
#       strip.text = element_text(size = 50, face = "bold")
#     ) 
#   
# }
# 
# PD <- partial_plot(BRT.SES.RQEP,
#                 var = c("PC1"))
# 
# ggsave("__Submission/Figures/04.BRT-par-SES.PD.png", PD, 
#        height=15, width=20, units="in", dpi=600, bg = "white")
# 
# # df <- partial_dependence(BRT.SES.RQEP, "is.forest") %>%
# #   mutate(f = ifelse(value == 1, "Non-Forest", "Forest")) %>%
# #   mutate(variable = "Vegetation")
# # 
# # df_mean <- df %>% 
# #   dplyr::group_by(variable) %>% 
# #   dplyr::summarise(mean = mean(fitted_function))
# # 
# # Fo <- ggplot(data = df, aes(x = f, y = fitted_function)) + 
# #   geom_point() +
# #   facet_wrap(~variable, ncol = 2, scales = "free_x") + 
# #   geom_hline(data = df_mean,
# #                aes(yintercept = mean), 
# #                colour = "red", linetype = "dashed",
# #                alpha = 0.75)  +
# #   labs(x = "Variable Values",
# #        y = "")
# 
# 
# FD <- partial_plot(x = BRT.SES.RQEF,
#                 vars = c(
#                   "stable.clim",
#                   "PC2",
#                   "PC3",
#                   "PC5",
#                   "plot.size")) + 
#   theme(
#     axis.title.x = element_text(size = 70),
#     axis.text.x = element_text(size = 50, color = "black"),
#     axis.title.y = element_text(size = 70),
#     axis.text.y = element_text(size = 50), 
#     strip.text = element_text(size = 50, face = "bold")
#   )
# 
# ggsave("__Submission/Figures/04.BRT-par-SES.FD.png", FD, 
#        height=45, width=30, units="in", dpi=600, bg = "white")
# 
# FDPD <- partial_plot(x = BRT.log,
#                    vars = c(
#                      "PC1")) + 
#   theme(
#     axis.title.x = element_text(size = 70),
#     axis.text.x = element_text(size = 50, color = "black"),
#     axis.title.y = element_text(size = 70),
#     axis.text.y = element_text(size = 50), 
#     strip.text = element_text(size = 50, face = "bold")
#   )
# 
# ggsave("__Submission/Figures/04.BRT-par-SES.FDPD.log.png", FDPD, 
#        height=15, width=20, units="in", dpi=600, bg = "white")


#### FD - PD by Biome weighted

BRT.BW <- readRDS("02.data/04.BRT-results-SES.BW.RDS")
BRT.BW.rel <- readRDS("02.data/04.BRT-results-SES.BW.rel.log.RDS")

BRT.sum <- summary(BRT.BW[[1]]) %>%
  as.data.frame() %>% 
  mutate(Response = "SES.FD",
         Direction = ifelse(var == "stable.clim", "+",
                            ifelse(var == "PC2", "+/-",
                                   ifelse(var == "PC5", "-", NA)))) %>% 
  bind_rows(., summary(BRT.BW[[2]]) %>%
              as.data.frame() %>% 
              mutate(Response = "SES.PD",
                     Direction = ifelse(var == "PC1", "+",
                                        ifelse(var == "is.forest", "-", NA)))) %>% 
  bind_rows(., summary(BRT.BW.rel) %>%
              as.data.frame() %>% 
              mutate(Response = "log",
                     Direction = ifelse(var == "PC1", "-/+",
                                        ifelse(var == "is.forest", "+", NA)))) %>%
  mutate(var = gsub("Plant_recorded", "Plants recorded", var)) %>% 
  mutate(var = gsub("is.forest", "Forest or Non-Forest", var)) %>% 
  mutate(var = gsub("stable.clim", "Climate variability after LGM", var)) %>%
  mutate(var = gsub("PC1", "PC1 - Annual precipitation", var)) %>%
  mutate(var = gsub("PC2", "PC2 - Temp. coldest quarter/month", var)) %>%
  mutate(var = gsub("PC3", "PC3 - Annual range temperature", var)) %>%
  mutate(var = gsub("PC4", "PC4 - Isothermality", var)) %>%
  mutate(var = gsub("PC5", "PC5 - Precipitation seasonality", var)) %>%
  mutate(var = gsub("`Plants recorded`", "Plants recorded", var)) 


(p3 <- ggplot(data = BRT.sum %>% filter(Response == "log"),
           aes(x=reorder(var, desc(var)), y=rel.inf, fill="log")) +
    geom_bar(stat = "identity", position = position_dodge2()) +
    geom_text(aes(label = Direction), size = 12, vjust = 0.5, hjust = -0.5) +
    # geom_errorbar(data = . %>% filter(!is.na(phylogenetic)) %>% arrange(phyl.group),
    #               aes(x=var, ymin=phylogenetic-sd, ymax=phylogenetic+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    # geom_bar(data = . %>% filter(Response != "log"),
    #          aes(x=reorder(var, desc(var)), y=rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
    # geom_errorbar(data = . %>% filter(!is.na(functional)) %>% arrange(funct.group),
    #               aes(x=var, y = -functional, ymin=-functional-sd, ymax=-functional+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    scale_fill_manual(labels = c(expression(paste("log(SES.FD"[Q]/"SES.PD"[Q],")")),
                                 expression(SES.FD[Q]), 
                                 expression(SES.PD[Q])),
                      values = c("darkcyan", "olivedrab3", "darkorchid4")) +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    labs(y = "Relative influence [%]",
         tag = "C") +
    geom_abline(slope = 0, intercept = 12.5, linetype = 2, color = "black", size = 1.2) +
    #geom_abline(slope = 0, intercept = -10, linetype = 2, color = "black", size = 1.2) +
    coord_flip() +
    guides(fill=guide_legend(title="Response")) +
    scale_y_continuous(breaks=c(0,25,50),
                       labels=c("0", "25", "50"),
                       limits=c(0,50)) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 40),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title=element_blank(), 
      legend.text=element_text(size=40),
      legend.position = "top",
      plot.tag = element_text(size = 40),
      plot.tag.position = "topleft"
    )
)

#ggsave(filename = "__Submission/Figures/04.BRT-summary.rel.png", width = 30, height = 20, units = "in", dpi = 600, p3)
(p3.0 <- ggplot(data = BRT.sum %>% filter(Response == "SES.FD"),
           aes(x=reorder(var, desc(var)), y=rel.inf, fill="SES.FD")) +
    geom_bar(stat = "identity", position = position_dodge2()) +
    geom_text(aes(label = Direction), size = 12, vjust = 0.5, hjust = -0.5) +
    # geom_errorbar(data = . %>% filter(!is.na(phylogenetic)) %>% arrange(phyl.group),
    #               aes(x=var, ymin=phylogenetic-sd, ymax=phylogenetic+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    # geom_bar(data = . %>% filter(Response != "log"),
    #          aes(x=reorder(var, desc(var)), y=rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
    # geom_errorbar(data = . %>% filter(!is.na(functional)) %>% arrange(funct.group),
    #               aes(x=var, y = -functional, ymin=-functional-sd, ymax=-functional+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    scale_fill_manual(labels = expression(SES.FD[Q]),
                      values = c("olivedrab3", "darkorchid4")) +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    labs(y = "Relative influence [%]",
         tag = "A") +
    geom_abline(slope = 0, intercept = 12.5, linetype = 2, color = "black", size = 1.2) +
    #geom_abline(slope = 0, intercept = -10, linetype = 2, color = "black", size = 1.2) +
    coord_flip() +
    guides(fill=guide_legend(title="Response")) +
    scale_y_continuous(breaks=c(0,25,50),
                       labels=c("0", "25", "50"),
                       limits=c(0,50)) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 40),
      #axis.text.y = element_text(size = 44, color = "black"),
      legend.title=element_blank(), 
      legend.text=element_text(size=40),
      legend.position = "top",
      plot.tag = element_text(size = 40),
      plot.tag.position = "topleft",
      axis.text.y = element_blank()
    )
)

(p3.01 <- ggplot(data = BRT.sum %>% filter(Response == "SES.PD"),
           aes(x=reorder(var, desc(var)), y=rel.inf, fill="SES.PD")) +
    geom_bar(stat = "identity", position = position_dodge2()) +
    geom_text(aes(label = Direction), size = 12, vjust = 0.5, hjust = -0.5) +
    # geom_errorbar(data = . %>% filter(!is.na(phylogenetic)) %>% arrange(phyl.group),
    #               aes(x=var, ymin=phylogenetic-sd, ymax=phylogenetic+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    # geom_bar(data = . %>% filter(Response != "log"),
    #          aes(x=reorder(var, desc(var)), y=rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
    # geom_errorbar(data = . %>% filter(!is.na(functional)) %>% arrange(funct.group),
    #               aes(x=var, y = -functional, ymin=-functional-sd, ymax=-functional+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    scale_fill_manual(labels = expression(SES.PD[Q]),
                      values = c("darkorchid4")) +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    labs(y = "Relative influence [%]",
         tag = "B") +
    geom_abline(slope = 0, intercept = 12.5, linetype = 2, color = "black", size = 1.2) +
    #geom_abline(slope = 0, intercept = -10, linetype = 2, color = "black", size = 1.2) +
    coord_flip() +
    guides(fill=guide_legend(title="Response")) +
    scale_y_continuous(breaks=c(0,25,50),
                       labels=c("0", "25", "50"),
                       limits=c(0,50)) +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title=element_blank(), 
      legend.text=element_text(size=40),
      legend.position = "top",
      plot.tag = element_text(size = 40),
      plot.tag.position = "topleft"
    )
)
(p3.1 <- BRT.sum %>% 
    filter(Response != "log") %>%
    ggplot() +
    geom_bar(data = . %>% filter(Response == "SES.FD"),
             aes(x=reorder(var, desc(var)), y=-rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
    # geom_errorbar(data = . %>% filter(!is.na(phylogenetic)) %>% arrange(phyl.group),
    #               aes(x=var, ymin=phylogenetic-sd, ymax=phylogenetic+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    geom_bar(data = . %>% filter(Response == "SES.PD"),
             aes(x=reorder(var, desc(var)), y=rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
    # geom_errorbar(data = . %>% filter(!is.na(functional)) %>% arrange(funct.group),
    #               aes(x=var, y = -functional, ymin=-functional-sd, ymax=-functional+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    scale_fill_manual(labels = c(expression(SES.FD[Q]), 
                                 expression(SES.PD[Q])),
                      values = c("olivedrab3", "darkorchid4")) +
    labs(tag = "A") +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    labs(y = "Relative influence [%]") +
    geom_abline(slope = 0, intercept = 12.5, linetype = 2, color = "black", size = 1.2) +
    geom_abline(slope = 0, intercept = -12.5, linetype = 2, color = "black", size = 1.2) +
    coord_flip() +
    guides(fill=guide_legend(title="Response")) +
    scale_y_continuous(breaks=c(-50,-25,0,25,50),
                       labels=c("50","25", "0", "25", "50"),
                       limits=c(-50,50)) +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      #axis.text.y = element_text(size = 44, color = "black"),
      legend.title=element_blank(), 
      legend.text=element_text(size=40),
      legend.position = "top",
      plot.tag = element_text(size = 40),
      plot.tag.position = "topleft",
      axis.text.y = element_blank()
    )
)

(p4 <- BRT.sum %>% 
    filter(Response != "log") %>%
    ggplot() +
    geom_bar(data = . %>% filter(Response == "SES.FD"),
             aes(x=reorder(var, desc(var)), y=-rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
    # geom_errorbar(data = . %>% filter(!is.na(phylogenetic)) %>% arrange(phyl.group),
    #               aes(x=var, ymin=phylogenetic-sd, ymax=phylogenetic+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    geom_bar(data = . %>% filter(Response == "SES.PD"),
             aes(x=reorder(var, desc(var)), y=rel.inf, fill=Response), stat = "identity", position = position_dodge2()) +
    # geom_errorbar(data = . %>% filter(!is.na(functional)) %>% arrange(funct.group),
    #               aes(x=var, y = -functional, ymin=-functional-sd, ymax=-functional+sd), 
    #               position = position_dodge2(width = .5, padding = .5), colour="black", alpha=0.9, size=1) +
    scale_fill_manual(labels = c(expression(SES.FD[Q]), 
                                 expression(SES.PD[Q])),
                      values = c("olivedrab3", "darkorchid4")) +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    geom_abline(slope = 0, intercept = 12.5, linetype = 2, color = "black", size = 1.2) +
    geom_abline(slope = 0, intercept = -12.5, linetype = 2, color = "black", size = 1.2) +
    coord_flip() +
    guides(fill=guide_legend(title="Response")) +
    scale_y_continuous(breaks=c(0),
                       labels=c(""),
                       limits=c(0,0)) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      axis.text.x = element_text(size = 40),
      axis.text.y = element_text(size = 44, color = "black"),
      legend.title = element_blank(), 
      legend.text = element_text(size=40),
      legend.position = "none",
      plot.margin = margin(120, 10, 50, 10)
    )
)

ps <- p4 + plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer()

pp <- p3.0 + p3.01 + p3 + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 40),
        legend.position = "top")

ppp <- p4 + inset_element(pp, left = 0.01, bottom = -0.2, right = 1, top = 1.24) 

ggsave(filename = "__Submission/Figures/05.BRT-summary-BW.png", width = 33, height = 10, 
       units = "in", dpi = 600, ppp)

ggsave(filename = "__Submission/Figures/Figure3.pdf",  width = 33, height = 10, 
       units = "in", dpi = 600, ppp)

#save source data
write.xlsx(BRT.sum, "__Submission/Source_Figure3.xlsx")

#remotes::install_github("njtierney/treezy")

library("treezy")

partial_plot <- function(x, vars){
  if(length(vars) > 1) {
    df_box <- list("vector", length(vars))
    for (i in (1:length(vars))) {
      df_box[[i]] <- partial_dependence(x, vars[[i]])
    }
    
    df <- dplyr::bind_rows(df_box) %>%
      mutate(variable = ifelse(variable == "plot.size", "Plot size", 
                               ifelse(variable == "stable.clim", "Climate variability\nafter LGM", variable)))
  } else(df <- partial_dependence(x, vars[[1]]))
  df_mean <- df %>%  
    mutate(variable = ifelse(variable == "is.forest", "Vegetation type", variable)) %>% 
    dplyr::group_by(variable) %>% 
    dplyr::summarise(mean = mean(fitted_function))
  
  ggplot2::ggplot(data = df %>%
                    filter(variable != "is.forest") %>% 
                    dplyr::group_by(variable) %>%
                    dplyr::filter(value < quantile(value, 0.975),
                                  value > quantile(value, 0.025)) %>% 
                    rbind(., df %>% 
                            filter(variable == "is.forest") %>% 
                            mutate(value = ifelse(value == 1, FALSE, TRUE),
                                   variable = "Vegetation type")) %>% 
                    filter(variable != "is.forest"), 
                  ggplot2::aes(x = value, y = fitted_function)) + 
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~variable, ncol = 2, scales = "free_x") + 
    ggplot2::geom_hline(data = df_mean,
                        ggplot2::aes(yintercept = mean), 
                        colour = "red", linetype = "dashed",
                        alpha = 0.75)  +
    ggplot2::labs(x = "",
                  y = "Model Predicted Values") +
    theme(
      axis.title.x = element_text(size = 40),
      axis.text.x = element_text(size = 40, color = "black"),
      axis.title.y = element_text(size = 40),
      axis.text.y = element_text(size = 36),
      legend.text = element_text(size = 40), 
      strip.text = element_text(size = 50, face = "bold")
    ) +
    scale_x_continuous(breaks = breaks_pretty(n = 1))
  
}

(PD <- (partial_plot(BRT.BW[[2]],
                   vars = c("PC1")) +
          scale_y_continuous(limits = c(-0.53, 2.6)) ) +
      (partial_plot(BRT.BW[[2]],
                    vars = c("is.forest")) +
         geom_line(color = "white", linetype = 6) +
         geom_point(size = 4) +
         scale_x_continuous(limits = c(-0.2, 1.2),
                            breaks = c(0, 1),
                            labels = c("Non-Forest", "Forest") ) +
         scale_y_continuous(limits = c(-0.53, 2.6)) +
         theme(
           axis.title.y = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks = element_blank()
         )))

ggsave("__Submission/Figures/04.BRT-par-SES.PD.BW.png", PD, 
       height=10, width=20, units="in", dpi=600, bg = "white")

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


FD <- partial_plot(x = BRT.BW[[1]],
                   vars = c(
                     "stable.clim",
                     "PC2",
                     "PC5"
                     )) + 
  scale_x_continuous(breaks = breaks_pretty(n = 4)) +
  theme(
    strip.text = element_text(size = 30, face = "bold")
  ) +
  facet_wrap(~variable, nrow = 1, scales = "free_x")

ggsave("__Submission/Figures/04.BRT-par-SES.FD.BW.png", FD, 
       height=10, width=30, units="in", dpi=600, bg = "white")

(FDPD <- (partial_plot(BRT.BW.rel,
                     vars = c("PC1")) +
          scale_y_continuous(limits = c(-0, 0.5)) ) +
    (partial_plot(BRT.BW.rel,
                  vars = c("is.forest")) +
       geom_line(color = "white", linetype = 6) +
       geom_point(size = 4) +
       scale_x_continuous(limits = c(-0.2, 1.2),
                          breaks = c(0, 1),
                          labels = c("Non-Forest", "Forest") ) +
       scale_y_continuous(limits = c(-0, 0.5)) +
       theme(
         axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks = element_blank()
       )))

ggsave("__Submission/Figures/04.BRT-par-SES.FDPD.log.BW.png", FDPD, 
       height=10, width=20, units="in", dpi=600, bg = "white")

ffout <- PD / FDPD +
  plot_annotation(tag_levels = list(c("A", "", "B", "")))  & 
  theme(plot.tag = element_text(size = 52))

ggsave("__Submission/Figures/04.BRT-par-SES.FD-FDPD.log.BW.png", ffout, 
       height=20, width=20, units="in", dpi=600, bg = "white")

ggsave("__Submission/Figures/04.BRT-par-SES.FD-FDPD.log.BW.pdf", ffout, 
       height=20, width=20, units="in", dpi=600, bg = "white")
