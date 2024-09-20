####recommended machine: 01#####

library("dplyr")
library("sp")
library("tidyr")
library("mgcv")
library("ggplot2")
library("ncdf4")
library("sf")
library("dggridR")
library("marginaleffects")
library("raster")
library("sf")
library("rnaturalearth")
library("viridis")
library("ggforce")
library("grid")

sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-BW-data.RDS")

#### histogram of species per plot

hist(log(sPlot.data$SR))

srp <- sPlot.data %>% 
   ggplot() +
   geom_histogram(aes(x = SR), binwidth = 4) +
   labs(x = "Species richness", y = "Count") +
   theme(
     axis.title.x = element_text(size = 40),
     axis.text.x = element_text(size = 36),
     axis.title.y = element_text(size = 40),
     axis.text.y = element_text(size = 36), 
     legend.title = element_text(size=36), 
     legend.text  =element_text(size=36),
     #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
     legend.key.height = unit(1.1, "cm"), 
     legend.key.width = unit(1.1, "cm"),
     legend.position = "none"
   ) +
   scale_y_log10()

ggsave("__Submission/Figures/00.srp.png", srp, 
       height=10, width=15, units="in", dpi=600, bg = "white")

ggsave("__Submission/Figures/00.srp.pdf", srp, 
       height=10, width=15, units="in", dpi=600, bg = "white")


###### RQEP ----
mod.RQEP <- readRDS("02.data/05.GAM_SES.RQEP.BWexp_1.RDS")
summary(mod.RQEP)

#pred.orig.RQEP <- predictions(mod.RQEP)
#saveRDS(pred.orig.RQEP, file = "02.data/05.GAM-prediction.RQEP.RDS")
#pred.orig.RQEP <- readRDS("02.data/05.GAM-prediction.RQEP.RDS")

sPlot.data <- sPlot.data %>% 
  filter(!is.na(SES.RQEP.BW),
         !is.na(SES.RQEF.BW),
         !is.na(is.forest),
         !is.na(PC1),
         !is.na(Longitude))

## PC1 ----

RQEP.PC1.preds <- predictions(mod.RQEP, newdata = datagrid(
  PC1 = seq(
    min(sPlot.data$PC1), 
    max(sPlot.data$PC1), 
    .1)
))


mod.P.PC1 <- readRDS("02.data/05.GAM_SES.RQEP.BWexp_2.RDS")
summary(mod.P.PC1)

res.P.PC1 <- residuals.gam(mod.P.PC1, type = "response")

(p.phyl.PC1 <- ggplot() +
  geom_hex(aes(y = res.P.PC1, x = sPlot.data$PC1),
           bins = 40) +
  #geom_ribbon(data = RQEP.PC1.preds, aes(ymin = conf.low, ymax = conf.high, x = PC1), outline.type = "full") +
  #geom_point(aes(y = pred.orig.RQEP$SES.RQEP-pred.orig.RQEP$predicted, x = pred.orig.RQEP$PC1), color = "grey80", alpha = 0.5) +
  geom_line(data = RQEP.PC1.preds, aes(x = PC1, y = estimate), linewidth = 2, lty = 1) +
  labs(y = expression(paste("Residuals SES.PD"[Q])), 
       x = "PC1 - Annual precipitation\n(dry to wet)") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 44),
    axis.text.x = element_text(size = 40),
    axis.title.y = element_text(size = 44),
    axis.text.y = element_text(size = 40),
    legend.title=element_text(size=36), 
    legend.text=element_text(size=36),
    legend.position = "right",
    #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
    legend.key.height = unit(1.1, "cm"), 
    legend.key.width = unit(1.1, "cm"),
    plot.title = element_text(hjust = .5, size = 70)
  ) +
    # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
    #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
    #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 100000),
                        breaks = c(10, 100, 1000, 10000)) +
  scale_x_continuous(limits = c(quantile(sPlot.data$PC1, 0.01), 
                                23)) +
  scale_y_continuous(limits = c(-9,9))
  )
  #geom_text(aes(x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, label = "A"), size = 20)

## spatial smooth ----

dggs <- dgconstruct(area=50*50, metric=T, resround='down')

sPlot.data$cell <- dgGEO_to_SEQNUM(dggs, sPlot.data$Longitude, sPlot.data$Latitude)$seqnum

mod.P.S <- readRDS("02.data/05.GAM_SES.RQEP.BWexp_4.RDS")
summary(mod.P.S)

sPlot.data$res.P.S <- residuals.gam(mod.P.S, type = "response")

d.phyl <- sPlot.data %>% 
  group_by(cell) %>% 
  summarise(smooth = mean(res.P.S, na.rm = T))

grid <- dgcellstogrid(dggs, d.phyl$cell) %>%
  st_as_sf() %>% 
  mutate(cell = d.phyl$cell) %>% 
  mutate(smooth = d.phyl$smooth) %>% 
  st_transform("+proj=eck4") %>% 
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))

world <- ne_countries(returnclass = "sf")%>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()

(p.phyl.geo <- grid %>% 
    mutate(area = as.numeric(st_area(.))) %>%
    filter(area <= quantile(area, c(0.999))) %>% 
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    # geom_sf(data = bb, col = "grey20", fill = NA) +
    coord_sf(crs = "+proj=eck4") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=36),
          legend.text=element_text(size=36),
          #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          legend.position = "bottom", 
          plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    geom_sf(aes(color = smooth, fill = smooth), lwd = 0) +
    scale_fill_viridis_c(option = "plasma",
                         #breaks = c(-4,0,4),
                         name = expression(paste("SES.PD"[Q]))) +
    scale_color_viridis_c(option = "plasma",
                          #breaks = c(-4,0,4),
                          name = expression(paste("SES.PD"[Q]))) +
    guides(
      fill = "none",
      color = guide_colourbar(
        #title = expression(paste("Higher\nSES.PD"[Q], symbol("\256"), "Higher\nSES.FD"[Q])),
        title.position="top", title.hjust = 0.5
      )))
    # guides(fill = guide_legend(title=expression(paste("Smooth SES.PD"[Q]))),
    #        color = guide_legend(title=expression(paste("Smooth SES.PD"[Q]))))
    
    # scale_fill_gradient2(low = "#91bfdb",
    #                      mid = "#ffffbf",
    #                      high = "#fc8d59",
    #                      # limits = c(-5,5),
    #                      breaks = c(-4,0,4)
    # )+
    # # midpoint = 0) +
    # scale_color_gradient2(low = "#91bfdb",
    #                       mid = "#ffffbf",
    #                       high = "#fc8d59",
    #                       breaks = c(-4,0,4)
    # )
  # guides(color = guide_colourbar(title.position="top", title.hjust = 0.5),
  #        size = guide_legend(title.position="top", title.hjust = 0.5))


## is forest ----

pred.orig.RQEP.new <- sPlot.data %>% 
  mutate(is.forest = ifelse(is.forest == "TRUE", "Forest", "Non-Forest"))

forest.data <- sPlot.data %>% 
  group_by(is.forest) %>% 
  split(.$is.forest)
        
pred.RQEP.for <- lapply(forest.data, function(data){ predictions(mod.RQEP, newdata = datagrid(
  PC1 = sample(size = 10, x = seq(min(sPlot.data$PC1), max(sPlot.data$PC1), .1), replace = T),
  is.forest = data$is.forest
))
})

pred.RQEP.for.out <- bind_rows(pred.RQEP.for)

pred.RQEP.for.out <- pred.RQEP.for.out %>% 
  mutate(is.forest = ifelse(is.forest == TRUE, "Forest", "Non-Forest"))

pred.RQEP.for1 <- predictions(mod.RQEP, newdata = datagrid(
  is.forest = rep(c(TRUE, FALSE), 100)
))

pred.RQEP.for.out <- pred.RQEP.for1 %>% 
  mutate(is.forest = ifelse(is.forest == TRUE, "Forest", "Non-Forest"))

mod.P.FF <- readRDS("02.data/05.GAM_SES.RQEP.BWexp_3.RDS")
summary(mod.P.FF)

res.P.FF <- residuals.gam(mod.P.FF, type = "response")

p.phyl.for <- ggplot() +
  geom_violin(aes(y = res.P.FF, 
                x = pred.orig.RQEP.new$is.forest), fill = "grey", color = "grey") +
  geom_boxplot(aes(x = as.factor(pred.RQEP.for.out$is.forest), y = pred.RQEP.for.out$estimate)) +
  #geom_violin(aes(x = pred.orig.RQEP$is.forest, y = pred.orig.RQEP$SES.RQEP, color = "SES.RQEP"), size = 1.5) +
  labs(y = expression(paste("Residuals SES.PD"[Q])), x = "", color = "") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 44),
    axis.text.x = element_text(size = 44, color = "black", vjust = -3.8),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 40),
    legend.text = element_text(size = 44),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(-9,9))
  #geom_text(aes(x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, label = "C"), size = 20)

library("patchwork")

layout <- "AABB"


p.out.RQEP <- p.phyl.PC1 + p.phyl.for + 
  plot_layout(design = layout, guides = "collect") +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 40))

ggsave("__Submission/Figures/05.GAM.expl.SES.RQEP.BW.png", p.out.RQEP, 
       height=10, width=25, units="in", dpi=600, bg = "white", limitsize = FALSE)

ggsave("__Submission/Figures/05.GAM.expl.SES.RQEP.BW.pdf", p.out.RQEP, 
       height=10, width=25, units="in", dpi=600, bg = "white", limitsize = FALSE)

# rm(list = ls())
# gc()
#### functional diversity ~ stable clim + recent climate + longlat ----

mod.RQEF <- readRDS("02.data/05.GAM_SES.RQEF.BWexp_5.RDS")
summary(mod.RQEF)

mod.F.SC <- readRDS("02.data/05.GAM_SES.RQEF.BWexp_6.RDS")
summary(mod.F.SC)

res.F.SC <- residuals.gam(mod.F.SC, type = "response")

# pred.orig.RQEF <- predictions(mod.RQEF)
# 
# saveRDS(pred.orig.RQEF, file = "02.data/05.GAM-prediction.RQEF.RDS")

#pred.orig.RQEF <- readRDS("02.data/05.GAM-prediction.RQEF.RDS")

#sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-PCA-data.RDS")

## stable climate ----

RQEF.sc.preds <- predictions(mod.RQEF, newdata = datagrid(
  stable.clim = seq(
    min(sPlot.data$stable.clim), 
    max(sPlot.data$stable.clim), 
    .1)
))

(p.RQEF.sc <- ggplot() +
 # geom_point(aes(y = pred.orig.RQEF$SES.RQEF-pred.orig.RQEF$predicted, x = pred.orig.RQEF$stable.clim), color = "grey80", alpha = 0.5) +
  geom_hex(aes(y = res.F.SC, 
               x = sPlot.data$stable.clim), 
           bins = 40) +
  geom_line(data = RQEF.sc.preds, aes(x = stable.clim, y = estimate), linewidth = 2, lty = 1) +
  labs(y = expression(paste("Residuals SES.FD"[Q])), 
       x = "Climate variability\nafter LGM [°C p.a.]") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(size = 36),
    axis.title.y = element_text(size = 40),
    legend.title=element_text(size=36), 
    legend.text=element_text(size=36),
    axis.text.y = element_text(size = 36),
    #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
    legend.key.height = unit(1.1, "cm"), 
    legend.key.width = unit(1.1, "cm"),
    legend.position = "right"
  ) +
  scale_x_continuous(breaks = c(0,1,2)) + 
    
    scale_y_continuous(limits = c(-6, 6)) +
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 100000),
                        breaks = c(10, 100, 1000, 10000)))

PC2.preds <- predictions(mod.RQEF, newdata = datagrid(
  PC2 = seq(
    min(sPlot.data$PC2), 
    max(sPlot.data$PC2), 
    .1)
))

mod.F.PC2 <- readRDS("02.data/05.GAM_SES.RQEF.BWexp_7.RDS")
summary(mod.F.PC2)

res.F.PC2 <- residuals.gam(mod.F.PC2, type = "response")

(p.RQEF.PC2 <- ggplot() +
  #geom_point(aes(x = pred.orig.FDis$annual.range.air.temp, y = pred.orig.FDis$SES.FDis), color = "purple4", alpha = 0.01) +
  #geom_point(aes(y = pred.orig.RQEF$SES.RQEF-pred.orig.RQEF$predicted, x = pred.orig.RQEF$PC2), color = "grey80", alpha = 0.5) +
  geom_hex(aes(y = res.F.PC2, x = sPlot.data$PC2),
           bins = 40) +
  geom_line(data = PC2.preds, aes(x = PC2, y = estimate), linewidth = 2) +
  labs(y = expression(paste("Residuals SES.FD"[Q])),
       x = "PC2 - Temp. coldest\nquarter/month (cold to hot)") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(size = 36),
    axis.title.y = element_text(size = 40),
    axis.text.y = element_text(size = 36),
    legend.title=element_text(size=36), 
    legend.text=element_text(size=36),
    legend.position = "right",
    #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
    legend.key.height = unit(1.1, "cm"), 
    legend.key.width = unit(1.1, "cm"),
    plot.title = element_text(hjust = .5, size = 70)
  ) +
  # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC2, 0.025), 
  #                               quantile(pred.orig.RQEF$PC2, 0.9775))) +
  #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
    
    scale_y_continuous(limits = c(-6, 6)) +
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 100000),
                        breaks = c(10, 100, 1000, 10000)) )

# PC3.preds <- predictions(mod.RQEF, newdata = datagrid(
#   PC3 = seq(
#     min(sPlot.data$PC3), 
#     max(sPlot.data$PC3), 
#     .05)
# ))
# 
# mod.F.PC3 <- readRDS("02.data/05.GAM_SES.RQEF.BWexp_8.RDS")
# summary(mod.F.PC3)
# 
# res.F.PC3 <- residuals.gam(mod.F.PC3, type = "response")
# (p.RQEF.PC3 <- ggplot() +
#   #geom_point(aes(x = pred.orig.FDis$annual.range.air.temp, y = pred.orig.FDis$SES.FDis), color = "purple4", alpha = 0.01) +
#   #geom_point(aes(y = pred.orig.RQEF$SES.RQEF-pred.orig.RQEF$predicted, x = pred.orig.RQEF$PC3), color = "grey80", alpha = 0.5) +
#   geom_hex(aes(y = res.F.PC3, x = pred.orig.RQEF$PC3),
#            bins = 40) +
#   geom_line(data = PC3.preds, aes(x = PC3, y = estimate), linewidth = 2) +
#   labs(y = expression(paste("Residuals SES.FD"[Q])), 
#        x = "PC3 - Annual\nrange temperature") +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 40),
#     axis.text.x = element_text(size = 36),
#     axis.title.y = element_text(size = 40),
#     axis.text.y = element_text(size = 36),
#     legend.title=element_text(size=36), 
#     legend.text=element_text(size=36),
#     legend.position = "right",
#     #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
#     legend.key.height = unit(1.1, "cm"), 
#     legend.key.width = unit(1.1, "cm"),
#     plot.title = element_text(hjust = .5, size = 70)
#   ) +
#   scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC3, 0.01), 
#                                 10)) +
#   #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
#     
#     scale_y_continuous(limits = c(-6, 6)) +
#     scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
#                         limits = c(1, 100000),
#                         breaks = c(10, 100, 1000, 10000)) )

PC5.preds <- predictions(mod.RQEF, newdata = datagrid(
  PC5 = seq(
    min(sPlot.data$PC5), 
    max(sPlot.data$PC5), 
    .1)
))

mod.F.PC5 <- readRDS("02.data/05.GAM_SES.RQEF.BWexp_8.RDS")
summary(mod.F.PC5)

res.F.PC5 <- residuals.gam(mod.F.PC5, type = "response")
(p.RQEF.PC5 <- ggplot() +
  #geom_point(aes(x = pred.orig.FDis$annual.range.air.temp, y = pred.orig.FDis$SES.FDis), color = "purple4", alpha = 0.01) +
  #geom_point(aes(y = pred.orig.RQEF$SES.RQEF-pred.orig.RQEF$predicted, x = pred.orig.RQEF$PC5), color = "grey80", alpha = 0.5) +
  geom_hex(aes(y = res.F.PC5, x = sPlot.data$PC5),
           bins = 40) +
  geom_line(data = PC5.preds, aes(x = PC5, y = estimate), linewidth = 2) +
  labs(y = expression(paste("Residuals SES.FD"[Q])), x = "PC5 - Precipitation\nseasonality (low to high)") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(size = 36),
    axis.title.y = element_text(size = 40),
    axis.text.y = element_text(size = 36),
    legend.title=element_text(size=36), 
    legend.text=element_text(size=36),
    legend.position = "right",
    #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
    legend.key.height = unit(1.1, "cm"), 
    legend.key.width = unit(1.1, "cm"),
    plot.title = element_text(hjust = .5, size = 70)
  ) +
  # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
  #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
  #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
    
    scale_y_continuous(limits = c(-6, 6)) +
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 100000),
                        breaks = c(10, 100, 1000, 10000)) )

# plot.preds <- predictions(mod.RQEF, newdata = datagrid(
#   plot.size = seq(
#     1, 
#     max(sPlot.data$plot.size, na.rm = TRUE), 
#     100)
# ))
# 
# mod.F.PP <- readRDS("02.data/05.GAM_SES.RQEF-exp_5.RDS")
# summary(mod.F.PP)
# 
# res.F.PP <- residuals.gam(mod.F.PP, type = "response")
# (p.RQEF.plot <- ggplot() +
#   #geom_point(aes(x = pred.orig.FDis$annual.range.air.temp, y = pred.orig.FDis$SES.FDis), color = "purple4", alpha = 0.01) +
#   #geom_point(aes(y = pred.orig.RQEF$SES.RQEF-pred.orig.RQEF$predicted, x = pred.orig.RQEF$plot.size), color = "grey80", alpha = 0.5) +
#   geom_hex(aes(y = res.F.PP, x = pred.orig.RQEF$plot.size),
#            bins = 40) +
#     geom_line(data = plot.preds, aes(x = plot.size, y = estimate), linewidth = 2) +
#   labs(y = expression(paste("Residuals SES.FD"[Q])), 
#        x = expression(paste("Plot size [m"^2,"]"))) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 40),
#     axis.text.x = element_text(size = 36),
#     axis.title.y = element_text(size = 40),
#     axis.text.y = element_text(size = 36),
#     legend.title=element_text(size=36), 
#     legend.text=element_text(size=36),
#     legend.position = "right",
#     #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
#     legend.key.height = unit(1.1, "cm"), 
#     legend.key.width = unit(1.1, "cm"),
#     plot.title = element_text(hjust = .5, size = 70)
#   ) +
#   scale_x_continuous(limits = c(1,
#                                 10001),
#                      breaks = c(1, 4000, 8000)) +
#   #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
#     scale_y_continuous(limits = c(-6, 6)) +
#     scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
#                         limits = c(1, 100000),
#                         breaks = c(10, 100, 1000, 10000))
#     # scale_fill_gradient(low = "blue", high = "red", 
#     #                     trans = "log"
#     # ) 
#   )

## forest, non-forest ----

# pred.orig.RQEF.new <- pred.orig.RQEF %>% 
#   mutate(is.forest = ifelse(is.forest == "TRUE", "Forest", "Non-Forest"))
# 
# forest.data <- sPlot.data %>% 
#   group_by(is.forest) %>% 
#   split(.$is.forest)

# pred.RQEF.for <- lapply(forest.data, function(data){ predictions(mod.RQEF, newdata = datagrid(
#   stable.clim = sample(size = 10, 
#                        x = seq(min(data$stable.clim), max(data$stable.clim), .1), replace = T),
#   PC2 = sample(size = 10, 
#                                        x= seq(min(data$PC2), 
#                                               max(data$PC2), 
#                                               .1), replace = T),
#   PC3 = sample(size = 10, 
#                                        x= seq(min(data$PC3), 
#                                               max(data$PC3), 
#                                               .1), replace = T),
#   PC5 = sample(size = 10, 
#                x= seq(min(data$PC5), 
#                       max(data$PC5), 
#                       .1), replace = T),
#   plot.size = sample(size = 10, 
#                x= seq(1, 
#                       max(data$plot.size, na.rm = TRUE), 
#                       100), replace = T),
#   is.forest = data$is.forest
# ))
# })
# 
# pred.RQEF.for1 <- predictions(mod.RQEF, newdata = datagrid(
#   is.forest = rep(c(TRUE, FALSE), 100)
# ))
# 
# pred.RQEF.for.out <- bind_rows(pred.RQEF.for)

# pred.RQEF.for.out <- pred.RQEF.for1 %>% 
#   mutate(is.forest = ifelse(is.forest == TRUE, "Forest", "Non-Forest"))

# saveRDS(pred.RQEF.for.out, "02.data/_temp-05.GAM-expl.for.RDS")
# 
# pred.RQEF.for.out <- readRDS("02.data/_temp-05.GAM-expl.for.RDS")

# mod.F.FF <- readRDS("02.data/05.GAM_SES.RQEF-exp_6.RDS")
# summary(mod.F.FF)
# 
# res.F.FF <- residuals.gam(mod.F.FF, type = "response")
# 
# 
# p.RQEF.for <- ggplot() +
#   geom_sina(aes(y = res.F.FF, x = pred.orig.RQEF.new$is.forest), color = "grey", alpha = 0.3) +
#   geom_boxplot(aes(x = as.factor(pred.RQEF.for.out$is.forest), y = pred.RQEF.for.out$estimate)) +
#   #geom_violin(aes(x = pred.orig.RQEP$is.forest, y = pred.orig.RQEP$SES.RQEP, color = "SES.RQEP"), size = 1.5) +
#   labs(y = "", x = "", color = "") +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 40),
#     axis.text.x = element_text(size = 40, color = "black", vjust = -4.5),
#     axis.title.y = element_text(size = 40),
#     axis.text.y = element_text(size = 36), 
#     legend.title=element_text(size=36), 
#     legend.text=element_text(size=36),
#     #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
#     legend.key.height = unit(1.1, "cm"), 
#     legend.key.width = unit(1.1, "cm"),
#     legend.position = "none"
#   ) +
#   scale_y_continuous(limits = c(-6, 6))

## spatial smooth ----

dggs <- dgconstruct(area=50*50, metric=T, resround='down')

world <- ne_countries(returnclass = "sf")%>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()

#pred.orig.RQEF$cell <- dgGEO_to_SEQNUM(dggs, pred.orig.RQEF$Longitude, pred.orig.RQEF$Latitude)$seqnum

mod.F.S <- readRDS("02.data/05.GAM_SES.RQEF.BWexp_9.RDS")
summary(mod.F.S)

sPlot.data$res.F.FF <- residuals.gam(mod.F.S, type = "response")

d.RQEF <- sPlot.data %>% 
  group_by(cell) %>% 
  summarise(smooth = mean(res.F.FF, na.rm = T))

grid.RQEF <- dgcellstogrid(dggs, d.RQEF$cell) %>%
  st_as_sf() %>% 
  mutate(cell = d.RQEF$cell) %>% 
  mutate(smooth = d.RQEF$smooth) %>%
  st_transform("+proj=eck4") %>% 
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))


(p.RQEF.geo <- grid.RQEF %>% 
    mutate(area = as.numeric(st_area(.))) %>%
    filter(area <= quantile(area, c(0.999))) %>%  
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
   # geom_sf(data = bb, col = "grey20", fill = NA) +
    coord_sf(crs = "+proj=eck4") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=36),
          legend.text=element_text(size=36),
          #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          legend.position = "bottom", 
          plot.margin=grid::unit(c(0,0,0,0), "mm")
          ) +
    geom_sf(aes(color = smooth, fill = smooth), lwd = 0) +
    scale_fill_viridis_c(option = "plasma",
                         breaks = c(-6,0,6),
                         name = expression(paste("SES.FD"[Q]))) +
    scale_color_viridis_c(option = "plasma",
                          breaks = c(-6,0,6),
                          name = expression(paste("SES.FD"[Q]))) +
    guides(
      fill = "none",
      color = guide_colourbar(
        #title = expression(paste("Higher\nSES.PD"[Q], symbol("\256"), "Higher\nSES.FD"[Q])),
        title.position="top", title.hjust = 0.5
      )) )
    # scale_fill_gradient2(low = "#91bfdb",
    #                      mid = "#ffffbf",
    #                      high = "#fc8d59",
    #                      # limits = c(-5,5),
    #                      breaks = c(-3,0,3)
    #                      )+
    #                      # midpoint = 0) +
    # scale_color_gradient2(low = "#91bfdb",
    #                       mid = "#ffffbf",
    #                       high = "#fc8d59",
    #                       breaks = c(-3,0,3)
    #                       )
    # guides(color = guide_colourbar(title.position="top", title.hjust = 0.5),
    #        size = guide_legend(title.position="top", title.hjust = 0.5))


#save.image(file = "02.data/-temp-05.GAM-expl.RData")

library("patchwork")

layout = "A
          B
          C"

p.out.RQEF <- p.RQEF.PC2 + #p.RQEF.PC3 + 
  p.RQEF.PC5 +
  p.RQEF.sc + 
  plot_layout(design = layout, guides = "collect") +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 40))

ggsave("__Submission/Figures/05.GAM.expl.SES.RQEF.BW.png", p.out.RQEF, 
       height=25, width=12, units="in", dpi=600, bg = "white", limitsize = FALSE)

ggsave("__Submission/Figures/05.GAM.expl.SES.RQEF.BW.pdf", p.out.RQEF, 
       height=25, width=12, units="in", dpi=600, bg = "white", limitsize = FALSE)

## rel log

mod.rel <- readRDS("02.data/05.GAM_rel.log.BW exp_10.RDS")
summary(mod.rel)

# PC1

rel.PC1.preds <- predictions(mod.rel, newdata = datagrid(
  PC1 = seq(
    min(sPlot.data$PC1), 
    max(sPlot.data$PC1), 
    .1)
))


mod.rel.PC1 <- readRDS("02.data/05.GAM_rel.log.BW exp_11.RDS")
summary(mod.rel.PC1)

res.rel.PC1 <- residuals.gam(mod.rel.PC1, type = "response")

(p.rel.PC1 <- ggplot() +
    geom_hex(aes(y = res.rel.PC1, x = sPlot.data$PC1),
             bins = 40) +
    #geom_ribbon(data = RQEP.PC1.preds, aes(ymin = conf.low, ymax = conf.high, x = PC1), outline.type = "full") +
    #geom_point(aes(y = pred.orig.RQEP$SES.RQEP-pred.orig.RQEP$predicted, x = pred.orig.RQEP$PC1), color = "grey80", alpha = 0.5) +
    geom_line(data = rel.PC1.preds, aes(x = PC1, y = estimate), linewidth = 2, lty = 1) +
    labs(y = expression(paste("log(SES.FD"[Q]/"SES.PD"[Q],")")), 
         x = "PC1 - Annual precipitation\n(dry to wet)") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40),
      legend.title=element_text(size=36), 
      legend.text=element_text(size=36),
      legend.position = "right",
      #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
      legend.key.height = unit(1.1, "cm"), 
      legend.key.width = unit(1.1, "cm"),
      plot.margin = margin(10, 10, 10, 40),
      plot.title = element_text(hjust = .5, size = 70)
    ) +
    # scale_x_continuous(limits = c(quantile(pred.orig.RQEF$PC5, 0.025), 
    #                               quantile(pred.orig.RQEF$PC5, 0.9775))) +
    #scale_fill_gradient(low = grey(0.8), high = grey(0.3), trans = "log")
    scale_fill_gradient(high = "#0DE8EF", low = "#373C3D", trans = "log",
                        limits = c(1, 140000),
                        breaks = c(10, 100, 1000, 10000)) +
    scale_x_continuous(limits = c(quantile(sPlot.data$PC1, 0.01), 
                                  23)) +
    scale_y_continuous(limits = c(-2.5,2.5))
)

# forest

pred.orig.rel.new <- sPlot.data %>% 
  mutate(is.forest = ifelse(is.forest == "TRUE", "Forest", "Non-Forest"))

forest.data <- sPlot.data %>% 
  group_by(is.forest) %>% 
  split(.$is.forest)

pred.rel.for <- lapply(forest.data, function(data){ predictions(mod.rel, newdata = datagrid(
  PC1 = sample(size = 10, x = seq(min(sPlot.data$PC1), max(sPlot.data$PC1), .1), replace = T),
  is.forest = data$is.forest
))
})

pred.rel.for.out <- bind_rows(pred.rel.for)

pred.rel.for.out <- pred.rel.for.out %>% 
  mutate(is.forest = ifelse(is.forest == TRUE, "Forest", "Non-Forest"))

pred.rel.for1 <- predictions(mod.rel, newdata = datagrid(
  is.forest = rep(c(TRUE, FALSE), 100)
))

pred.rel.for.out <- pred.rel.for1 %>% 
  mutate(is.forest = ifelse(is.forest == TRUE, "Forest", "Non-Forest"))

mod.rel.FF <- readRDS("02.data/05.GAM_rel.log.BW exp_12.RDS")
summary(mod.rel.FF)

res.rel.FF <- residuals.gam(mod.rel.FF, type = "response")

p.rel.for <- ggplot() +
  geom_violin(aes(y = res.rel.FF, 
                x = pred.orig.rel.new$is.forest), fill = "grey", color = "grey") +
  geom_boxplot(aes(x = as.factor(pred.rel.for.out$is.forest), y = pred.rel.for.out$estimate)) +
  #geom_violin(aes(x = pred.orig.RQEP$is.forest, y = pred.orig.RQEP$SES.RQEP, color = "SES.RQEP"), size = 1.5) +
  labs(y = expression(paste("Residuals log(SES.FD"[Q]/"SES.PD"[Q],")")), x = "", color = "") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 44),
    axis.text.x = element_text(size = 44, color = "black", vjust = -3.8),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 40),
    legend.text = element_text(size = 44),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(-2.5,2.5))
#geom_text(aes(x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, label = "C"), size = 20)

library("patchwork")

layout <- "AABB"


p.out.rel <- p.rel.PC1 + p.rel.for + 
  plot_layout(design = layout, guides = "collect") +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 40))

pp.out <- wrap_elements(p.out.rel) +
  labs(tag = "Residuals") +
  theme(
    plot.tag = element_text(size = 47, angle = 90, vjust = -6#, hjust = 20
                            ),
    plot.tag.position = "left"
  )

ggsave("__Submission/Figures/05.GAM.expl.SES.rel.BW.png", pp.out, 
       height=10, width=25, units="in", dpi=600, bg = "white", limitsize = FALSE)

ggsave("__Submission/Figures/05.GAM.expl.SES.rel.BW.pdf", pp.out, 
       height=10, width=25, units="in", dpi=600, bg = "white", limitsize = FALSE)

##smooth
mod.rel.S <- readRDS("02.data/05.GAM_rel.log.BW exp_13.RDS")
summary(mod.rel.S)

sPlot.data$res.rel.S <- residuals.gam(mod.rel.S, type = "response")

d.rel <- sPlot.data %>% 
  group_by(cell) %>% 
  summarise(smooth = mean(res.rel.S, na.rm = T))

grid.rel <- dgcellstogrid(dggs, d.rel$cell) %>%
  st_as_sf() %>% 
  mutate(cell = d.rel$cell) %>% 
  mutate(smooth = d.rel$smooth) %>%
  st_transform("+proj=eck4") %>% 
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))


(p.rel.geo <- grid.rel %>% 
    mutate(area = as.numeric(st_area(.))) %>%
    filter(area <= quantile(area, c(0.999))) %>%  
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    # geom_sf(data = bb, col = "grey20", fill = NA) +
    coord_sf(crs = "+proj=eck4") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=36),
          legend.text=element_text(size=36),
          #legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          legend.position = "bottom", 
          plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) +
    geom_sf(aes(color = smooth, fill = smooth), lwd = 0) +
    scale_fill_viridis_c(option = "plasma",
                         breaks = c(-1,0,1),
                         name = expression(paste("log(SES.FD"[Q]/"SES.PD"[Q],")"))) +
    scale_color_viridis_c(option = "plasma",
                          breaks = c(-1,0,1),
                          name = expression(paste("log(SES.FD"[Q]/"SES.PD"[Q],")"))) +
    guides(
      fill = "none",
      color = guide_colourbar(
        #title = expression(paste("Higher\nSES.PD"[Q], symbol("\256"), "Higher\nSES.FD"[Q])),
        title.position="top", title.hjust = 0.5
      )) )

## smooths

p.smooth <- p.RQEF.geo + p.phyl.geo + p.rel.geo +
  plot_layout() +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 40))

ggsave("__Submission/Figures/05.GAM.expl.smooth.png", p.smooth, 
       height=10, width=25, units="in", dpi=600, bg = "white", limitsize = FALSE)

# ggsave("__Submission/Figures/05.GAM.expl.smooth.pdf", p.smooth, 
#        height=10, width=25, units="in", dpi=600, bg = "white", limitsize = FALSE)


#### Plots by category status ----

# ll3 <- readRDS("02.Data/04.GAM-status-ll3.RDS")
# 
# 
# p1 <- ll3[["PC1"]] |>
#   filter(fitted < quantile(fitted, 0.95),
#          fitted > quantile(fitted, 0.05)) %>% 
#   ggplot(aes(x = PC1, y = fitted, group = rep)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
#   geom_line(alpha = 0.05) +
#   geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 44),
#     axis.text.x = element_text(size = 40),
#     axis.title.y = element_text(size = 44),
#     axis.text.y = element_text(size = 40) ) +
#   labs(x = "PC1 - Annual precipitation\n(dry to wet)") +
#   annotate("text", x = -3, y = -2.5, label = "Coupling", vjust = -0.5, size = 10) +
#   annotate("text", x = -3, y = -2.5, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
#   
#   annotate("text", x = -3, y = -0.1, label = "Decoupling", vjust = -0.5, size = 10) +
#   annotate("text", x = -3, y = -0.1, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7) +
#   
#   annotate("text", x = -3, y = 2, label = "Decoupling", vjust = -0.5, size = 10) +
#   annotate("text", x = -3, y = 2, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7)
# 
# 
# p2 <- ll3[["PC2"]] |>
#   filter(fitted < quantile(fitted, 0.95),
#          fitted > quantile(fitted, 0.05)) %>% 
#   ggplot(aes(x = PC2, y = fitted, group = rep)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
#   geom_line(alpha = 0.05) +
#   geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 44),
#     axis.text.x = element_text(size = 40),
#     axis.title.y = element_blank(),
#     axis.text.y = element_text(size = 40) ) +
#   labs(x = "PC2 - Temp. coldest\nquarter/month (cold to hot)") +
#   annotate("text", x = -12.5, y = -2.5, label = "Coupling", vjust = -0.5, size = 10) +
#   annotate("text", x = -12.5, y = -2.5, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
#   
#   annotate("text", x = -12.5, y = -0, label = "Decoupling", vjust = -0.5, size = 10) +
#   annotate("text", x = -12.5, y = -0, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7) +
#   
#   annotate("text", x = -12.5, y = 1.5, label = "Decoupling", vjust = -0.5, size = 10) +
#   annotate("text", x = -12.5, y = 1.5, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7)
# 
# 
# 
# p5 <- ll3[["PC5"]] |>
#   filter(fitted < quantile(fitted, 0.95),
#          fitted > quantile(fitted, 0.05)) %>% 
#   ggplot(aes(x = PC5, y = fitted, group = rep)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
#   geom_line(alpha = 0.05) +
#   geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 44),
#     axis.text.x = element_text(size = 40),
#     axis.title.y = element_text(size = 44),
#     axis.text.y = element_text(size = 40) ) +
#   labs(x = "PC5 - Precipitation\nseasonality (low to high)")+
#   annotate("text", x = -4, y = -2, label = "Coupling", vjust = -0.5, size = 10) +
#   annotate("text", x = -4, y = -2, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
#   
#   annotate("text", x = -4, y = -0, label = "Decoupling", vjust = -0.5, size = 10) +
#   annotate("text", x = -4, y = -0, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7) +
#   
#   annotate("text", x = -4, y = 1, label = "Decoupling", vjust = -0.5, size = 10) +
#   annotate("text", x = -4, y = 1, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7)
# 
# 
# ps <- ll3[["Stable.clim"]] |>
#   filter(fitted < quantile(fitted, 0.95),
#          fitted > quantile(fitted, 0.05)) %>% 
#   ggplot(aes(x = stable.clim, y = fitted, group = rep)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
#   geom_line(alpha = 0.05) +
#   geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 44),
#     axis.text.x = element_text(size = 40),
#     axis.title.y = element_blank(),
#     axis.text.y = element_text(size = 40) ) +
#   labs(x = "Climate variability\nafter LGM [°C p.a.]")+
#   annotate("text", x = 0.5, y = -2.5, label = "Coupling", vjust = -0.5, size = 10) +
#   annotate("text", x = 0.5, y = -2.5, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
#   
#   annotate("text", x = 0.5, y = -0, label = "Decoupling", vjust = -0.5, size = 10) +
#   annotate("text", x = 0.5, y = -0, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7) +
#   
#   annotate("text", x = 0.5, y = 1.5, label = "Decoupling", vjust = -0.5, size = 10) +
#   annotate("text", x = 0.5, y = 1.5, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7)
# 
# pf <- ll3[["Forest"]] |>
#   mutate(is.forest = ifelse(is.forest == TRUE, "Forest", "Non-Forest")) %>% 
#   ggplot(aes(x = is.forest, y = fitted, group = rep)) +
#   #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
#   geom_point(alpha = 0.05) +
#   geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 40, color = "black", vjust = -0.5),
#     axis.title.y = element_text(size = 44),
#     axis.text.y = element_text(size = 40) ) +
#   scale_y_continuous(limits = c(-2, 1.5)) +
#   annotate("text", x = "Forest", y = -1.5, label = "Coupling", vjust = -0.5, size = 10, hjust = 1) +
#   annotate("text", x = "Forest", y = -1.5, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7, hjust = 1) +
#   
#   annotate("text", x = "Forest", y = -0, label = "Decoupling", vjust = -0.5, size = 10, hjust = 1) +
#   annotate("text", x = "Forest", y = -0, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7, hjust = 1) +
#   
#   annotate("text", x = "Forest", y = 1, label = "Decoupling", vjust = -0.5, size = 10, hjust = 1) +
#   annotate("text", x = "Forest", y = 1, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7, hjust = 1)
# 
# pp <- ll3[["pv"]] %>% 
#   filter(Expl != "Dev. expl.") %>% 
#   ggplot() +
#   geom_point(aes(x = Expl, y = V1, group = rep), alpha = 0.2) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 40, color = "black", vjust = -0.5),
#     axis.title.y = element_text(size = 44),
#     axis.text.y = element_text(size = 40) ) +
#   labs(y = "p-value")
# 
# devv <- ll3[["pv"]] %>% 
#   filter(Expl == "Dev. expl.") %>% 
#   ggplot() +
#   geom_point(aes(x = Expl, y = V1, group = rep), alpha = 0.2) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.title.y = element_text(size = 44),
#     axis.text.y = element_text(size = 40) ) +
#   labs(y = "Explained deviance", x = "") +
#   scale_y_continuous(breaks = scales::breaks_pretty())
# 
# 
# library(patchwork)
# 
# style <- "AAAAB"
# 
# pd <- pp + devv +
#   plot_layout(design = style)
# 
# ppp <- (p1 + p2) / (p5 + ps) / (pf + pd) +
#   plot_annotation(tag_levels = "A")  & 
#   theme(plot.tag = element_text(size = 40))
# 
# ggsave("__Submission/Figures/05.GAM-cat.png", ppp, 
#        height=30, width=22, units="in", dpi=600, bg = "white", limitsize = FALSE)
# 
# ### 2nd order
# 
# ll2 <- readRDS("02.Data/04.GAM-status-ll2.RDS")
# 
# (p1 <- ll2[["PC1"]] |>
#     filter(fitted < quantile(fitted, 0.95),
#            fitted > quantile(fitted, 0.05)) %>% 
#     ggplot(aes(x = PC1, y = fitted, group = rep)) +
#     geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
#     geom_line(alpha = 0.05) +
#     geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
#     theme_bw() +
#     theme(
#       axis.title.x = element_text(size = 44),
#       axis.text.x = element_text(size = 40),
#       axis.title.y = element_text(size = 44),
#       axis.text.y = element_text(size = 40) ) +
#     labs(x = "PC1 - Annual precipitation\n(dry to wet)") +
#     annotate("text", x = -3, y = -0.1, label = "Coupling", vjust = -0.5, size = 10) +
#     annotate("text", x = -3, y = -0.1, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
#     
#     annotate("text", x = -3, y = -2.5, label = "Decoupling", vjust = -0.5, size = 10) +
#     annotate("text", x = -3, y = -2.5, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7) +
#     
#     annotate("text", x = -3, y = 2, label = "Decoupling", vjust = -0.5, size = 10) +
#     annotate("text", x = -3, y = 2, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7)
#   
# )
# 
# (p2 <- ll2[["PC2"]] |>
#     filter(fitted < quantile(fitted, 0.95),
#            fitted > quantile(fitted, 0.05)) %>% 
#     ggplot(aes(x = PC2, y = fitted, group = rep)) +
#     geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
#     geom_line(alpha = 0.05) +
#     geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
#     theme_bw() +
#     theme(
#       axis.title.x = element_text(size = 44),
#       axis.text.x = element_text(size = 40),
#       axis.title.y = element_blank(),
#       axis.text.y = element_text(size = 40) ) +
#     labs(x = "PC2 - Temp. coldest\nquarter/month (cold to hot)") +
#     annotate("text", x = -12.5, y = -0, label = "Coupling", vjust = -0.5, size = 10) +
#     annotate("text", x = -12.5, y = -0, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
#     
#     annotate("text", x = -12.5, y = -2.5, label = "Decoupling", vjust = -0.5, size = 10) +
#     annotate("text", x = -12.5, y = -2.5, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7) +
#     
#     annotate("text", x = -12.5, y = 1.5, label = "Decoupling", vjust = -0.5, size = 10) +
#     annotate("text", x = -12.5, y = 1.5, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7)
# )
# 
# 
# (p5 <- ll2[["PC5"]] |>
#     filter(fitted < quantile(fitted, 0.95),
#            fitted > quantile(fitted, 0.05)) %>% 
#     ggplot(aes(x = PC5, y = fitted, group = rep)) +
#     geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
#     geom_line(alpha = 0.05) +
#     geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
#     theme_bw() +
#     theme(
#       axis.title.x = element_text(size = 44),
#       axis.text.x = element_text(size = 40),
#       axis.title.y = element_text(size = 44),
#       axis.text.y = element_text(size = 40) ) +
#     labs(x = "PC5 - Precipitation\nseasonality (low to high)")+
#     annotate("text", x = -4, y = -0, label = "Coupling", vjust = -0.5, size = 10) +
#     annotate("text", x = -4, y = -0, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
#     
#     annotate("text", x = -4, y = -2, label = "Decoupling", vjust = -0.5, size = 10) +
#     annotate("text", x = -4, y = -2, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7) +
#     
#     annotate("text", x = -4, y = 1, label = "Decoupling", vjust = -0.5, size = 10) +
#     annotate("text", x = -4, y = 1, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7) )
# 
# 
# (ps <- ll2[["Stable.clim"]] |>
#     filter(fitted < quantile(fitted, 0.95),
#            fitted > quantile(fitted, 0.05)) %>% 
#     ggplot(aes(x = stable.clim, y = fitted, group = rep)) +
#     geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
#     geom_line(alpha = 0.05) +
#     geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
#     theme_bw() +
#     theme(
#       axis.title.x = element_text(size = 44),
#       axis.text.x = element_text(size = 40),
#       axis.title.y = element_blank(),
#       axis.text.y = element_text(size = 40) ) +
#     labs(x = "Climate variability\nafter LGM [°C p.a.]")+
#     annotate("text", x = 0.5, y = -0, label = "Coupling", vjust = -0.5, size = 10) +
#     annotate("text", x = 0.5, y = -0, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
#     
#     annotate("text", x = 0.5, y = -2.5, label = "Decoupling", vjust = -0.5, size = 10) +
#     annotate("text", x = 0.5, y = -2.5, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7) +
#     
#     annotate("text", x = 0.5, y = 1.5, label = "Decoupling", vjust = -0.5, size = 10) +
#     annotate("text", x = 0.5, y = 1.5, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7)
# )
# 
# 
# pf <- ll2[["Forest"]] |>
#   mutate(is.forest = ifelse(is.forest == TRUE, "Forest", "Non-Forest")) %>% 
#   ggplot(aes(x = is.forest, y = fitted, group = rep)) +
#   #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
#   geom_point(alpha = 0.05) +
#   geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 40, color = "black", vjust = -0.5),
#     axis.title.y = element_text(size = 44),
#     axis.text.y = element_text(size = 40) ) +
#   scale_y_continuous(limits = c(-2, 1.5)) +
#   annotate("text", x = "Forest", y = -1.5, label = "Coupling", vjust = -0.5, size = 10, hjust = 1) +
#   annotate("text", x = "Forest", y = -1.5, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7, hjust = 1) +
#   
#   annotate("text", x = "Forest", y = -0, label = "Decoupling", vjust = -0.5, size = 10, hjust = 1) +
#   annotate("text", x = "Forest", y = -0, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7, hjust = 1) +
#   
#   annotate("text", x = "Forest", y = 1, label = "Decoupling", vjust = -0.5, size = 10, hjust = 1) +
#   annotate("text", x = "Forest", y = 1, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7, hjust = 1)
# 
# pp <- ll2[["pv"]] %>% 
#   filter(Expl == "Dev. expl.") %>% 
#   ggplot() +
#   geom_point(aes(x = Expl, y = V1, group = rep), alpha = 0.2) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 40, color = "black", vjust = -0.5),
#     axis.title.y = element_text(size = 44),
#     axis.text.y = element_text(size = 40) ) +
#   labs(y = "p-value")
# 
# devv <- ll2[["pv"]] %>% 
#   filter(Expl == "Dev. expl.") %>% 
#   ggplot() +
#   geom_point(aes(x = Expl, y = V1, group = rep), alpha = 0.2) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.title.y = element_text(size = 44),
#     axis.text.y = element_text(size = 40) ) +
#   labs(y = "Explained deviance", x = "") +
#   scale_y_continuous(breaks = scales::breaks_pretty())
# 
# 
# style <- "AAAAB"
# 
# pd <- pp + devv +
#   plot_layout(design = style)
# 
# ppp <- (p1 + p2) / (p5 + ps) / (pf + pd) +
#   plot_annotation(tag_levels = "A")  & 
#   theme(plot.tag = element_text(size = 40))
# 
# ggsave("__Submission/Figures/05.GAM-cat1.png", ppp, 
#        height=30, width=22, units="in", dpi=600, bg = "white", limitsize = FALSE)

### change again

ll3 <- readRDS("02.Data/04.GAM-status-ll3.RDS")

(p1 <- ll3[["PC1"]] |>
    filter(fitted < quantile(fitted, 0.95),
           fitted > quantile(fitted, 0.05)) %>% 
    ggplot(aes(x = PC1, y = fitted, group = rep)) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = min(ll3[["Thresh"]]$threshold), 
              color = "darkorchid4", fill = NA, linewidth = 3) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = min(ll3[["Thresh"]]$threshold)+0.01, 
              ymax = max(ll3[["Thresh"]]$threshold), 
              color = "darkcyan", fill = NA, linewidth = 3) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = max(ll3[["Thresh"]]$threshold)+0.01, 
              ymax = Inf, color = "olivedrab3", fill = NA, linewidth = 3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
    geom_line(alpha = 0.05) +
    geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40) ) +
    labs(x = "PC1 - Annual precipitation\n(dry to wet)", y = "Fit") +
    annotate("text", x = -3, y = -0.1, label = "Coupling", vjust = -0.5, size = 11, color = "darkcyan")  +
    annotate("text", x = -3, y = -0.1, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7, color = "darkcyan") +
    
    annotate("text", x = -3, y = -2.5, label = "Decoupling", vjust = -0.5, size = 11, color = "darkorchid4") +
    annotate("text", x = -3, y = -2.5, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7, color = "darkorchid4") +
    
    annotate("text", x = -3, y = 2, label = "Decoupling", vjust = -0.5, size = 11, color = "olivedrab3") +
    annotate("text", x = -3, y = 2, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7, color = "olivedrab3")
  
)

(p2 <- ll3[["PC2"]] |>
    filter(fitted < quantile(fitted, 0.95),
           fitted > quantile(fitted, 0.05)) %>% 
    ggplot(aes(x = PC2, y = fitted, group = rep)) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = min(ll3[["Thresh"]]$threshold), 
              color = "darkorchid4", fill = NA, linewidth = 3) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = min(ll3[["Thresh"]]$threshold)+0.01, 
              ymax = max(ll3[["Thresh"]]$threshold), 
              color = "darkcyan", fill = NA, linewidth = 3) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = max(ll3[["Thresh"]]$threshold)+0.01, 
              ymax = Inf, color = "olivedrab3", fill = NA, linewidth = 3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
    geom_line(alpha = 0.05) +
    geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 40) ) +
    labs(x = "PC2 - Temp. coldest\nquarter/month (cold to hot)") 
    # annotate("text", x = -12.5, y = -0, label = "Coupling", vjust = -0.5, size = 10) +
    # annotate("text", x = -12.5, y = -0, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
    # 
    # annotate("text", x = -12.5, y = -2.5, label = "Decoupling", vjust = -0.5, size = 10) +
    # annotate("text", x = -12.5, y = -2.5, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7) +
    # 
    # annotate("text", x = -12.5, y = 1.5, label = "Decoupling", vjust = -0.5, size = 10) +
    # annotate("text", x = -12.5, y = 1.5, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7)
)


(p5 <- ll3[["PC5"]] |>
    filter(fitted < quantile(fitted, 0.95),
           fitted > quantile(fitted, 0.05)) %>% 
    ggplot(aes(x = PC5, y = fitted, group = rep)) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = min(ll3[["Thresh"]]$threshold), 
              color = "darkorchid4", fill = NA, linewidth = 3) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = min(ll3[["Thresh"]]$threshold)+0.01, 
              ymax = max(ll3[["Thresh"]]$threshold), 
              color = "darkcyan", fill = NA, linewidth = 3) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = max(ll3[["Thresh"]]$threshold)+0.01, 
              ymax = Inf, color = "olivedrab3", fill = NA, linewidth = 3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
    geom_line(alpha = 0.05) +
    geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40) ) +
    labs(x = "PC5 - Precipitation\nseasonality (low to high)", y = "Fit")
    # annotate("text", x = -4, y = -0, label = "Coupling", vjust = -0.5, size = 10) +
    # annotate("text", x = -4, y = -0, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
    # 
    # annotate("text", x = -4, y = -2, label = "Decoupling", vjust = -0.5, size = 10) +
    # annotate("text", x = -4, y = -2, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7) +
    # 
    # annotate("text", x = -4, y = 1, label = "Decoupling", vjust = -0.5, size = 10) +
    # annotate("text", x = -4, y = 1, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7) 
  )


(ps <- ll3[["Stable.clim"]] |>
    filter(fitted < quantile(fitted, 0.95),
           fitted > quantile(fitted, 0.05)) %>% 
    ggplot(aes(x = stable.clim, y = fitted, group = rep)) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = min(ll3[["Thresh"]]$threshold), 
              color = "darkorchid4", fill = NA, linewidth = 3) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = min(ll3[["Thresh"]]$threshold)+0.01, 
              ymax = max(ll3[["Thresh"]]$threshold), 
              color = "darkcyan", fill = NA, linewidth = 3) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = max(ll3[["Thresh"]]$threshold)+0.01, 
              ymax = Inf, color = "olivedrab3", fill = NA, linewidth = 3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
    geom_line(alpha = 0.05) +
    geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
    theme_bw() +
    scale_x_continuous(limits = c(0, 2)) +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 40) ) +
    labs(x = "Climate variability\nafter LGM [°C p.a.]")
    # annotate("text", x = 0.5, y = -0, label = "Coupling", vjust = -0.5, size = 10) +
    # annotate("text", x = 0.5, y = -0, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
    # 
    # annotate("text", x = 0.5, y = -2.5, label = "Decoupling", vjust = -0.5, size = 10) +
    # annotate("text", x = 0.5, y = -2.5, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7) +
    # 
    # annotate("text", x = 0.5, y = 1.5, label = "Decoupling", vjust = -0.5, size = 10) +
    # annotate("text", x = 0.5, y = 1.5, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7)
)


(pf <- ll3[["Forest"]] |>
  mutate(is.forest = ifelse(is.forest == TRUE, "Forest", "Non-Forest")) %>% 
  ggplot(aes(x = is.forest, y = fitted, group = rep)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = min(ll3[["Thresh"]]$threshold), 
            color = "darkorchid4", fill = NA, linewidth = 3) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = min(ll3[["Thresh"]]$threshold)+0.01, 
            ymax = max(ll3[["Thresh"]]$threshold), 
            color = "darkcyan", fill = NA, linewidth = 3) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = max(ll3[["Thresh"]]$threshold)+0.01, 
            ymax = Inf, color = "olivedrab3", fill = NA, linewidth = 3) +
  #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01) +
  geom_point(alpha = 0.05) +
  geom_hline(data = ll3[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
  theme_bw() +
  labs(y = "Fit") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 40, color = "black", vjust = -0.5),
    axis.title.y = element_text(size = 44),
    axis.text.y = element_text(size = 40) ) +
  scale_y_continuous(limits = c(-2, 1.5)) 
  )
  # annotate("text", x = "Forest", y = -1.5, label = "Coupling", vjust = -0.5, size = 10, hjust = 1) +
  # annotate("text", x = "Forest", y = -1.5, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7, hjust = 1) +
  # 
  # annotate("text", x = "Forest", y = -0, label = "Decoupling", vjust = -0.5, size = 10, hjust = 1) +
  # annotate("text", x = "Forest", y = -0, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7, hjust = 1) +
  # 
  # annotate("text", x = "Forest", y = 1, label = "Decoupling", vjust = -0.5, size = 10, hjust = 1) +
  # annotate("text", x = "Forest", y = 1, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7, hjust = 1)

(pp <- ll3[["pv"]] %>% 
  filter(Expl != "Dev. expl.") %>% 
  ggplot() +
  geom_point(aes(x = Expl, y = V1, group = rep), alpha = 0.2) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 40, color = "black", vjust = -0.5),
    axis.title.y = element_text(size = 44),
    axis.text.y = element_text(size = 40) ) +
  labs(y = "p-value") +
  scale_y_continuous(limits = c(0, 2e-5),
                     breaks = scales::breaks_pretty())
)

(devv <- ll3[["pv"]] %>% 
  filter(Expl == "Dev. expl.") %>% 
  ggplot() +
  geom_violin(aes(x = Expl, y = V1), alpha = 0.2) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 44),
    axis.text.y = element_text(size = 40) ) +
  labs(y = "Explained deviance", x = "") +
  scale_y_continuous(limits = c(0.13, 0.15), breaks = c(0.13, 0.14, 0.15)) 
 # scale_y_continuous(breaks = scales::breaks_pretty())
)


style <- "AAAAB"

pd <- pp + devv +
  plot_layout(design = style)

ppp <- (p1 + p2) / (p5 + ps) / (pf + devv) +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 40))

ggsave("__Submission/Figures/05.GAM-cat2.png", ppp, 
       height=30, width=22, units="in", dpi=600, bg = "white", limitsize = FALSE)



## linear expl


ll4 <- readRDS("02.Data/04.GAM-status-ll4.RDS")

ll5 <- readRDS("02.Data/04.GAM-status-ll5.RDS")

(p1 <- ll4[["PC1"]] |>
    filter(fitted < quantile(fitted, 0.95),
           fitted > quantile(fitted, 0.05)) %>% 
    ggplot() +
    geom_rect(data = data.frame(
      ymin = c(-Inf, min(ll4[["Thresh"]]$threshold), max(ll4[["Thresh"]]$threshold) ), 
      ymax = c(min(ll4[["Thresh"]]$threshold), max(ll4[["Thresh"]]$threshold), Inf), 
      fill = factor(c("Decoupling with FD < PD", "Coupling with FD ~ PD", "Decoupling with FD > PD"), 
                    levels = c("Decoupling with FD > PD", "Coupling with FD ~ PD", "Decoupling with FD < PD"))),
      aes(xmin = -Inf, xmax = Inf, 
          ymin = ymin, ymax = ymax,
          fill = fill),
      alpha = 0.3) +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = min(ll4[["Thresh"]]$threshold)+0.01, 
    #           ymax = max(ll4[["Thresh"]]$threshold), 
    #           color = "darkcyan", fill = NA, linewidth = 3) +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = max(ll4[["Thresh"]]$threshold)+0.01, 
    #           ymax = Inf, color = "olivedrab3", fill = NA, linewidth = 3) +
    geom_ribbon(aes(x = PC1, y = fitted, group = rep, ymin = lower, ymax = upper), alpha = 0.01) +
    geom_line(aes(x = PC1, y = fitted, group = rep), alpha = 0.05) +
    geom_ribbon(data = ll5[["PC1"]], aes(x = PC1, y = fitted, ymin = lower, ymax = upper), alpha = 0.3, color = "darkblue") +
    geom_line(data = ll5[["PC1"]], aes(x = PC1, y = fitted), color = "darkblue", linewidth = 1.333) +
    geom_hline(data = ll4[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
    geom_hline(data = ll5[["Thresh"]], aes(yintercept = threshold), color = "darkblue", linewidth = 1.3, linetype = 2) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40),
      legend.title = element_blank(),
      legend.text = element_text(size = 40)) +
    scale_fill_viridis_d(direction = -1) +
    labs(x = "PC1 - Annual precipitation\n(dry to wet)", y = "Fitted function")
    # annotate("text", x = -3, y = -0.3, label = "Coupling", vjust = -0.5, size = 11, color = "darkcyan")  +
    # annotate("text", x = -3, y = -0.3, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7, color = "darkcyan") +
    # 
    # annotate("text", x = -3, y = -2.2, label = "Decoupling", vjust = -0.5, size = 11, color = "darkorchid4") +
    # annotate("text", x = -3, y = -2.2, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7, color = "darkorchid4") +
    # 
    # annotate("text", x = -3, y = 1.3, label = "Decoupling", vjust = -0.5, size = 11, color = "olivedrab3") +
    # annotate("text", x = -3, y = 1.3, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7, color = "olivedrab3")
  
)

(p2 <- ll4[["PC2"]] |>
    filter(fitted < quantile(fitted, 0.95),
           fitted > quantile(fitted, 0.05)) %>% 
    ggplot() +
    geom_rect(data = data.frame(
      ymin = c(-Inf, min(ll4[["Thresh"]]$threshold), max(ll4[["Thresh"]]$threshold) ), 
      ymax = c(min(ll4[["Thresh"]]$threshold), max(ll4[["Thresh"]]$threshold), Inf), 
      fill = factor(c("Decoupling with FD < PD", "Coupling with FD ~ PD", "Decoupling with FD > PD"), 
                    levels = c("Decoupling with FD > PD", "Coupling with FD ~ PD", "Decoupling with FD < PD"))),
      aes(xmin = -Inf, xmax = Inf, 
          ymin = ymin, ymax = ymax,
          fill = fill),
      alpha = 0.3) +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = min(ll4[["Thresh"]]$threshold)+0.01, 
    #           ymax = max(ll4[["Thresh"]]$threshold), 
    #           color = "darkcyan", fill = NA, linewidth = 3) +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = max(ll4[["Thresh"]]$threshold)+0.01, 
    #           ymax = Inf, color = "olivedrab3", fill = NA, linewidth = 3) +
    geom_ribbon(aes(x = PC2, y = fitted, group = rep, ymin = lower, ymax = upper), alpha = 0.01) +
    geom_line(aes(x = PC2, y = fitted, group = rep), alpha = 0.05) +
    geom_ribbon(data = ll5[["PC2"]], aes(x = PC2, y = fitted, ymin = lower, ymax = upper), alpha = 0.3, color = "darkblue") +
    geom_line(data = ll5[["PC2"]], aes(x = PC2, y = fitted), color = "darkblue", linewidth = 1.333) +
    geom_hline(data = ll4[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
    geom_hline(data = ll5[["Thresh"]], aes(yintercept = threshold), color = "darkblue", linewidth = 1.3, linetype = 2) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      #axis.title.y = element_blank(),
      axis.text.y = element_text(size = 40),
      legend.title = element_blank(),
      legend.text = element_text(size = 40)) +
    scale_fill_viridis_d(direction = -1) +
    labs(x = "PC2 - Temp. coldest\nquarter/month (cold to hot)", y = "") 
  # annotate("text", x = -12.5, y = -0, label = "Coupling", vjust = -0.5, size = 10) +
  # annotate("text", x = -12.5, y = -0, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
  # 
  # annotate("text", x = -12.5, y = -2.5, label = "Decoupling", vjust = -0.5, size = 10) +
  # annotate("text", x = -12.5, y = -2.5, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7) +
  # 
  # annotate("text", x = -12.5, y = 1.5, label = "Decoupling", vjust = -0.5, size = 10) +
  # annotate("text", x = -12.5, y = 1.5, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7)
)


(p5 <- ll4[["PC5"]] |>
    filter(fitted < quantile(fitted, 0.95),
           fitted > quantile(fitted, 0.05)) %>% 
    ggplot() +
    geom_rect(data = data.frame(
      ymin = c(-Inf, min(ll4[["Thresh"]]$threshold), max(ll4[["Thresh"]]$threshold) ), 
      ymax = c(min(ll4[["Thresh"]]$threshold), max(ll4[["Thresh"]]$threshold), Inf), 
      fill = factor(c("Decoupling with FD < PD", "Coupling with FD ~ PD", "Decoupling with FD > PD"), 
                    levels = c("Decoupling with FD > PD", "Coupling with FD ~ PD", "Decoupling with FD < PD"))),
      aes(xmin = -Inf, xmax = Inf, 
          ymin = ymin, ymax = ymax,
          fill = fill),
      alpha = 0.3) +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = min(ll4[["Thresh"]]$threshold)+0.01, 
    #           ymax = max(ll4[["Thresh"]]$threshold), 
    #           color = "darkcyan", fill = NA, linewidth = 3) +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = max(ll4[["Thresh"]]$threshold)+0.01, 
    #           ymax = Inf, color = "olivedrab3", fill = NA, linewidth = 3) +
    geom_ribbon(aes(x = PC5, y = fitted, group = rep, ymin = lower, ymax = upper), alpha = 0.01) +
    geom_line(aes(x = PC5, y = fitted, group = rep), alpha = 0.05) +
    geom_ribbon(data = ll5[["PC5"]], aes(x = PC5, y = fitted, ymin = lower, ymax = upper), alpha = 0.3, color = "darkblue") +
    geom_line(data = ll5[["PC5"]], aes(x = PC5, y = fitted), color = "darkblue", linewidth = 1.333) +
    geom_hline(data = ll4[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
    geom_hline(data = ll5[["Thresh"]], aes(yintercept = threshold), color = "darkblue", linewidth = 1.3, linetype = 2) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40),
      legend.title = element_blank(),
      legend.text = element_text(size = 40)) +
    scale_fill_viridis_d(direction = -1) +
    labs(x = "PC5 - Precipitation\nseasonality (low to high)", y = "Fitted function")
  # annotate("text", x = -4, y = -0, label = "Coupling", vjust = -0.5, size = 10) +
  # annotate("text", x = -4, y = -0, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
  # 
  # annotate("text", x = -4, y = -2, label = "Decoupling", vjust = -0.5, size = 10) +
  # annotate("text", x = -4, y = -2, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7) +
  # 
  # annotate("text", x = -4, y = 1, label = "Decoupling", vjust = -0.5, size = 10) +
  # annotate("text", x = -4, y = 1, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7) 
)


(ps <- ll4[["Stable.clim"]] |>
    filter(fitted < quantile(fitted, 0.95),
           fitted > quantile(fitted, 0.05)) %>% 
    ggplot() +
    geom_rect(data = data.frame(
      ymin = c(-Inf, min(ll4[["Thresh"]]$threshold), max(ll4[["Thresh"]]$threshold) ), 
      ymax = c(min(ll4[["Thresh"]]$threshold), max(ll4[["Thresh"]]$threshold), Inf), 
      fill = factor(c("Decoupling with FD < PD", "Coupling with FD ~ PD", "Decoupling with FD > PD"), 
                    levels = c("Decoupling with FD > PD", "Coupling with FD ~ PD", "Decoupling with FD < PD"))),
      aes(xmin = -Inf, xmax = Inf, 
          ymin = ymin, ymax = ymax,
          fill = fill),
      alpha = 0.3) +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = min(ll4[["Thresh"]]$threshold)+0.01, 
    #           ymax = max(ll4[["Thresh"]]$threshold), 
    #           color = "darkcyan", fill = NA, linewidth = 3) +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = max(ll4[["Thresh"]]$threshold)+0.01, 
    #           ymax = Inf, color = "olivedrab3", fill = NA, linewidth = 3) +
    geom_ribbon(aes(x = stable.clim, y = fitted, group = rep, ymin = lower, ymax = upper), alpha = 0.01) +
    geom_line(aes(x = stable.clim, y = fitted, group = rep), alpha = 0.05) +
    geom_ribbon(data = ll5[["Stable.clim"]], aes(x = stable.clim, y = fitted, ymin = lower, ymax = upper), alpha = 0.3, color = "darkblue") +
    geom_line(data = ll5[["Stable.clim"]], aes(x = stable.clim, y = fitted), color = "darkblue", linewidth = 1.333) +
    geom_hline(data = ll4[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
    geom_hline(data = ll5[["Thresh"]], aes(yintercept = threshold), color = "darkblue", linewidth = 1.3, linetype = 2) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 44),
      axis.text.x = element_text(size = 40),
      #axis.title.y = element_blank(),
      axis.text.y = element_text(size = 40),
      legend.title = element_blank(),
      legend.text = element_text(size = 40)) +
    scale_fill_viridis_d(direction = -1) +
    labs(x = "Climate variability\nafter LGM [°C p.a.]", y = "")
  # annotate("text", x = 0.5, y = -0, label = "Coupling", vjust = -0.5, size = 10) +
  # annotate("text", x = 0.5, y = -0, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7) +
  # 
  # annotate("text", x = 0.5, y = -2.5, label = "Decoupling", vjust = -0.5, size = 10) +
  # annotate("text", x = 0.5, y = -2.5, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7) +
  # 
  # annotate("text", x = 0.5, y = 1.5, label = "Decoupling", vjust = -0.5, size = 10) +
  # annotate("text", x = 0.5, y = 1.5, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7)
)


(pf <- ll4[["Forest"]]  %>% 
    mutate(is.forest = ifelse(is.forest == TRUE, "Forest", "Non-Forest")) %>% 
    ggplot() +
    geom_rect(data = data.frame(
      ymin = c(-Inf, min(ll4[["Thresh"]]$threshold), max(ll4[["Thresh"]]$threshold) ), 
      ymax = c(min(ll4[["Thresh"]]$threshold), max(ll4[["Thresh"]]$threshold), Inf), 
      fill = factor(c("Decoupling with FD < PD", "Coupling with FD ~ PD", "Decoupling with FD > PD"), 
                    levels = c("Decoupling with FD > PD", "Coupling with FD ~ PD", "Decoupling with FD < PD"))),
      aes(xmin = -Inf, xmax = Inf, 
          ymin = ymin, ymax = ymax,
          fill = fill),
      alpha = 0.3) +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = min(ll4[["Thresh"]]$threshold)+0.01, 
    #           ymax = max(ll4[["Thresh"]]$threshold), 
    #           color = "darkcyan", fill = NA, linewidth = 3) +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = max(ll4[["Thresh"]]$threshold)+0.01, 
    #           ymax = Inf, color = "olivedrab3", fill = NA, linewidth = 3) +
    #geom_ribbon(aes(x = PC1, y = fitted, group = rep, ymin = lower, ymax = upper), alpha = 0.01) +
    #geom_line(aes(x = PC1, y = fitted, group = rep), alpha = 0.05) +
    geom_point(aes(x = is.forest, y = fitted, group = rep)) +
    geom_point(data = ll5[["Forest"]] %>% 
                 mutate(is.forest = ifelse(is.forest == TRUE, "Forest", "Non-Forest")), 
               aes(x = is.forest, y = fitted), color = "darkblue", size = 3.666) +
    geom_hline(data = ll4[["Thresh"]], aes(yintercept = threshold, group = rep), alpha = 0.01) +
    geom_hline(data = ll5[["Thresh"]], aes(yintercept = threshold), color = "darkblue", linewidth = 1.3, linetype = 2) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 30, color = "black", vjust = 0.2),
      axis.title.y = element_text(size = 44),
      axis.text.y = element_text(size = 40),
      legend.title = element_blank(),
      legend.text = element_text(size = 40)) +
    scale_fill_viridis_d(direction = -1) +
    labs(y = "Fitted function") +
    scale_y_continuous(limits = c(-2, 1.5)) 
)
# annotate("text", x = "Forest", y = -1.5, label = "Coupling", vjust = -0.5, size = 10, hjust = 1) +
# annotate("text", x = "Forest", y = -1.5, label = expression(paste("with FD" %~~% "PD")), vjust = 0.8, size = 7, hjust = 1) +
# 
# annotate("text", x = "Forest", y = -0, label = "Decoupling", vjust = -0.5, size = 10, hjust = 1) +
# annotate("text", x = "Forest", y = -0, label =  expression(paste("with FD"<"PD")), vjust = 0.8, size = 7, hjust = 1) +
# 
# annotate("text", x = "Forest", y = 1, label = "Decoupling", vjust = -0.5, size = 10, hjust = 1) +
# annotate("text", x = "Forest", y = 1, label =  expression(paste("with FD">"PD")), vjust = 0.8, size = 7, hjust = 1)

(pp <- ll4[["pv"]] %>% 
    filter(Expl != "Dev. expl.") %>% 
    ggplot() +
    geom_point(aes(x = Expl, y = V1, group = rep), alpha = 0.2) +
    geom_point(data = ll5[["pv"]] %>% 
                 filter(Expl != "Dev. expl."), 
               aes(x = Expl, y = V1), color = "darkblue", size = 3.666) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 40, color = "black", vjust = -0.5),
      axis.title.y = element_text(size = 34),
      axis.text.y = element_text(size = 30) ) +
    labs(y = "p-value") +
    scale_y_continuous(breaks = scales::breaks_pretty())
)

(devv <- ll4[["pv"]] %>% 
    filter(Expl == "Dev. expl.") %>% 
    ggplot() +
    geom_violin(aes(x = Expl, y = V1), alpha = 0.2) +
    geom_point(data = ll5[["pv"]] %>% 
                 filter(Expl == "Dev. expl."), 
               aes(x = Expl, y = V1), color = "darkblue", size = 3.666) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 34),
      axis.text.y = element_text(size = 30) ) +
    labs(y = "Expl. deviance", x = "") +
    scale_y_continuous(breaks = c(0.11, 0.12, 0.13, 0.14, 0.15)) 
  # scale_y_continuous(breaks = scales::breaks_pretty())
)


style <- "AAAABBB \n AAAABBB \n CCCCDDD \n CCCCDDD
          EEEEFGG \n EEEEHHH"

ppp <- p1 + p2 + p5 + ps + pf + pp + devv + guide_area() +
  plot_layout(guides = "collect",
              design = style) &
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 40))

ggsave("__Submission/Figures/05.GAM-cat-lin.png", ppp, 
       height=30, width=22, units="in", dpi=600, bg = "white", limitsize = FALSE)

ggsave("__Submission/Figures/05.GAM-cat-lin.pdf", ppp, 
       height=30, width=22, units="in", dpi=600, bg = "white", limitsize = FALSE)
