####packages ----
library("tidyverse")
library("tidyr")
library("dplyr")
library("data.table")
library("raster")
library("ncdf4")
library("ggplot2")
library("hexbin")
library("sf")
library("rnaturalearth")
library("viridis")
library("dggridR")
library("scales")

sPlot.data <- readRDS("02.data/03.sPlot.PD.FD.CD-data.RDS")

#### Rao'S QE ~ phyl + funct ----

## Hotspots

dggs <- dgconstruct(area=50*50, metric=T, resround='down')

sPlot.data$cell <- dgGEO_to_SEQNUM(dggs, sPlot.data$Longitude, sPlot.data$Latitude)$seqnum

d.RQE <- sPlot.data %>% 
  filter(SES.RQEP>quantile(SES.RQEP, 0.95) &
           SES.RQEF>quantile(SES.RQEF, 0.95) ) %>%
  group_by(cell) %>% 
  summarise(RQEP = mean(SES.RQEP, na.rm = T),
            RQEF = mean(SES.RQEF, na.rm = T),
            n = n())

grid <- dgcellstogrid(dggs, d.RQE$cell) %>%
  st_as_sf() %>% 
  mutate(cell = d.RQE$cell) %>% 
  mutate(RQEP = d.RQE$RQEP) %>% 
  mutate(RQEF = d.RQE$RQEF) %>% 
  mutate(n = d.RQE$n) %>%
  st_transform("+proj=eck4") %>% 
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))

world <- ne_countries(returnclass = "sf")%>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()


d<-expand.grid(x=1:4,y=1:4)
RQEP.RQEF.legend <- ggplotGrob(
  ggplot(d, aes(x,y,fill=atan(y/x), alpha = atan(x+y)))+
    geom_tile()+
    scale_fill_viridis() +
    labs(y = expression(paste("SES Functional RQ Entropy",
                              symbol("\256"))),
         x = expression(paste("SES Phylogenetic RQ Entropy",
                              symbol("\256")))) +
    theme_minimal() +
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(colour = "black"),
          legend.position="none") +
    theme(
      axis.title = element_text(size = 13)
    ) +
    coord_fixed()
)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

(p.hot <- grid %>% 
    mutate(area = as.numeric(st_area(.))) %>%
    filter(area <= quantile(area, c(0.99))) %>% 
    filter(!is.na(RQEP),
           !is.na(RQEF)) %>% 
    mutate(RQEP = range01(RQEP)) %>% 
    mutate(RQEF = range01(RQEF)) %>% 
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=12), 
          legend.text=element_text(size=12),
          legend.background = element_rect(size=1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm")) +
    geom_sf(aes(fill = atan(RQEF/RQEP), color = atan(RQEF+RQEP)), 
            lwd = 0) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    theme(legend.position = "none") +
    annotation_custom(grob = RQEP.RQEF.legend,
                      #xmin =  -45e6,
                      xmin = -51e6,
                      ymax =  -0.5e6,
                      ymin =  -8.5e6
    )
)
#Construct a global grid
dggs <- dgconstruct(area=50*50, metric=T, resround='down')

sPlot.data$cell <- dgGEO_to_SEQNUM(dggs, sPlot.data$Longitude, sPlot.data$Latitude)$seqnum

d.RQE <- sPlot.data %>% 
  group_by(cell) %>% 
  summarise(RQEP = mean(SES.RQEP, na.rm = T),
            RQEF = mean(SES.RQEF, na.rm = T),
            n = n())

grid <- dgcellstogrid(dggs, d.RQE$cell) %>%
  st_as_sf() %>% 
  mutate(cell = d.RQE$cell) %>% 
  mutate(RQEP = d.RQE$RQEP) %>% 
  mutate(RQEF = d.RQE$RQEF) %>% 
  mutate(n = d.RQE$n) %>%
  st_transform("+proj=eck4") %>% 
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))

world <- ne_countries(returnclass = "sf")%>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()


d<-expand.grid(x=1:4,y=1:4)
RQEP.RQEF.legend <- ggplotGrob(
  ggplot(d, aes(x,y,fill=atan(y/x), alpha = atan(x+y)))+
    geom_tile()+
    scale_fill_viridis() +
    labs(y = expression(paste("SES Functional RQ Entropy",
                    symbol("\256"))),
         x = expression(paste("SES Phylogenetic RQ Entropy",
                              symbol("\256")))) +
    theme_minimal() +
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(colour = "black"),
          legend.position="none") +
    theme(
      axis.title = element_text(size = 13)
    ) +
    coord_fixed()
)

(p.n <- grid %>% 
    mutate(area = as.numeric(st_area(.))) %>%
    filter(area <= quantile(area, c(0.99))) %>% 
    filter(!is.na(RQEP),
           !is.na(RQEF)) %>% 
    mutate(n = scale(n, center = .5)) %>%
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    coord_sf(crs = "+proj=eck4") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=12), 
          legend.text=element_text(size=12),
          legend.background = element_rect(size=1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm")) +
    geom_sf(aes(fill = n, color = n), lwd = 0) +
    scale_fill_gradient(low = grey(.8), high = grey(.05), trans = "log") +
    scale_color_gradient(low = grey(.8), high = grey(.05), trans = "log") +
    theme(legend.position = "none") )

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

p.RQEP.RQEF <- sPlot.data %>%
  filter(!is.na(Longitude),
         !is.na(Latitude)) %>%
  mutate(RQEF = SES.RQEF) %>% 
  mutate(RQEP = SES.RQEP) %>% 
  filter(!is.na(RQEP),
         !is.na(RQEF)) %>% 
  mutate(RQEP = range01(RQEP)) %>% 
  mutate(RQEF = range01(RQEF)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "WGS84") %>% 
  st_transform(., crs = "+proj=eck4") %>%
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=12), 
          legend.text=element_text(size=12),
          legend.background = element_rect(size=1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm")) +
    geom_sf(aes(fill = atan(RQEF/RQEP), color = atan(RQEF/RQEP) , alpha = atan(RQEF+RQEP)), 
                lwd = 0) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    theme(legend.position = "none") +
    annotation_custom(grob = RQEP.RQEF.legend,
                      xmin =  -45e6,
                      ymax =  -1.1e6,
                      ymin =  -8.5e6
                      ) 

g <- ggplot_build(p.RQEP.RQEF)
color.legend <- ggplot_build(ggplot(d, aes(x,y,fill=atan(y/x),alpha = atan(x+y)))+
                            geom_tile()+
                            scale_fill_viridis())
out <- sPlot.data %>%
  mutate(color = g$data[[2]]$fill,
         alpha = g$data[[2]]$alpha) %>%
  dplyr::select(color, alpha)

color_distance <- function(color1, color2, alpha1, alpha2) {
  sum((round(alpha1*col2rgb(color1)+(1-alpha1)*255, digits = 0) - #transform to RGB and add alpha
         round(alpha2*col2rgb(color2)+(1-alpha2)*255, digits = 0))^2)
}

df.legend <- color.legend$data[[1]] %>%
  as.data.frame() %>%
  dplyr::select(fill, alpha, x, y) %>%
  mutate(fill.alpha = paste0(fill, round(alpha, digits = 3)))

closest_colors <- mapply(color1 = out$color, alpha1 = out$alpha, function(color1, alpha1) {
  distances <- mapply(color2 = df.legend$fill, alpha2 = df.legend$alpha,
                      function(color2, alpha2) color_distance(color1, color2, alpha1, alpha2))
  df.legend$fill.alpha[which.min(distances)]
})

rm(g)

col.out <- closest_colors
table(col.out)

col.freq <- round(100 * table(col.out) / length(col.out), 5) %>% as.data.frame()

df.legend <- df.legend %>%
  as.data.frame() %>%
  dplyr::select(fill, x, y, fill.alpha) %>%
  left_join(col.freq, by = c("fill.alpha" = "col.out"))

new.legend <- ggplotGrob(
  ggplot(df.legend %>%
          mutate(Freq = ifelse(is.na(Freq) | Freq < 0.01, "<0.01", round(Freq, digits = 2)))
       , aes(x,y,fill=atan(y/x),alpha =atan(x+y)))+
  geom_tile()+
  scale_fill_viridis() +
  labs(y = expression(paste("SES.FD"[Q],
                            symbol("\256"))),
       x = expression(paste("SES.PD"[Q],
                            symbol("\256")))) +
  theme_minimal() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(colour = "black"),
        legend.position="none",
        axis.title = element_text(size = 36)
  ) +
  coord_fixed() #+
  # geom_text(aes(label = Freq), alpha = 1,
  #               size = 6.5)
)

(p.new <- grid %>% 
    mutate(area = as.numeric(st_area(.))) %>%
    filter(area <= quantile(area, c(0.99))) %>% 
    filter(!is.na(RQEP),
           !is.na(RQEF)) %>% 
    mutate(RQEP = range01(RQEP)) %>% 
    mutate(RQEF = range01(RQEF)) %>% 
    ggplot() +
    geom_sf(data = world, fill = "grey90", col = NA, lwd = 0.3) +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.title=element_text(size=12), 
          legend.text=element_text(size=12),
          legend.background = element_rect(size=1, linetype="solid", colour = 1), 
          legend.key.height = unit(1.1, "cm"), 
          legend.key.width = unit(1.1, "cm"),
          plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    geom_sf(aes(fill = atan(RQEF/RQEP), color = atan(RQEF/RQEP)), 
            lwd = 0) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    theme(legend.position = "none") +
    annotation_custom(grob = new.legend,
                      #xmin =  -45e6,
                      xmin = -51e6,
                      ymax =  -0.5e6,
                      ymin =  -8.5e6
    )
)

# 
# (lat.RQEF <-  
#     ggplot(sPlot.data ,aes(x=Latitude, 
#                                y = SES.RQEF))+
#     geom_hex(na.rm = T, 
#              alpha=.8, show.legend=FALSE)+
#     geom_smooth(method = "loess",  se=F, color = "red")+
#     labs(title = '', 
#          x = '', 
#          y = "SES Functional RQ Entropy")+
#     scale_x_continuous(breaks = c(-60, -40, -20, 0, 20, 
#                                   40, 60, 80),
#                        limits = c(-65,85),
#                        label = c("60°S", "40°S", "20°S", 
#                                  "0", "20°N", "40°N", 
#                                  "60°N", "80°N")) + 
#     # scale_y_continuous(breaks = c(500, 1000),
#     #                    limits = c(100,1400),
#     #                    label = c("500", "1000")) + 
#     scale_fill_gradientn(colors = c("grey90","grey10"),
#                          trans = "log")+
#     theme_classic() +
#     theme(text = element_text(size = 18),
#           axis.title = element_text(size = 14),
#           plot.margin = unit(c(0,0,0,0), "mm")) +
#     coord_flip()
# )




# ggsave("03.results/04.grid-cell-plot-number.png", p.n, height=10, width=15, units="in", dpi=600, bg = "white")
# 

load("03.results/_temp.04.FD-PD-GAM-SES.RData")


library("patchwork")


p.out <- (out +
  labs(y=expression(paste("SES.FD"[Q])))) /p.new + 
  plot_layout(widths = c(1,1), heights = c(1,1)) + 
  plot_annotation(tag_levels = list(c("A", "", "", "", "B")))  & 
  theme(plot.tag = element_text(size = 40)) 

ggsave("__Submission/Figures/04.map-gam.png", p.out, height=30, width=20, units="in", dpi=600, bg = "white")

