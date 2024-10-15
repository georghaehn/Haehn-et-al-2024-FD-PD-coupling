#### recommended machine: 01#####

library("tidyverse")
library("ggdendro")
library("viridis")
library("patchwork")

hc <- hclust(dist(mtcars))

ggdendrogram(hc, theme_dendro = TRUE, rotate = TRUE)

(top.left <- mtcars %>%
  as.data.frame() %>%
  mutate(car = rownames(mtcars)) %>%
  filter(car == "Merc 280C" | car == "Mazda RX4" |
    car == "Merc 240D" | car == "Merc 230") %>%
  dplyr::select(-car) %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendro_data(., type = "rectangle") %>%
  segment() %>%
  ggplot() +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1.5) +
  geom_segment(aes(x = 2.5, y = 150, xend = 2.5, yend = 65), linewidth = 1.5) +
  # geom_point(aes(x = 2, y = -4), shape = 17, colour = "red", size = 10) +
  # geom_point(aes(x = 1, y = -4), shape = 19, colour = "darkgreen", size = 10) +
  labs(
    title = "Decoupling",
    subtitle = expression(paste("with FD" > "PD"))
  ) +
  coord_flip() +
  scale_y_reverse(
    expand = c(0, 0),
    limits = c(150, -150 / 4)
  ) +
  scale_x_continuous(limits = c(0.9, 4.1)) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 40),
    plot.subtitle = element_text(hjust = 0.5, size = 30),
    # panel.border = element_rect(colour = "darkgrey", fill = NA, linewidth = 2),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
)

(top.right <- mtcars %>%
  as.data.frame() %>%
  mutate(car = rownames(mtcars)) %>%
  filter(car == "Merc 280C" | car == "Mazda RX4" |
    car == "Merc 240D" | car == "Merc 230") %>%
  dplyr::select(-car) %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendro_data(., type = "rectangle") %>%
  segment() %>%
  ggplot() +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1.5) +
  geom_segment(aes(x = 2.5, y = 70, xend = 2.5, yend = 65), linewidth = 1.5) +
  # geom_point(aes(x = 2, y = -4), shape = 17, colour = "red", size = 10) +
  # geom_point(aes(x = 1, y = -4), shape = 19, colour = "darkgreen", size = 10) +
  labs(
    title = "Coupling",
    subtitle = expression(paste("with FD" %~~% "PD"))
  ) +
  coord_flip() +
  scale_y_reverse(
    expand = c(0, 0),
    limits = c(70, -70 / 4)
  ) +
  scale_x_continuous(limits = c(0.9, 4.1)) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 40),
    plot.subtitle = element_text(hjust = 0.5, size = 30),
    # panel.border = element_rect(colour = "darkgrey", fill = NA, linewidth = 2),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
)

(mid <- mtcars %>%
  as.data.frame() %>%
  mutate(car = rownames(mtcars)) %>%
  filter(car == "Merc 280C" | car == "Mazda RX4" |
    car == "Merc 240D" | car == "Merc 230") %>%
  dplyr::select(-car) %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendro_data(., type = "rectangle") %>%
  segment() %>%
  ggplot() +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1.5) +
  geom_segment(aes(x = 2.5, y = 110, xend = 2.5, yend = 65), linewidth = 1.5) +
  # geom_point(aes(x = 2, y = -4), shape = 17, colour = "red", size = 10) +
  # geom_point(aes(x = 1, y = -4), shape = 19, colour = "darkgreen", size = 10) +
  labs(
    title = "Coupling",
    subtitle = expression(paste("with FD" %~~% "PD"))
  ) +
  coord_flip() +
  scale_y_reverse(
    expand = c(0, 0),
    limits = c(110, -110 / 4)
  ) +
  scale_x_continuous(limits = c(0.9, 4.1)) +
  # scale_x_continuous(limits = c(0,5)) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 40),
    plot.subtitle = element_text(hjust = 0.5, size = 30),
    # panel.border = element_rect(colour = "darkgrey", fill = NA, linewidth = 2),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
)


(bottom.left <- mtcars %>%
  as.data.frame() %>%
  mutate(car = rownames(mtcars)) %>%
  filter(car == "Merc 280C" | car == "Mazda RX4" |
    car == "Merc 240D" | car == "Merc 230") %>%
  dplyr::select(-car) %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendro_data(., type = "rectangle") %>%
  segment() %>%
  ggplot() +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1.5) +
  geom_segment(aes(x = 2.5, y = 150, xend = 2.5, yend = 65), linewidth = 1.5) +
  # geom_point(aes(x = 2, y = -4), shape = 17, colour = "red", size = 10) +
  # geom_point(aes(x = 1, y = -4), shape = 19, colour = "darkgreen", size = 10) +
  labs(
    title = "Coupling",
    subtitle = expression(paste("with FD" %~~% "PD"))
  ) +
  coord_flip() +
  scale_y_reverse(
    expand = c(0, 0),
    limits = c(150, -150 / 4)
  ) +
  scale_x_continuous(limits = c(0.9, 4.1)) +
  # scale_x_continuous(limits = c(-2,7)) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 40),
    plot.subtitle = element_text(hjust = 0.5, size = 30),
    # panel.border = element_rect(colour = "darkgrey", fill = NA, linewidth = 2),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
)


(bottom.right <- mtcars %>%
  as.data.frame() %>%
  mutate(car = rownames(mtcars)) %>%
  filter(car == "Merc 280C" | car == "Mazda RX4" |
    car == "Merc 240D" | car == "Merc 230") %>%
  dplyr::select(-car) %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendro_data(., type = "rectangle") %>%
  segment() %>%
  ggplot() +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1.5) +
  geom_segment(aes(x = 2.5, y = 70, xend = 2.5, yend = 65), linewidth = 1.5) +
  # geom_point(aes(x = 2, y = -4), shape = 17, colour = "red", size = 10) +
  # geom_point(aes(x = 1, y = -4), shape = 19, colour = "darkgreen", size = 10) +
  labs(
    title = "Decoupling",
    subtitle = expression(paste("with FD" < "PD"))
  ) +
  coord_flip() +
  scale_y_reverse(
    expand = c(0, 0),
    limits = c(70, -70 / 4)
  ) +
  scale_x_continuous(limits = c(0.9, 4.1)) +
  # scale_x_continuous(limits = c(-2,7)) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 40),
    plot.subtitle = element_text(hjust = 0.5, size = 30),
    # panel.border = element_rect(colour = "darkgrey", fill = NA, linewidth = 2),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
)


bg <- expand.grid(
  x = c(1:1000),
  y = c(1:1000)
) %>%
  ggplot(aes(x, y, fill = atan(y / x), alpha = atan(y + x))) +
  geom_abline(slope = 1, color = "grey60", intercept = seq(-100, 100, 0.01), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.5, lty = 3) +
  geom_raster(alpha = 0.4) +
  scale_fill_viridis() +
  labs(
    y = expression(paste(
      "Functional diversity",
      symbol("\256")
    )),
    x = expression(paste(
      "Phylogenetic diversity",
      symbol("\256")
    ))
  ) +
  theme_minimal() +
  scale_x_discrete(
    breaks = c(10, 500, 1000),
    labels = c(0, 0.5, 1)
  ) +
  scale_y_discrete(
    breaks = c(10, 500, 1000),
    labels = c(0, 0.5, 1)
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(colour = "darkgrey", linewidth = 2),
    legend.position = "none",
    axis.title = element_text(size = 50),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30)
  ) +
  coord_fixed()

concept <- bg +
  inset_element(top.left, left = 0.05, bottom = 0.65, right = 0.35, top = 0.95) + # top left
  inset_element(bottom.left, left = 0.05, bottom = 0.05, right = 0.35, top = 0.35) + # bottom left
  inset_element(mid, left = 0.35, bottom = 0.35, right = 0.65, top = 0.65) + # mid
  inset_element(top.right, left = 0.65, bottom = 0.65, right = 0.95, top = 0.95) + # top right
  inset_element(bottom.right, left = 0.65, bottom = 0.05, right = 0.95, top = 0.35) # bottom right


ggsave(file = "__Submission/Figures/00.concept-coupling.png", concept, height = 20, width = 20, units = "in", dpi = 600)
