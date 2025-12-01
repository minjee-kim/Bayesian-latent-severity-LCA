
###############################################
library(ggplot2)
library(grid)

# ---- inputs ----
gamma_thr <- 1.2
mu        <- 1.6
sd_latent <- 1

# ---- derived metrics ----
spec <- pnorm(gamma_thr, mean = 0, sd = sd_latent)
sens <- 1 - pnorm(gamma_thr, mean = mu, sd = sd_latent)

# ---- data grid ----
xmax <- max(gamma_thr, mu) + 4
xmin <- -4
x <- seq(xmin, xmax, length.out = 2000)
dh <- dnorm(x, mean = 0,  sd = sd_latent)
dd <- dnorm(x, mean = mu, sd = sd_latent)
ymax <- max(dh, dd)

# groups
df <- rbind(
  data.frame(x, y = dh, group = "healthy"),
  data.frame(x, y = dd, group = "diseased")
)
df$group <- factor(df$group, levels = c("healthy", "diseased"))

# shaded regions
shade_spec <- data.frame(x = x[x <= gamma_thr], y = dh[x <= gamma_thr], what = "spec")
shade_sens <- data.frame(x = x[x >= gamma_thr], y = dd[x >= gamma_thr], what = "sens")
shade_df   <- rbind(shade_spec, shade_sens)

# colors
col_healthy  <- "#0072B2"
col_diseased <- "#D55E00"

p <- ggplot(df, aes(x, y)) +
  geom_area(data = shade_df, aes(x, y, fill = what), alpha = 0.25, color = NA) +
  geom_line(aes(color = group), linewidth = 1.1, lineend = "round") +
  geom_vline(xintercept = gamma_thr, linetype = 2, linewidth = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey25") +
  # ---- clearer γ_j label: above the curves with a light backdrop ----
annotate("label",
         x = gamma_thr, y = ymax * 0.92,
         label = "gamma[j]", parse = TRUE,
         size = 4.6, label.size = 0, fill = "white", alpha = 0.85) +
  # callouts with subscripts
  annotate("label",
           x = xmin + 0.5, y = ymax * 0.96,
           label = "Specificity == Phi(gamma[j])",
           parse = TRUE, hjust = 0, size = 3.6, label.size = 0) +
  annotate("label",
           x = xmax - 0.5, y = ymax * 0.88,
           label = "Sensitivity == 1 - Phi(gamma[j] - beta[j]*S[i])",
           parse = TRUE, hjust = 1, size = 3.6, label.size = 0) +
  scale_color_manual(
    values = c(healthy = col_healthy, diseased = col_diseased),
    breaks  = c("healthy", "diseased"),
    labels  = c(
      expression(Healthy ~ N(0, 1)),
      expression(Diseased ~ N(beta[j] * S[i], 1))
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(spec = col_healthy, sens = col_diseased),
    breaks = c("spec","sens"),
    labels = c(
      expression("Specificity:  " * Pr(V[ij] <= gamma[j] ~ "|" ~ Healthy)),
      expression("Sensitivity: " * Pr(V[ij] >  gamma[j] ~ "|" ~ Diseased))
    ),
    name = NULL
  ) +
  labs(
    x = expression("Latent test score  " * V[ij]),
    y = "Density"
  ) +
  coord_cartesian(xlim = c(xmin, xmax), ylim = c(0, ymax * 1.12)) +
  theme_classic(base_size = 14) +
  theme(
    plot.title.position = "plot",
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "vertical",
    axis.ticks.length = grid::unit(3, "pt"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(8, 12, 8, 12)
  )

print(p)

ggsave("Vij_plot.png", p, width = 6, height = 6.5, dpi = 1000)


# ==== Publication-ready panel: Gamma(4.5, sqrt(4.5)) vs Normal(0, 1) ====
library(ggplot2)

col_gamma  <- "#E69F00"
col_normal <- "#999999"

shape <- 4.5
rate  <- sqrt(4.5)

# Full Normal(0,1), one-sided Gamma
x_n <- seq(-4, 4, length.out = 2400)
x_g <- seq(0, 10, length.out = 2400)

df2 <- rbind(
  data.frame(
    x = x_g,
    y = dgamma(x_g, shape = shape, rate = rate),
    dist = "Gamma(4.5, sqrt(4.5))"
  ),
  data.frame(
    x = x_n,
    y = dnorm(x_n, mean = 0, sd = 1),
    dist = "Normal(0, 1)"
  )
)

y_max <- max(df2$y) * 1.1

p2 <- ggplot(df2, aes(x, y, color = dist, fill = dist)) +
  geom_area(alpha = 0.15, position = "identity") +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6, color = "grey40") +
  scale_color_manual(values = c("Gamma(4.5, sqrt(4.5))" = col_gamma,
                                "Normal(0, 1)" = col_normal)) +
  scale_fill_manual(values = c("Gamma(4.5, sqrt(4.5))" = col_gamma,
                               "Normal(0, 1)" = col_normal)) +
  labs(
    x = expression(S[i]~"(latent severity)"),
    y = "Density",
    color = NULL,
    fill  = NULL
  ) +
  coord_cartesian(xlim = c(-4, 8), ylim = c(0, y_max)) +
  theme_classic(base_size = 14) +   # match base_size too if you want
  theme(
    axis.ticks.length = grid::unit(3, "pt"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box        = "vertical",
    plot.margin       = margin(8, 12, 8, 12)
  )

ggsave("Figure_gamma.png", p2, width = 6, height = 6.5, dpi = 800)





