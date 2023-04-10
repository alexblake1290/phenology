
#### Setup ####

library(tidyverse)
library(brms)
library(tidybayes)
library(ecodatamisc)

# Issues installing from instructions on R 4.2. Use:
# https://discourse.mc-stan.org/t/rstan-fails-to-build-example-model/27248
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# test with
# example(stan_model, package = "rstan", run.dontrun = TRUE)

dat = read.csv("./aphids_2023/osf/Aggregated Aphid Data/Idaho pan trap with degree days.csv") %>%
  mutate(
    Site = str_remove(SiteName, "[:digit:]{4}-")
  )

vdat = read.csv("./aphids_2023/osf/Aggregated Aphid Data/Vetch 2019 and 2020 with degree days.csv")


#### Model A: Pan trap, Aphid DD ####
m1 <- brm(AphidCount ~ s(Aphids_DD),
          family=negbinomial,
          data = dat,
          chains = 4, cores = 4,
          iter = 5000, warmup = 1000,
          seed = 42,
          control = list(adapt_delta = 0.95,
                         max_treedepth = 15)
          )
# 1 divergent transition after warmup.

#write_rds(m1, "./aphids_2023/fit_models/pooled_gam_cov.rds")

# Extract peak
# Get curve and 95% CIs for plotting
m1_p <- data.frame(Aphids_DD = seq(from = 0, to = max(dat$Aphids_DD, na.rm = T), by = 1)) %>%
  add_fitted_draws(m1) %>%
  median_hdci()

# Get curve and 95% CI for the peak
m1_peak_dist <- data.frame(Aphids_DD = seq(from = 0, to = max(dat$Aphids_DD, na.rm = T), by = 1)) %>%
  add_fitted_draws(m1) %>%
  ungroup() %>%
  group_by(.draw) %>%
  slice_max(.value, n = 1) %>%
  ungroup() %>%
  median_hdci()

# 1204, 1338, 1452
# 1201, 1338, 1452 with hdis rather hdcis
cdd_med <- m1_peak_dist$Aphids_DD
cdd_lo <- m1_peak_dist$Aphids_DD.lower
cdd_hi <- m1_peak_dist$Aphids_DD.upper

# Remove what we don't need anymore
# rm(m1_peak_dist)
# rm(m1)

# Plot
p1 <- m1_p %>%
  ggplot(aes(x = Aphids_DD, y = .value, ymin = .lower, ymax = .upper)) +
  geom_vline(xintercept = cdd_med, color = "#FF0000", linetype = 2) +
  annotate(geom = 'rect', xmin = cdd_lo, xmax = cdd_hi,
           ymin = -Inf, ymax = Inf,
           fill = "#FF0000", alpha = 0.1) +
  annotate(geom = 'text', label = paste("CDD: ", cdd_med),
           x = cdd_med + 600, y = 7.5, color = "#FF0000",
           family = formals(theme_ecodata)$base_family) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  xlim(0,4500) +
  labs(y = "Estimated Sampling Abundance",
       title = "Estimated Aphid densities",
       subtitle = "Data aggregated across sites and years from Idaho Pan Trapping.\nRed line is median peak, pink band is 95% credible interval.") +
  theme()
ggsave(filename="agg_gam.png",path="./aphids_2023/figures",
       width=15,height=11,units="cm",dpi=300,device=ragg::agg_png())


#### Model B: Pan trap, Vetch DD, currently unused ####
mb <- brm(AphidCount ~ s(Vetch_DD),
          family=negbinomial,
          data = dat,
          chains = 4, cores = 4,
          iter = 5000, warmup = 1000,
          control = list(adapt_delta = 0.95,
                         max_treedepth = 15)
)

write_rds(mb, "./aphids_2023/fit_models/pooled_gam_B.rds")

# extract peaks

mb_p <- data.frame(Vetch_DD = seq(from = 0, to = max(dat$Vetch_DD, na.rm = T), by = 1)) %>%
  add_fitted_draws(mb) %>%
  median_hdci()

mb_peak_dist <- data.frame(Vetch_DD = seq(from = 0, to = max(dat$Vetch_DD, na.rm = T), by = 1)) %>%
  add_fitted_draws(mb) %>%
  ungroup() %>%
  group_by(.draw) %>%
  slice_max(.value, n = 1) %>%
  ungroup() %>%
  median_hdci()

cdd_med <- mb_peak_dist$Vetch_DD
cdd_lo <- mb_peak_dist$Vetch_DD.lower
cdd_hi <- mb_peak_dist$Vetch_DD.upper

rm(mb_peak_dist)
rm(mb)

pb <- mb_p %>%
  ggplot(aes(x = Vetch_DD, y = .value, ymin = .lower, ymax = .upper)) +
  geom_vline(xintercept = cdd_med, color = "#FF0000", linetype = 2) +
  annotate(geom = 'rect', xmin = cdd_lo, xmax = cdd_hi,
           ymin = -Inf, ymax = Inf,
           fill = "#FF0000", alpha = 0.1) +
  annotate(geom = 'text', label = paste("CDD: ", cdd_med),
           x = cdd_med + 600, y = 7.5, color = "#FF0000",
           family = formals(theme_ecodata)$base_family) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  xlim(0,4500) +
  labs(y = "Estimated Sampling Abundance",
       title = "Estimated Aphid densities",
       subtitle = "Data aggregated across sites and years from Idaho Pan Trapping.\nRed line is median peak, pink band is 95% credible interval.") +
  theme()
ggsave(filename="agg_gam_b.png",path="./aphids_2023/figures",
       width=15,height=11,units="cm",dpi=300,device=ragg::agg_png())

rm(mb_p)

#### Model C: Sweep, Aphid DD ####
mc <- brm(AphidCount ~ s(Aphids_DD),
          family=negbinomial,
          data = vdat,
          chains = 4, cores = 4,
          iter = 5000, warmup = 1000,
          seed = 42,
          control = list(adapt_delta = 0.99,
                         max_treedepth = 15)
)
# 1 divergent transition after warmup. If you only get a few divergences and a good Rhat and ESS, resulting posterior is often good enough to move forward.

#write_rds(mc, "./aphids_2023/fit_models/pooled_gam_c.rds")

# Extract peak
# Get curve and 95% CIs for plotting
mc_p <- data.frame(Aphids_DD = seq(from = 0, to = max(vdat$Aphids_DD, na.rm = T), by = 1)) %>%
  add_fitted_draws(mc) %>%
  median_hdci()

# Get curve and 95% CI for the peak
# try using max of 2300 to get sensible upper CI (the upper CI goes berserk above that point, right near the end of the data)
# subbed out max(vdat$Aphids_DD, na.rm = T)
mc_peak_dist <- data.frame(Aphids_DD = seq(from = 0, to = 2300, by = 1)) %>%
  add_fitted_draws(mc) %>%
  ungroup() %>%
  group_by(.draw) %>%
  slice_max(.value, n = 1) %>%
  ungroup() %>%
  median_hdci()

# 1602, 1660 2615
# Sparse data at the upper end of the x values causes CIs to blow out. Might try better priors to tighten them.
# Simply truncating the series at 2300 CDD yields the following (used in predictions, for now):
# 1578, 1654, 1771
cdd_med <- mc_peak_dist$Aphids_DD
cdd_lo <- mc_peak_dist$Aphids_DD.lower
cdd_hi <- mc_peak_dist$Aphids_DD.upper

# Remove what we don't need anymore
# rm(mc_peak_dist)
# rm(mc)

# Plot
pc <- mc_p %>%
  ggplot(aes(x = Aphids_DD, y = .value, ymin = .lower, ymax = .upper)) +
  geom_vline(xintercept = cdd_med, color = "#FF0000", linetype = 2) +
  annotate(geom = 'rect', xmin = cdd_lo, xmax = cdd_hi,
           ymin = -Inf, ymax = Inf,
           fill = "#FF0000", alpha = 0.1) +
  annotate(geom = 'text', label = paste("CDD: ", cdd_med),
           x = cdd_med + 200, y = 1000, color = "#FF0000",
           family = formals(theme_ecodata)$base_family) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  labs(y = "Estimated Sampling Abundance",
       title = "Estimated Aphid densities",
       subtitle = "Data aggregated across sites and years from Vetch Sweeping.\nRed line is median peak, pink band is 95% credible interval.") +
  theme()
ggsave(filename="agg_gam_c.png",path="./aphids_2023/figures",
       width=15,height=11,units="cm",dpi=300,device=ragg::agg_png())

rm(mc_p)


#### Model D: Sweep, Vetch DD, currently unused ####
md <- brm(AphidCount ~ s(Vetch_DD),
          family=negbinomial,
          data = vdat,
          chains = 4, cores = 4,
          iter = 5000, warmup = 1000,
          control = list(adapt_delta = 0.95,
                         max_treedepth = 15)
)
# 1 divergent transition after warmup. If you only get a few divergences and a good Rhat and ESS, resulting posterior is often good enough to move forward.

write_rds(md, "./aphids_2023/fit_models/pooled_gam_d.rds")

# extract peaks

md_p <- data.frame(Vetch_DD = seq(from = 0, to = max(vdat$Vetch_DD, na.rm = T), by = 1)) %>%
  add_fitted_draws(md) %>%
  median_hdci()

md_peak_dist <- data.frame(Vetch_DD = seq(from = 0, to = max(vdat$Vetch_DD, na.rm = T), by = 1)) %>%
  add_fitted_draws(md) %>%
  ungroup() %>%
  group_by(.draw) %>%
  slice_max(.value, n = 1) %>%
  ungroup() %>%
  median_hdci()

cdd_med <- md_peak_dist$Vetch_DD
cdd_lo <- md_peak_dist$Vetch_DD.lower
cdd_hi <- md_peak_dist$Vetch_DD.upper

rm(md_peak_dist)
rm(md)

pd <- md_p %>%
  ggplot(aes(x = Vetch_DD, y = .value, ymin = .lower, ymax = .upper)) +
  geom_vline(xintercept = cdd_med, color = "#FF0000", linetype = 2) +
  annotate(geom = 'rect', xmin = cdd_lo, xmax = cdd_hi,
           ymin = -Inf, ymax = Inf,
           fill = "#FF0000", alpha = 0.1) +
  annotate(geom = 'text', label = paste("CDD: ", cdd_med),
           x = cdd_med + 200, y = 1000, color = "#FF0000",
           family = formals(theme_ecodata)$base_family) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
#  xlim(0,4500) +
  labs(y = "Estimated Sampling Abundance",
       title = "Estimated Aphid densities",
       subtitle = "Data aggregated across sites and years from Idaho Pan Trapping.\nRed line is median peak, pink band is 95% credible interval.") +
  theme()
ggsave(filename="agg_gam_d.png",path="./aphids_2023/figures",
       width=15,height=11,units="cm",dpi=300,device=ragg::agg_png())

rm(md_p)
