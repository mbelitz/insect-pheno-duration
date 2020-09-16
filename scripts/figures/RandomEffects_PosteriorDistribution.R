library(phyr)
library(ape)
library(tidyverse)
library(stringr)

# read in pglmm data
insect_tree <- read.tree("PGLMM_data/insect_tree4.tre")

mdf <- readr::read_csv("PGLMM_data/model_df_5cells_treespp.csv")

insect_tree$tip.label <- word(insect_tree$tip.label, start = 1, end = 2, sep = "_") %>% 
  sub(pattern = "_", replacement = " ")

# Onset PGLMM

pglmm_onset <- pglmm(onset ~ temp + pop + prec + temp_seas + 
                       temp:prec + temp:pop + 
                       diapause.stage + immature.habitat +
                       temp:immature.habitat +
                       pop:immature.habitat +
                       (1|id_cells) + (1|scientificName__) +
                       (0 + temp | scientificName__) + 
                       (0 + prec | scientificName__),
                     data = mdf, 
                     cov_ranef = list(scientificName = insect_tree), 
                     bayes = TRUE)

# Offset PGLMM

pglmm_offset <- pglmm(offset ~ prec + temp_seas + 
                        seas + diapause.stage + immature.habitat + larval.diet +
                        temp_seas:seas + temp_seas:immature.habitat + 
                        prec:immature.habitat +
                        (1|id_cells) + (1|scientificName__) +
                        (0 + temp_seas | scientificName__) + 
                        (0 + prec | scientificName__),
                      data = mdf, 
                      cov_ranef = list(scientificName = insect_tree), 
                      bayes = TRUE)

## Duration PGLMM

pglmm_dur <- pglmm(duration ~ temp + pop + temp_seas + temp:pop +
                     seas + diapause.stage + flights + immature.habitat + larval.diet +
                     (1|id_cells) + (1|scientificName__) +
                     (0 + temp | scientificName__) +
                     (0 + temp_seas | scientificName__) + 
                     pop:seas + temp_seas:seas + 
                     temp:diapause.stage +
                     temp:flights + temp_seas:flights + 
                     pop:immature.habitat + 
                     temp:larval.diet,
                   data = mdf, 
                   cov_ranef = list(scientificName = insect_tree), 
                   bayes = TRUE)

# Onset Posterior distributions

x <- pglmm_onset
n_samp <- 1000

random_samps <- lapply(x$inla.model$marginals.hyperpar, 
                       function(x) INLA::inla.rmarginal(n_samp, INLA::inla.tmarginal(function(x) sqrt(1 / x), x))) %>% 
  setNames(names(x$random.effects))

random_samps2 <- random_samps[-length(random_samps)] %>% 
  setNames(names(x$random.effects)) %>%
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Random Effects")

fixed_samps <- lapply(x$inla.model$marginals.fixed, function(x) INLA::inla.rmarginal(n_samp, x)) %>% 
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Fixed Effects")

samps <- dplyr::bind_rows(random_samps2, fixed_samps) %>%
  dplyr::mutate(effect_type = factor(effect_type, 
                                     levels = c("Random Effects", "Fixed Effects")))

ci <- samps %>%
  dplyr::group_by(var, effect_type) %>%
  dplyr::summarise(lower = quantile(val, 0.025),
                   upper = quantile(val, 0.975),
                   mean = mean(val),
                   .groups = "drop_last")
sig_vars <- ci %>%
  dplyr::mutate(sig = ifelse(effect_type == "Random Effects",
                             "CI no overlap with zero",
                             ifelse(sign(lower) == sign(upper),
                                    "CI no overlap with zero",
                                    "CI overlaps zero"))) %>%
  dplyr::select(var, sig)

samps <- samps %>%
  dplyr::left_join(sig_vars, by = "var") %>%
  dplyr::group_by(var) %>%
  dplyr::filter(abs(val - mean(val)) < (10 * sd(val))) %>% 
  dplyr::ungroup()

pal <- c("#fc8d62", "#8da0cb")

re_on <- filter(samps, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Onset")
re_ci_on <- filter(ci, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Onset")

on_re_pd <- ggplot(re_on, aes(val, var, height = ..density..)) +
  ggridges::geom_density_ridges(aes(val, var, alpha = sig, fill = sig), 
                                 color = "gray70") +
  geom_point(aes(x = mean, y = var), data = re_ci_on, inherit.aes = FALSE) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = var), data = re_ci_on,
                 inherit.aes = FALSE, height = 0.1) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  scale_fill_manual(values = rev(pal)) +
  ylab("") +
  xlab("Estimate") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) + 
  ggtitle("Onset")

# Offset Posterior distributions

x <- pglmm_offset
n_samp <- 1000

random_samps <- lapply(x$inla.model$marginals.hyperpar, 
                       function(x) INLA::inla.rmarginal(n_samp, INLA::inla.tmarginal(function(x) sqrt(1 / x), x))) %>% 
  setNames(names(x$random.effects))

random_samps2 <- random_samps[-length(random_samps)] %>% 
  setNames(names(x$random.effects)) %>%
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Random Effects")

fixed_samps <- lapply(x$inla.model$marginals.fixed, function(x) INLA::inla.rmarginal(n_samp, x)) %>% 
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Fixed Effects")

samps <- dplyr::bind_rows(random_samps2, fixed_samps) %>%
  dplyr::mutate(effect_type = factor(effect_type, 
                                     levels = c("Random Effects", "Fixed Effects")))

ci <- samps %>%
  dplyr::group_by(var, effect_type) %>%
  dplyr::summarise(lower = quantile(val, 0.025),
                   upper = quantile(val, 0.975),
                   mean = mean(val),
                   .groups = "drop_last")
sig_vars <- ci %>%
  dplyr::mutate(sig = ifelse(effect_type == "Random Effects",
                             "CI no overlap with zero",
                             ifelse(sign(lower) == sign(upper),
                                    "CI no overlap with zero",
                                    "CI overlaps zero"))) %>%
  dplyr::select(var, sig)

samps <- samps %>%
  dplyr::left_join(sig_vars, by = "var") %>%
  dplyr::group_by(var) %>%
  dplyr::filter(abs(val - mean(val)) < (10 * sd(val))) %>% 
  dplyr::ungroup()

pal <- c("#fc8d62", "#8da0cb")

re_off <- filter(samps, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Offset")
re_ci_off <- filter(ci, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Offset")

off_re_pd <- ggplot(re_off, aes(val, var, height = ..density..)) +
  ggridges::geom_density_ridges(aes(val, var, alpha = sig, fill = sig), 
                                color = "gray70") +
  geom_point(aes(x = mean, y = var), data = re_ci_off, inherit.aes = FALSE) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = var), data = re_ci_off,
                 inherit.aes = FALSE, height = 0.1) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  scale_fill_manual(values = rev(pal)) +
  ylab("") +
  xlab("Estimate") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) + 
  ggtitle("Offset")

# Duration Posterior distributions
x <- pglmm_dur
n_samp <- 1000

random_samps <- lapply(x$inla.model$marginals.hyperpar, 
                       function(x) INLA::inla.rmarginal(n_samp, INLA::inla.tmarginal(function(x) sqrt(1 / x), x))) %>% 
  setNames(names(x$random.effects))

random_samps2 <- random_samps[-length(random_samps)] %>% 
  setNames(names(x$random.effects)) %>%
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Random Effects")

fixed_samps <- lapply(x$inla.model$marginals.fixed, function(x) INLA::inla.rmarginal(n_samp, x)) %>% 
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Fixed Effects")

samps <- dplyr::bind_rows(random_samps2, fixed_samps) %>%
  dplyr::mutate(effect_type = factor(effect_type, 
                                     levels = c("Random Effects", "Fixed Effects")))

ci <- samps %>%
  dplyr::group_by(var, effect_type) %>%
  dplyr::summarise(lower = quantile(val, 0.025),
                   upper = quantile(val, 0.975),
                   mean = mean(val),
                   .groups = "drop_last")
sig_vars <- ci %>%
  dplyr::mutate(sig = ifelse(effect_type == "Random Effects",
                             "CI no overlap with zero",
                             ifelse(sign(lower) == sign(upper),
                                    "CI no overlap with zero",
                                    "CI overlaps zero"))) %>%
  dplyr::select(var, sig)

samps <- samps %>%
  dplyr::left_join(sig_vars, by = "var") %>%
  dplyr::group_by(var) %>%
  dplyr::filter(abs(val - mean(val)) < (10 * sd(val))) %>% 
  dplyr::ungroup()

pal <- c("SLATEBLUE2", "tomato1", "springgreen2")

re_dur <- filter(samps, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Duration")
re_ci_dur <- filter(ci, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Duration")

dur_re_pd <- ggplot(re_dur, aes(val, var, height = ..density..)) +
  ggridges::geom_density_ridges(aes(val, var, alpha = sig, fill = sig), 
                                color = "gray70") +
  geom_point(aes(x = mean, y = var), data = re_ci_dur, inherit.aes = FALSE) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = var), data = re_ci_dur,
                 inherit.aes = FALSE, height = 0.1) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  scale_fill_manual(values = rev(pal)) +
  ylab("") +
  xlab("Estimate") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) + 
  ggtitle("Duration")

re_comb <- rbind(re_on, re_off, re_dur)
re_ci_comb <- rbind(re_ci_on, re_ci_off, re_ci_dur)

re_comb$Phenometric <- factor(re_comb$Phenometric, levels = c("Onset", "Offset", "Duration"))

comb_re_pd <- ggplot(re_comb, aes(val, var, height = ..density..)) +
  ggridges::geom_density_ridges(aes(val, var, fill = Phenometric), 
                                alpha = 0.35, color = "gray70") +
  geom_point(aes(x = mean, y = var), data = re_ci_comb, inherit.aes = FALSE,
             alpha = 0.75) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = var), data = re_ci_comb,
                 inherit.aes = FALSE, height = 0.1, alpha = 0.75) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  scale_fill_manual(values = rev(pal)) +
  ylab("") +
  xlab("Estimate") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) 

ggsave(filename = "Tables&Figures/RandomEffects_PosteriorDistribution.png", 
       dpi = 400, height = 8, width = 8)

