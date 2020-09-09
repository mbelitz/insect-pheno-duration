
## try the ROTL package
library(rotl)
nmz <- tnrs_match_names(names = unique(model_df4$scientificName), context_name = "Animals")
nmz2 <- filter(nmz, ott_id != 7146697)
insect_tree <- tol_induced_subtree(ott_ids = nmz2$ott_id)
insect_tree$tip.label <- word(insect_tree$tip.label, start = 1, end = 2, sep = "_") %>% 
  sub(pattern = "_", replacement = " ")

plot(insect_tree, type = "fan")
plot(insect_tree)


insect_tree_bl <- ape::compute.brlen(insect_tree)

insect_tree_bl

on_coef
physig_col = function(inte){
  col_to = names(inte)[-1]
  out = vector("list", length = length(col_to))
  names(out) = col_to
  for(i in col_to){
    x = inte[[i]]
    names(x) = inte$scientificName
    xp1 = phytools::phylosig(insect_tree_bl, x, method = "K", test = T)
    xp2 = phytools::phylosig(insect_tree_bl, x, method = "lambda", test = T)
    out[[i]] = tibble::tibble(statistic = c(xp1$K, xp2$lambda),
                              P = c(ifelse(xp1$P < 0.5, xp1$P, 1 - xp1$P), xp2$P),
                              test = c("K", "Lambda"))
  }
  bind_rows(out, .id = "terms")
}

inte = left_join(dplyr::select(on_coef, scientificName, intercept_ave_onset),
                 dplyr::select(off_coef, scientificName, intercept_ave_offset)) %>% 
  left_join(dplyr::select(dur_coef, scientificName, intercept_ave_duration))

physig_intercept = physig_col(inte) %>% 
  arrange(test) %>% 
  mutate(statistic = round(statistic, 4),
         P = round(P, 3)) 
physig_intercept

physig_slopes = left_join(dplyr::select(on_coef, scientificName, temp_onset = temp, prec_onset = prec, 
                                        temp_seas_onset = temp_seas, tempprecint_onset = temp:prec),
                          dplyr::select(off_coef, scientificName, prec_offset = prec, temp_seas_offset = temp_seas,
                                        )) %>% 
  left_join(dplyr::select(dur_coef, scientificName, temp_dur = temp, prec_dur = prec, 
                          temp_seas_dur = temp_seas, tempprecint_dur = temp:prec,
                          )) %>% 
  physig_col() %>% 
  arrange(test) %>% 
  mutate(statistic = round(statistic, 4),
         P = round(P, 3)) 
as.data.frame(physig_slopes)



####### PGLMM

# onset 
library(phyr)
library(INLA)

# Burnsius communis not in phylogeny, remove

model_df5 <- filter(model_df4, scientificName %in% insect_tree_bl$tip.label)

a <- insect_tree_bl$tip.label %in% model_df5$scientificName
a[140]
a[32]

insect_tree_bl$tip.label[32]
insect_tree_bl$tip.label[140]


insect_tree_bl2 <- ape::drop.tip(insect_tree_bl, c("Enodia anthedon", "Aeshna cyanea")) # = (C:1,(D:1,E:1):1);

testm <- pglmm(onset ~ temp + (1 | scientificName__),
               data = model_df5,   
               cov_ranef = list(scientificName = insect_tree_bl2), 
               bayes = T)



pm_onset <- pglmm(onset ~ temp + prec + temp_seas + temp:prec +
                            diapause.stage + immature.habitat +
                            (1 | id_cells) + (1 | scientificName__) + 
                             (0 + temp | scientificName__) + 
                             (0 + prec | scientificName__) +
                             (0 + temp_seas | scientificName__)  +    
                            temp_seas:diapause.stage + 
                             temp:immature.habitat + 
                             prec:immature.habitat,
                          data = model_df5, 
                          cov_ranef = list(scientificName = insect_tree_bl2), 
                          bayes = T)


ranef(pm_onset)
fixef(pm_onset) %>% knitr::kable()
 inla_onset = pm_onset$inla.model
 summary(inla_onset)
 inla_onset$summary.fixed # kld are small, good
 plot(inla_onset$marginals.fixed$`(Intercept)`)
 plot(inla_onset$marginals.fixed$temp)
 plot(inla_onset$marginals.fixed$`temp:prec`)
# 
 names(inla_onset$marginals.fixed)
 names(inla_onset$marginals.hyperpar)
 length(inla_onset$summary.random)
 inla_onset$marginals.random
 
 install.packages("remotes")
 remotes::install_github("julianfaraway/brinla")
 library(brinla)
 
# bri.hyperpar.plot(inla_onset, F)
# bri.lmresid.plot(inla_onset)
 inla_onset$marginals.fitted.values
 invsqrt <- function(x) 1/sqrt(x)
 invsqrt(inla_onset$summary.hyperpar[, -2]) # SD of random terms
# bri.hyperpar.summary(inla_onset)
 inla_onset$marginals.random # species level random term
# bri.random.plot(inla_onset)
# 
# # phylo and spatial
 
 pm_onset2 <- pglmm(onset ~ temp + prec + temp_seas + temp:prec +
                     diapause.stage + immature.habitat +
                     (1 | id_cells) + (1 | scientificName__) + 
                     (0 + temp | scientificName__) + 
                     (0 + prec | scientificName__) +
                     (0 + temp_seas | scientificName__)  +    
                     temp_seas:diapause.stage + 
                     temp:immature.habitat + 
                     prec:immature.habitat,
                   data = model_df5, 
                   cov_ranef = list(id_cells = V.space, scientificName = insect_tree_bl2), 
                   bayes = T)
