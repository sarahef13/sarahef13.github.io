library(tidyverse)
library(stats)
library(easystats)
library(modelr)

# NOTE: data set is too small to split into train vs. test subsets for cross-validation

datafile <- list.files(pattern = 'all_norm_wt_data.csv', recursive = TRUE)
df <- read_csv(datafile)
df <- df %>% 
      mutate(species = as.factor(species))

#### RANDOM FOREST PREDICTORS ####

mir_315a_mod <- glm(data = df, formula = adjusted_devel ~ mir_315a * species)
report(mir_315a_mod)
# R2 = 0.75

mir_305_mod <- glm(data = df, formula = adjusted_devel ~ mir_305 * species)
report(mir_305_mod)
# R2 = 0.71

mir_375_mod <- glm(data = df, formula = adjusted_devel ~ mir_375 * species)
report(mir_375_mod)
# R2 = 0.72

mir_9a_mod <- glm(data = df, formula = adjusted_devel ~ mir_9a * species)
report(mir_9a_mod)
# R2 = 0.55; "The effect of mir 9a is statistically non-significant"

all_rf_mir_mod <- glm(data = df, 
                      formula = adjusted_devel ~ (mir_315a + mir_305 + mir_375 + mir_9a) * species)
report(all_rf_mir_mod)
# R2 = 0.92

mir_315a_375_mod <- glm(data = df, formula = adjusted_devel ~ (mir_315a + mir_375) * species)
report(mir_315a_375_mod)
# R2 = 0.87

mir_315a_305_mod <- glm(data = df, formula = adjusted_devel ~ (mir_315a + mir_305) * species)
report(mir_315a_305_mod)
# R2 = 0.84

comparison <- compare_performance(mir_315a_mod, all_rf_mir_mod, mir_315a_375_mod,
                                  rank=TRUE)
plot(comparison)
# mir_315a_mod is poor all around
# all_rf_mir_mod has better result stats, but worse complexity stats
# mir_315a_375_mod has okay result stats, but good complexity stats

mir_315a_375_noluci <- glm(data = (df %>% filter(species != 'Lucilia')),
                           formula = adjusted_devel ~ (mir_315a + mir_375) * species)
report(mir_315a_375_noluci)
# R2 = 0.88

mir_315a_375_noluci_simple <- glm(data = (df %>% filter(species != 'Lucilia')),
                           formula = adjusted_devel ~ mir_315a + mir_375 + species)
report(mir_315a_375_noluci_simple)
# R2 = 0.63

mir_all_rf_noluci_nospecies <- glm(data = (df %>% filter(species != 'Lucilia')),
                                  formula = adjusted_devel ~ mir_315a + mir_375 + mir_305 + mir_9a)
report(mir_all_rf_noluci_nospecies)
# R2 = 0.58

all_mirs_noluci_nospecies <- glm(data = (df %>% filter(species != 'Lucilia')),
                                  formula = adjusted_devel ~ mir_315a + mir_375 + mir_305 + mir_9a + mir_bft + mir_let_7 + mir_100 + mir_125)
report(all_mirs_noluci_nospecies)
# R2 = 0.82

comparison2 <- compare_performance(all_rf_mir_mod, mir_315a_375_mod, mir_315a_375_noluci,
                                   all_mirs_noluci_nospecies,
                                   rank=TRUE)
plot(comparison2)
# noluci has better result stats than its counterpart, but still worse than all_rf_mir_mod
# combining 8 miRs with no species data is worse than using 2 miRs with all species


model_results_noluci <- 
  df %>% 
  filter(species != 'Lucilia') %>% 
  gather_predictions(all_rf_mir_mod, mir_315a_375_noluci, mir_315a_mod,
                     type='response') %>% 
  ggplot(aes(x=adjusted_devel, y=pred, color=model)) +
  geom_segment(aes(x=-0.25, y=-0.25, xend=1.25, yend=1.25),linetype=2, color="hotpink",alpha=.5) +
  geom_smooth(method = "lm", se=FALSE, alpha=.5) +
  geom_point(alpha=0.5) +
  scale_color_viridis_d(end=0.8) +
  labs(title = "Actual vs Predicted Adjusted Puparial Development, without Lucilia Data",
       subtitle = "(dashed line indicates perfect overlap between actual & predicted values)")

model_results_luci <- 
  df %>% 
  gather_predictions(all_rf_mir_mod, mir_315a_375_mod, mir_315a_mod,
                     type='response') %>% 
  ggplot(aes(x=adjusted_devel, y=pred, color=model)) +
  geom_segment(aes(x=-0.25, y=-0.25, xend=1.25, yend=1.25),linetype=2, color="hotpink",alpha=.5) +
  geom_smooth(method = "lm", se=FALSE, alpha=.5) +
  geom_point(alpha=0.5) +
  scale_color_viridis_d(end=0.8) +
  labs(title = "Actual vs Predicted Adjusted Puparial Development, with Lucilia Data",
       subtitle = "(dashed line indicates perfect overlap between actual & predicted values)")

model_results_faceted <- 
  df %>% 
  gather_predictions(all_rf_mir_mod, mir_315a_375_mod, mir_315a_mod,
                     type='response') %>% 
  ggplot(aes(x=adjusted_devel, y=pred, color=model)) +
  geom_segment(aes(x=-0.25, y=-0.25, xend=1.25, yend=1.25),linetype=2, color="hotpink",alpha=.5) +
  geom_smooth(method = "lm", se=FALSE, alpha=.5) +
  geom_point(alpha=0.5) +
  scale_color_viridis_d(end=0.8) +
  labs(title = "Actual vs Predicted Adjusted Puparial Development, Faceted by Species",
       subtitle = "(dashed line indicates perfect overlap between actual & predicted values)") +
  facet_wrap(~species)

pdf(paste0(getwd(), '/glm_prediction_luci_impact.pdf'))
model_results_noluci
model_results_luci
model_results_faceted
dev.off()



#### DESEQ PREDICTORS ####

mir_bft_mod <- glm(data = df, formula = adjusted_devel ~ mir_bft * species)
report(mir_bft_mod)
# R2 = 0.07

mir_bft_mod_parab <- glm(data = df, formula = adjusted_devel ~ poly(mir_bft, 2) * species)
report(mir_bft_mod_parab)
# R2 = 0.17


mir_let7_mod <- glm(data = df, formula = adjusted_devel ~ mir_let_7 * species)
report(mir_let7_mod)
# R2 = 0.34

mir_let7_mod_parab <- glm(data = df, formula = adjusted_devel ~ poly(mir_let_7, 2) * species)
report(mir_let7_mod_parab)
# R2 = 0.69


mir_100_mod <- glm(data = df, formula = adjusted_devel ~ mir_100 * species)
report(mir_100_mod)
# R2 = 0.20
# note: mir_100 didn't seem resistant to environment changes

mir_100_mod_parab <- glm(data = df, formula = adjusted_devel ~ poly(mir_100, 2) * species)
report(mir_100_mod_parab)
# R2 = 0.45
# note: mir_100 didn't seem resistant to environment changes


mir_125_mod <- glm(data = df, formula = adjusted_devel ~ mir_125 * species)
report(mir_125_mod)
# R2 = 0.53
# note: mir_125 didn't seem resistant to environment changes

mir_125_mod_parab <- glm(data = df, formula = adjusted_devel ~ poly(mir_125, 2) * species)
report(mir_125_mod_parab)
# R2 = 0.61
# note: mir_125 didn't seem resistant to environment changes


mir_bft_let7_mod_parab <- glm(data = df, 
                              formula = adjusted_devel ~ (poly(mir_bft, 2) + poly(mir_let_7, 2)) * species)
report(mir_bft_let7_mod_parab)
# R2 = 0.84


all_deseq_mir_mod <- glm(data = df, 
                      formula = adjusted_devel ~ (mir_bft + mir_let_7 + mir_100 + mir_125) * species)
report(all_deseq_mir_mod)
# R2 = 0.75


comparison3 <- compare_performance(mir_315a_mod, mir_let7_mod_parab, mir_bft_let7_mod_parab,
                                   mir_125_mod_parab, all_deseq_mir_mod, rank=TRUE)
plot(comparison3)
# DEDeq miR mods need to be complex & involve parabolas to have 
# a predictive ability like the simplest Random Forest miR mod's
# Also, combining all four miRs made the model worse, somehow



#### RANDOM FOREST PREDICTORS, NOW WITH PARABOLAS/CUBES ####

mir_315a_mod_parab <- glm(data = df, formula = adjusted_devel ~ poly(mir_315a, 2) * species)
report(mir_315a_mod_parab)
# R2 = 0.78; linear version was 0.75

mir_315a_mod_cube <- glm(data = df, formula = adjusted_devel ~ poly(mir_315a, 3) * species)
report(mir_315a_mod_cube)
# R2 = 0.79; linear version was 0.75

# actual expression graph for miR 315a has a different line shape per species


mir_305_mod_parab <- glm(data = df, formula = adjusted_devel ~ poly(mir_305, 2) * species)
report(mir_305_mod_parab)
# R2 = 0.78; linear version was 0.71

mir_305_mod_cube <- glm(data = df, formula = adjusted_devel ~ poly(mir_305, 3) * species)
report(mir_305_mod_cube)
# R2 = 0.82; linear version was 0.71

# actual expression graph for miR 305 looks like shallow cube curves


mir_375_mod_parab <- glm(data = df, formula = adjusted_devel ~ poly(mir_375, 2) * species)
report(mir_375_mod_parab)
# R2 = 0.76; linear version was 0.72

mir_375_mod_cube <- glm(data = df, formula = adjusted_devel ~ poly(mir_375, 3) * species)
report(mir_375_mod_cube)
# R2 = 0.79; linear version was 0.72

# actual expression graph for miR 375 looks very cubic


mir_9a_mod_parab <- glm(data = df, formula = adjusted_devel ~ poly(mir_9a, 2) * species)
report(mir_9a_mod_parab)
# R2 = 0.68; linear version was 0.55

mir_9a_mod_cube <- glm(data = df, formula = adjusted_devel ~ poly(mir_9a, 3) * species)
report(mir_9a_mod_cube)
# R2 = 0.70; linear version was 0.55

# actual expression graph for miR 9a looks mostly cubic


comparison4 <- compare_performance(mir_315a_mod, mir_315a_mod_parab, mir_305_mod_cube,
                                   mir_375_mod_cube, mir_9a_mod_cube, rank=TRUE)
plot(comparison4)
# the linear model is still by far the simplest, but the parabolas/cubes predict better
# mir 9a improved the most from linear to cubic, but it was still the worst model
# mir 305 had a large improvement AND became the best predictor among the single-miR models


mir_315a_line_305_cube <- glm(data = df, 
                              formula = adjusted_devel ~ (mir_315a + poly(mir_305, 3)) * species)
report(mir_315a_line_305_cube)
# R2 = 0.90

mir_315a_line_375_cube <- glm(data = df, 
                              formula = adjusted_devel ~ (mir_315a + poly(mir_375, 3)) * species)
report(mir_315a_line_375_cube)
# R2 = 0.90

mir_315a_305_375_9a_complex <- glm(data = df, 
                                     formula = adjusted_devel ~ (mir_315a + poly(mir_305, 3) + poly(mir_375, 3) + poly(mir_9a, 3)) * species)
report(mir_315a_305_375_9a_complex)
# R2 = 0.98

comparison5 <- compare_performance(mir_305_mod_cube, mir_315a_375_mod, 
                                   mir_315a_line_305_cube, mir_315a_line_375_cube, 
                                   all_rf_mir_mod, mir_bft_let7_mod_parab, 
                                   mir_315a_305_375_9a_complex,
                                   rank=TRUE)
plot(comparison5)
# no apparent difference between combining linear 315a with cubic 305 vs. with cubic 375
# all_rf_mir_mod is still the best simple-ish predictor, but 315a with 305^3 or with 375^3 comes close
# 315a + 375 (both linear) has the best complexity score by far and is the 4th best predictor
# combining all 4 miRs with complex curves is INCREDIBLE predictor, but probably way overfit


model_results_cubic <- 
  df %>% 
  gather_predictions(all_rf_mir_mod, mir_315a_375_mod, 
                     mir_315a_line_375_cube, #mir_315a_line_305_cube, # indistinguishable
                     type='response') %>% 
  ggplot(aes(x=adjusted_devel, y=pred, color=model)) +
  geom_segment(aes(x=-0.25, y=-0.25, xend=1.25, yend=1.25),linetype=2, color="hotpink",alpha=.5) +
  geom_smooth(method = "lm", se=FALSE, alpha=.5) +
  geom_point(alpha=0.5) +
  scale_color_viridis_d(end=0.8) +
  labs(title = "Actual vs Predicted Adjusted Puparial Development, with Varying Exponents",
       subtitle = "(dashed line indicates perfect overlap between actual & predicted values)")

pdf(paste0(getwd(), '/glm_prediction_exponent_impact.pdf'))
model_results_cubic
dev.off()



#### COMBINING DESEQ + RANDOM FOREST PREDICTORS ####

mir_315a_375_bft_let7_complex <- glm(data = df, 
                              formula = adjusted_devel ~ (mir_315a + poly(mir_375, 3) + poly(mir_bft, 2) + poly(mir_let_7, 2)) * species)
report(mir_315a_375_bft_let7_complex)
# R2 = 0.97

mir_315a_375_let7_complex <- glm(data = df, 
                                 formula = adjusted_devel ~ (mir_315a + poly(mir_375, 3) + poly(mir_let_7, 2)) * species)
report(mir_315a_375_let7_complex)
# R2 = 0.94

mir_315a_375_bft_let7_simple <- glm(data = df, 
                              formula = adjusted_devel ~ (mir_315a + mir_375 + mir_bft + mir_let_7) * species)
report(mir_315a_375_bft_let7_simple)
# R2 = 0.90

mir_315a_375_let7_simple <- glm(data = df, 
                                    formula = adjusted_devel ~ (mir_315a + mir_375 + mir_let_7) * species)
report(mir_315a_375_let7_simple)
# R2 = 0.88


comparison6 <- compare_performance(mir_315a_line_375_cube, mir_315a_375_let7_complex,
                                   mir_315a_375_bft_let7_complex, mir_315a_375_let7_simple,
                                   all_rf_mir_mod, rank=TRUE)
plot(comparison6)


model_results_deseq_vs_rf <- 
  df %>% 
  gather_predictions(all_rf_mir_mod, all_deseq_mir_mod,
                     mir_315a_375_bft_let7_complex,
                     mir_315a_305_375_9a_complex,
                     type='response') %>% 
  ggplot(aes(x=adjusted_devel, y=pred, color=model)) +
  geom_segment(aes(x=-0.25, y=-0.25, xend=1.25, yend=1.25),linetype=2, color="hotpink",alpha=.5) +
  geom_smooth(method = "lm", se=FALSE, alpha=.5) +
  geom_point(alpha=0.5) +
  scale_color_viridis_d(end=0.8) +
  labs(title = "Actual vs Predicted Adjusted Puparial Development, Comparing miR Sources",
       subtitle = "(dashed line indicates perfect overlap between actual & predicted values)")

pdf(paste0(getwd(), '/glm_prediction_mir_source_impact.pdf'))
model_results_deseq_vs_rf
dev.off()




