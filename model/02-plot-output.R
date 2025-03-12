# Plot model fit and parameters

library(coda)
library(broom.mixed)
library(udunits2)
library(tidyverse)

# Load data
wrc <- read_csv("data_clean/rwc-swp.csv") |> 
  mutate(psi_hm = ud.convert(swp, "MPa", "cmH2O") |> 
           ud.convert("cm", "m"),
         population = factor(population, levels = c("CCR", "JLA", "NRV", "TSZ")))



# Load codas
load(file = "model/coda/jm_coda.Rdata")
load(file = "model/coda/jm_rep.Rdata")

# Tidy parameters
sum_param <- tidyMCMC(jm_coda, 
                      conf.int =  TRUE, 
                      conf.method = "HPDinterval",
                      conf.level = 0.95)

# Tidy replicated and add to original data
sum_rep <- tidyMCMC(jm_rep, 
                    conf.int =  TRUE, 
                    conf.method = "HPDinterval")

pred <- cbind(wrc, sum_rep)


#### Plot model fit ####
m1 <- lm(rwc ~ estimate, data = pred)
summary(m1) # R2 = 0.931

# Observed vs. fitted
pred %>%
  ggplot(aes(x = rwc, y = estimate)) +
  geom_abline(slope = 1, intercept = 0, lty = 1,
              color = "red") +
  geom_pointrange(aes(ymin = conf.low, 
                      ymax = conf.high,
                      color = population)) +
  facet_wrap(~population) +
  theme_bw()

# VG curve observed and fitted
pred %>% 
  ggplot(aes(x = -psi_hm)) +
  geom_point(aes(y = rwc, 
                 color = "Observed")) +
  geom_pointrange(aes(y = estimate, 
                      ymin = conf.low,
                      ymax = conf.high,
                      color = "Predicted")) +
  facet_wrap(~population) +
  # scale_color_manual() +
  scale_x_log10(name = "Pressure Head (-m)",
                breaks = c(0.01, 0.1, 1, 10, 100, 1000),
                labels = c(0.01, 0.1, 1, 10, 100, 1000)) +
  theme_bw()

#### Calculate R2, coverage, and bias for each TRT
pred2 <- pred %>%
  mutate(cover = if_else(rwc >= conf.low & rwc <= conf.high, 1, 0))

# total data frame
tot.df <- data.frame(population = "all",
                     R2 = summary(m1)$adj.r.squared,
                     coverage = mean(pred2$cover),
                     bias = summary(m1)$coef[2,1])

# for each TRT
trt.df <- data.frame(population = levels(pred2$population),
                     R2 = pred2 %>%
                       group_by(population) %>%
                       group_map(~ summary(lm(estimate ~ rwc, data = .x))$adj.r.squared) %>%
                       purrr:::simplify(),
                     coverage = pred2 %>%
                       group_by(population) %>%
                       summarize(coverage = mean(cover)) %>%
                       pull(coverage),
                     bias = pred2 %>%
                       group_by(population) %>%
                       group_map(~ summary(lm(estimate ~ rwc, data = .x))$coef[2,1]) %>%
                       purrr:::simplify())

all.df <- bind_rows(trt.df, tot.df)

# write out model fits
write_csv(all.df, "model/output/model_fit.csv")

# Replot observed vs fitted

trt.df.labels <- trt.df |> 
  mutate(lab = paste0("R^2==", round(R2, 3)))
pred %>%
  ggplot(aes(x = rwc, y = estimate)) +
  geom_abline(slope = 1, intercept = 0, lty = 1,
              color = "red") +
  geom_pointrange(aes(ymin = conf.low, 
                      ymax = conf.high,
                      color = population)) +
  geom_text(data = trt.df.labels, 
            aes(x = 0.22, y = 0,
            label = lab),
            parse = TRUE,
            hjust = 0,
            vjust = 0) +
  facet_wrap(~population) +
  theme_bw()

ggsave(filename = "model/output/fit_by_pop.png",
       height = 4,
       width = 6, 
       units = "in")

#### Plot model parameters ####

pop <- sum_param %>%
  filter(grepl("E\\.", term)) %>%
  mutate(parameter = sub("\\[[0-9]\\]", "", term),
         parameter = sub("E\\.", "", parameter),
         population = "all")

trt <- sum_param %>%
  filter(grepl("^alpha\\.cm", term) |
           grepl("^n\\[", term) |
           grepl("^theta.r" , term)) %>%
  mutate(parameter = sub("\\[[0-9]\\]", "", term),
         population = case_when(grepl("1", term) ~ "CCR",
                               grepl("2", term) ~ "JLA",
                               grepl("3", term) ~ "NRV",
                               grepl("4", term) ~ "TSZ"),
         parameter = factor(parameter, levels = c("theta.r", "alpha.cm", "n")))

ggplot() +
  geom_pointrange(data = trt,
                  aes(x = parameter,
                      y = estimate,
                      ymin = conf.low, 
                      ymax = conf.high,
                      color = population),
                  position = position_dodge(width = 1)) +
  geom_pointrange(data = pop,
                  aes(x = parameter,
                      y = estimate,
                      ymin = conf.low, 
                      ymax = conf.high)) +
  facet_wrap(~parameter, scales = "free") +
  theme_bw(base_size = 14) 

ggsave(filename = "model/output/params_by_pop.png",
       height = 4,
       width = 6, 
       units = "in")

sum_param %>%
  mutate(parameter = sub("\\[[0-9]\\]", "", term),
         TRT = case_when(grepl("1", term) ~ "CCR",
                         grepl("2", term) ~ "JLA",
                         grepl("3", term) ~ "NRV",
                         grepl("4", term) ~ "TSZ"),
         parameter = factor(parameter, levels = c("theta.r", "alpha.cm", "n"))) %>%
  filter(parameter %in% c("alpha.cm", "n", "theta.r")) %>%
  ggplot(aes(x = TRT, y = estimate)) +
  geom_pointrange(aes(ymin = conf.low, 
                      ymax = conf.high)) +
  scale_y_continuous("Posterior estimate") +
  facet_wrap(~parameter, scales = "free") +
  theme_bw(base_size = 14)  

all_params <- bind_rows(pop, trt)
# write out model fits
write_csv(all_params, "model/output/model_params.csv")
