
library(readr)     # load data
library(survival)  # base survival (coxph, etc)
library(survminer) # better plotting of KM curves
library(tidyverse) # data manipulation

# load data
data_url <- 'https://raw.githubusercontent.com/hammerlab/paper-aocs-chemo-neoantigens/master/additional-files/Additional%20File%201.csv'
d <- readr::read_csv(data_url) 

## modified or mutated data stored in md
above_quantile <- function(x, q, ...) {
  x > quantile(x, probs = q, ...)
}

above_quantile_str <- function(x, q, ...) {
  threshold <- as.integer(quantile(x, probs = q, ...))
  ifelse(x > threshold, stringr::str_c('>', threshold), stringr::str_c('<=', threshold))
}

md <- d %>%
  dplyr::mutate(event = donor_vital_status == 'deceased',
                mutations_excl = ifelse(mutations >= 30000, NA, mutations)) %>%
  dplyr::filter(specific_treatment == 'primary/untreated' & tissue_type == 'solid') %>%
  tidyr::drop_na(event, donor_survival_time, tumour_stage, mutations) %>%
  dplyr::mutate_each(funs = funs(above_median = above_quantile(., q = 0.5, na.rm = T),
                                 above_median_str = above_quantile_str(., q = 0.5, na.rm = T),
                                 above_p75 = above_quantile(., q = 0.75, na.rm = T),
                                 above_p75_str = above_quantile_str(., q = 0.75, na.rm = T)
                                 ),
                     mutations, peptides, mutations_excl)

md_all <- d %>%
  dplyr::mutate(event = donor_vital_status == 'deceased') %>%
  tidyr::drop_na(event, donor_survival_time, tumour_stage, mutations) %>%
  dplyr::mutate_each(funs = funs(above_median = above_quantile(., q = 0.5),
                                 above_median_str = above_quantile_str(., q = 0.5),
                                 above_p75 = above_quantile(., q = 0.75),
                                 above_p75_str = above_quantile_str(., q = 0.75)
  ),
  mutations, peptides)

# ---- graphical analysis ---- 

# null model
m0 <- survival::survfit(Surv(time = donor_survival_time, event = event) ~ 1, data = md)
# plot KM curves for this cohort
survminer::ggsurvplot(m0, risk.table = TRUE) 

# plot by tumor stage
m1 <- update(m0, . ~ . + tumour_stage, data = md)
survminer::ggsurvplot(m1, risk.table = TRUE, conf.int = TRUE) 

# plot by mutation count <> median
m2 <- update(m0, . ~ . + mutations_above_median, data = md)
survminer::ggsurvplot(m2, risk.table = TRUE, conf.int = TRUE)

# plo by mutation count <> p75
m3 <- update(m0, . ~ . + mutations_above_p75, data = md)
survminer::ggsurvplot(m3, risk.table = TRUE, conf.int = TRUE)

# ---- review association results ---- 

# null model 
mph0 <- survival::coxph(Surv(time = donor_survival_time, event = event) ~ 1, data = md)

# ---- look at tumour stage ---- 
mph1 <- update(mph0, . ~ . + tumour_stage)
summary(mph1)

# do we have non-proportional hazards with tumour stage?
cox.zph(mph1)

# look at residuals - are any observations overly influential?
mph1.dfbeta <- residuals(mph1, type="dfbeta")
plot(mph1.dfbeta)

# is there heterogeneity in tumor stage?
md %>% 
  dplyr::group_by(tumour_stage, tumour_stage_system, tumour_grade, tumour_histological_type, metastasis) %>% 
  dplyr::tally()

md %>% 
  dplyr::group_by(tumour_stage, tumour_stage_system, tumour_grade, tumour_histological_type, metastasis) %>% 
  dplyr::summarise_each(funs = funs(mean(. , na.rm = T), pct = length(na.omit(.))/n(), n()),
                        percentage_cellularity)


# ---- also adjust for grade? --- 
mph2 <- update(mph1, . ~ . + tumour_grade)

cox.zph(mph2, transform = log)
mph2.dfbeta <- residuals(mph2, type="dfbeta")
par(mfrow=c(1, 2))
for (j in 1:2) {
  plot(mph2.dfbeta[, j], ylab=names(coef(mph2))[j])
  abline(h=0, lty=2)
  }

# ---- also adjust for cellularity? ---- 

mph3 <- update(mph1, . ~ . + percentage_cellularity)

# non-prop hazards
cox.zph(mph3, transform = log)

# overly influential?
mph3.dfbeta <- residuals(mph3, type="dfbeta")
par(mfrow=c(1, 2))
for (j in 1:2) {
  plot(mph3.dfbeta[, j], ylab=names(coef(mph3))[j])
  abline(h=0, lty=2)
}
dev.off() ## reset plot area

# non-linearity?
mph3.martin <- residuals(mph3, type="martingale")
mph3.X <- as.matrix(md %>% tidyr::drop_na(percentage_cellularity) %>% dplyr::select(percentage_cellularity))
plot(mph3.X[, 1], mph3.martin, xlab="percentage_cellularity", ylab="residuals")
abline(h=0, lty=2)
lines(lowess(mph3.X[, 1], mph3.martin, iter=0))


# ---- start analysis of mutations ---- 

mph4 <- update(mph0, . ~ . + tumour_stage + percentage_cellularity + mutations)

# non-prop hazards
cox.zph(mph4, transform = log)

# overly influential?
mph4.dfbeta <- residuals(mph4, type="dfbeta")
par(mfrow=c(1, 3))
for (j in 1:3) {
  plot(mph4.dfbeta[, j], ylab=names(coef(mph4))[j])
  abline(h=0, lty=2)
}

# non-linearity?
mph4.martin <- residuals(mph4, type='martingale')
mph4.X <- as.matrix(md %>% tidyr::drop_na(percentage_cellularity) %>% dplyr::select(percentage_cellularity, mutations))
par(mfrow=c(1, 2))
for (j in 1:2) {
  plot(mph4.X[, j], mph4.martin, xlab=c('percentage_cellularity', 'mutations')[j], ylab="residuals")
  abline(h=0, lty=2)
  lines(lowess(mph4.X[, j], mph4.martin, iter=0))
}
# this confirms that we should probably drop this observation with very high mutation count
dev.off() ## reset plot area

## ---- ... drop v high mutation count ---- 

mph4e <- update(mph0, . ~ . + tumour_stage + percentage_cellularity + mutations_excl)

# non-prop hazards
cox.zph(mph4e, transform = log)

# overly influential?
mph4e.dfbeta <- residuals(mph4e, type="dfbeta")
par(mfrow=c(1, 3))
for (j in 1:3) {
  plot(mph4e.dfbeta[, j], ylab=names(coef(mph4e))[j])
  abline(h=0, lty=2)
}

# non-linearity?
mph4e.martin <- residuals(mph4e, type='martingale')
mph4e.X <- as.matrix(md %>% tidyr::drop_na(percentage_cellularity, mutations_excl) %>% dplyr::select(percentage_cellularity, mutations_excl))
par(mfrow=c(1, 2))
for (j in 1:2) {
  plot(mph4e.X[, j], mph4e.martin, xlab=c('percentage_cellularity', 'mutations_excl')[j], ylab="residuals")
  abline(h=0, lty=2)
  lines(lowess(mph4e.X[, j], mph4e.martin, iter=0))
}
# this confirms that we should probably drop this observation with very high mutation count
dev.off() ## reset plot area

## ---- stratify by tumour stage ---- 

mph4s <- update(mph0, . ~ . + strata(tumour_stage) + percentage_cellularity + mutations_excl)

# non-prop hazards
cox.zph(mph4s, transform = log)

# overly influential?
mph4s.dfbeta <- residuals(mph4s, type="dfbeta")
par(mfrow=c(1, 2))
for (j in 1:2) {
  plot(mph4s.dfbeta[, j], ylab=names(coef(mph4s))[j])
  abline(h=0, lty=2)
}

# non-linearity?
mph4s.martin <- residuals(mph4s, type='martingale')
mph4s.X <- as.matrix(md %>% tidyr::drop_na(percentage_cellularity, mutations_excl) %>% dplyr::select(percentage_cellularity, mutations_excl))
par(mfrow=c(1, 2))
for (j in 1:2) {
  plot(mph4s.X[, j], mph4s.martin, xlab=c('percentage_cellularity', 'mutations_excl')[j], ylab="residuals")
  abline(h=0, lty=2)
  lines(lowess(mph4s.X[, j], mph4s.martin, iter=0))
}
# this confirms that we should probably drop this observation with very high mutation count
dev.off() ## reset plot area

## ---- use buckets of cellularity ---- 

mph5 <- update(mph0, . ~ . + strata(tumour_stage) + level_of_cellularity + mutations_excl)

# non-prop hazards
cox.zph(mph5, transform = log)

# overly influential?
mph5.dfbeta <- residuals(mph5, type="dfbeta")
par(mfrow=c(1, length(coef(mph5))))
for (j in 1:length(coef(mph5))) {
  plot(mph5.dfbeta[, j], ylab=names(coef(mph5))[j])
  abline(h=0, lty=2)
}

# non-linearity?
mph5.martin <- residuals(mph5, type='martingale')
mph5.X <- as.matrix(md %>% tidyr::drop_na(percentage_cellularity, mutations_excl) %>% dplyr::select(percentage_cellularity, mutations_excl))
par(mfrow=c(1, 2))
for (j in 1:2) {
  plot(mph5.X[, j], mph5.martin, xlab=c('percentage_cellularity', 'mutations_excl')[j], ylab="residuals")
  abline(h=0, lty=2)
  lines(lowess(mph5.X[, j], mph5.martin, iter=0))
}
# this confirms that we should probably drop this observation with very high mutation count
dev.off() ## reset plot area

## ---- look at neoantigens (peptides) ---- 

mph6 <- update(mph0, . ~ . + strata(tumour_stage) + level_of_cellularity + peptides)

# non-prop hazards
cox.zph(mph6, transform = log)

# overly influential?
mph6.dfbeta <- residuals(mph6, type="dfbeta")
par(mfrow=c(1, length(coef(mph6))))
for (j in 1:length(coef(mph6))) {
  plot(mph6.dfbeta[, j], ylab=names(coef(mph6))[j])
  abline(h=0, lty=2)
}

# non-linearity?
mph6.martin <- residuals(mph6, type='martingale')
mph6.X <- as.matrix(md %>% tidyr::drop_na(percentage_cellularity, peptides) %>% dplyr::select(percentage_cellularity, peptides))
par(mfrow=c(1, 2))
for (j in 1:2) {
  plot(mph6.X[, j], mph6.martin, xlab=c('percentage_cellularity', 'peptides')[j], ylab="residuals")
  abline(h=0, lty=2)
  lines(lowess(mph6.X[, j], mph6.martin, iter=0))
}
# this confirms that we should probably drop this observation with very high mutation count
dev.off() ## reset plot area

