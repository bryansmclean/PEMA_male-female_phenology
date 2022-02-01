##############################################
## combined (male + female) phenology models##
## for Peromyscus maniculatus             	##
## Bryan S. McLean    5 April 2021          ##
##############################################

#library(standardize)
library(lme4)
library(sjPlot)

phen.m <- read.csv("./data/_final_w_clim_/male_data_ClimNA.csv")
phen.f <- read.csv("/Users/mclean/Box/Projects/PEMA_phenology/data/phen.final.data.csv")

##############################################
## data merge
##############################################

## put in a combined data frame, harmonize data fields
phens <- data.frame(
  source = c(phen.m$source, phen.f$source),
  sex = c(rep('M', nrow(phen.m)), rep("F", nrow(phen.f))),
  jday = c(phen.m$julianday, phen.f$jday.corr),
  #month = phen$month,
  #year = phen$year,
  state = c(phen.m$t_state, phen.f$state),
  ecoregion = as.factor(c(phen.m$ecoregion2, phen.f$ecoregion.corr)),
  photoperiod = c(phen.m$photoperiod, phen.f$photoperiod.corr),
  Temp = c(phen.m$Temp_BM, phen.f$Tave.weighted.corr),
  Prec = c(phen.m$Prec_BM, phen.f$PPT.weighted.corr)
  #massing = scale(phen$massing)
)
phens$source[which(phens$source == 'nacsm')] <- 'NACSM'
phens$source[which(phens$source == 'neon')] <- 'NEON'
phens$source[which(phens$source == 'VN' | phens$source == 'vertnet')] <- 'VERTNET'
phens$source[which(phens$source == 'SPEC')] <- 'VERTNET'

## clean incomplete records for GLMMs
no.data <- which(is.na(phens$Temp) | is.na(phens$Prec) | is.na(phens$photoperiod) == TRUE)
phens.glmm <- phens[-no.data,]; rm(no.data)

##############################################
## mixed models, two ways
##############################################

# combined M+F GLMM, sex as fixed effect
mod.quad.full.sex <- standardize::standardize(state ~ 
  (poly(Temp, 2, raw = F) + poly(Prec, 2, raw = F) + poly(photoperiod, 2, raw = F) + sex)^2 + 
  (sex * source) +
  (1 | ecoregion) + 
  (0 + poly(Temp, 2, raw = F) | ecoregion) +
  (0 + poly(Prec, 2, raw = F) | ecoregion) +
  (0 + poly(photoperiod, 2, raw = F) | ecoregion),
  data = phens.glmm,
  family = binomial(link='logit')
  )

# alternative (longhand) syntax to the above, to allow easier effect plotting later
# first specifying polys
Temp_poly <- poly(phens.glmm$Temp, 2, raw = F, simple = T); colnames(Temp_poly) <- c("Temp_1","Temp_2")
Prec_poly <- poly(phens.glmm$Prec, 2, raw = F, simple = T); colnames(Prec_poly) <- c("Prec_1","Prec_2")
photoperiod_poly <- poly(phens.glmm$photoperiod, 2, raw = F, simple = T); colnames(photoperiod_poly) <- c("photoperiod_1","photoperiod_2")
phens.glmm <- cbind(phens.glmm, Temp_poly, Prec_poly, photoperiod_poly)

mod.quad.full.sex <- standardize::standardize(state ~ 
  Temp_1 + Temp_2 +
  Prec_1 + Prec_2 +
  photoperiod_1 + photoperiod_2 + 
  (Temp_1 + Prec_1 + photoperiod_1 + sex)^2 + 
  (Temp_1 + Prec_2 + photoperiod_1 + sex)^2 + 	
  (Temp_1 + Prec_1 + photoperiod_2 + sex)^2 + 
  (Temp_2 + Prec_1 + photoperiod_1 + sex)^2 + 
  (Temp_2 + Prec_2 + photoperiod_1 + sex)^2 + 
  (Temp_2 + Prec_2 + photoperiod_2 + sex)^2 + 
  (sex * source) +
  (1 | ecoregion) + 
  (0 + (Temp_1 + Temp_2) | ecoregion) +
  (0 + (Prec_1 + Prec_2) | ecoregion) +
  (0 + (photoperiod_1 + photoperiod_2) | ecoregion),
  data = phens.glmm,
  family = binomial(link='logit')
  )

# save the R model objects
saveRDS(
  glmer(
    mod.quad.full.sex$formula,
    family = binomial(),
    data = mod.quad.full.sex$data,
    control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))
    ),
  file="/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/glmm_male_female_full_bobyqaopt.rds"
  )
write.table(
  broom::tidy(full.mf),
  file = "/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/glmm_male_female_full_bobyqaopt.csv",
  sep = ",",
  row.names = F,
  quote = F
  )

##############################################
## some effects plots
##############################################
mod <- readRDS(file = "/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/glmm_male_female_full_bobyqaopt.rds")

sjPlot::plot_model(mod, type = "pred", terms = c("Prec_1 [all]", "sex"), colors = c('goldenrod3','dodgerblue4')) + theme_sjplot2()
sjPlot::plot_model(mod, type = "pred", terms = c("Temp_1 [all]", "sex"), colors = c('goldenrod3','dodgerblue4')) + theme_sjplot2()
sjPlot::plot_model(mod, type = "pred", terms = c("photoperiod_1 [all]", "sex"), colors = c('goldenrod3','dodgerblue4')) + theme_sjplot2()
sjPlot::plot_model(mod, type = "pred", terms = c("source", "sex"), colors = c('goldenrod3','dodgerblue4')) + theme_sjplot2()






