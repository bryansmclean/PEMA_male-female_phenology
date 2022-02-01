##############################################
## male phenology models (and comparison)   ##
## for Peromyscus maniculatus             	##
## Bryan S. McLean    5 April 2021          ##
##############################################

phen.m <- read.csv("/users/mclean/Box/Projects/PEMA_male-female_phenology/data/_final_w_clim_/male_data_ClimNA.csv")
phen.m$source <- as.factor(phen.m$source)

##############################################
## testes state mixed models
##############################################

library(standardize)
library(lme4)

# (very longhand) syntax, to allow easier effect plotting later
# first specify polys and add to data frame
Temp_poly <- poly(phen.m$Temp_BM, 2, raw = F, simple = T); colnames(Temp_poly) <- c("Temp_1","Temp_2")
Prec_poly <- poly(phen.m$Prec_BM, 2, raw = F, simple = T); colnames(Prec_poly) <- c("Prec_1","Prec_2")
photoperiod_poly <- poly(phen.m$photoperiod, 2, raw = F, simple = T); colnames(photoperiod_poly) <- c("photoperiod_1","photoperiod_2")
phen.m <- cbind(phen.m, Temp_poly, Prec_poly, photoperiod_poly)

## full model, all 2-way interactions, ridiculous longhand version
## set up a ::standardize:: object
mod.quad.full <- standardize::standardize(t_state ~ 
  Temp_1 + Temp_2 + 
  Prec_1 + Prec_2 + 
  photoperiod_1 + photoperiod_2 + 
  (Temp_1 + Prec_1 + photoperiod_1)^2 + 
  (Temp_1 + Prec_2 + photoperiod_1)^2 + 
  (Temp_1 + Prec_1 + photoperiod_2)^2 + 
  (Temp_2 + Prec_1 + photoperiod_1)^2 + 
  (Temp_2 + Prec_2 + photoperiod_1)^2 + 
  (Temp_2 + Prec_2 + photoperiod_2)^2 + 
  massing + 
  source + 
  (1 | ecoregion2) + 
  (0 + (Temp_1 + Temp_2) | ecoregion2) + 
  (0 + (Prec_1 + Prec_2) | ecoregion2) + 
  (0 + (photoperiod_1 + photoperiod_2) | ecoregion2),
  data = phen.m,
  family = binomial(link='logit')	
	)

saveRDS(
  m1 <- glmer(
    mod.quad.full$formula,
    family = binomial(),
    data = mod.quad.full$data,
    control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))
    ),
  file="/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/glmm_male_tstate_full2_bobyqaopt.rds"
  )

write.table(
	broom::tidy(full2),
	file = "/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/glmm_male_tstate_full2_bobyqaopt.csv",
	sep = ",",
	row.names = F,
	quote = F
	)

## GLMM with no interactions btw fixed terms

mod.quad.noint <- standardize::standardize(t_state ~ 
  poly(Temp_BM, 2, raw = F) +
  poly(Prec_BM, 2, raw = F) +
  poly(photoperiod, 2, raw = F) +
  massing +
  (1 | source) +
  (1 | ecoregion2) + 
  (0 + poly(Temp_BM, 2, raw = F) | ecoregion2) +
  (0 + poly(Prec_BM, 2, raw = F) | ecoregion2) +
  (0 + poly(photoperiod, 2, raw = F) | ecoregion2),
  data = phen.m,
  family = binomial(link='logit')
  )

saveRDS(
  glmer(
    mod.quad.noint$formula,
    family = binomial(),
    data = mod.quad.noint$data,
    control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))
  ),
  file="/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final/glmm_male_tstate_nointer_bobyqaopt.rds"
)

## GLMM with no polynomials on fixed terms

mod.quad.nopoly <- standardize::standardize(t_state ~ 
  (Temp_BM +
  Prec_BM +
  photoperiod)^2 +
  massing +
  (1 | source) +
  (1 | ecoregion2) + 
  (0 + Temp_BM | ecoregion2) +
  (0 + Prec_BM | ecoregion2) +
  (0 + photoperiod | ecoregion2),
  data = phen,
  family = binomial(link='logit')
  )

saveRDS(
  glmer(
    mod.quad.nopoly$formula,
    family = binomial(),
    data = mod.quad.nopoly$data,
    control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))
    ),
  file="/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/glmm_male_tstate_nopoly_bobyqaopt.rds"
  )


##############################################
## GLMM comparisons
##############################################

aics <- AIC(
  readRDS("/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/glmm_male_tstate_full2_bobyqaopt.rds"),
	readRDS("/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/glmm_male_tstate_nointer_bobyqaopt.rds"),
	readRDS("/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/glmm_male_tstate_nopoly_bobyqaopt.rds"),
	readRDS("/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/glmm_male_tstate_nopolynoint_bobyqaopt.rds")
	)
aics <- cbind(aics, aics$AIC - min(aics$AIC))
write.table(aics, "/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/glmm_male_tstate_aic.txt", quote = F)


