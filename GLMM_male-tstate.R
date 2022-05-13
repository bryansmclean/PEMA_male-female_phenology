#######################################################################
### MALE TESTES STATE GLMMs (generalized linear mixed models)		###
### used for testing environmental drivers of breeding			 	###
### Bryan S. McLean													###
### 3 September 2021, last modified 12 May 2022						###
#######################################################################

library(standardize)
library(lme4)

############
## DATA ##
############

phen.m <- read.csv("PEMA_male_data.csv")
phen.m$source <- as.factor(phen.m$source)

############
## MODELS ##
############

# full GLMM with all interactions and poly terms

mod.full <- standardize::standardize(
	testes_state ~ 
	(poly(temp, 2, raw = F) +
	poly(prec, 2, raw = F) +
	poly(photoperiod, 2, raw = F))^2 +
	massing +
	(1 | source) +
	(1 | ecoregion) + 
	(0 + poly(temp, 2, raw = F) | ecoregion) +
	(0 + poly(prec, 2, raw = F) | ecoregion) +
	(0 + poly(photoperiod, 2, raw = F) | ecoregion),
	data = phen.m,
	family = binomial(link='logit')	
	)
m1 <- glmer(
	mod.full$formula,
	family = binomial(),
	data = mod.full$data,
	control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))
	)
	
# reduced GLMM, lacking interactions among environmental terms

mod.noint <- standardize::standardize(
	testes_state ~ 
	poly(temp, 2, raw = F) +
	poly(prec, 2, raw = F) +
	poly(photoperiod, 2, raw = F) +
	massing +
	(1 | source) +
	(1 | ecoregion) + 
	(0 + poly(temp, 2, raw = F) | ecoregion) +
	(0 + poly(prec, 2, raw = F) | ecoregion) +
	(0 + poly(photoperiod, 2, raw = F) | ecoregion),
	data = phen.m,
	family = binomial(link='logit')	
	)
m2 <- glmer(
	mod.noint$formula,
	family = binomial(),
	data = mod.noint$data,
	control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))
	)
	
# reduced GLMM, lacking second-order environmental terms

mod.nopoly <- standardize::standardize(
	testes_state ~ 
	(temp +
	prec +
	photoperiod)^2 +
	massing +
	(1 | source) +
	(1 | ecoregion) + 
	(0 + temp | ecoregion) +
	(0 + prec | ecoregion) +
	(0 + photoperiod | ecoregion),
	data = phen.m,
	family = binomial(link='logit')	
	)
m3 <- glmer(
	mod.nopoly$formula,
	family = binomial(),
	data = mod.nopoly$data,
	control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))
	)
	
# reduced GLMM, lacking second-order environmental terms and interactions among environmental terms

mod.noint.nopoly <- standardize::standardize(
	testes_state ~ 
	temp +
	prec +
	photoperiod +
	massing +
	(1 | source) +
	(1 | ecoregion) + 
	(0 + temp | ecoregion) +
	(0 + prec | ecoregion) +
	(0 + photoperiod | ecoregion),
	data = phen.m,
	family = binomial(link='logit')	
	)
m4 <- glmer(
	mod.noint.nopoly$formula,
	family = binomial(),
	data = mod.noint.nopoly$data,
	control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))
	)
