#######################################################################
### MALE - FEMALE GLMMs (generalized linear mixed models)			###
### used for testing sex-specific environmental drivers of breeding	###
### Bryan S. McLean													###
### 3 September 2021, last modified 12 May 2022						###
#######################################################################

library(standardize)
library(lme4)

############
## DATA ##
############

# read the data
phen.m <- read.csv("PEMA_male_data.csv")
phen.f <- read.csv("PEMA_phenology_final_dataframe.csv") #also at: https://datadryad.org/stash/dataset/doi:10.5061/dryad.zgmsbcc8w

# combine the data sets
phens <- data.frame(
		source = c(phen.m$source, phen.f$source),
		sex = c(rep('M', nrow(phen.m)), rep("F", nrow(phen.f))),
		jday = c(phen.m$julianday, phen.f$julianday.corrected),
		state = c(phen.m$testes_state, phen.f$state),
		ecoregion = as.factor(c(phen.m$ecoregion, phen.f$ecoregion)),
		photoperiod = c(phen.m$photoperiod, phen.f$photoperiod),
		Temp = c(phen.m$temp, phen.f$Tave),
		Prec = c(phen.m$prec, phen.f$PPT)
		)
		
# factorize some fields, standardize other fields
phens$source[which(phens$source == 'nacsm')] <- 'NACSM'
phens$source[which(phens$source == 'neon')] <- 'NEON'
phens$source[which(phens$source == 'VN' | phens$source == 'vertnet')] <- 'VERTNET'
phens$source[which(phens$source == 'SPEC')] <- 'VERTNET'
phens$source <- as.factor(phens$source)
phens$sex <- as.factor(phens$sex)

############
## MODELS ##
############

mod.full.mf <- standardize::standardize(
	state ~ 
	(poly(Temp, 2, raw = F) + poly(Prec, 2, raw = F) + poly(photoperiod, 2, raw = F) + sex)^2 + 
	(sex * source) +
	(1 | ecoregion) + 
	(0 + poly(Temp, 2, raw = F) | ecoregion) +
	(0 + poly(Prec, 2, raw = F) | ecoregion) +
	(0 + poly(photoperiod, 2, raw = F) | ecoregion),
	data = phens,
	family = binomial(link='logit')	
	)
m1 <- glmer(
	mod.full.mf$formula,
	family = binomial(),
	data = mod.full.mf$data,
	control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))
	)
	
	
	
	
	
	
	
	
	
	