#######################################################################
### MALE + FEMALE HGAMs	(heirarchical generalized additive models)	###
### used for reconstructing + visualizing sex-specific phenologies 	###
### Bryan S. McLean													###
### 3 September 2021, last modified 12 May 2022						###
#######################################################################

library(mgcv)
library(ggplot2)

############
## DATA ##
############

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

# full HGAM 

mf.gam <- gam(
	state ~ 
	(sex * source) + 
	s(jday, by = ecoregion:sex, k = 5, bs = "cc", m = 2) + 
	s(ecoregion, sex, bs = "re"),
	data = phens, 
	knots = list(jday = c(0, 364)),
	family = binomial(), 
	method = 'REML',
	drop.unused.levels = FALSE 
	)

# "reduced" HGAM lacking interactions between parametric terms 

mf.gam.reduced1 <- gam(
	state ~ 
	(sex + source) + 
	s(jday, by = ecoregion:sex, k = 5, bs = "cc", m = 2) + 
    	s(ecoregion, sex, bs = "re"),
    	data = phens, 
	knots = list(jday = c(0, 364)),
	family = binomial(), 
	method = 'REML',
	drop.unused.levels = FALSE 
	)

# "reduced" HGAM lacking random effect of sex 

mf.gam.reduced2 <- gam(
	state ~ 
	(sex * source) + 
	s(jday, by = ecoregion:sex, k = 5, bs = "cc", m = 2) + 
    	s(ecoregion, bs = "re"),
    	data = phens, 
	knots = list(jday = c(0, 364)),
	family=binomial(), 
	method = 'REML',
	drop.unused.levels = FALSE 
	)

# "reduced" HGAM lacking interactions between parametric terms and random effect of sex 

mf.gam.reduced3 <- gam(
	state ~ 
	(sex + source) + 
	s(jday, by = ecoregion:sex, k = 5, bs = "cc", m = 2) + 
    	s(ecoregion, bs = "re"),
    	data = phens, 
	knots = list(jday = c(0, 364)),
	family = binomial(), 
	method = 'REML',
	drop.unused.levels = FALSE 
	)



############
## PLOT   ##
############

# prediction (from full HGAM only)
pred.i <- with(
	phens, 
	expand.grid(
		jday = seq(0, 364, length = 365),
		sex = unique(phens$sex),
		ecoregion = unique(phens$ecoregion),
		source = "VERTNET" # vertnet and nacsm had higher (and not differet) avgs -> predict to this level
		)
	)
pred.i <- cbind(pred.i, predict(mf.gam, pred.i, se.fit = T, type = 'response'))

# some lines for equinox/solstice
solst <- c(
	strptime(paste(c(20,3,2019),collapse="-"),"%d-%m-%Y")$yday,
	strptime(paste(c(21,6,2019),collapse="-"),"%d-%m-%Y")$yday,
	strptime(paste(c(23,9,2019),collapse="-"),"%d-%m-%Y")$yday
	)

# cols for sexes
MF.cols <- c("M" = "dodgerblue4", "F" = "goldenrod3")

# plot
ggplot(data = phens, aes(x = jday, y = state, group = ecoregion)) +
	geom_vline(
		alpha = 0.15,
		colour = 'grey10',
		xintercept = solst) +
	geom_ribbon(
		aes(x = jday, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)),
		data = pred.i[which(pred.i$sex == "M"),],
		fill = 'dodgerblue4',
		alpha = 0.6,
		inherit.aes = F,
		show.legend = F 
		) + 
	geom_line(
		aes(y = fit), 
		data = pred.i[which(pred.i$sex == "M"),], 
		#linetype='dashed',
		size = 1,
		colour = "dodgerblue4"
		) +
	geom_ribbon(
		aes(x = jday, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)),
		data = pred.i[which(pred.i$sex == "F"),],
		fill = 'goldenrod1',
		alpha = 0.6,
		inherit.aes = F,
		show.legend = F
		) + 
	geom_line(
		aes(y = fit), 
		data = pred.i[which(pred.i$sex == "F"),], 
		#linetype='dashed',
		size = 1,
		colour = "goldenrod3"
		) +
	geom_rug(
		data = subset(phens, sex == "M"),
		aes(x = jday),
		colour = 'dodgerblue4',
		#mapping = aes(color = sex), 
		alpha = 0.1, 
		sides = c('t'), 
		length = unit(0.03, "npc")
		) + 
	geom_rug(
		data = subset(phens, sex == "F"),
		aes(x = jday),
		colour = 'goldenrod3',
		#mapping = aes(color = sex), 
		alpha = 0.1, 
		sides = c('b'), 
		length = unit(0.03, "npc")
		) + 
	scale_color_manual(values = MF.cols) +
	facet_wrap(~ecoregion, scales = 'fixed') + 
	scale_y_continuous(
		position = 'left', 
		breaks = c(0,0.25,0.5,0.75,1), 
		labels = c(0,0.25,0.5,0.75,1), 
		expand = expansion(add = 0.05)
		) +
	theme(strip.text.x = element_text(size=0)) + 
	theme_bw() +
	theme(
		panel.grid.major.y = element_line(color='transparent'),
		panel.grid.minor.y = element_line(color='transparent'),
		panel.grid.major.x = element_line(color='transparent'),
		panel.grid.minor.x = element_line(color='transparent')
		)


