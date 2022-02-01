##############################################
## combined (male + female) phenoplots      ##
## for Peromyscus maniculatus             	##
## Bryan S. McLean    5 April 2021          ##
##############################################

library(mgcv)
library(ggplot2)
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
no.data <- which(is.na(phens$state) | is.na(phens$ecoregion) | is.na(phens$sex) == TRUE)
phens.gam <- phens[-no.data,]; rm(no.data)
phens.gam$sex <- as.factor(phens.gam$sex)

##############################################
## hierarchical GAMs
##############################################

# full HGAM with ecoregion:sex random effect
mf.gam <- gam(state ~ 
  (sex * source) + 
  s(jday, by = ecoregion:sex, k = 5, bs = "cc", m = 2) + 
  s(ecoregion, sex, bs = "re"),
  data = phens.gam, 
  knots = list(jday=c(0,364)),
  family=binomial(), 
  method='REML',
  drop.unused.levels=FALSE
  )

saveRDS(mf.gam, file="/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/hgam_breedingstate_malefemalecombined.rds")

# reduced HGAM lacking interactions btw parametric terms (reduced1)
mf.gam.reduced1 <- gam(state ~ 
  (sex + source) + 
  s(jday, by = ecoregion:sex, k = 5, bs = "cc", m = 2) + 
  s(ecoregion, sex, bs = "re"),
  data = phens.gam, 
  knots = list(jday=c(0,364)),
  family=binomial(), 
  method='REML',
  drop.unused.levels=FALSE
  )

# reduced HGAM with no random effect for sex (reduced2)
mf.gam.reduced2 <- gam(state ~ 
  (sex * source) + 
  s(jday, by = ecoregion:sex, k = 5, bs = "cc", m = 2) + 
  s(ecoregion, bs = "re"),
  data = phens.gam, 
  knots = list(jday=c(0,364)),
  family=binomial(), 
  method='REML',
  drop.unused.levels=FALSE
  )

# reduced HGAM lacking interactions btw parametric terms and with no random effect for sex (reduced3)
mf.gam.reduced3 <- gam(state ~ 
  (sex + source) + 
  s(jday, by = ecoregion:sex, k = 5, bs = "cc", m = 2) + 
  s(ecoregion, bs = "re"),
  data = phens.gam, 
  knots = list(jday=c(0,364)),
  family=binomial(), 
  method='REML',
  drop.unused.levels=FALSE
  )


##############################################
## HGAM comparisons
##############################################

mf.scores <- AIC(mf.gam, mf.gam.reduced1, mf.gam.reduced2, mf.gam.reduced3)
mf.calls <- c(formula(mf.gam), formula(mf.gam.reduced1), formula(mf.gam.reduced2), formula(mf.gam.reduced3))
delta.AIC <- mf.scores$AIC - min(mf.scores$AIC)
mf.scores <- cbind(mf.scores, delta.AIC, as.character(mf.calls))

# write.table(mf.scores, file="/users/mclean/Box/Projects/PEMA testes/models-final-bsm/hgam_breedingstate_malefemalecombined_comparison.txt", sep = ' ', row.names=T, quote=F)


##############################################
## phenoplots
##############################################

# working from full model
# read the Rdata file
mf.gam <- readRDS(file="/users/mclean/Box/Projects/PEMA_male-female_phenology/data/_final_w_clim_/mod.gam.MFcombined.rds")
# generate matrix of values from which to predict
pred.i <- with(phens.gam,
               expand.grid(
                jday = seq(0,364,length=365),
                sex = unique(phens.gam$sex),
                ecoregion = unique(phens.gam$ecoregion),
                source = "VERTNET" # vertnet and nacsm had higher (and not differet) avgs, so predict to this level
                )
              )
# HGAM predict
pred.i <- cbind(pred.i, predict(mf.gam, pred.i, se.fit = T, type = 'response'))

# building some lines to represent equinox/solstices
solst <- c(
  strptime(paste(c(20,3,2019),collapse="-"),"%d-%m-%Y")$yday,
  strptime(paste(c(21,6,2019),collapse="-"),"%d-%m-%Y")$yday,
  strptime(paste(c(23,9,2019),collapse="-"),"%d-%m-%Y")$yday
  )

# cols for sexes (if needed)
MF.cols <- c("M" = "dodgerblue4", "F" = "goldenrod3")

# PHENOPLOTS!
ggplot(data = phens.gam, aes(x = jday, y = state, group = ecoregion)) +
  geom_vline(
    alpha=0.15,
    colour='grey10',
    xintercept=solst) +
  geom_ribbon(
    aes(x = jday, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)),
    data=pred.i[which(pred.i$sex == "M"),],
    fill = 'dodgerblue4',
    alpha=0.6,
    inherit.aes=F,
    show.legend=F
  ) + 
  geom_line(
    aes(y=fit), 
    data=pred.i[which(pred.i$sex == "M"),], 
    #linetype='dashed',
    size=1,
    colour = "dodgerblue4"
  ) +
  geom_ribbon(
    aes(x = jday, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)),
    data=pred.i[which(pred.i$sex == "F"),],
    fill = 'goldenrod1',
    alpha=0.6,
    inherit.aes=F,
    show.legend=F
  ) + 
  geom_line(
    aes(y=fit), 
    data=pred.i[which(pred.i$sex == "F"),], 
    #linetype='dashed',
    size=1,
    colour = "goldenrod3"
  ) +
  geom_rug(
    data = subset(phens.gam, sex == "M"),
    aes(x = jday),
    colour = 'dodgerblue4',
    #mapping = aes(color = sex), 
    alpha = 0.1, 
    sides = c('t'), 
    length=unit(0.03, "npc")
  ) + 
  geom_rug(
    data = subset(phens.gam, sex == "F"),
    aes(x = jday),
    colour = 'goldenrod3',
    #mapping = aes(color = sex), 
    alpha = 0.1, 
    sides = c('b'), 
    length=unit(0.03, "npc")
  ) + 
  scale_color_manual(values = MF.cols) +
  facet_wrap(~ecoregion ,scales='fixed') + 
  scale_y_continuous(position = 'left', breaks = c(0,0.25,0.5,0.75,1), labels = c(0,0.25,0.5,0.75,1), expand = expansion(add = 0.05)) +
  theme(strip.text.x = element_text(size=0)) + 
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color='transparent'),
    panel.grid.minor.y = element_line(color='transparent'),
    panel.grid.major.x = element_line(color='transparent'),
    panel.grid.minor.x = element_line(color='transparent')
  )

