#######################################################################
### MALE HGAMs	(heirarchical generalized additive models)			###
### used for reconstructing phenology of testis length			 	###
### Bryan S. McLean													###
### 5 April 2021, last modified 12 May 2022							###
#######################################################################

library(mgcv)
library(ggplot2)

############
## DATA ##
############

phen.m <- read.csv("PEMA_male_data.csv")
no.data <- which(is.na(phen.m$testes_max_length) == TRUE)
phen.m.gam <- phen.m[-no.data,]; rm(no.data)

# assign testes size observations in the filtered data frame to NFFD quantiles
# using 5 quantiles
nffd.bin5 <- as.factor(.bincode(phen.m.gam$nffd, breaks = quantile(phen.m.gam$nffd, seq(0,1,(1/5)), na.rm=T), include.lowest = TRUE))
phen.m.gam <- data.frame(phen.m.gam, nffd.bin5)	

# re-cleaning
no.data <- which(is.na(phen.m.gam$nffd.bin5) == TRUE)
phen.m.gam <- phen.m.gam[-no.data,]; rm(no.data)


############
## MODELS ##
############

# HGAM with random slope and intercept on NFFD quantile
tsize.gam <- gam(testes_max_length ~ 
  s(julianday, by = nffd.bin5, k = 5, bs = "cc", m = 2) + 
  s(nffd.bin5, bs = "re"),
  data = phen.m.gam, 
  knots = list(jday=c(0,364)), 
  method = 'REML',
  drop.unused.levels = FALSE
  )


############
## PLOT   ##
############

# generate matrix of values from which to predict
pred.t <- with(phen.m.gam,
			expand.grid(
                julianday = seq(0,364,length=365),
                nffd.bin5 = unique(phen.m.gam$nffd.bin5)
                )
              )


# HGAM predict
pred.t <- cbind(pred.t, predict(tsize.gam, pred.t, se.fit = T, type = 'response'))

# Phenoplots of testes length 
ggplot(data = phen.m.gam, aes(x = julianday, y = testes_max_length, group = nffd.bin5)) +
  geom_vline(
    alpha=0.15,
    colour='grey10',
    xintercept=solst) +
  geom_ribbon(
    aes(x = julianday, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)),
    data = pred.t,
    #fill = 'dodgerblue4',
    alpha = 0.6,
    inherit.aes = F,
    show.legend = F
  	) + 
  geom_line(
    aes(y = fit), 
    data = pred.t, 
    #linetype='dashed',
    size = 1,
    #colour = "dodgerblue4"
  	) +
  geom_rug(
    data = phen.m.gam,
    aes(x = julianday),
    #colour = 'dodgerblue4',
    #mapping = aes(color = sex), 
    alpha = 0.1, 
    sides = c('b'), 
    length = unit(0.05, "npc")
  	) +
  #scale_color_manual(values = MF.cols) +
  facet_wrap(~nffd.bin5, scales = 'fixed', nrow=1) +
  scale_y_continuous(
    position = 'left', 
    limits = c(1.5,11.5),
    breaks = c(1,3,5,7,9,11), 
    labels = c(1,3,5,7,9,11)
  	) +
  theme(strip.text.x = element_text(size=0)) + 
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color='transparent'),
    panel.grid.minor.y = element_line(color='transparent'),
    panel.grid.major.x = element_line(color='transparent'),
    panel.grid.minor.x = element_line(color='transparent')
  )



