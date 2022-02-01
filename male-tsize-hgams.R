##############################################
## Phenoplots of male testes size           ##
## for Peromyscus maniculatus             	##
## Bryan S. McLean    5 April 2021          ##
##############################################

phen.m <- read.csv("./data/_final_w_clim_/male_data_ClimNA.csv")

# create new testes size-specific data frame
phen.m2 <- data.frame(
  source = c(phen.m$source),
  jday = c(phen.m$julianday),
  state = c(phen.m$t_state),
  tsize = c(phen.m$testes_max_length),
  ecoregion = as.factor(c(phen.m$ecoregion2)),
  nffd = c(phen.m$NFFD),
  massing = phen.m$massing
)

# harmonizing for future female comparisons...
phen.m2$source[which(phen.m2$source == 'nacsm')] <- 'NACSM'
phen.m2$source[which(phen.m2$source == 'neon')] <- 'NEON'
phen.m2$source[which(phen.m2$source == 'VN' | phen.m2$source == 'vertnet')] <- 'VERTNET'
phen.m2$source[which(phen.m2$source == 'SPEC')] <- 'VERTNET'
no.data <- which(is.na(phen.m2$tsize) == TRUE)
phen.m2.gam <- phen.m2[-no.data,]; rm(no.data)

# assign testes size observations in the filtered data frame to NFFD quantiles
# using 5 quantiles
nffd.bin5 <- as.factor(.bincode(phen.m2.gam$nffd, breaks = quantile(phen.m2.gam$nffd,seq(0,1,(1/5)),na.rm=T), include.lowest = TRUE))
phen.m2.gam <- data.frame(phen.m2.gam,nffd.bin5)	

# re-cleaning
no.data <- which(is.na(phen.m2.gam$nffd.bin5) == TRUE)
phen.m2.gam <- phen.m2.gam[-no.data,]; rm(no.data)


##############################################
## hierarchical GAM on testes size
##############################################
# random slope and intercept for NFFD quantile
tsize.gam <- gam(tsize ~ 
  s(jday, by = nffd.bin5, k = 5, bs = "cc", m = 2) + 
  s(nffd.bin5, bs = "re"),
  data = phen.m2.gam, 
  knots = list(jday=c(0,364)), 
  method = 'REML',
  drop.unused.levels = FALSE
  )
saveRDS(tsize.gam, file="/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/hgam_male_tsize.rds")


##############################################
## phenoplots
##############################################

# read the Rdata file
tsize.gam <- readRDS(file="/users/mclean/Box/Projects/PEMA_male-female_phenology/models-final-bsm/hgam_male_tsize.rds")
# generate matrix of values from which to predict
pred.t <- with(phen.m2.gam,
               expand.grid(
                jday = seq(0,364,length=365),
                nffd.bin5 = unique(phen.m2.gam$nffd.bin5)
                )
               )
# HGAM predict
pred.t <- cbind(pred.t, predict(tsize.gam, pred.t, se.fit = T, type = 'response'))

# TESTES!
ggplot(data = phen.m2.gam, aes(x = jday, y = tsize, group = nffd.bin5)) +
  geom_vline(
    alpha=0.15,
    colour='grey10',
    xintercept=solst) +
  geom_ribbon(
    aes(x = jday, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)),
    data=pred.t,
    #fill = 'dodgerblue4',
    alpha=0.6,
    inherit.aes=F,
    show.legend=F
  ) + 
  geom_line(
    aes(y=fit), 
    data=pred.t, 
    #linetype='dashed',
    size=1,
    #colour = "dodgerblue4"
  ) +
  geom_rug(
    data = phen.m2.gam,
    aes(x = jday),
    #colour = 'dodgerblue4',
    #mapping = aes(color = sex), 
    alpha = 0.1, 
    sides = c('b'), 
    length=unit(0.05, "npc")
  ) +
  #scale_color_manual(values = MF.cols) +
  facet_wrap(~nffd.bin5, scales='fixed', nrow=1) +
  scale_y_continuous(
    position = 'left', 
    limits = c(1.5,11.5),
    breaks = c(1,3,5,7,9,11), 
    labels = c(1,3,5,7,9,11)
    # , expand = expansion(add = 0.05)
  ) +
  theme(strip.text.x = element_text(size=0)) + 
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color='transparent'),
    panel.grid.minor.y = element_line(color='transparent'),
    panel.grid.major.x = element_line(color='transparent'),
    panel.grid.minor.x = element_line(color='transparent')
  )


##############################################
## additional analyses of testes size phenology
## in relation to female PEMA breeding phenology
##############################################

# regressing testes size on female breeding season length
ecoregions <- unique(phens.gam$ecoregion)
tsize_fbreed <- data.frame(matrix(NA, nrow = length(ecoregions), ncol = 4, dimnames = list(ecoregions, c("female_breeding_season_length","tsize", "tsize_scaled","tnobs"))))

for(i in rownames(tsize_fbreed)){
  
  tmp <- pred.i[which(pred.i$sex == "F" & pred.i$ecoregion == i),]
  tsize_fbreed[i,1] <- length(which(tmp$fit > (1/3))) # 25% of females breeding
  tmp <- phen.m2[which(phen.m2$ecoregion == i),]
  tsize_fbreed[i,2] <- max(tmp$tsize, na.rm = T)
  tsize_fbreed[i,3] <- max(tmp$tsize/(tmp$massing^(1/3)), na.rm = T) # cube root of body mass as scalar
  tsize_fbreed[i,4] <- length(na.omit(tmp$tsize/(tmp$massing^(1/3))))
  
}
plot(tsize_fbreed[which(tsize_fbreed$tnobs > 10), c(1,3)]), ylim = c(4,8))
summary(lm(tsize_scaled ~ female_breeding_season_length, tsize_fbreed))


# regressions of testes size on nffd

tmp <- phen.m2.gam
tfxn <- max
tsize_nffd5 <- data.frame(
  nffd = aggregate(nffd ~ nffd.bin5, FUN = mean, data = tmp)[,2],
  tsize = aggregate(tsize ~ nffd.bin5, FUN = tfxn, data = tmp)[,2],
  tsize_scaled = aggregate(tmp$tsize/(tmp$massing^(1/3)) ~ tmp$nffd.bin5, FUN = tfxn)[,2],
  tnobs = as.vector(table(tmp$nffd.bin5))
)
tsize_nffd8 <- data.frame(
  nffd = aggregate(nffd ~ nffd.bin8, FUN = mean, data = tmp)[,2],
  tsize = aggregate(tsize ~ nffd.bin8, FUN = tfxn, data = tmp)[,2],
  tsize_scaled = aggregate(tmp$tsize/(tmp$massing^(1/3)) ~ tmp$nffd.bin8, FUN = tfxn)[,2],
  tnobs = as.vector(table(tmp$nffd.bin8))
)
tsize_nffd11 <- data.frame(
  nffd = aggregate(nffd ~ nffd.bin11, FUN = mean, data = tmp)[,2],
  tsize = aggregate(tsize ~ nffd.bin11, FUN = tfxn, data = tmp)[,2],
  tsize_scaled = aggregate(tmp$tsize/(tmp$massing^(1/3)) ~ tmp$nffd.bin11, FUN = tfxn)[,2],
  tnobs = as.vector(table(tmp$nffd.bin11))
)
tsize_ecoregion <- data.frame(
  row.names = aggregate(nffd ~ ecoregion, FUN = mean, data = tmp)[,1],
  nffd = aggregate(nffd ~ ecoregion, FUN = mean, data = tmp)[,2],
  tsize = aggregate(tsize ~ ecoregion, FUN = tfxn, data = tmp)[,2],
  tsize_scaled = aggregate(tmp$tsize/(tmp$massing^(1/3)) ~ tmp$ecoregion, FUN = tfxn)[,2],
  tnobs = as.vector(table(tmp$ecoregion))
)

# compare and plot

summary(lm(tsize_scaled ~ female_breeding_season_length, tsize_fbreed))
summary(lm(tsize_scaled ~ nffd, tsize_ecoregion))

g <- ggplot(
  data = tsize_ecoregion, 
  aes(x = nffd, y = tsize_scaled)
) + 
  geom_smooth(	
    method = "lm", 
    formula = y~x,
    se = T,
    colour='grey40',
    alpha=0.3
  ) +
  geom_point() + 
  theme_bw()

g2 <- ggplot(
  data = tsize_fbreed, 
  aes(x = female_breeding_season_length, y = tsize_scaled)
) + 
  geom_smooth(	
    method = "lm", 
    formula = y~x,
    se = T,
    colour='grey40',
    alpha=0.3
  ) +
  geom_point() + 
  theme_bw()


