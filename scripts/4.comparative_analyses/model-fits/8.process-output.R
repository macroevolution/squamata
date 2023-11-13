basedir <- '~/Dropbox/Oz_Crown_Ages/dataArchive'

setwd(basedir)

setwd('./scripts/4.comparative_analyses/model-fits')

library(ape)

skull <- read.csv("fit_skull.csv")
mass  <- read.csv("fit_mass.csv")
vert  <- read.csv("fit_vert.csv")
elong <- read.csv("fit_elong.csv")
diet  <- read.csv("fit_diet.csv")
climate <- read.csv("fit_climate.csv")

source("1.load-data.R")


#------------------# 
# Summary model fits:

plotSetup <- function(x, y, trait){
	
	bb <- 1
	if (max(x) > 5){
		bb <- 2
	}

	plot.new()
	plot.window(xlim=c(-0.5, max(x)+0.5), ylim=c(0, max(y)))
	
	axis(1, at=c(-1, seq(0, max(x)+2, by=bb)))
	axis(2, at=c(round(-0.5*max(y)), round(seq(0, max(y), length.out=5))), las=1)
	
	mtext("Rate partitions", side=1, line=2.7, cex=0.8)
	mtext("AIC (relative to best)", side=2, line=3.5, cex=0.9)
	mtext(trait, side=3, line=0, font=2, cex=1)
}

base <- "white"
col_shift <- "cornflowerblue"
col_cens  <- "darkorange1"
pcex <- 1.3


quartz.options(height=6, width=10)
# pdf('~/Downloads/figS16.pdf', height=6, width=10)
par(mfrow=c(2,3), mar=c(4,4,3,1), oma=c(2,1,1,1))

letterX <- -0.1
letterY <- -1

#---------------------#

DD <- mass

plotSetup(DD$nshifts, DD$aic, trait="Mass")

#jx <- jitter(DD$nshifts, factor=1) 
colvec <- rep(col_shift, nrow(DD))
colvec[DD$type == "base"] <- "white"
colvec[DD$type == "censored"] <- col_cens
points(x = DD$nshifts, y = DD$aic, cex=pcex, pch=21, bg=colvec)
mtext('A', adj = letterX, padj = letterY, font = 2, xpd = NA, cex = 1.2)

legend(x=2.5, y=8200,  fill = c(base, col_shift, col_cens), legend = c("Base", "Rate shift only", "Fully decoupled"))


#---------------------#

DD <- elong

plotSetup(DD$nshifts, DD$aic, "Elongation")

#jx <- jitter(DD$nshifts, factor=1) 
colvec <- rep(col_shift, nrow(DD))
colvec[DD$type == "base"] <- "white"
colvec[DD$type == "censored"] <- col_cens
points(x = DD$nshifts, y = DD$aic, cex=pcex, pch=21, bg=colvec)
mtext('B', adj = letterX, padj = letterY, font = 2, xpd = NA, cex = 1.2)


#---------------------#

DD <- vert

plotSetup(DD$nshifts, DD$aic, "Presacral vertebrae")

#jx <- jitter(DD$nshifts, factor=1) 
colvec <- rep(col_shift, nrow(DD))
colvec[DD$type == "base"] <- "white"
colvec[DD$type == "censored"] <- col_cens
points(x = DD$nshifts, y = DD$aic, cex=pcex, pch=21, bg=colvec)
mtext('C', adj = letterX, padj = letterY, font = 2, xpd = NA, cex = 1.2)

#------------------------# 

DD <- skull
plotSetup(DD$nshifts, DD$aic, "Multivariate skull")

jx <- jitter(DD$nshifts, factor=1)
 
colvec <- rep(col_shift, nrow(DD))
colvec[DD$type == "base"] <- "white"
colvec[DD$type == "censored"] <- col_cens
points(x = jx, y = DD$aic, cex=1.3, pch=21, bg=colvec)
mtext('D', adj = letterX, padj = letterY, font = 2, xpd = NA, cex = 1.2)


#---------------------#

DD <- climate

plotSetup(DD$nshifts, DD$aic, "Multivariate climate")

#jx <- jitter(DD$nshifts, factor=1) 
colvec <- rep(col_shift, nrow(DD))
colvec[DD$type == "base"] <- "white"
colvec[DD$type == "censored"] <- col_cens
points(x = DD$nshifts, y = DD$aic, cex=pcex, pch=21, bg=colvec)
mtext('E', adj = letterX, padj = letterY, font = 2, xpd = NA, cex = 1.2)

#---------------------#

DD <- diet

plotSetup(DD$nshifts, DD$aic, "Multivariate diet")

#jx <- jitter(DD$nshifts, factor=1) 
colvec <- rep(col_shift, nrow(DD))
colvec[DD$type == "base"] <- "white"
colvec[DD$type == "censored"] <- col_cens
points(x = DD$nshifts, y = DD$aic, cex=pcex, pch=21, bg=colvec)
mtext('F', adj = letterX, padj = letterY, font = 2, xpd = NA, cex = 1.2)

dev.off()
#--------------------------------------# 






