# vertebral evolution vs elongation index

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

require(ape)
require(scales)

alldat <- read.csv('./data/alldat.csv')

snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']


dat <- alldat[, c('numberPresacralVert', 'elongationIndex')]
rownames(dat) <- alldat$treename
dat <- dat[complete.cases(dat),]


pdf('~/Downloads/figS12.pdf', width = 7, height = 7)

plot.new()
plot.window(xlim = range(dat[,1]), ylim = range(dat[,2]))
axis(1, lwd = 0, lwd.ticks = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0))
axis(2, lwd = 0, lwd.ticks = 1, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0))
box(which = "plot", bty = "l")
mtext('presacral vertebral count', side = 1, line = 2.5, cex = 1)
mtext('elongation index', side = 2, line = 2.5, cex = 1)

points(dat[rownames(dat) %in% snakes,], bg = adjustcolor('coral2', alpha.f = 0.75), lwd = 0.5, pch = 21)
points(dat[rownames(dat) %in% lizards,], bg = adjustcolor('slateblue1', alpha.f = 0.75), lwd = 0.5, pch = 21)

legend('bottomright', legend = c('lizards', 'snakes'), pt.bg = adjustcolor(c('slateblue1', 'coral2'), alpha.f = 0.75), bty = 'n', pch = 21, pt.cex = 1.5)

dev.off()







