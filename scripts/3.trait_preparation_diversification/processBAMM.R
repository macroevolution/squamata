# BAMM analyses
## check convergence
## compare outputs of 3 independent runs
## export for analysis

# expectedNumberOfShifts = 10

library(BAMMtools)

setwd('~/Dropbox/Oz_Crown_Ages')

tree <- read.tree('./final-trees/best_ultrametric_fulltree_ddBD_revision.tre')
master <- read.csv('./taxon-attributes/squamatesCoreTaxonomy_treeTaxa.csv')

round(setBAMMpriors(tree, total.taxa = nrow(master), outfile = NULL), 3)

devtools::source_url('https://gist.githubusercontent.com/ptitle/4c98cdcfb4ff507b2aec4ba8b9a914ef/raw/3b24b5ffe99f5892321bbe0e82e401f6089e5af3/bammCheck.R')

setwd('~/Dropbox/Oz_Crown_Ages/diversification/bamm')

bammCheck(expectedNumberOfShifts = 10, burnin = 0.1, mcmcFile = './run1/mcmc_out.txt')

bammCheck(expectedNumberOfShifts = 10, burnin = 0.1, mcmcFile = './run2/mcmc_out.txt')

bammCheck(expectedNumberOfShifts = 10, burnin = 0.1, mcmcFile = './run3/mcmc_out.txt')


# read in BAMM analyses

ed1 <- getEventData(tree, eventdata = './run1/event_data.txt', nsamples = 1000, burnin = 0.1)

ed2 <- getEventData(tree, eventdata = './run2/event_data.txt', nsamples = 1000, burnin = 0.1)

ed3 <- getEventData(tree, eventdata = './run3/event_data.txt', nsamples = 1000, burnin = 0.1)

tipRates1 <- getTipRates(ed1)$lambda.avg
tipRates2 <- getTipRates(ed2)$lambda.avg
tipRates3 <- getTipRates(ed3)$lambda.avg

par(mfrow = c(1,3))
plot(tipRates1, tipRates2); abline(0, 1)
plot(tipRates1, tipRates3); abline(0, 1)
plot(tipRates2, tipRates3); abline(0, 1)

cor(tipRates1, tipRates2)
cor(tipRates1, tipRates3)
cor(tipRates2, tipRates3)

# We will use run2 since it has the best convergence properties (best ESS).

# write tip rates to file
write.csv(as.matrix(tipRates2, rownames.force = T), file = 'bammTipRates.csv', row.names = TRUE)


