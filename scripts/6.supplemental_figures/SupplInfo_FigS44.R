library(bm)
library(LSD)

setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")


# Randomly sample a clade from a larger phylogeny such that all clade
# sizes larger than min.size have a uniform probability of being selected
#
# @param phy The simulation phylogeny.
# @param min.size The minimum clade size to consider.
# @param replace Sample with replacement?
# @param include.root Should the root clade be included in sampling?
rsubtree = function(phy, min.size=0, replace=FALSE, include.root=FALSE)
{    
    ntip = phylo::Ntip(phy)
    nnode = phylo::Nnode(phy)
    ancestors = phylo::root(phy):nnode
    if (!include.root)
        ancestors = ancestors[-1]
    ndesc = sapply(ancestors, function(n) length(phylo::tips(n, phy)))
    keep = which(ndesc > min.size)
    ndesc = ndesc[keep]
    ancestors = ancestors[keep]
    w = 1 / table(ndesc)
    sample_weights = unname(w[match(ndesc, as.integer(names(w)))])
    function(n)
    {
        stopifnot(n >= 1)
        subtrees = sample(ancestors, size=n, replace=replace, prob=sample_weights)
        return (subtrees)
    }
}

# @param phy The simulation phylogeny.
# @param d The number of traits.
# @param nsim The number of simulations.
# @param es The expected number of shifts in each simulation.
simulate.bm = function(phy, d, nsim, es, minSize = floor(0.1*phylo::Ntip(phy)))
{
    ret = array(dim=c(phylo::Ntip(phy), 3, nsim))

    for (sim in seq_len(nsim))
    {
        sample_shifts = rsubtree(phy, min.size = minSize)

        tmax = max(phylo::brlens(phy))
        node.state = matrix(0, phylo::Nnode(phy), d)

#       n_shift = rpois(1, es)
        n_shift = es
        rate = vector("list", n_shift+1)
        branch.regime = numeric(phylo::Nnode(phy))
        branch.rate = numeric(phylo::Nnode(phy))

        preorder = vector("list", n_shift+1)

        shift_nodes = integer(0)

        preorder[[1]] = phylo::descendants(phylo::root(phy), phy)

        if (n_shift)
        {
            shift_nodes = sample_shifts(n_shift)
            shift_nodes = shift_nodes[match(preorder[[1]], shift_nodes, nomatch=0)]
            for (s in seq_along(shift_nodes))
            {
                preorder[[s+1]] = phylo::descendants(shift_nodes[s], phy)
                preorder[[1]] = setdiff(preorder[[1]], preorder[[s+1]])
            }
        }

        node.state[phylo::root(phy), ] = rnorm(d)

        for (i in 1:(n_shift+1))
        {
            if (i > 1)
                tmax = max(phylo::brlens(phy)[phylo::descendants(shift_nodes[i-1], phy)])

            # phylogenetic signal defined as the height of the transition 
            # kernel at 0. the higher this height, the more probability is
            # concentrated on close to zero change away from the ancestral
            # value.
            delta = runif(1)
            
            tau = runif(1, 0, 1 / (2*pi*tmax*delta*delta))
            tau = rWishart(1, d, diag(tau, d, d))[,,1]

            if (i == 1)
            {
                branch.regime[phylo::root(phy)] = i
                branch.rate[phylo::root(phy)] = sum(diag(tau))
            }

            branch.regime[preorder[[i]]] = i

            for (nd in preorder[[i]])
            {
                anc = phylo::ancestors(phy)[[nd]][1]
                node.state[nd,] = mvtnorm::rmvnorm(1, node.state[anc,], tau*phylo::brlens(phy)[nd])
                branch.rate[nd] = sum(diag(tau))
            }
        }

        rownames(node.state) = c(phylo::tiplabels(phy), phylo::root(phy):phylo::Nnode(phy))

        tiprate = bm::bm.mvtiprate(node.state[1:phylo::Ntip(phy),], phy)
        tiprate = sapply(tiprate, function(v) sum(diag(v)))
        ret[,,sim] = cbind(branch.rate[1:phylo::Ntip(phy)], tiprate, 
            branch.regime[1:phylo::Ntip(phy)])
    }

    ret
}


fulltree = phylo::read.newick("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")

nshifts <- seq(0, 10, by = 1)

ntrait = 5
nsim = 10

resList <- vector('list', length(nshifts))
for (i in 1:length(nshifts)) {
	
	message('\t', i)
	resList[[i]] <-  simulate.bm(fulltree, ntrait, nsim, nshifts[i], minSize = 5)
}

# For each simulation, returns a n dimensional array (where n is number of sims): column 1 true tip rates, column 2 estimated tip rates, column 3 regime id.


# reorganize
mat <- lapply(resList, function(x) apply(x, 3, function(y) y, simplify = FALSE))
mat <- lapply(mat, function(x) do.call(rbind, x))
for (i in 1:length(mat)) {
	mat[[i]] <- cbind(mat[[i]], rep(nshifts[i], nrow(mat[[i]])), rep(1:nsim, each = Ntip(fulltree)))
	colnames(mat[[i]]) <- c('trueRates', 'estimatedRates', 'regimeID', 'nshifts', 'sim')
}
mat <- do.call(rbind, mat)
mat <- mat[order(mat[, 'nshifts'], mat[, 'sim'], mat[, 'regimeID']),]
mat <- as.data.frame(mat)
dim(mat)

# add:
# unique ID to regimes within sims within nshifts
# regime size
# mean estimated tip rate per regime
# mean proportional error
mat <- cbind(mat, uniqueID = NA, regimeSize = NA, meanTipRate = NA, meanPropError = NA, RMSE = NA)
for (i in 1:length(nshifts)) {
	for (j in 1:nsim) {
		ind <- which(mat[, 'nshifts'] == nshifts[i] & mat[, 'sim'] == j)
		regimes <- sort(unique(mat[ind, 'regimeID']))
		# table(mat[ind, 'regimeID'])
		for (k in 1:length(regimes)) {
			kind <- intersect(which(mat[, 'regimeID'] == regimes[k]), ind)
			uniqueID <- paste(nshifts[i], (1:nsim)[j], regimes[k], sep = '_')
			mat[kind, 'uniqueID'] <- uniqueID
			mat[kind, 'regimeSize'] <- length(kind)
			mat[kind, 'meanTipRate'] <- mean(mat[kind, 'estimatedRates'])
			mat[kind, 'meanPropError'] <- sum(((mat[kind, 'estimatedRates'] - mat[kind, 'trueRates']) / mat[kind, 'trueRates'])) / length(kind)
			mat[kind, 'RMSE'] <- sqrt(mean((mat[kind, 'estimatedRates'] - mat[kind, 'trueRates']) ^ 2))
		}
	}
}
head(mat)

# double-check: there should be nshifts + 1 regimes in each set.

# get per-regime mean tip rates


#############
png('~/Downloads/mvTipRateSummary.png', width = 10, height = 8, units = 'in', res = 400)

par(mfrow = c(2,2), oma = rep(0.5, 4), mar = c(4, 4, 2, 2))

### plot pooled true vs estimated tip rates

plot.new()
plot.window(xlim = range(log(mat[, c('trueRates', 'estimatedRates')])), ylim = range(log(mat[, c('trueRates', 'estimatedRates')])))

axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")

mtext('pooled true rates (log)', side = 1, line = 2.5)
mtext('pooled estimated tip rates (log)', side = 2, line = 2.5)

heatscatterpoints(log(mat[, 'trueRates']), log(mat[, 'estimatedRates']), grid = 200)
abline(a = 0, b = 1, lty = 2)

text(x = grconvertX(-0.1, from = 'npc', to = 'user'), y = grconvertY(1.1, from = 'npc', to = 'user'), 'A', xpd = NA, cex = 1.2, font = 2)

# plot true vs estimated regime means
xx <- split(mat, f = mat$uniqueID)

plot.new()
plot.window(xlim = range(log(mat[, c('trueRates', 'meanTipRate')])), ylim = range(log(mat[, c('trueRates', 'meanTipRate')])))

axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")

mtext('by-regime mean true rates (log)', side = 1, line = 2.5)
mtext('by-regime mean estimated tip rates (log)', side = 2, line = 2.5)

heatscatterpoints(log(sapply(xx, function(x) x[1, 'trueRates'])), log(sapply(xx, function(x) x[1, 'meanTipRate'])), grid = 200)
abline(a = 0, b = 1, lty = 2)

text(x = grconvertX(-0.1, from = 'npc', to = 'user'), y = grconvertY(1.1, from = 'npc', to = 'user'), 'B', xpd = NA, cex = 1.2, font = 2)


## mean proportional error as a function of regime size
plot.new()

xvals <- log(sapply(xx, function(x) x[1, 'regimeSize']))
yvals <- sapply(xx, function(x) x[1, 'meanPropError'])

# there are outliers between 5 and 15, 20 and 85, one at ~ 500, one at 145k.

plot.window(xlim = range(log(sapply(xx, function(x) x[1, 'regimeSize']))), ylim = c(-0.5, 5))

axis(1, at = axTicks(1), labels = round(exp(axTicks(1)), 0), lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")

mtext('log regime size', side = 1, line = 2.5)
mtext('mean proportional error', side = 2, line = 2.5)

points(xvals, yvals, col = adjustcolor('black', alpha.f = 0.5))
abline(h = 0, lty = 2)

text(x = grconvertX(-0.1, from = 'npc', to = 'user'), y = grconvertY(1.1, from = 'npc', to = 'user'), 'C', xpd = NA, cex = 1.2, font = 2)

### RMSE by nshifts
xx <- split(mat, mat$nshifts)
plot.new()
plot.window(xlim = range(nshifts), ylim = c(0, 3))

axis(1, at = nshifts, labels = nshifts, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")

mtext('n shifts', side = 1, line = 2.5)
mtext('root mean square error', side = 2, line = 2.5)

width <- 0.5

for (i in nshifts) {
	
	#violin(xx[[i + 1]][, 'RMSE'], center = i, width = 0.25, col = 'red', border = NA)
	qStats <- quantile(xx[[i + 1]][, 'RMSE'], c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE)
	iqr <- qStats[4] - qStats[2]
	iqr1.5 <- iqr * 1.5
	lowerWhisker <- qStats[2] - iqr1.5
	upperWhisker <- qStats[4] + iqr1.5
	if (lowerWhisker < min(xx[[i + 1]][, 'RMSE'])) {
		lowerWhisker <- min(xx[[i + 1]][, 'RMSE'])
	}
	if (upperWhisker > max(xx[[i + 1]][, 'RMSE'])) {
		upperWhisker <- max(xx[[i + 1]][, 'RMSE'])
	}

	rect(i - width/2, qStats[2], i + width/2, qStats[4], col=gray(0.75))
	segments(i, lowerWhisker, i, qStats[2], lty=2, lend=1)
	segments(i, qStats[4], i, upperWhisker, lty=2, lend=1)
	segments(i - width/3, lowerWhisker, i + width/3, lowerWhisker, lend=1)
	segments(i - width/3, upperWhisker, i + width/3, upperWhisker, lend=1)
	segments(i - width/3, qStats[3], i + width/3, qStats[3], lwd=2, lend=1)

}

text(x = grconvertX(-0.1, from = 'npc', to = 'user'), y = grconvertY(1.1, from = 'npc', to = 'user'), 'D', xpd = NA, cex = 1.2, font = 2)

dev.off()


