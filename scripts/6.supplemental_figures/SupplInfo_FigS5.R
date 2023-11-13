
require(ape)
require(MCMCtreeR)
require(coda)
require(pbapply)

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

mcmcTreeDir <- './data/2.time_calibration_imputation/MCMCTreeAnalyses'
outdir <- '~/Downloads'


# getSpanningTips
#	returns a pair of tips that span a given node. 
#	if the node is terminal, includes "NA"
getSpanningTips <- function(phy, node){	
    if (node <= length(phy$tip.label)){
        return(c(phy$tip.label[node], 'NA'));
    }else{
        dnodes <- phy$edge[,2][phy$edge[,1] == node];
        
        while (dnodes[1] > length(phy$tip.label)){
            dnodes[1] <- phy$edge[,2][phy$edge[,1] == dnodes[1]][1];
        }
        while (dnodes[2] > length(phy$tip.label)){
            dnodes[2] <- phy$edge[,2][phy$edge[,1] == dnodes[2]][1];
        }		
        dset <- phy$tip.label[dnodes];
        return(dset);
    }	
}

# fossil calibration data
fossildat <- read.csv('./data/2.time_calibration_imputation/fossildat.csv', stringsAsFactors=FALSE)

# check convergence properties
lnlList <- list()
counter <- 1

reps <- c('iterA', 'iterB', 'iterC')

for (i in 1:length(reps)) {
	setwd(paste0(mcmcTreeDir, '/mcmctree_useData2/', reps[i]))
	message('\t', reps[i])
            
	mcmc <- data.table::fread('mcmc.txt', data.table=FALSE)
	divtimeES <- pbapply(mcmc[, grep('t_n', colnames(mcmc))], 2, effectiveSize)
	sort(divtimeES)
	divtimeES <- round(divtimeES, 2)
	message('\t\tlnL ES: ', round(effectiveSize(mcmc$lnL), 2))
	message('\t\tdivtime ES range: ', min(divtimeES), '-', max(divtimeES))
	message('\t\tmedian lnL: ', median(mcmc[, 'lnL']))
	setwd('../../../../../')
}

## Iteration A has best properties -> this will be the run we use.


# retrieve prior and posterior information by fossil


tr1 <- readMCMCtree(paste0(mcmcTreeDir, '/mcmctree_useData2/', reps[1], '/FigTree.tre'))
fossildat <- fossildat[which(sapply(fossildat$spanning.taxa, function(x) all(strsplit(x, '\\|')[[1]] %in% tr1[[1]]$tip.label)) == TRUE),]

# receiving list
tmp <- replicate(length(reps), list(), simplify = FALSE)
names(tmp) <- reps
tmp <- replicate(1, tmp, simplify = FALSE)
names(tmp) <- 'AR'
tmp <- replicate(1, tmp, simplify = FALSE)
names(tmp) <- 'calibS3'
tmp <- replicate(nrow(fossildat)+1, tmp, simplify = FALSE)
calibPlotData <- tmp

trFile <- paste0(mcmcTreeDir, '/fossilCalibrations_S3.tre')
mcmcPrior <-data.table::fread(paste0(mcmcTreeDir, '/mcmctree_useData0/mcmc.txt'), data.table=F)
mcmcPost <- data.table::fread(paste0(mcmcTreeDir, '/mcmctree_useData2/', reps[1], '/mcmc.txt'), data.table=F)
            
qq <- priorPosterior(MCMCPrior = mcmcPrior, MCMCPosterior = mcmcPost, inputTree = trFile, return.density = TRUE, rootCalibration = NULL)
            
for (i in 1:length(qq[[3]])) {
                
	message('\tfossil ', i)
                
	nodeNum <- strsplit(names(qq[[3]])[i], '_')[[1]][1]
	priorType <- strsplit(names(qq[[3]])[i], '_')[[1]][2]
                
	userPrior <- qq[[3]][[i]] # mcmcTree parameters
	effectivePrior <- qq$prior[[i]] # marginal prior density
	posterior <- qq$posterior[[i]] # posterior density
                
	if (!grepl('=', priorType)) {
		priorType <- c(ST='skewT', SN='skewNormal', L='cauchy', G='gamma', B= 'bound')[priorType]
                    
		userPrior <- plotMCMCtree(userPrior, method= priorType, lowerTime = min(effectivePrior$x)*0.95, upperTime = max(effectivePrior$x)*1.05, plotMCMCtreeData=FALSE)
                    
		calibPlotData[[i]][[1]][[1]][[1]] <- list(userPrior = userPrior,
												effectivePrior = effectivePrior,
												posterior = posterior,
												spanningTaxa = getSpanningTips(tr1[[1]], nodeNum)
											)
	} else {
                    
		# fixed root
		calibPlotData[[i]][[1]][[1]][[1]] <- list(userPrior = userPrior,
												effectivePrior = effectivePrior,
												posterior = posterior,
												spanningTaxa = getSpanningTips(tr1[[1]], nodeNum)
											)
	}
}
            
# add crown squamata root
nodeNum <- getMRCA(tr1[[1]], c('Coleonyx_variegatus', 'Hydrophis_platurus'))
            
effectivePrior <- density(mcmcPrior[, paste0('t_n', nodeNum)])
posterior <- density(mcmcPost[, paste0('t_n', nodeNum)])
            
calibPlotData[[length(qq[[3]])+1]][[1]][[1]][[1]] <- list(userPrior = NA,
														effectivePrior = effectivePrior,
														posterior = posterior,
														spanningTaxa = getSpanningTips(tr1[[1]], nodeNum))



# colors: user prior, effective prior, posterior
cols <- c('#66c2a5', '#fc8d62', '#8da0cb')

# order chronologically
## use effective prior
chronoOrder <- numeric(length(calibPlotData))
for (i in 1:length(calibPlotData)) {
    dat <- calibPlotData[[i]]
    chronoOrder[i] <- min(dat[[1]][[1]][[1]]$effectivePrior$x)
}	
chronoOrder <- order(chronoOrder, decreasing = TRUE)

timeRange <- c(3.50, 0)

    
pdf(paste0(outdir, '/calibrationDensities_ddBD_part1.pdf'), width=15, height=20)
        
layout(matrix(1:16, ncol = 2))
        
for (i in 1:16) {
            
	dat <- calibPlotData[[chronoOrder[i]]]
	names(dat)
	names(dat[[1]])
	names(dat[[1]][[1]])
            
	allnodes <- c()
	for (l in 1:length(fossildat$spanning.taxa)) {
		if (all(strsplit(fossildat$spanning.taxa[l], '\\|')[[1]] %in% tr1[[1]]$tip.label)) {
			allnodes[l] <- getMRCA(tr1[[1]], strsplit(fossildat$spanning.taxa[l], '\\|')[[1]])
		}
	}
            
	fossilInd <- which(allnodes == getMRCA(tr1[[1]], dat[[1]][[1]][[1]]$spanningTaxa))
            
	if (all(is.na(dat[[1]][[1]][[1]]$userPrior))) {
		y1 <- 0
	} else {
		y1 <- unlist(lapply(dat[[1]][[1]], function(y) y$userPrior[,2]))
	}
	y2 <- unlist(lapply(dat[[1]][[1]], function(y) y$effectivePrior$y))
	y3 <- unlist(lapply(dat[[1]][[1]], function(y) y$posterior$y))
	yrange <- c(0, max(c(y1, y2, y3)))
            
	plot.new()
	plot.window(xlim = timeRange, ylim = yrange)
	axis(1, at = seq(from = timeRange[1], to = timeRange[2], by = -0.25), labels = seq(from = timeRange[1], to = timeRange[2], by = -0.25)*100, cex.axis = 0.85, lwd = 0, lwd.ticks = 1)
	axis(2, cex.axis = 0.85, lwd = 0, lwd.ticks = 1, las = 1)
	box(which = 'plot', bty = 'l')
	mtext('mya before present', side = 1, line = 2.5)
	mtext('density', side=2, line = 2.5)			
            
	# User prior is the same regardless of treatment
	if (all(dat[[1]][[1]][[1]]$userPrior == 1) | anyNA(dat[[1]][[1]][[1]]$userPrior)) {
	} else {
		lines(dat[[1]][[1]][[1]]$userPrior, col = cols[1], lwd = 2)
	}
            
	# Effective prior
	## effective priors seem to be essentially the same regardless of clock
	for (l in 1:length(reps)) {
		lines(dat[[1]][[1]][[l]]$effectivePrior, col = cols[2], lwd = 2, xpd = NA)
	}
            
	# Posterior
	for (l in 1:length(reps)) {
		lines(dat[[1]][[1]][[l]]$posterior, col = cols[3], lwd = 2, xpd = NA)
	}
            
	if (getMRCA(tr1[[1]], dat[[1]][[1]][[1]]$spanningTaxa) == getMRCA(tr1[[1]], c('Coleonyx_variegatus', 'Hydrophis_platurus'))) {
		nodeName <- 'crown squamata'
	} else {
		nodeName <- fossildat[fossilInd, 'name']
		abline(v = fossildat[fossilInd, c('min.age', 'max.age')]/100, lty = 2)
	}
            
	if (i %in% 4:5) {
		mtext(paste(nodeName, '\n', fossildat[fossilInd, 'group']), side = 3, font = 2)	
	} else {
		mtext(bquote(atop(bolditalic(.(nodeName)), bold(.(fossildat[fossilInd, 'group'])))), side = 3, font = 2)
	}        
            
	legend('topright', legend = c('prior', 'effective prior', 'posterior'), fill = cols, cex = 1, bty = 'n')
            
}
dev.off()



# Part2
      
pdf(paste0(outdir, '/calibrationDensities_ddBD_part2.pdf'), width=15, height=20)

layout(matrix(1:16, ncol = 2))
        
for (i in 17:length(calibPlotData)) { 
 
	dat <- calibPlotData[[chronoOrder[i]]]
	names(dat)
	names(dat[[1]])
	names(dat[[1]][[1]])
            
	allnodes <- c()
	for (l in 1:length(fossildat$spanning.taxa)) {
		if (all(strsplit(fossildat$spanning.taxa[l], '\\|')[[1]] %in% tr1[[1]]$tip.label)) {
			allnodes[l] <- getMRCA(tr1[[1]], strsplit(fossildat$spanning.taxa[l], '\\|')[[1]])
		}
	}
            
	fossilInd <- which(allnodes == getMRCA(tr1[[1]], dat[[1]][[1]][[1]]$spanningTaxa))
            
	if (all(is.na(dat[[1]][[1]][[1]]$userPrior))) {
		y1 <- 0
	} else {
		y1 <- unlist(lapply(dat[[1]][[1]], function(y) y$userPrior[,2]))
	}
	y2 <- unlist(lapply(dat[[1]][[1]], function(y) y$effectivePrior$y))
	y3 <- unlist(lapply(dat[[1]][[1]], function(y) y$posterior$y))
	yrange <- c(0, max(c(y1, y2, y3)))
            
	plot.new()
	plot.window(xlim = timeRange, ylim = yrange)
	axis(1, at = seq(from = timeRange[1], to = timeRange[2], by = -0.25), labels = seq(from = timeRange[1], to = timeRange[2], by = -0.25)*100, cex.axis = 0.85, lwd = 0, lwd.ticks = 1)
	axis(2, cex.axis = 0.85, lwd = 0, lwd.ticks = 1, las = 1)
	box(which = 'plot', bty = 'l')
	mtext('mya before present', side = 1, line = 2.5)
	mtext('density', side=2, line = 2.5)			
            
	# User prior is the same regardless of treatment
	if (all(dat[[1]][[1]][[1]]$userPrior == 1) | anyNA(dat[[1]][[1]][[1]]$userPrior)) {
	} else {
		lines(dat[[1]][[1]][[1]]$userPrior, col = cols[1], lwd = 2)
	}
            
	# Effective prior
	## effective priors seem to be essentially the same regardless of clock
	for (l in 1:length(reps)) {
		lines(dat[[1]][[1]][[l]]$effectivePrior, col = cols[2], lwd = 2, xpd = NA)
	}
            
	# Posterior
	for (l in 1:length(reps)) {
		lines(dat[[1]][[1]][[l]]$posterior, col = cols[3], lwd = 2, xpd = NA)
	}
            
	if (getMRCA(tr1[[1]], dat[[1]][[1]][[1]]$spanningTaxa) == getMRCA(tr1[[1]], c('Coleonyx_variegatus', 'Hydrophis_platurus'))) {
		nodeName <- 'crown squamata'
	} else {
		nodeName <- fossildat[fossilInd, 'name']
		abline(v = fossildat[fossilInd, c('min.age', 'max.age')]/100, lty = 2)
	}
            
	if (i %in% 4:5) {
		mtext(paste(nodeName, '\n', fossildat[fossilInd, 'group']), side = 3, font = 2)	
	} else {
		mtext(bquote(atop(bolditalic(.(nodeName)), bold(.(fossildat[fossilInd, 'group'])))), side = 3, font = 2)
	}        
            
	legend('topleft', legend = c('prior', 'effective prior', 'posterior'), fill = cols, cex = 1, bty = 'n')
 
}
dev.off()
         
        
        
       