library(ape)

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

setwd('./scripts/4.comparative_analyses/model-fits')

source('1.load-data.R')

alldat <- read.csv('../../../data/alldat.csv')

subfamTree <- read.tree('subfamily-pairs-squams-r1.tre')

getNodeShiftEvidence <- function(dat, subfamTree, alldat, trait) {
        
    # best model type
    bestModelType <- dat[which.min(dat$aic), 'type']
    
    # for model type, get the shift locations (in terms of the subfam tree), which will be the same with added partitions, and the improvement in AIC
    modInd <- which(dat$type == bestModelType)
    
    xx <- vector('list', length(modInd))
    for (i in 1:length(modInd)) {
        
        shifts <- which(!is.na(dat[modInd[i], paste0('shiftnodes.', 1:10)]))
        nn <- rep(NA, length(shifts))
        for (k in 1:length(shifts)) {
            tax <- dat[modInd[i], paste0('span.', c(seq(1,20, by = 2)[k], seq(2,20, by = 2)[k]))]
            tax <- as.character(tax)    
                            
            # if (!all(tax %in% subfamTree$tip.label)) stop()
            subfams <- alldat[alldat$treename %in% tax, 'subfamily']
            subfamInt <- intersect(alldat[alldat$subfamily %in% subfams, 'treename'], subfamTree$tip.label)
            # if (getMRCA(subfamTree, subfamInt) != nn1) stop()
            nn[k] <- getMRCA(subfamTree, subfamInt)
        }
        xx[[i]] <- sort(nn)
    }
    
    # reorder such that common nodes are first
    ord <- sort(table(unlist(xx)), decreasing = TRUE)
    for (i in 1:length(xx)) {
        xx[[i]] <- xx[[i]][match(names(ord), as.character(xx[[i]]))]
        xx[[i]] <- xx[[i]][!is.na(xx[[i]])]
    }
    
    
    # vector of node number that is new to partition set, and
    # vector of AIC improvements, starting with comparison to base, and then relative to 1 less partition
    
    aic <- dat[modInd, 'aic']
    nodes <- xx[[which.min(aic)]]
    aicdiffs <- abs(diff(c(dat[1,'aic'], aic) - dat[1, 'aic']))
    return(list(nodes = nodes, aicdiffs = aicdiffs))

}



subfamTreeLabeled <- subfamTree
subfamLabels <- sapply(subfamTree$tip.label, function(x) paste0(alldat[alldat$treename == x, c('subfamily', 'genus')], collapse = '_/_'))
subfamTreeLabeled$tip.label <- subfamLabels


files <- c("fit_diet_noTree.csv", "fit_diet_PC_noTree.csv", "fit_diet_noTree_randomSubset.csv")

traitLabels <- c("diet full", "diet PC", "diet random subsets")

# Make more compact version for manuscript

aicRange <- list()
for (i in 1:length(files)) {
    
    resdat <- read.csv(files[i])
    trait <- gsub('(fit_)(.+)(\\.csv)', '\\2', files[i])
    
    if (!grepl('randomSubset', files[i])) {
        
        xx <- getNodeShiftEvidence(resdat, subfamTree, alldat, trait)
        nn <- xx$nodes
        nnAic <- xx$aicdiffs
            
    } else {
        
        reslist <- split(resdat, resdat$sim)
        table(sapply(reslist, function(x) x[which.min(x$aic), 'type']))
        xx <- lapply(reslist, function(x) getNodeShiftEvidence(x, subfamTree, alldat, trait))
        
        # alternative approach, showing median AIC improvement per node
        uniqueNodes <- unique(unlist(lapply(xx, function(x) x$nodes)))
        nodeAic <- matrix(nrow = length(uniqueNodes) * length(xx), ncol = 3)
        colnames(nodeAic) <- c('sim', 'node', 'aicDiff')
        nodeAic <- as.data.frame(nodeAic)
        counter <- 1
        for (j in 1:length(xx)) {
                
            for (k in 1:length(uniqueNodes)) {
                nodeAic[counter, 'sim'] <- j
                nodeAic[counter, 'node'] <- uniqueNodes[k]
                if (uniqueNodes[k] %in% xx[[j]]$nodes) {
                    nodeAic[counter, 'aicDiff'] <- xx[[j]]$aicdiffs[which(xx[[j]]$nodes == uniqueNodes[k])]
                } else {
                    nodeAic[counter, 'aicDiff'] <- 0
                }
                counter <- counter + 1
            }
        }
        
        sort(sapply(split(nodeAic$aicDiff, nodeAic$node), median))
        nn <- as.numeric(names(sort(sapply(split(nodeAic$aicDiff, nodeAic$node), median))))
        nnAic <- sort(sapply(split(nodeAic$aicDiff, nodeAic$node), median))
        nn <- nn[nnAic > 0]
        nnAic <- nnAic[nnAic > 0]
        
    }
    
    aicRange[[i]] <- nnAic
}
lapply(aicRange, range)
globalAicRange <- range(unlist(aicRange))

#maxAIC <- 3
#aicRange[2] <- maxAIC
ptSizeRange <- c(1.2,5)

pdf('~/Downloads/figS17.pdf', width = 10, height = 6)

mat <- matrix(1:8, nrow = 2, ncol = 4, byrow = TRUE)
mat[1,] <- c(1,2,3,4)
mat[2,] <- c(1,2,3,4)

layout(mat, widths = c(rep(0.9/3, 3), 0.1))
par(mar = c(1,1,1,1))

for (i in 1:length(files)) {
    
    resdat <- read.csv(files[i])
    trait <- gsub('(fit_)(.+)(\\.csv)', '\\2', files[i])
    
    if (!grepl('randomSubset', files[i])) {
        
        xx <- getNodeShiftEvidence(resdat, subfamTree, alldat, trait)
        nn <- xx$nodes
        nnAic <- xx$aicdiffs
            
    } else {
        
        reslist <- split(resdat, resdat$sim)
        table(sapply(reslist, function(x) x[which.min(x$aic), 'type']))
        xx <- lapply(reslist, function(x) getNodeShiftEvidence(x, subfamTree, alldat, trait))
        
        # alternative approach, showing median AIC improvement per node
        uniqueNodes <- unique(unlist(lapply(xx, function(x) x$nodes)))
        nodeAic <- matrix(nrow = length(uniqueNodes) * length(xx), ncol = 3)
        colnames(nodeAic) <- c('sim', 'node', 'aicDiff')
        nodeAic <- as.data.frame(nodeAic)
        counter <- 1
        for (j in 1:length(xx)) {
                
            for (k in 1:length(uniqueNodes)) {
                nodeAic[counter, 'sim'] <- j
                nodeAic[counter, 'node'] <- uniqueNodes[k]
                if (uniqueNodes[k] %in% xx[[j]]$nodes) {
                    nodeAic[counter, 'aicDiff'] <- xx[[j]]$aicdiffs[which(xx[[j]]$nodes == uniqueNodes[k])]
                } else {
                    nodeAic[counter, 'aicDiff'] <- 0
                }
                counter <- counter + 1
            }
        }
        
        sort(sapply(split(nodeAic$aicDiff, nodeAic$node), median))
        nn <- as.numeric(names(sort(sapply(split(nodeAic$aicDiff, nodeAic$node), median))))
        nnAic <- sort(sapply(split(nodeAic$aicDiff, nodeAic$node), median))
        nn <- nn[nnAic > 0]
        nnAic <- nnAic[nnAic > 0]
        
    }
    
    message('aic range: ', min(nnAic), ' - ', max(nnAic))
    # nnAic[nnAic > maxAIC] <- maxAIC
    aicRange <- range(nnAic)
    aicRange[1] <- globalAicRange[1]

    snakeNode <- intersect(subfamTree$tip.label, alldat[alldat$clade_Serpentes == 1, 'treename'])
    snakeNode <- getMRCA(subfamTree, snakeNode)

    ColubriformesNode <- intersect(subfamTree$tip.label, alldat[alldat$clade_Colubriformes == 1, 'treename'])
    ColubriformesNode <- getMRCA(subfamTree, ColubriformesNode)

    edge.color = rep(gray(0.65), nrow(subfamTree$edge))
    edge.color[which(subfamTree$edge[,1] %in% c(snakeNode, phangorn::Descendants(subfamTree, snakeNode, type = 'all')))] = '#d95f02'
    edge.color[which(subfamTree$edge[,1] %in% c(ColubriformesNode, phangorn::Descendants(subfamTree, ColubriformesNode, type = 'all')))] = '#386cb0'
    table(edge.color)
    
    plot.phylo(subfamTreeLabeled, show.tip.label = FALSE, type = 'phylogram', open.angle = 10, lwd = 0.4, edge.color = edge.color, lend = 1)

    nodeSize <- scales::rescale(nnAic, from = aicRange, to = ptSizeRange)
    
    nodelabels(node = nn, frame = 'circle', pch = 21, bg = adjustcolor('coral', alpha.f = 0.75), col = adjustcolor('black', alpha.f = 0.75), cex = nodeSize, lwd = 0.5)
    
    
    title(main = traitLabels[i], line = 0)
}

plot.new()
plot.window(xlim = c(0,10), ylim = c(0,10))

ptY <- c(0.1, 0.14, 0.2)

points(x = grconvertX(rep(0.2, 3), from = 'npc', to = 'user'), y = grconvertY(ptY, from = 'npc', to = 'user'), cex = seq(ptSizeRange[1], ptSizeRange[2], length.out = 3), bg = adjustcolor('coral', alpha.f = 0.75), col = adjustcolor('black', alpha.f = 0.75), pch = 21, xpd = NA)
# text(x = grconvertX(rep(0.37, 3), from = 'npc', to = 'user'), y = grconvertY(ptY, from = 'npc', to = 'user'), labels = round(seq(min(aicRange), max(aicRange), length.out = 3), 3), pos = 4)
text(grconvertX(0.2, from = 'npc', to = 'user'), y = grconvertY(0.24, from = 'npc', to = 'user'), labels = 'improvement\nin AIC', pos = 3, xpd = NA, font = 2)

points(x = rep(-0.5, 3), y = seq(8, 9, length.out = 3), pch = 22, bg = c('#386cb0', '#d95f02', gray(0.65)), lwd = 0.5, cex = 3, xpd = NA)
text(x = rep(0, 3), y = seq(8, 9, length.out = 3), labels = c('Colubriformes', 'snakes', 'squamates'), pos = 4, xpd = NA)




dev.off()

