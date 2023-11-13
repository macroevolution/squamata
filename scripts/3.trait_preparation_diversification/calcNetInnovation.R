
calcNetInnovation <- function(x, phy) {
	
	library(bm)
	
	# x should be a named vector or a matrix/dataframe with rownames that match the tree
	if (is.vector(x)) {
		x <- x[complete.cases(x)]

		common <- intersect(names(x), tiplabels(phy))
		if (length(common) == 0) stop('no names shared between tree and dataset.')
		phy <- keep.tip(phy, common)
		x <- x[names(x) %in% tiplabels(phy)]
		
		tree.height <- sum(brlens(phy)[c(1, ancestors(phy)[[1]])])

		p = bm::bm.pic(x, phy)
		aTrait <- p$ace[1]
		vTrait <- attr(p$pic, "rate")

	    e1 <- sqrt((x - aTrait) ^ 2)
	    e2 <- sqrt((x - aTrait) ^ 2 / (vTrait * tree.height))

	} else {
	
		x <- x[complete.cases(x), ]
	
		common = intersect(rownames(x), tiplabels(phy))
		if (length(common) == 0) stop('no names shared between tree and dataset.')
		phy <- keep.tip(phy, common)
		x <- x[rownames(x) %in% tiplabels(phy), ]
		
		y <- data.matrix(x)
		y <- y[tiplabels(phy), ]
		
		tree.height <- sum(brlens(phy)[c(1, ancestors(phy)[[1]])])
	
		aTrait <- structure(numeric(ncol(y)))
		vTrait <- structure(numeric(ncol(y)))
		
		for (i in 1:ncol(y)) {
		    p <- bm.pic(y[,i], phy)
		    aTrait[i] <- p$ace[1]
		    vTrait[i] <- attr(p$pic, "rate")
		}
		
		e1 <- structure(numeric(nrow(x)), names = rownames(y))
		e2 <- structure(numeric(nrow(x)), names = rownames(y))
		for (i in 1:nrow(y)) {
		    e1[i] <- sqrt(sum((y[i, ] - aTrait) ^ 2))
		    e2[i] <- sqrt(mahalanobis(y[i, ], aTrait, diag(tree.height * vTrait)))
		}
	}
	
	ancdist <- cbind(euclidean = unname(e1), mahalanobis = unname(e2))
	rownames(ancdist) <- names(e1)

	return(ancdist)
}
