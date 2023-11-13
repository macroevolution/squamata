violin <- function(dat, center, width, truncate = TRUE, barlwd = 1, ptCex = 1, ...) {
	
	dat <- dat[complete.cases(dat)]
	
	if (truncate) {
		truncateRange <- range(dat, na.rm = TRUE)
	}
	
	dens <- density(dat)
	dens <- cbind(dens$x, dens$y)
	
	densR <- dens
	densR[,1] <- dens[,2]
	densR[,2] <- dens[,1]
	densR[,1] <- scales::rescale(densR[,1], to = c(center, center + width))
	if (truncate) {
		densR[densR[,2] > truncateRange[2], 2] <- truncateRange[2]
		densR[densR[,2] < truncateRange[1], 2] <- truncateRange[1]
	}
	# polygon(densR, ...)

	densL <- dens
	densL[,1] <- dens[,2]
	densL[,2] <- dens[,1]
	densL[,1] <- scales::rescale(densL[,1], to = c(center, center - width))
	if (truncate) {
		densL[densL[,2] > truncateRange[2], 2] <- truncateRange[2]
		densL[densL[,2] < truncateRange[1], 2] <- truncateRange[1]
	}
	# polygon(densL, ...)
	
	tol <- 1e-12
	densR <- densR[abs(densR[,1] - center) > tol,]
	densL <- densL[abs(densL[,1] - center) > tol,]
	densComb <- rbind(densL, densR[nrow(densR):1, ])

	polygon(densComb, ...)
	
	# 5-95% percentile bar
	qq <- quantile(dat, c(0.05, 0.95))
	segments(x0 = center, x1 = center, y0 = qq[1], y1 = qq[2], lwd = barlwd, lend = 1)
	
	# median
	points(x = center, y = median(dat), pch = 19, cex = ptCex)
	
}
