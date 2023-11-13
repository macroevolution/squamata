


 

isShiftConfigValid <- function( dmat, bt, shifts , min.per.regime = 3){
	
	# Requires 3 tips in each shift regime
	bt <- sort(bt[as.character(shifts)], decreasing=TRUE)
	nodes <- as.numeric(names(bt))
	
	tips <- numeric(nrow(dmat))
	for (ii in 1:length(nodes)){
		tips[dmat[ ,nodes[ii] ] == 1] <- ii
	}
	tx <- table(factor(tips, levels=as.character(0:length(nodes))))
	if (sum(tx < min.per.regime) > 0){
		return(FALSE)
	}else{
		return(TRUE)
	}
}


getShiftNodes <- function(obj){
	
	if (length(obj$rate) == 1){
		return(NULL)
	}else{
		return(as.integer(names(obj$rate)[2:length(obj$rate)]))
	}
 
}



NA_shift_obj <- function(){
 
 	structure(list(
        rate=NA
        , logL=NA
        , aic=NA
        , npar=NA
        , converged=NA
        , msg= "likelihood-error"
        , auto.dropped=NA
        , phy.scaled=NA
        , phy.original=NA
    ), class="fit-error")
	
}

return_all_bad <- function(KMAX){
	
	fit0 <- NA_shift_obj()
	
	fit_NC <- vector(mode = "list", length = KMAX)
	fit_C <- vector(mode = "list", length = KMAX)
	
	for (ii in 1:KMAX){
		fit_NC[[ii]] <- NA_shift_obj()
		fit_C[[ii]]  <- NA_shift_obj()
	}
	ll <- list(fit0, fit_NC , fit_C)
	return(ll)
}

DEBUG <- FALSE
VERBOSE <- FALSE

if (DEBUG){
 x <- D$y[,c(1,3,5, 7 , 9)]
 phy <- D$phy
 shifts <- possible_shifts
 KMAX <- 5
 ii <- 1
 ignore.cov=F
 dmat <- dmat	
 run.shift = F
 run.censored=T
}
 
 
shift_search_step <- function(x, phy, shifts, KMAX = 5, ignore.cov = T, dmat, print.status=T, run.shift = T, run.censored = T ){
 
	bt <- ape::branching.times(phy)
 
 	fit0 <- fit_shift(x, phy, ignore.cov=ignore.cov, auto.drop=F) 
	
	if (inherits(fit0, "error")){
		cat("Initial likelihood (fit0) evaluation failed\n")
		return(return_all_bad(KMAX) )
	}
	
	fit_C <- vector(mode = "list", length = KMAX)
	fit_NC <- fit_C
	
	best_C <- fit0
	best_NC <- fit0
		
	shifts_C <- shifts
	shifts_NC <- shifts
	
	for (ii in 1:KMAX){
		
		if (print.status){
			cat("Fitting: Shift with ", ii, "shifts\n")			
		}

		AA <- getShiftNodes(best_C)
		BB <- getShiftNodes(best_NC) 
 		shifts_C <- setdiff(shifts_C, AA)
		shifts_NC <- setdiff(shifts_NC, BB)

		# Track likelihood for each node partition
		#  this can also index how many failures are 
		#  occurring (-Inf if so).

		llvec_C <- rep(NA, length(shifts_C))
		llvec_NC <- rep(NA, length(shifts_NC))
 		names(llvec_C) <- shifts_C
 		names(llvec_NC) <- shifts_NC
 
		# uncensored fits:
 		if (run.shift){
 			
 		for (kk in 1:length(shifts_NC)){
			
			if (VERBOSE) cat("NC: ", kk, "\t", shifts_NC[kk], "\n")
								
			tmp <- fit_shift(x, phy, c(BB, shifts_NC[kk]), 
						ignore.cov = ignore.cov, auto.drop=F )
			
			if (inherits(tmp, "error")){
				llvec_NC[ kk ] <- -Inf
 
			}else{
 				llvec_NC[ kk ] <- tmp$logL
 				if (kk == 1){
					best_NC <- tmp
				}
			
				if (tmp$logL > best_NC$logL){
					best_NC <- tmp
				}					
			}
		} } # IF run.shift
		
		best_NC$llvec <- llvec_NC 
		fit_NC[[ii]] <- best_NC
		
		if (print.status){
			cat("Fitting: Censored with ", ii+1, "partitions\n")			
		}		
		
		
		if (run.censored){
		for (kk in 1:length(shifts_C)){
			
			if (VERBOSE) cat("C: ", kk, "\t", shifts_C[kk], "\n")	
			
			tset <- c(AA, shifts_C[kk])
			isGood <- isShiftConfigValid(dmat, bt, tset)
			
			if (isGood){
				
				tmp <- fit_censored(x, phy, tset, 
							ignore.cov=ignore.cov)
				
				if (inherits(tmp, "error")){
					llvec_C[ kk ] <- -Inf
				}else{
					llvec_C[ kk ] <- tmp$logL
					if (kk == 1){
						best_C <- tmp
					}
					if (tmp$logL > best_C$logL){
						best_C <- tmp
					}					
				}
 
			}
 
		} } # IF run.censored
		best_C$llvec <- llvec_C
		fit_C[[ii]] <- best_C
	}
	
	ll <- list(fit0, fit_NC , fit_C)
	return(ll)

}



summarizeStepwiseFits <- function(x, max_shifts = 10){
	
	nn <- length(x[[2]])
	
	tt <- c("base", rep("shift", nn), rep("censored", nn))
	aicvec <- rep(NA, length(tt))
	likvec <- rep(NA, length(tt))
	nshifts <- c(0, 1:nn, 1:nn)
	npar <- rep(NA, length(tt))
	
	shiftnodemat <- matrix(NA, nrow=length(tt), ncol=max_shifts)
	spanmat      <- matrix("", nrow=length(tt), ncol=2*max_shifts)
	
	aicvec[1] <- x[[1]]$aic
	likvec[1] <- x[[1]]$logL
	npar[1] <- x[[1]]$npar
	for (ii in 1:nn){
		
		tmp <- getShiftNodes(x[[2]][[ii]])
		
		for (kk in 1:length(tmp)){
			m1 <- 2*(kk - 1) + 1
			spanmat[ii+1,  m1:(m1+1)] <- getSpanningTips(x[[2]][[ii]]$phy.original, tmp[kk])
		}
		
		shiftnodemat[ii+1, 1:length(tmp)] <- tmp 
		aicvec[ii+1] <- x[[2]][[ii]]$aic
		likvec[ii+1] <- x[[2]][[ii]]$logL
		npar[ii + 1] <- x[[2]][[ii]]$npar
		
	}
	
	for (ii in 1:nn){
		tmp <- getShiftNodes(x[[3]][[ii]])
		
		for (kk in 1:length(tmp)){
			m1 <- 2*(kk - 1) + 1
			spanmat[ii+ nn + 1,  m1:(m1+1)] <- getSpanningTips(x[[3]][[ii]]$phy.original, tmp[kk])
		}		
		
		
		shiftnodemat[ii+1 + nn, 1:length(tmp)] <- tmp 
		aicvec[ii + nn + 1] <- x[[3]][[ii]]$aic
		likvec[ii+1 + nn] <- x[[3]][[ii]]$logL
		npar[[ii + nn + 1]] <- x[[3]][[ii]]$npar
	}
	
	aicvec <- aicvec - min(aicvec)
	
	return(data.frame(
			   aic=aicvec, 
			   logL = likvec, 
			   type=tt, 
			   npar, 
			   nshifts, 
			   shiftnodes=shiftnodemat,
			   span = spanmat
			   ))
	
}


subsetDataRandom <- function(x, ntips){
	
	keep <- sample(x$phy$tip.label, size=ntips)
	p2 <- ape::drop.tip(x$phy, setdiff(x$phy$tip.label, keep))
	y2 <- x$y[p2$tip.label, 1, drop=FALSE]
 
	return(list(y=y2, phy=p2))
}

 






