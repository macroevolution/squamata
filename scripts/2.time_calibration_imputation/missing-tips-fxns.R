



# Get elements of x that are not in y,
#   when x and y are numeric vectors.
#   Things are "same" when abs(x - y) < tol

setdiff_numeric <- function(x, y, tol = 1e-07){
	
	noMatch <- logical(length(x))
	for (ii in 1:length(x)){
		mm <- min(abs(x[ii] - y))
		if (mm > tol){
			noMatch[ii] <- TRUE		 
		}
	}
	return(x[noMatch])
}

 

addMissingTips <- function(phy, c1, c2, prefix = NULL, missing, in_tree, stem=TRUE, lambda=NULL){
	
 
		
	tree_age <- as.numeric(max(branching.times(phy)))
	tipnames <- paste("xzx", 1:missing, sep="")
	if (!is.null(prefix)){
		 tipnames <- paste(prefix, "_", tipnames, sep="")
	}

	# Get Stem age for the clade, whether crown or stem
    #   CRBD estimator applied to this
    stem_age <- NA
    if (is.na(c2)){
    	nn <- which(phy$tip.label == c1)
    	stem_age <- phy$edge.length[phy$edge[,2] == nn]
    }else{
    	nn <- getMRCA(phy, tip=c(c1, c2))
    	nn <- phy$edge[,1][phy$edge[,2] == nn]
    	stem_age <- as.numeric(branching.times(phy)[as.character(nn)])
    }
	if (is.null(lambda)){
		lambda <- bd.ms(time = stem_age, n = missing + in_tree, crown=F, epsilon=0)		
	}

	
	# case 1: stem addition to clade with >= 2 sampled lineages
	if (!is.na(c2) & stem){
		
		anc <- getMRCA(phy, tip = c(c1, c2))
		subtree <- extract.clade(phy, node = anc)
		subtimes <- as.numeric(branching.times(subtree))
		clade_age <- max(subtimes)
		csbt      <- corsim(subtimes, lambda=lambda, mu=0, missing = missing, told=stem_age)
		stimes    <- tree_age - setdiff_numeric(csbt, subtimes)

		if (length(stimes) != missing){
			stop("mismatched speciation times / case1")
		}

		for (ii in 1:length(tipnames)){
			
			phy <- BAMMtools:::getStartStopTimes(phy)
			
			in_tree_tip_set <- extract.clade(phy, node = getMRCA(phy, c(c1, c2)))$tip.label
			in_tree_tip_set <- c(in_tree_tip_set, intersect(tipnames, phy$tip.label))
			
		 
			goodset <- getDescNodes(phy, anc_node = getMRCA(phy, tip= in_tree_tip_set)) 
			
			good2   <- phy$edge[,2][phy$begin <= stimes[ii] & phy$end > stimes[ii]]
			
			goodset <- intersect(goodset, good2)
			
			focal_node <- NA
			if (length(goodset) == 1){
				focal_node <- goodset
			}else{
				focal_node <- sample(goodset, 1)
			}
			 
			edge_length <- tree_age - stimes[ii]
			
			pos <- NA
			
			if (focal_node <= length(phy$tip.label)){
				pos  <- edge_length
			}else{
				pos <- as.numeric(edge_length - branching.times(phy)[as.character(focal_node)])
			}
		 
			phy <- phytools:::bind.tip(phy, tip.label=tipnames[ii], where=focal_node, 
							    edge.length=edge_length, position=pos)
			
		}
		
		
	}
	
	
	# case 2: crown addition to clade with >= 2 samples lineages
	if (!is.na(c2) & !stem){
		
		anc <- getMRCA(phy, tip = c(c1, c2))
		subtree <- extract.clade(phy, node = anc)
		subtimes <- as.numeric(branching.times(subtree))
		clade_age <- max(subtimes)
		csbt      <- corsim(subtimes, lambda=lambda, mu=0, missing = missing, told=clade_age)
		stimes    <- tree_age - setdiff_numeric(csbt, subtimes, tol=1e-07)

		if (length(stimes) != missing){
			stop("mismatched speciation times / case2")
		}

		for (ii in 1:length(tipnames)){
			
			phy <- BAMMtools:::getStartStopTimes(phy)
			
			# next line not safe for stem addition as MRCA spanning set can change
			goodset <- getDescNodes(phy, anc_node = getMRCA(phy, tip=c(c1, c2))) 
			
			good2   <- phy$edge[,2][phy$begin <= stimes[ii] & phy$end > stimes[ii]]
			
			goodset <- intersect(goodset, good2)
			
			focal_node <- NA
			if (length(goodset) == 1){
				focal_node <- goodset
			}else{
				focal_node <- sample(goodset, 1)
			}
			
			edge_length <- tree_age - stimes[ii]
			
			pos <- NA
			
			if (focal_node <= length(phy$tip.label)){
				pos  <- edge_length
			}else{
				pos <- as.numeric(edge_length - branching.times(phy)[as.character(focal_node)])
			}
		 
			phy <- phytools:::bind.tip(phy, tip.label=tipnames[ii], where=focal_node, 
							    edge.length=edge_length, position=pos)
			
		}
				
		
		
	}
	
	
	# case 3: addition to singleton lineage 
	if (is.na(c2)){
		
		nd <- which(phy$tip.label == c1)
 		clade_age <- phy$edge.length[phy$edge[,2] == nd]
		csbt      <- corsim(clade_age, lambda=lambda, mu=0, missing = missing, told=clade_age)
		stimes    <- tree_age - setdiff_numeric(csbt, clade_age)

		edge_length <- tree_age - stimes[1]

		## add the first tip:
		phy <- phytools:::bind.tip(phy, tip.label=tipnames[1], where=nd, 
							    edge.length=edge_length, position=edge_length)
		
		if (missing >= 2){
			
			for (ii in 2:missing){
				taxset <- c(c1, tipnames[1:(ii-1)])		
				phy <- BAMMtools:::getStartStopTimes(phy)	
				goodset <- getDescNodes(phy, anc_node = getMRCA(phy, tip=taxset)) 
				good2   <- phy$edge[,2][phy$begin <= stimes[ii] & phy$end > stimes[ii]]
				goodset <- intersect(goodset, good2)

				focal_node <- NA
				
				if (length(goodset) == 1){
					focal_node <- goodset
				}else{
					focal_node <- sample(goodset, 1)
				}

				edge_length <- tree_age - stimes[ii]
				pos <- NA
				if (focal_node <= length(phy$tip.label)){
					pos  <- edge_length
				}else{
					pos <- as.numeric(edge_length - branching.times(phy)[as.character(focal_node)])
				}
				phy <- phytools:::bind.tip(phy, tip.label=tipnames[ii], where=focal_node, 
							    edge.length=edge_length, position=pos)
			
			} # for ii in 2:missing
			
		} # if missing >= 2

	} # case 3
	
	return(phy)
}


getDescNodes <- function(phy, anc_node){
	
	if (anc_node <= length(phy$tip.label)){
		return(anc_node)
	}else{
		dset <- phy$edge[,2][phy$edge[,1] == anc_node]
		return(c(anc_node, getDescNodes(phy, dset[1]), getDescNodes(phy, dset[2])))
	}
	
} 


## cc = dataframe with inclusion/exclusion information

addMissingTips_wExclusions <- function(phy, cc, prefix = NULL, lambda=NULL){
 
	
	# We will estimate CRBD rates for the crown subtree after dropping
	#   the exclusions
	
	cc_i <- cc[cc$constraint_type == "include", ]
	cc_e <- cc[cc$constraint_type == "exclude", ]
	
	missing <- cc_i$n_missing
	prefix  <- cc_i$genus
	c1 <- cc_i$X1
	c2 <- cc_i$X2
	
	
	tree_age <- as.numeric(max(branching.times(phy)))
	tipnames <- paste("xzx", 1:missing, sep="")
	if (!is.null(prefix)){
		 tipnames <- paste(prefix, "_", tipnames, sep="")
	}

	tip_sets <- list()
	
	nn <- getMRCA(phy, tip=c(c1, c2))
	xtmp <- extract.clade(phy, node = nn)
	
	tip_sets[[1]] <- xtmp$tip.label
	
	for (ii in 1:nrow(cc_e)){
		tmp <- extract.clade(phy, node= getMRCA(phy, c(cc_e$X1[ii], cc_e$X2[ii])))
		tip_sets[[ii+1]] <- tmp$tip.label
		tip_sets[[1]] <- setdiff(tip_sets[[1]], tip_sets[[ii+1]])
	}
	
	# tip_sets[[1]] is the full set of taxa for the "include" set after eliminating all the excludes
	# tip_sets[[2]] and higher are all of the exclusions


	# Get Stem age for the clade, whether crown or stem
    #   CRBD estimator applied to this
 
 	in_tree <- length(tip_sets[[1]])
 	
    nn <- phy$edge[,1][phy$edge[,2] == nn]
    
    stem_age <- as.numeric(branching.times(phy)[as.character(nn)])
   
	if (is.null(lambda)){
		lambda <- bd.ms(time = stem_age, n = missing + in_tree, crown=F, epsilon=0)		
	}
 
 	## Simulate the missing times:
 	
 	anc <- getMRCA(phy, tip = c(c1, c2))
	subtree <- extract.clade(phy, node = anc)
	
	# Clip out the exclusions
	for (ii in 2:length(tip_sets)){
		subtree <- drop.tip(subtree, tip = tip_sets[[ii]])
	}
	subtimes <- as.numeric(branching.times(subtree))
	clade_age <- max(subtimes)
	csbt      <- corsim(subtimes, lambda=lambda, mu=0, missing = missing, told=stem_age)
 
	stimes    <- tree_age - setdiff_numeric(csbt, subtimes)

	if (length(stimes) != missing){
		stop("mismatched speciation times / case1")
	}

	# Taxon additions:
	for (ii in 1:length(tipnames)){
		
		phy <- BAMMtools:::getStartStopTimes(phy)
 
		in_tree_tip_set <- c(tip_sets[[1]], intersect(tipnames, phy$tip.label))
		
		# Full set of nodes before removing exclusions:
		pre_set <- getDescNodes(phy, anc_node = getMRCA(phy, tip= in_tree_tip_set)) 
		
		# Eliminate the nodes from clade exclusions
		for (kk in 2:length(tip_sets)){
			tmp <- getDescNodes(phy, anc_node = getMRCA(phy, tip = tip_sets[[kk]]))
			pre_set <- setdiff(pre_set, tmp)
		}
		
		# Nodes that have matching time criteria
		good_times   <- phy$edge[,2][phy$begin <= stimes[ii] & phy$end > stimes[ii]]
		
		goodset <- intersect(pre_set, good_times)
		
		focal_node <- NA
		if (length(goodset) == 1){
			focal_node <- goodset
		}else{
			focal_node <- sample(goodset, 1)
		}
		 
		edge_length <- tree_age - stimes[ii]
		
		pos <- NA
		
		if (focal_node <= length(phy$tip.label)){
			pos  <- edge_length
		}else{
			pos <- as.numeric(edge_length - branching.times(phy)[as.character(focal_node)])
		}
	 
		phy <- phytools:::bind.tip(phy, tip.label=tipnames[ii], where=focal_node, 
						    edge.length=edge_length, position=pos)
		
	}
		
 	return(phy)
}








