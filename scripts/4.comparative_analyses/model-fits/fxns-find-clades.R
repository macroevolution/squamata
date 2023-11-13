find_clades = function(phy, fulltree, clades)
{
    # spanning pairs in fulltree. modify to add clades of interest, which
    # can then be references in the `clades` argument
    CLADE_LOOKUP = list(
        Serpentes=c("Namibiana_labialis","Coluber_constrictor")
        , Macrostomata=c("Python_regius","Coluber_constrictor")
        , Colubroides=c("Fimbrios_klossi","Coluber_constrictor")
        , Viperidae=c("Proatheris_superciliaris","Bothrops_marajoensis")
        , Elapidae=c("Calliophis_melanurus","Hydrophis_lamberti")
        , Dipsadinae=c("Heterodon_platirhinos","Alsophis_sibonius")
    )

    stopifnot(all(clades %in% names(CLADE_LOOKUP)))

    sapply(clades, function(g) {

        sp = match(CLADE_LOOKUP[[g]], fulltree$tip.label)

        tips = ape::extract.clade(fulltree, ape::getMRCA(fulltree, sp))$tip.label

        tips = tips[ tips %in% phy$tip.label ]

        if (length(tips) > 1L)
            return (ape::getMRCA(phy, tips))
        else
            return (NA_integer_)
    })
}


# Get spanning tips for a given node
getSpanningTips <- function(phy, node = "root"){
	if (node == "root"){
		node <- length(phy$tip.label) + 1
	}
	
	dset <- phy$edge[,2][phy$edge[,1] == node]
	
	d1 <- dset[1]
	d2 <- dset[2]
	while (d1 > length(phy$tip.label)){
		d1 <- phy$edge[,2][phy$edge[,1] == d1][1]
	}
	while (d2 > length(phy$tip.label)){
		d2 <- phy$edge[,2][phy$edge[,1] == d2][1]
	}
	
	return(phy$tip.label[c(d1, d2)])
}


buildSpanningSet_DF <- function(phy){
	# all nodes except the root:
	nodes <- (length(phy$tip.label)+2):(2*length(phy$tip.label) - 1)
	mm    <- matrix("", nrow=2, ncol=length(nodes))
	colnames(mm) <- paste("a", nodes, sep="")
	for (ii in 1:length(nodes)){
		mm[ , ii] <- getSpanningTips(phy, node=nodes[ii])
	}
	return(as.data.frame(mm))
}




# pset: dataframe of spanning tips for set of clades of interest
#       Number of columns = number of clades
#			Column names  = clade names
#           2 rows: each the spanning pair 

find_clades_pairs = function(phy, fulltree, pset, clip.NA = TRUE, mintips = 3)
{
 
    tipset <- unique(as.vector(unlist(pset)))
	stopifnot(all(tipset %in% fulltree$tip.label))

	# build clade lookup table:

	clades <- colnames(pset)

    zz <- sapply(clades, function(g) {

        sp = match(pset[[g]], fulltree$tip.label)

        tips = ape::extract.clade(fulltree, ape::getMRCA(fulltree, sp))$tip.label

        tips = tips[ tips %in% phy$tip.label ]

        if (length(tips) > mintips)
            return (ape::getMRCA(phy, tips))
        else
            return (NA_integer_)
    })
    if (clip.NA){
    		return(zz[!is.na(zz)])
    }else{
    		return(zz)
    }
    
}


# very inefficient as it redoes the full rootwards traversal for each node
getDescendantMatrix <- function(phy){
	
	rootnode <- length(phy$tip.label) + 1
	ntip <- length(phy$tip.label)
	dmat <- matrix(0, nrow=length(phy$tip.label), ncol= 2*ntip - 1)
	
	for (ii in 1:ntip){
		nd <- ii 
		while (nd != rootnode){
			nd <- phy$edge[,1][phy$edge[,2] == nd]
			dmat[ii, nd] <- 1
		}
	}
	
	diag(dmat[,1:ntip]) <- 1
	dmat[,rootnode] <- 1
	return(dmat)
}

 


























