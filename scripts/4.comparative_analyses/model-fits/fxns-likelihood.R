 
##' Compute the marginal log likelihood of a dataset under Brownian motion
##'
##' @param x A matrix of trait data with rownames corresponding
##' to tip labels.
##' @param phy A phylogeny of class "phylo". The phylogeny must be
##' in "pruningwise" order.
##' @param rate A numeric vector with the relative rate of trait
##' evolution governing each branch. rate[k] is the relative rate
##' of evolution on the branch leading to the node with index k.
##' @param ignore.cov Logical. Should evolutionary covariances be
##' set to zero?
##' @param auto.drop Logical. If evolutionary covariances are not
##' ignored should collinear variables that render the covariance matrix
##' singular be dropped?
downpass = function(x, phy, rate, ignore.cov, auto.drop)
{
    stopifnot(inherits(phy, "phylo"))
    if (ape::Ntip(phy) == 1)
        return (structure(0, rate=0,np.cov=0,np.mu=0,auto.dropped=integer(0)))
    if (ape::Ntip(phy) > 2)
        stopifnot(attr(phy, "order") == "pruningwise")
    stopifnot(is.matrix(x))
    stopifnot(!is.null(rownames(x)))

    edge = phy$edge
    ntip = ape::Ntip(phy)
    nnode = ape::Nnode(phy) + ntip
    x = x[phy$tip.label, , drop=FALSE]
    v = numeric(nnode)
    v[edge[,2]] = phy$edge.length
    v = v * rate

    d = ncol(x)
    m = rbind(x, matrix(0, nnode-ntip, d))
    u = matrix(0, nnode, d)
    vu = numeric(nnode)

    for (i in seq.int(1, nnode-1, 2))
    {
        node = edge[i, 1]
        lf = edge[i, 2]
        rt = edge[i+1, 2]
        v[node] = v[node] + v[lf]*v[rt] / (v[lf] + v[rt])
        m[node,] = (v[rt]*m[lf,] + v[lf]*m[rt,]) / (v[lf] + v[rt])
        vu[node] = v[lf] + v[rt]
        u[node,] = (m[lf,] - m[rt,]) / sqrt(vu[node])
    }

    u = u[-(1:ntip), , drop=FALSE]
    vu = vu[-(1:ntip)]

    # be aware: for high-dimensional data the sample covariance matrix of
    # the contrasts (MLE of brownian rate) can sometimes be non positive 
    # definite due to high colinearity of some variables. i've seen this
    # crop up with both the skull data and the diet data. this results in 
    # non-finite likelihoods. not a problem when ignore.cov=TRUE, and this 
    # is the reason it is the default.
    sigma = (t(u) %*% u) / (ntip - 1)
    auto.dropped = integer(0)

    if (ignore.cov)
    {
        sigma = diag(diag(sigma), d, d)
        # number of variance parameters
        np.cov = d
        # number of mean parameters
        np.mu = d
    }
    else
    {
        np.cov = d*(d+1)/2
        np.mu = d
        if (auto.drop)
        {
            if (inherits(try(chol(sigma), silent=TRUE), "try-error"))
            {
                i = rep(FALSE, nrow(sigma))
                i[1] = TRUE
                for (j in 2:nrow(sigma))
                {
                    i[j] = TRUE
                    if (inherits(try(chol(sigma[i,i,drop=FALSE]), silent=TRUE), "try-error"))
                        i[j] = FALSE
                }
                auto.dropped = which(!i)
                msg = paste0(c(
                    "sample covariance is not positive definite,"
                    , "dropping collinear variables:\n\t%s"), collapse=" ")
                msg = sprintf(msg, paste0(auto.dropped, collapse=", "))
                warning(msg)
                sigma = sigma[i, i, drop=FALSE]
                u = u[, i, drop=FALSE]
                d = sum(i)
                np.cov = d*(d+1)/2
                np.mu = d
            }
        }
    }

    l = sum(mvtnorm::dmvnorm(u, sigma=sigma, log=TRUE)) - d*sum(log(sqrt(vu)))
    r = sum(diag(sigma))

    structure(l, rate=r, np.cov=np.cov, np.mu=np.mu, auto.dropped=auto.dropped)
}

##' Optimize the marginal log likelihood of a multi-rate or
##' constant-rate Brownian motion model
##'
##' @param x A matrix of trait data with rownames corresponding
##' to tip labels.
##' @param phy A phylogeny of class "phylo".
##' @param shift An integer vector of node indices corresponding to
##' locations of rate shifts. If missing or of zero-length then
##' just fit a constant rate model.
##' @param ignore.cov Logical. Should evolutionary covariances be
##' set to zero?
##' @param auto.drop Logical. If evolutionary covariances are not
##' ignored should collinear variables that render the covariance matrix
##' singular be dropped?
##' @return A list with the following components:
##'  $rate The rate of each shift. If the trait is multivariate this is the
##'        sum of the diagonal rate matrix.
##'  $logL The marginal log likelihood. NA on failure.
##'  $aic The Akaike Information Criterion (uncorrected). NA on failure.
##'  $npar The number of estimated parameters
##'  $converged A logical indicating if the optimization converged.
##'  $auto.dropped A vector of variable indices that were dropped to make
##'                the sample covariance matrix of contrasts non-singular.
##'  $phy.scaled The phylogeny with branch lengths equal to their rate
##'              of evolution.
fit_shift = function(x, phy, shift, ignore.cov=TRUE, auto.drop=FALSE, ...)
{
    stopifnot(inherits(phy, "phylo"))
    stopifnot(attr(phy, "order") == "cladewise")
    if (missing(shift) || !length(shift))
    {
        phy = ape::reorder.phylo(phy, "pruningwise")
        rate = rep(1, ape::Ntip(phy) + ape::Nnode(phy))
        ans = downpass(x, phy, rate, ignore.cov, auto.drop)
        
        if (is.nan(ans) || !is.finite(ans))
        {
            return (errorCondition(
                "sample covariance matrix for contrasts is singular", 
                class="singular-covariance-error"))
        }

        npar = attr(ans, "np.cov") + attr(ans, "np.mu")
        ll = ans
        aic = -2*ll + 2*npar
        
        phy.scaled = ape::reorder.phylo(phy, "cladewise")
        phy.scaled$edge.length = rep(ans[2], nrow(phy.scaled$edge))
        return(list(
            rate=setNames(attr(ans, "rate"), "background")
            , logL=ll
            , aic=aic
            , npar=npar
            , converged=TRUE
            , auto.dropped=attr(ans, "auto.dropped")
            , phy.scaled=phy.scaled
            , phy.original=phy
        ))
    }

    preorder_edge = phy$edge
    brlen = numeric(ape::Ntip(phy) + ape::Nnode(phy))
    brlen[preorder_edge[,2]] = phy$edge.length
    phy = ape::reorder.phylo(phy, "pruningwise")
    
    treelen = sum(brlen)

    lik = function(pars)
    {
        if (any(pars < 0)) # L-BFGS-B doesn't always respect bounds.
            return (-1e9)
        ratemultipl = rep(1, ape::Ntip(phy)+ape::Nnode(phy))
        ratemultipl[shift] = pars
        rltvrate = rep(1, ape::Ntip(phy)+ape::Nnode(phy))
        for (i in 1:nrow(preorder_edge))
        {
            p = preorder_edge[i, 1]
            nd = preorder_edge[i, 2]
            rltvrate[nd] = ratemultipl[nd] * rltvrate[p]
        }
        treescal = sum(brlen * rltvrate)
        norm = treelen / treescal
        downpass(x, phy, norm * rltvrate, ignore.cov, auto.drop)
    }

    l = log(1e-3)
    u = log(1e3)

    init = exp(rnorm(length(shift), mean=0, sd=0.25))

    # if ignore.cov=FALSE, auto.drop=FALSE and the covariance matrix is singular
    # optim will fail immediately, so we do an early check here
    if (!ignore.cov && !auto.drop)
    {
        ans = lik(init)
        if (is.nan(ans) || !is.finite(ans))
        {
            return (errorCondition(
                "sample covariance matrix for contrasts is singular", 
                class="singular-covariance-error"))
        }
    }

    fit = try(optim(init, lik, method="L-BFGS-B", lower=l, upper=u, 
        control=list(fnscale=-1, ...)), silent=TRUE)

    if (inherits(fit, "try-error"))
    {
        return (errorCondition(
            attr(fit, "condition")$message, class="optim-error"))
    }

    if (is.nan(fit$value) || !is.finite(fit$value))
    {
        return (errorCondition(
            "non-finite likelihood returned by optim", 
            class="optim-likelihood-error"))
    }
   
    ratemultipl = rep(1, ape::Ntip(phy)+ape::Nnode(phy))
    ratemultipl[shift] = fit$par
    rltvrate = rep(1, ape::Ntip(phy)+ape::Nnode(phy))
    for (i in 1:nrow(preorder_edge))
    {
        p = preorder_edge[i, 1]
        nd = preorder_edge[i, 2]
        rltvrate[nd] = ratemultipl[nd] * rltvrate[p]
    }
    treescal = sum(brlen * rltvrate)
    norm = treelen / treescal
    ans = downpass(x, phy, norm * rltvrate, ignore.cov, auto.drop)

    rate = attr(ans, "rate")

    stopifnot(all.equal(sum(norm*rltvrate*brlen), treelen))

    # parameter counting:
    #
    #   1 rate multiplier for each shift + the mean vector + the covariances
    #
    # technically we don't estimate the mean vector but we count it for all
    # models for ease of comparison
    npar = length(shift) + attr(ans, "np.cov") + attr(ans, "np.mu")

    phy.scaled = ape::reorder.phylo(phy, "cladewise")
    phy.scaled$edge.length = norm*rltvrate[phy.scaled$edge[,2]]*rate

    ll = fit$value
    aic = -2*ll + 2*npar
 
    structure(list(
        rate=setNames(c(norm*rate, norm*rltvrate[shift]*rate), c("background", shift))
        , logL=ll
        , aic=aic
        , npar=npar
        , converged=fit$convergence == 0L
        , msg=fit$message
        , auto.dropped=attr(ans, "auto.dropped")
        , phy.scaled=phy.scaled
        , phy.original=phy
    ), class="fit-shift")
}

##' Optimize the marginal log likelihood of a multi-rate censored
##' Brownian motion model
##'
##' @param x A matrix of trait data with rownames corresponding
##' to tip labels.
##' @param phy A phylogeny of class "phylo".
##' @param shift An integer vector of node indices corresponding to
##' locations of shifts. These branches will be "censored" or cut
##' such that the tree is partitioned into independent subtrees having
##' their own rate of evolution and phylogenetic mean.
##' @param ignore.cov Logical. Should evolutionary covariances be
##' set to zero?
##' @param auto.drop Logical. If evolutionary covariances are not
##' ignored should collinear variables that render the covariance matrix
##' singular be dropped?
##' @return A list with the following components:
##'  $rate The rate of each subtree. If the trait is multivariate this is the
##'        sum of the diagonal rate matrix.
##'  $logL The marginal log likelihood. NA on failure.
##'  $aic The Akaike Information Criterion (uncorrected). NA on failure.
##'  $npar The number of estimated parameters
##'  $converged A logical indicating if the optimization converged.
##'  $auto.dropped A list. E.g., $autodropped[[i]] is a vector variable indices
##'                that were dropped to make the sample covariance matrix of
##'                contrasts non-singular for the i-th subtree in the partition. 
##'  $phy.scaled A list of the subtrees with branch lengths equal to their rate
##'              of evolution.
fit_censored = function(x, phy, shift, ignore.cov=TRUE, auto.drop=FALSE)
{
    stopifnot(inherits(phy, "phylo"))
    stopifnot(attr(phy, "order") == "cladewise")
    stopifnot(!missing(shift) && length(shift))
    bt = sort(ape::branching.times(phy), decreasing=TRUE)
    badshift = rep(FALSE, length(shift))
    for (i in 1:length(shift))
    {
        s = shift[i]
        if (all(phy$edge[which(s == phy$edge[, 1]), 2] %in% shift))
        {
            msg = sprintf("removing shift %d: both children in shift set", s)
            warning(msg)
            badshift[i] = TRUE
        }
    }
    if (all(phy$edge[which((ape::Ntip(phy)+1) == phy$edge[,1]), 2] %in% shift))
    {
        s = phy$edge[match(ape::Ntip(phy)+1, phy$edge[,1]), 2]
        i = match(s, shift)
        msg = sprintf("removing shift %d: both children of root in shift set", s)
        warning(msg)
        badshift[i] = TRUE
    }
    shift = shift[!badshift]
    ord = order(match(as.character(shift), names(bt)))
    shift = shift[ord] # preorder
    shift = rev(shift) # postorder
    keep = vector("list", length(shift)+1)
    phy.scaled = structure(vector("list", length(shift)+1), class="multiPhylo")
    names(phy.scaled) = c(shift, "background")
    for (i in seq_along(shift))
    {
        keep[[i]] = setdiff(
            ape::extract.clade(phy, shift[i])$tip.label, unlist(keep[1:i]))
    }
    keep[[length(shift)+1]] = setdiff(phy$tip.label, unlist(keep))

    r = setNames(numeric(length(shift)+1), c(shift, "background"))
    ll = 0
    # parameter counting:
    #
    #   for each shift and for the background process we have a mean vector
    #   and a covariance matrix, so just sum these up into np.cov and np.mu
    np.cov = 0
    np.mu = 0
    auto.dropped = vector("list", length(keep))
    for (i in seq_along(keep))
    {
        phy2 = ape::keep.tip(phy, keep[[i]])
        rate = rep(1, ape::Ntip(phy2)+ape::Nnode(phy2))
        phy2 = ape::reorder.phylo(phy2, "pruningwise")
        ans = downpass(x[phy2$tip.label, , drop=FALSE], phy2, rate, ignore.cov, 
            auto.drop)
        if (is.nan(ans) || !is.finite(ans))
        {
            return (errorCondition(
                "sample covariance matrix for contrasts is singular", 
                class="singular-covariance-error"))
        }
        ll = ll + ans
        r[i] = attr(ans, "rate")
        np.cov = np.cov + attr(ans, "np.cov")
        np.mu = np.mu + attr(ans, "np.mu")
        phy2 = ape::reorder.phylo(phy2, "cladewise")
        phy2$edge.length = rep(r[i], nrow(phy2$edge))
        phy.scaled[[i]] = phy2
        auto.dropped[[i]] = attr(ans, "auto.dropped")
    }

    npar = np.mu + np.cov
    
    aic = -2*ll + 2*npar
    
    structure(list(
        rate=rev(r)
        , logL=ll
        , aic=aic
        , npar=npar
        , converged=TRUE
        , auto.dropped=rev(auto.dropped)
        , phy.scaled=rev(phy.scaled)
        , phy.original = phy
    ), class="fit-censored")
}
