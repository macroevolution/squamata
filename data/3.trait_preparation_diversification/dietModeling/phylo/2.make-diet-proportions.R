library(macroevolution)

X = vector("list", 4)
w = numeric(4)
names(w) = names(X) = c(20,50,100,1000)

for (i in 1:4)
{
    mcmc = read.rcm.dmm(sprintf("mcmc-%s.out", names(X)[i]))

    D = xtabs(mcmc$dataset[,3] ~ mcmc$dataset[,1] + mcmc$dataset[,2])
    rownames(D) = tiplabels(mcmc$phy)[as.integer(rownames(D))+1L]
    colnames(D) = names(mcmc$dirichlet.prior)[as.integer(colnames(D))+1L]

    end = nrow(mcmc$pars)
    start = floor(end / 1.95)
    n = end-start+1

    w[i] = mean(rowSums(mcmc$pars[start:end,1:2]))

    post.states = mcmc$stateid[start:end,]

    X[[i]] = matrix(0, Ntip(mcmc$phy), 31, 
        dimnames=list(tiplabels(mcmc$phy), names(mcmc$dirichlet.prior)))

    invisible(apply(post.states, 1, function(p) {
        up = unique(p)
        op = order(up)
        p = setNames(match(p, up[op]), names(p))
        x = t(sapply(lapply(split(D, p), matrix, ncol=31), colSums))+1
        x = sweep(x, 1, rowSums(x), "/")
        X[[i]] <<- X[[i]] + x[p,] / n
    }))
}

w = exp(w - max(w)) / sum(exp(w - max(w)))

Y = matrix(0, Ntip(mcmc$phy), 31, 
        dimnames=list(tiplabels(mcmc$phy), names(mcmc$dirichlet.prior)))

for (i in 1:4)
    Y = Y + w[i] * X[[i]]

Y = Y[,sort(colnames(Y))]

write.csv(Y, "../diet-proportions-phylo.csv")



