cutree_bhc = function(res, alpha=0.5) {
    n = nrow(res$merge)
    if (res$height[n] > res$height[n-1])
        root = res$merge[n,]
    else
        root = n
    tips = function(node) {
        ans = c()
        foo = function(node) {
            if (node > 0) {
                lf = res$merge[node, 1]
                rt = res$merge[node, 2]
                foo(lf)
                foo(rt)
            } else {
                ans <<- c(ans, node)
            }
        }
        foo(node)
        -ans
    }

    bar = function(node) {
        if (node > 0) {
            if (exp(res$height[node]) < alpha) {
                lf = res$merge[node, 1]
                rt = res$merge[node, 2]
                bar(lf)
                bar(rt)
            } else {
                g <<- g + 1
                #ans[[g]] <<- tips(node)
                ans[tips(node)] <<- g
            }
        } else {
            g <<- g + 1
            #ans[[g]] <<- -node
            ans[tips(node)] <<- g
        }
    }
    #ans = list()
    ans = integer(n+1)
    g = 0
    for (node in root)
        bar(node)
    ans
}

DATAFILE = "../diet-matrix.csv"
TREEFILE = "../../../1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre"
NAMEMAP = "../dietname-to-treename.csv"

x = data.matrix(read.csv(DATAFILE, row.names=1))
fit = readRDS("bhc.rds")

diet.states = t(sapply(split(rownames(x), cutree_bhc(fit, 0.95)), 
    function(d)
    {
        colSums(x[d, , drop=FALSE])+1
    }
))

diet.states = sweep(diet.states, 1, rowSums(diet.states), "/")

y = structure(diet.states[cutree_bhc(fit, 0.95),], dimnames=list(rownames(x),colnames(diet.states)))

names_map = read.csv(NAMEMAP, header=TRUE)

diet_names = rownames(y)
tree_names = names_map[match(diet_names, names_map[,1]), 2]
rownames(y) = tree_names
diet_names = rownames(y)

phy = ape::read.tree(TREEFILE)

common_spp = intersect(phy$tip.label, diet_names)

y = y[diet_names %in% common_spp, ]

phy = ape::keep.tip(phy, rownames(y))

y = y[phy$tip.label, ]


write.csv(y, file="../diet-proportions-nonphylo.csv")
