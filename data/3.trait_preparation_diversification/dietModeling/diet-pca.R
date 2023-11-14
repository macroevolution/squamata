# diet PCA calculation

for (TYPE in c("phylo", "nonphylo")) {

    DATAFILE = sprintf("diet-proportions-%s.csv", TYPE)
    x = read.csv(DATAFILE, row.names=1)
    z = data.matrix(x)
    z = sweep(log(z), 1, rowMeans(log(z)))
    pca = prcomp(z)
    
}