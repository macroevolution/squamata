# diet breadth calculation

for (TYPE in c("phylo", "nonphylo")) {

    DATAFILE = sprintf("diet-proportions-%s.csv", TYPE)
    x = data.matrix(read.csv(DATAFILE, row.names=1))
    diet.breadth = apply(x, 1, function(p) 1 / sum(p*p))

}