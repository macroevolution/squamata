
setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

library(prtree) #remotes::install_github('blueraleigh/prtree')

x <- read.csv('./data/3.trait_preparation_diversification/dietModeling/diet-proportions-phylo.csv', row.names = 1)

z = data.matrix(x)

y = sweep(log(z), 1, rowMeans(log(z)))

ans = prtree(y, nrow(y), 200, .002, 100)

g = factor(c(rep(2, 1021), rep(4, 1314-1021)))

ans = plot(ans, cex=125, by=g, piebg=c(2,4), main="Trophic principal tree")
legend("topright", legend=c("snakes", "lizards"), pt.cex=1, pch=21, pt.bg=c(2,4), bty="n",cex=0.6)


####

z2=data.frame(lapply(
    list(
        worms_mollusks=c(2,21),
        decapods=12,
        arachnids_centipedes=c(9,27,29),
        insects=c(3,4,7,10,16,19,22,23,30), # 19 is misc_inverts
        plants=24,
        mammals_birds=c(5:6,18),
        amphibians=c(1,8,15,26),
        squamates=c(17,25,28),
        fishes=c(13,14),
        other_verts=c(11,20,31) # crocs, misc_verts, turtles
    ), function(idx)
    {
        rowSums(z[,idx,drop=FALSE])
    }
))

z2 = round(data.matrix(z2),2)
# assign the 'other_verts' category to 'squamates'. the only species with non-neglible
# values for this are a few lizards, which are probably eating other lizards
z2[,8] = z2[,8] + z2[,10]
z2 = z2[,-10]
# ignore minor stuff for now
z2[z2 < 0.05] = 0

pie = matrix(NA, nrow(ans$r), ncol(z2))

tmp = tapply(1:nrow(ans$r), apply(ans$r, 1, which.max), function(idx) {
    colMeans(z2[idx,,drop=FALSE])
})

for (i in 1:length(tmp))
    pie[as.integer(names(tmp)[i]),] = tmp[[i]]

colv=c(
"#8c510a", # worms_mollusks
"#3D3C3D", # decapods
"orange",  # arachnids_centipedes
"yellow",  # insects
"#f5f5f5", # plants
"#c7eae5", # mammals_birds
"#80A063", # amphibians
"#35978f", # squamates
"#014B5E"  # fishes
)
ans = plot(ans, cex=125, by=pie, piebg=colv)
legend("topright",legend=colnames(z2),pch=21,pt.bg=colv,bty="n", ncol=1,pt.cex=1,cex=0.6)


