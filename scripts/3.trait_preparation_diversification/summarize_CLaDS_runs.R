library(ape)

clads = list.files("/media/babs/brains2/best-full-trees/", pattern = "Rdata",
                    full.names = TRUE)
load(clads[1])
tips = sort(CladsOutput$tree$tip.label)

stats = data.frame(gelm = rep(NA, length(clads)), length = rep(NA, length(clads)), run = clads)
res = matrix(NA, nrow = length(tips),
             ncol = length(clads), dimnames = list(tips))
for (i in 1:length(clads)) {
  load(clads[i])
  d = CladsOutput
  cladsrates = d$lambdatip_map
  names(cladsrates) = d$tree$tip.label
  res[ ,i] = cladsrates[tips]
  stats[i, "gelm"] = d$gelm[2]
  stats[i, "length"] = length(d$rtt_chains[[1]])
  print(i)
}

write.csv(stats, "best_tree_100_imputations.convergence_stats.csv", row.names = F)
res2 = res[grep("xzx", dimnames(res)[[1]], invert = T), ]
write.table(res2, "best_tree_100_imputations.CLaDS_rates.txt", 
             row.names = T, col.names = F)
# write.table(res2, "~/Dropbox (Personal)/squamate_tree/diversification/pseudoposterior_100_imputations.CLaDS_# rates.txt", 
#            row.names = T, col.names = F)
