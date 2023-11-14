library(bhc) # github.com/blueraleigh/bhc

DATAFILE = "../diet_matrix.csv"

x = data.matrix(read.csv(DATAFILE, row.names=1))

fit = bhc.multinomial(x)

saveRDS(fit, "bhc.rds")
