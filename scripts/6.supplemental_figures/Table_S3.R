rm(list = ls())
setwd("~/Dropbox (Personal)/research_projects-ACTIVE/squamate_tree/")

c = readRDS("./diversification/pgls-model-fits_v5.rds")
b = readRDS("./diversification/pgls-model-fits_bamm_v5.rds")
c1 = readRDS("./diversification/pgls-model-fits_v5_with_brlen.rds")
a = readRDS("./diversification/pgls-model-fits_alethniophidia_v5.rds")

newnames <- c(M_TROPHIC1 = 'diet', M_TROPHIC2 = 'foraging mode', 
              M_TROPHIC3 = 'chemosensory', M_PARITY = 'parity mode', 
              M_KINESIS = 'skull kinesis', M_BIOCLIM = 'climate',
              M_MORPH = 'body size/shape', M_MERISTIC1 = 'digit/limb count',
              M_MERISTIC2 = 'vertebral count', M_SKULL = 'skull shape', 
              M_GEO = 'geography', M_BIOGEO = 'biogeog. theater', 
              M_CLADE = 'clade membership')

sub = sapply(c, function(x) paste(names(x$coefficients), collapse = ", "))
r_c <- round(sapply(c, function(x) summary(x)$adj.r.squared), 3)
ss_c = sapply(c, function(x) length(x$residuals))
r_b <- round(sapply(b, function(x) summary(x)$adj.r.squared), 3)
r_c1 <- round(sapply(c1, function(x) summary(x)$adj.r.squared), 3)
ss_a = sapply(a, function(x) length(x$residuals))
r_a <- round(sapply(a, function(x) summary(x)$adj.r.squared), 3)

# type, sample size, r-sq for bamm & clads
res = data.frame(cbind(ss_c, sub, r_c, r_b, r_c1))
res = cbind( var = newnames[ rownames(res) ], res)
names(res) = c("predictor", "n", "components",
               "r-squared CLaDS", "r-squared BAMM", 
               "r-squared CLaDS, br. len.")
res = res[order(res$`r-squared CLaDS`, decreasing = T), ]
write.csv(res, "./manuscript/figures/Table_S3.csv",
          row.names = F)

ivan = data.frame( var = res$predictor, 
                   r2 = res$`r-squared CLaDS`)
write.csv(ivan, 
          sprintf("ivan-figure-work/data/speciation_rate_phylo_r2.%s.csv", Sys.Date()))

# individual p-values of all coefficients
coefd = do.call("rbind", sapply(c, function(x) summary(x)$coefficients))


res = data.frame(cbind(ss_a, sub, r_a))
res = cbind( var = newnames[ rownames(res) ], res)
names(res) = c("predictor", "n", "components",
               "r-squared CLaDS")
res = res[order(res$`r-squared CLaDS`, decreasing = T), ]
write.csv(res, "~/Desktop/alethino.csv")

coefa = do.call("rbind", sapply(a, function(x) summary(x)$coefficients))
