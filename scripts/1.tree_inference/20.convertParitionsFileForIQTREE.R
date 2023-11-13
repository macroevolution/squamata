# Convert partition file to be compatible with IQTREE

# need to change naming to: DNA, mtDNA_c1 (DNA, single name)
# replace commas between seq ranges with spaces

setwd('~/Dropbox/Oz_Crown_Ages/finalAlignments_June2020/renamed/')

inputFile <- 'concatenated.partitions'
outputFile <- 'concatenated.partitions.iqtree'

qq <- scan(inputFile, what = 'character', sep = '\n')
qq

qq <- strsplit(qq, ' = ')

# collapse names
partNames <- sapply(qq, function(x) x[1])
partNames <- gsub('GTR\\+G\\,\\s', '', partNames)
partNames <- gsub(',\\s', '_', partNames)
partNames <- paste0('DNA, ', partNames)

# replace commas with spaces for seq ranges
seqRanges <- sapply(qq, function(x) x[2])
seqRanges <- gsub(',\\s?', ' ', seqRanges)

newPartitions <- paste(partNames, seqRanges, sep = ' = ')

write(newPartitions, outputFile)