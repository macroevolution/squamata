
Tree inference README
=====================

concatenatedAlignment: directory containing the concatenated alignment and partition files
fullTree_constrainedInference_mergeAll: directory containing output files from the tree inference stage where separate family-level trees were merged together into a single tree. This is the tree inference constrained by the genomic backbone.
fullTree_unconstrainedInference_mergeAll: directory containing output files from the tree inference stage where separate family-level trees were merged together into a single tree. This is the tree inference that did not include the genomic backbone constraint.
mainTrees: final constrained and unconstrained molecular trees, and final constrained ultrametric tree.
pseudoposterior: set of 100 pseudo-posterior trees that represent several sources of uncertainty.

In mainTrees directory:

squamates_Title_Science2024_ultrametric_constrained.tre: 6885-tip time-calibrated genomic-constrained phylogeny (main tree in paper)
squamates_Title_Science2024_molecular_constrained.tre: 6885-tip molecular genomic-constrained non-time-calibrated phylogeny
squamates_Title_Science2024_unconstrained.tre: 6885-tip molecular non-time-calibrated phylogeny, not constrained by genomic backbone