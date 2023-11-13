using DelimitedFiles
using PANDA

my_tree = load_tree("best_ultrametric_fulltree_ddBD_revision.tre")
miss = readdlm("family_samplingFracs.tip_order.txt")
miss2 = vec(miss)
out = infer_ClaDS(my_tree, f = miss2)
save_ClaDS_in_R(out, "best_ultrametric_fulltree_ddBD_revision.CLaDS_samplingfrac.Rdata")