using PANDA
my_tree = load_tree("/media/babs/brains2/best-full-trees/squams-full-imputed-xx1.tre")
out = infer_ClaDS(my_tree)
save_ClaDS_in_R(out, "/media/babs/brains2/best-full-trees/squams-full-imputed-xx1.Rdata")
