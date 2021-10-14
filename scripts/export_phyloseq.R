suppressMessages(library(optparse))
suppressMessages(library(phyloseq))
option_list = list(
    make_option(c("-b", "--biom"),
    type="character"),
    make_option(c("-t", "--tree"),
    type="character"),
    make_option(c("-o", "--outp"),
    type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#print(opt$biom)
ps <- phyloseq::import_biom(opt$biom, treefilename=opt$tree)


saveRDS(ps, opt$outp)