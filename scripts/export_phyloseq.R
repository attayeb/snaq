suppressMessages(library(optparse))
suppressMessages(library(phyloseq))
suppressMessages(library(dplyr))
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
df <- tidyr::separate(as.data.frame(tax_table(ps)), "Rank1", c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Rank8", "Rank9", "Rank10"), sep=";")
not_all_na <- function(x) any(!is.na(x))

tax_table(ps) <- apply(df%>%select(where(not_all_na)), 2, trimws)

saveRDS(ps, opt$outp)