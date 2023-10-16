suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))

source('scripts/helper_functions.R')
source('scripts/interactive_args.R')

minimum_buscos = 5


R_df <- read_buscos(args$busco1, 'R')
Q_df <- read_buscos(args$busco2, 'Q') # was Agrochola_circellaris.tsv

R_chromosomes <- read.table(args$chrom1, sep = '\t', header = TRUE)
R_chromosomes <- R_chromosomes %>% arrange(order)

# apply any filters
R_chromosomes <- R_chromosomes %>% filter(chr %in% alignments$chrR) # removing those with no / fwe BUSCO genes fromt he reference chr table

######## till here it's what is in the origin thingy

generate_auto_query_order(R_chromosomes, R_df, Q_df)


# write.table(chromosomal_correspondences, 'test_data/Vanessa_cardui_info_generated.tsv', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)


