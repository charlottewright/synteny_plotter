suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))

parser <- ArgumentParser()
parser$add_argument("-busco1",
                    help="Busco \"full_table.tsv\" of the reference (top) species")
parser$add_argument("-busco2", 
                    help="Busco \"full_table.tsv\" of the query (bottom) species")
parser$add_argument("-o", "--output_prefix", default = "synteny_plot",
                    help="Name pattern for the output")
parser$add_argument("-f", "-filter", type = "integer", default=5,
                    help="The minimal number of BUSCOs on a chromosome to include")

args <- parser$parse_args()
source('scripts/helper_functions.R')
#source('scripts/interactive_args.R')

minimum_buscos = args$f

print(paste('Running with a minimum number of markers of:', minimum_buscos))
# uncomment if running manually
#args$busco1 <-'test_data/Melitaea_cinxia.tsv'
#args$busco2 <-'test_data/Vanessa_cardui.tsv'
#args$output_prefix <- 'test_two_species'

print(args$busco1)
R_df <- read_buscos(args$busco1, 'R')
Q_df <- read_buscos(args$busco2, 'Q') # was Agrochola_circellaris.tsv
R_chromosome_chr <- unique(R_df$chr)
print(R_chromosome_chr)
R_chromosome_order <- seq(1, length(R_chromosome_chr)) # make a list of just 1 to n based on number of chr
R_chromosomes <- data.frame(chr = R_chromosome_chr, order = R_chromosome_order)
R_chromosomes$invert <- 'F' # set all to false at the start
#R_chromosomes <- read.table(args$chrom1, sep = '\t', header = TRUE)
R_chromosomes <- R_chromosomes %>% arrange(order)

# apply any filters
alignments <- merge(Q_df, R_df, by='busco')
# apply any filters
alignments <- alignments %>% group_by(chrR) %>% filter(n() > minimum_buscos) %>% ungroup()
alignments <- alignments %>% group_by(chrQ) %>% filter(n() > minimum_buscos) %>% ungroup()
R_chromosomes <- R_chromosomes %>% filter(chr %in% alignments$chrR) # removing those with no / few BUSCO genes from the reference chr table
print(head(Q_df))
Q_df <- Q_df %>% group_by(chrQ) %>% filter(n() > minimum_buscos) %>% ungroup()
#Q_chromosomes <- Q_chromosomes %>% filter(chr %in% alignments$chrQ) # removing those with no / few BUSCO genes from the query chr table

######## till here it's what is in the origin thingy

chromosomal_correspondences <- generate_auto_query_order(R_chromosomes, R_df, Q_df)


write.table(chromosomal_correspondences, paste0(args$output_prefix, '.auto_generated.tsv'), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
