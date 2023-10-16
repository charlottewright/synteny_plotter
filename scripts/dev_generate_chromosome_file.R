suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))

source('scripts/helper_functions.R')
source('scripts/interactive_args.R')

minimum_buscos = 5


R_df <- read_buscos(args$busco1, 'R')
Q_df <- read_buscos(args$busco2, 'Q') # was Agrochola_circellaris.tsv
alignments <- merge(Q_df, R_df, by='busco')
alignments <- alignments %>% group_by(chrR) %>% filter(n() > minimum_buscos) %>% ungroup() # filtering out those with very few BUSCO genes

R_chromosomes <- read.table(args$chrom1, sep = '\t', header = TRUE)
Q_chromosomes <- read.table(args$chrom2, sep = '\t', header = TRUE)
R_chromosomes <- R_chromosomes %>% arrange(order)
Q_chromosomes <- Q_chromosomes %>% arrange(order)
row.names(R_chromosomes) <- R_chromosomes$chr

# apply any filters
R_chromosomes <- R_chromosomes %>% filter(chr %in% alignments$chrR) # removing those with no / fwe BUSCO genes fromt he reference chr table
Q_chromosomes <- Q_chromosomes %>% filter(chr %in% alignments$chrQ) # removing those with no / fwe BUSCO genes fromt he quert chr table

get_most_frequent_ref_hit <- function(chr){
    hit_table <- alignments %>% filter(chrQ == chr) %>% select(chrR) %>% table
    names(hit_table)[which.max(hit_table)]
}
chromosomal_correspondences <- sapply(unique(alignments$chrQ), get_most_frequent_ref_hit)
chromosomal_correspondences <- data.frame(chrQ = names(chromosomal_correspondences), chrR = as.character(chromosomal_correspondences))
chromosomal_correspondences$R_order <- R_chromosomes[chromosomal_correspondences$chrR, 'order']
chromosomal_correspondences <- chromosomal_correspondences %>% arrange(R_order)  

test_invert <- function(chr_row){
    chr_aln <- alignments %>% filter(chrQ == as.character(chr_row['chrQ'])) %>%  filter(chrR == as.character(chr_row['chrR']))
    cortest <- cor.test(unlist(chr_aln[, 'Qstart']), unlist(chr_aln[, 'Rstart']))
    if (cortest$estimate < 0){
        return (TRUE)
    }
    return (FALSE)    
}

chromosomal_correspondences$order <- 1:nrow(chromosomal_correspondences)
chromosomal_correspondences$invert <- apply(chromosomal_correspondences, 1, test_invert)


chromosomal_correspondences <- chromosomal_correspondences %>% select(c(chrQ, order, invert)) 
colnames(chromosomal_correspondences)[1] <- 'chr'

write.table(chromosomal_correspondences, 'test_data/Vanessa_cardui_info_generated.tsv', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)


