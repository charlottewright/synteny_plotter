# cutoff of block in terms of same strand and stuff, and distance until next ortholog

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))

parser <- ArgumentParser()
parser$add_argument("-busco1",
                    help="Busco \"full_table.tsv\" of the top species")
parser$add_argument("-busco2", 
                    help="Busco \"full_table.tsv\" of the bottom species")
parser$add_argument("-chrom1", 
                    help="A .tsv table with chromosomes of the species 1 (4 columns expected chromosome, length, order, direction)")
parser$add_argument("-chrom2", 
                    help="A .tsv table with chromosomes of the species 2 (4 columns expected chromosome, length, order, direction)")
parser$add_argument("-o", "--output_prefix", default = "synteny_plot",
                    help="Name pattern for the output")
parser$add_argument("-g", "-gap", type = "integer",
                    help="Gap between two chromosomal sets")
parser$add_argument("-f", "-filter", type = "integer",
                    help="The minimal number of BUSCOs on a chromosome to include")
parser$add_argument("-alpha", type = "integer", default = 0,
                    help="Set transparency to colours [%]")

args <- parser$parse_args()

### specify arguments ###

busco1 <- args$busco1
busco2 <- args$busco2
# gap <- args$gap
chrom1 <- args$chrom1
chrom2 <- args$chrom2
# minimum_buscos <- args$filter
# output_prefix <- args$output

source('scripts/helper_functions.R') # import functions

### specify parameters
gap=6
show_outline = TRUE
chr_offset = 20000000 # TODO make this automatically a prop of chr length
minimum_buscos = 5

### read in data ###

# TODO: allow for algs in the arguments
# algs <- read.csv(args$alg_file, sep='\t', header=FALSE)[,c(1,3)]
# colnames(algs) <- c('busco', 'alg')
# alignments <- merge(alignments, algs, by='busco')

R_df <- read_buscos(args$busco1, 'R')
Q_df <- read_buscos(args$busco2, 'Q') # was Agrochola_circellaris.tsv
R_chromosomes <- read.table(args$chrom1, sep = '\t', header = TRUE)
Q_chromosomes <- read.table(args$chrom2, sep = '\t', header = TRUE)

R_chromosomes <- R_chromosomes %>% arrange(order)
Q_chromosomes <- Q_chromosomes %>% arrange(order)
alignments <- merge(Q_df, R_df, by='busco')

# apply any filters
alignments <- alignments %>% group_by(chrR) %>% filter(n() > minimum_buscos) %>% ungroup()
R_chromosomes <- R_chromosomes %>% filter(chr %in% alignments$chrR)  
Q_chromosomes <- Q_chromosomes %>% filter(chr %in% alignments$chrQ) 
chr_order_R <- R_chromosomes$chr # extract order of chr
chr_order_Q <- Q_chromosomes$chr #extract order of chr

alignments$alg <- NA
alignments <- manual_invert_chr(alignments, Q_chromosomes)

offset_alignments_Q <- offset_chr(alignments, 'Q', chr_offset, chr_order_Q)
offset_alignments_RQ <- offset_chr(offset_alignments_Q$df, 'R', chr_offset, chr_order_R)
alignments <- offset_alignments_RQ$df
offset_list_R <- offset_alignments_RQ$offset_list
offset_list_Q <- offset_alignments_Q$offset_list
#chr_order_R <- offset_alignments_RQ$chr_order
#chr_order_Q <- offset_alignments_Q$chr_order

rough_max_end <-  max(max(alignments$Qend), max(alignments$Rend))
plot_length = rough_max_end  # make this nicer
# if want to centre the query chr:
if (max(alignments$Rend) > max(alignments$Qend)) { # if total ref length is greater than total query length, then plot size & adjustment is dictated by ref
  adjustment_length_Q <- (max(alignments$Rend) - max(alignments$Qend)) / 2  # i.e. half the difference between the two
  adjustment_length_R <- 0
} else{
  adjustment_length_R <- (max(alignments$Qend) - max(alignments$Rend)) / 2
  adjustment_length_Q <-0
}
#adjustment_length <- (max(alignments$Rend) - max(alignments$Qend)) / 2 # i.e. half the difference between the two

print(head(alignments))

pdf(paste0(args$output_prefix, '.pdf'))

plot(0,cex = 0, xlim = c(1, plot_length), ylim = c(((gap+1)*-1),(gap+1.5)), xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")


if (length(chr_order_R) <= 6){
  col_list <- c("#ffc759","#FF7B9C", "#607196", "#BABFD1", '#BACDB0', '#C6E2E9', '#F3D8C7')
} else { col_list <- c('#577590', '#617A8B', '#6B8086', '#748581', '#7E8A7C', '#889077', '#929572', '#9B9A6D',
                       '#A5A068', '#AFA563', '#B9AA5E', '#C2AF59', '#CCB554', '#D6BA4F', '#E0BF4A', '#E9C545', 
                       '#F3CA40', '#F1C544', '#EFC148', '#EDBC4C', '#EBB84F', '#E9B353', '#E7AF57', '#E5AA5B',
                       '#E3A55F', '#E1A163', '#DF9C67', '#DD986A', '#DB936E', '#D98F72', '#D78A76')
}
counter <- 1

col_list <- sapply(col_list, t_col, args$alpha)

for (i in chr_order_R){
  temp <- alignments[alignments$chrR == i,]
  plot_one_ref_chr(temp, col_list[[counter]], "red", adjustment_length_R, adjustment_length_Q)
  counter <- counter + 1
}

## chr outlines for query:
counter <- 1
for (i in chr_order_Q){
  temp <- alignments[alignments$chrQ == i,]
  Qfirst <- min(temp$Qstart)
  Qlast <- max(temp$Qend)
  segments(Qfirst+adjustment_length_Q, 1-gap, Qlast+adjustment_length_Q, 1-gap, lwd = 5)
}

## chr outlines for ref:
counter <- 1
for (i in chr_order_R){
  temp <- alignments[alignments$chrR == i,]
  Rfirst <- min(temp$Rstart)
  Rlast <- max(temp$Rend)
  segments(Rfirst+adjustment_length_R, gap, Rlast+adjustment_length_R, gap, lwd = 5)
}

## text labels for ref:
counter <- 1
for (i in chr_order_R){
  temp <- alignments[alignments$chrR == i,]
  Rfirst <- min(temp$Rstart)
  Rlast <- max(temp$Rend)
  offset <- offset_list_R[[counter]]
  text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = (gap+1.5), label = i,
       srt = 90, cex=0.5) # Rotation
  counter <- counter + 1
}

## text labels for query:
counter <- 1
for (i in chr_order_Q){
  temp <- alignments[alignments$chrQ == i,]
  Qfirst <- min(temp$Qstart)
  Qlast <- max(temp$Qend)
  text(x = ((Qlast+Qfirst+1)/2)+adjustment_length_Q, y = (gap+0.5)*-1, label = i,
       srt = 90, cex=0.5) # Rotation
  counter <- counter + 1
}

dev.off()

# chrom1 - M. cinxia REF
# chrom2 - V. cardui  QUERY

