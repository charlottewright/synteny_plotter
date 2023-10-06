# cutoff of block in terms of same strand and stuff, and distance until next ortholog
library(dplyr)

setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Data/')
assignments <- read.csv('../../Chromosome_evolution_Lepidoptera_MS/sup_tables/TableS4_Merian_element_definitions.tsv', sep='\t', header=FALSE)[,c(1,3)]
colnames(assignments) <- c('busco', 'merian')

read_buscos <- function(file_path, file_name, prefix){
  df <- read.csv(paste(file_path, file_name, sep=''), sep='\t', comment.char = '#', header = FALSE)[,c(0:6)]
  colnames(df) <- c('busco', 'status', paste0('chr', prefix), paste0(prefix, 'start'), paste0(prefix, 'end'), paste0(prefix, 'strand'))
  df <- df[df$status == "Complete",]
  df <- subset(df, select=-c(status))
  df <- df[!grepl(':', df[[paste0('chr', prefix)]]),]
#  df <- df[!grepl(':', df$chrQ),]
  return(df)
}

#make a polygon from two lines - SIMON
lines.to.poly <- function(l1,l2, col=NULL, border=NULL, lwd=NULL){
  polygon(c(l1[,1],rev(l2[,1])), c(l1[,2],rev(l2[,2])), col=col, border=border)
}

sigmoid.connector <- function(x1,y1,x2,y2, curvature=10, steps=50, vertical=FALSE){
  if (vertical==TRUE) {
    vals <- c(x1,y1,x2,y2)
    x1 <- vals[2]
    x2 <- vals[4]
    y1 <- vals[1]
    y2 <- vals[3]
  }
  x <- seq(x1,x2,(x2-x1)/steps)
  x_norm <- (x-(x1+x2)/2)/((x2-x1)/2)
  y_norm <- 1/(1+exp(-x_norm*curvature))
  y <- y1+(y_norm*(y2-y1))
  if (vertical==TRUE) return(cbind(y,x))
  cbind(x,y)
}

## check order is matching in the two chr i.e. that one chr isn't flipped relative to the other
orientate_chr <- function(df, ref_chr){
  df <- df[df$chrR == ref_chr,]
  df <- df %>% arrange(Qstart) # sort
  num_rows <- nrow(df) /20  # 5% of rows
  print(num_rows)
  R_start_value <- mean(head(df, n=num_rows)$Rstart)
  R_end_value <- mean(tail(df, n=num_rows)$Rend)
  print(R_start_value)
  print(R_end_value)
  Q_start_value <- mean(head(df, n=num_rows)$Qstart)
  Q_end_value <- mean(tail(df, n=num_rows)$Qend)
  print(Q_start_value)
  print(Q_end_value)
  R_diff <- sign(R_end_value - R_start_value) # positive means (+) direction 
  Q_diff <- sign(Q_end_value - Q_start_value) # negative means (-) direction
  print(R_diff)
  print(Q_diff)
  # This code block aims to a ref chr if its in different direction to the query
  if (R_diff != Q_diff){
    df$Rstart <- (df$Rstart - R_end_value)*-1 # should technically be length of chr not R_end value
    df$Rend <- (df$Rend - R_end_value)*-1 # should technically be length of chr not R_end value
    df$Rstrand[df$Rstrand == '-'] <- '--'
    df$Rstrand[df$Rstrand == '+'] <- '-'
    df$Rstrand[df$Rstrand == '--'] <- '+'
    df$Rstart_temp <- df$Rstart
    df$Rstart <- df$Rend
    df$Rend <- df$Rstart_temp
    df <- subset(df, select = -c(Rstart_temp) )
  }
  return(df)
}

offset_chr <- function(df, q_or_r, chr_offset){
  counter = 0
  offset = 0
  offset_list <- list()
  offset_alignments <- data.frame(matrix(ncol = 10, nrow = 0))
  chr_prefix <- paste0('chr', q_or_r)
  chr_order <- unique(df[[chr_prefix]])
  chrStart <- paste0(q_or_r, 'start')
  chrEnd <- paste0(q_or_r, 'end')
  for (i in chr_order){
    chr_df <- df[df[[chr_prefix]] == i,]
    max_end <- max(chr_df[[chrStart]])
    print(max_end)
    if (counter != 0){ # only need to offset start/end if this is not the first chr
      chr_df <- chr_df %>% arrange(Rstart) # sort by ref regardless of it ref or query
      chr_df[[chrStart]] <- chr_df[[chrStart]] + offset  # allows for accumulative chr positions
      chr_df[[chrEnd]] <- chr_df[[chrEnd]] + offset # allows for accumulative chr positions
    }
    offset_alignments <- rbind(offset_alignments, chr_df)
    offset_list <- append(offset_list, offset)
    offset <- offset + max_end + chr_offset # accumulative offset
    counter <- counter + 1
  }
  output_list <- list('df' = offset_alignments, 'offset_list' = offset_list, 'chr_order'=chr_order)
  return(output_list)
}
    
plot_one_ref_chr <- function(df, col1, col2, adjustment_length_R, adjustment_length_Q){ 
  df$Qstart <- df$Qstart + adjustment_length_Q
  df$Qend <- df$Qend + adjustment_length_Q
  df$Rstart <- df$Rstart + adjustment_length_R
  df$Rend <- df$Rend + adjustment_length_R
  Qstarts <-df$Qstart
  Qends <- df$Qend
  Rstarts <- df$Rstart # was alignments
  Rends <- df$Rend # was alignments
  #cols =  c("blue", "red")
  cols = c(col1, col2)
  border = ifelse(sign(Rends - Rstarts) == sign(Qends - Qstarts), cols[1], cols[2]) # if R&Q are same sign, use blue, else use red
  col = border
  for (i in 1:nrow(df)){ # curved lines- sigmoid connector = x1,y1,x2,y2
    lines.to.poly(sigmoid.connector(Qstarts[i], 1-gap, Rstarts[i], 0+gap, vertical=T),
                  sigmoid.connector(Qends[i], 1-gap, Rends[i], 0+gap, vertical=T),
                  col = col[i], border=ifelse(show_outline==FALSE, NA, border[i]), lwd=lwd)
  }
}

### specify arguments ###
file_path = 'BUSCOs/All/'
gap =6
show_outline = TRUE
chr_offset = 20000000

#chr_offset = 5000000
R_ID <- 'Melitaea_cinxia.tsv'
Q_ID <- 'Lysandra_coridon.tsv'
### read in data ###
R_df <- read_buscos(file_path, R_ID, 'R')
Q_df <- read_buscos(file_path, Q_ID, 'Q') # was Agrochola_circellaris.tsv
#R_df <- read_buscos(file_path, "Lysandra_coridon.tsv", 'R') #

# OU611839.1 is fused in A. circellaris
alignments <- merge(Q_df, R_df, by='busco')
alignments <- merge(alignments, assignments, by='busco')
#alignments <- alignments[alignments$chrQ %in% c('HG992236.1','HG992211.1','HG992210.1', 'HG992209.1','HG992235.1','HG992237.1'),] # fused chr
alignments <- alignments %>% group_by(chrR) %>% filter(n() > 5) %>% ungroup()
#alignments <- merge(alignments, assignments)

offset_alignments_Q <- offset_chr(alignments, 'Q', chr_offset)
offset_alignments_RQ <- offset_chr(offset_alignments_Q$df, 'R', chr_offset)
alignments <- offset_alignments_RQ$df
offset_list_R <- offset_alignments_RQ$offset_list
offset_list_Q <- offset_alignments_Q$offset_list
chr_order_R <- offset_alignments_RQ$chr_order
chr_order_Q <- offset_alignments_Q$chr_order

rough_max_end <-  max(max(alignments$Qend), max(alignments$Rend))
plot_length = rough_max_end  # make this nicer
# if want to centre the query chr:
if (max(alignments$Rend) > max(alignments$Qend)) { # if total ref length is greater than total query length, then plot size & adjustment is dictated by ref
  adjustment_length_Q <- (max(alignments$Rend) - max(alignments$Qend)) / 2  # i.e. half the difference between the two
  adjustment_length_Q <- 0
} else{
  adjustment_length_R <- (max(alignments$Qend) - max(alignments$Rend)) / 2 
  adjustment_length_Q <-0
}
#adjustment_length <- (max(alignments$Rend) - max(alignments$Qend)) / 2 # i.e. half the difference between the two

plot(0,cex = 0, xlim = c(1, plot_length), ylim = c(((gap+0.2)*-1),(gap+1)), xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")




if (length(chr_order_R) <= 6){
  col_list <- c("#ffc759","#FF7B9C", "#607196", "#BABFD1", '#BACDB0', '#C6E2E9', '#F3D8C7')
} else { col_list <- c('#577590', '#617A8B', '#6B8086', '#748581', '#7E8A7C', '#889077', '#929572', '#9B9A6D',
                '#A5A068', '#AFA563', '#B9AA5E', '#C2AF59', '#CCB554', '#D6BA4F', '#E0BF4A', '#E9C545', 
                '#F3CA40', '#F1C544', '#EFC148', '#EDBC4C', '#EBB84F', '#E9B353', '#E7AF57', '#E5AA5B',
                '#E3A55F', '#E1A163', '#DF9C67', '#DD986A', '#DB936E', '#D98F72', '#D78A76')
}
counter <- 1

for (i in chr_order_R){
  temp <- alignments[alignments$chrR == i,]
  plot_one_ref_chr(temp, col_list[[counter]], "red", adjustment_length_R, adjustment_length_Q)
  counter <- counter + 1
}

## chr outlines for ref:
counter <- 1
for (i in chr_order_Q){
  temp <- alignments[alignments$chrQ == i,]
  Qfirst <- min(temp$Qstart)
  Qlast <- max(temp$Qend)
  segments(Qfirst+adjustment_length, 1-gap, Qlast+1+adjustment_length_Q, 1-gap, lwd = 5)
}

## chr outlines for query:
counter <- 1
for (i in chr_order_R){
  temp <- alignments[alignments$chrR == i,]
  Rfirst <- min(temp$Rstart)
  Rlast <- max(temp$Rend)
  segments(Rfirst+adjustment_length, gap, Rlast+1+adjustment_length_R, gap, lwd = 5)
}

## text labels for ref:
counter <- 1
for (i in chr_order_R){
  temp <- alignments[alignments$chrR == i,]
  Rfirst <- min(temp$Rstart)
  Rlast <- max(temp$Rend)
  offset <- offset_list_R[[counter]]
  text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = (gap+1), label = i,
       srt = 90, cex=0.5) # Rotation
  counter <- counter + 1
  }

## text labels for query:
counter <- 1
for (i in chr_order_Q){
  temp <- alignments[alignments$chrQ == i,]
  Qfirst <- min(temp$Qstart)
  Qlast <- max(temp$Qend)
  text(x = ((Qlast+Qfirst+1)/2)+adjustment_length_Q, y = (gap)*-1, label = i,
       srt = 90, cex=0.5) # Rotation
  counter <- counter + 1
}


