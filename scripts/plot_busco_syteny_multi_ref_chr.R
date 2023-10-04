# cutoff of block in terms of same strand and stuff, and distance until next ortholog
library(dplyr)

setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Data/')
assignments <- read.csv('../Chromosome_evolution_Lepidoptera_MS/sup_tables/TableS4_Merian_element_definitions.tsv', sep='\t', header=FALSE)[,c(1,3)]
colnames(assignments) <- c('busco', 'merian')

read_buscos <- function(file_path, file_name, prefix){
  df <- read.csv(paste(file_path, file_name, sep=''), sep='\t', comment.char = '#', header = FALSE)[,c(0:6)]
  colnames(df) <- c('busco', 'status', paste('chr', prefix, sep=''), paste(prefix, 'start', sep=''), paste(prefix, 'end', sep=''), paste(prefix, 'strand', sep=''))
  df <- df[df$status == "Complete",]
  df <- subset(df, select=-c(status))
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


file_path = 'Data/BUSCOs/All/'
Q_df <- read_buscos(file_path, "Melitaea_cinxia.tsv", 'Q') # was Agrochola_circellaris.tsv
R_df <- read_buscos(file_path, "Micropterix_aruncella.tsv", 'R')
#R_df <- read_buscos(file_path, "Lysandra_coridon.tsv", 'R') #

gap =6
show_outline = TRUE


# OU611839.1 is fused in A. circellaris
alignments <- merge(Q_df, R_df, by='busco')
alignments <- merge(alignments, assignments, by='busco')
alignments <- alignments[alignments$chrQ == 'HG992211.1',] # fused chr
alignments <- alignments %>% group_by(chrR) %>% filter(n() > 5) %>% ungroup()
#alignments <- merge(alignments, assignments)

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

# in this case, using one query chr - could have >1 ref chr
# find which ref chr is "first" in terms of synteny to the query
temp <- alignments %>% arrange(Qstart) 
R_chr_order <- unique(temp$chrR)

reorientated_alignents <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(reorientated_alignents) <- colnames(alignments)
offset_list <- list()
offset <- 0
chr_offset <- 5000000 # distance between chr - make this a prop of chr length
counter <- 0

for (i in R_chr_order){
  print(i)
  R_chr_df <- orientate_chr(alignments, i)
  max_Rend <- max(R_chr_df$Rend)
  if (counter != 0){ # only need to offset start/end if this is not the first chr
    R_chr_df$Rstart <- R_chr_df$Rstart + offset  # allows for accumulative chr positions
    R_chr_df$Rend <- R_chr_df$Rend + offset # allows for accumulative chr positions
  }
  reorientated_alignents <- rbind(reorientated_alignents, R_chr_df)
  offset_list <- append(offset_list, offset)
  offset <- offset + max_Rend + chr_offset# accumulative offset
  counter <- counter + 1
}

#first_ref_chr <- (alignments %>% arrange(Qstart) %>% head(n=1))$chrR
#second_ref_chr <- (alignments %>% arrange(Qstart) %>% tail(n=1))$chrR
#R1_df <- orientate_chr(alignments, first_ref_chr)
#R2_df <- orientate_chr(alignments, second_ref_chr)
#alignments <- rbind(R1_df, R2_df)
alignments <- reorientated_alignents

plot_one_ref_chr <- function(df, col1, col2, adjustment=0){ # by default dont adjust positio of query chr
  df$Qstart <- df$Qstart + adjustment_length
  df$Qend <- df$Qend + adjustment_length
  Qstarts <-df$Qstart
  Qends <- df$Qend
  Rstarts <- df$Rstart # was alignments
  Rends <- df$Rend # was alignments
  
  #cols =  c("blue", "red")
  cols = c(col1, col2)
  border = ifelse(sign(Rends - Rstarts) == sign(Qends - Qstarts), cols[1], cols[2]) # if R&Q are same sign, use blue, else use red
 # col = paste0(border,"50")
  col = border
  # sigmoid connector = x1,y1,x2,y2
  # x1 = qstart
  # x2 = rstart
  for (i in 1:nrow(df)){ # curved lines
    lines.to.poly(sigmoid.connector(Qstarts[i], 1-gap, Rstarts[i], 0+gap, vertical=T),
                  sigmoid.connector(Qends[i], 1-gap, Rends[i], 0+gap, vertical=T),
                  col = col[i], border=ifelse(show_outline==FALSE, NA, border[i]), lwd=lwd)
  }
}

# plot_length = max(c(Qlast, Rlast) + 1)
rough_max_end <- max(alignments$Rend)
plot_length = rough_max_end + 1  # make this nicer
# if want to centre the query chr:
adjustment_length <- (max(alignments$Rend) - max(alignments$Qend)) / 2 # i.e. half the difference between the two

plot(0,cex = 0, xlim = c(1, plot_length), ylim = c(((gap+0.2)*-1),(gap+0.2)), xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")

col_list <- c("#ffc759","#FF7B9C", "#607196", "#BABFD1")
counter <- 1

for (i in R_chr_order){
  temp <- alignments[alignments$chrR == i,]
  plot_one_ref_chr(temp, col_list[[counter]], "red", adjustment =adjustment_length)
  counter <- counter + 1
}

max(alignments$Rend)/2
# query is easy - just plot one chr:
Qfirst <- min(alignments$Qstart)
Qlast <- max(alignments$Qend)

segments(Qfirst+adjustment_length, 1-gap, Qlast+1+adjustment_length, 1-gap, lwd = 5) 
#segments(Qfirst, 0, Qlast+1, 0, lwd = 5)  # if don't want to adapt the query chr (keep start at 0)
counter <- 1
for (i in R_chr_order){
  temp <- alignments[alignments$chrR == i,]
  segments(min(temp$Rstart), gap, max(temp$Rend), gap, lwd = 5) # I swapped (Rlast-Rfirst+1) & (Qlast-Qfirst+1) around
}

# query is easy to label:
mtext(text=alignments$chrQ[1], side=1, at=((Qlast+Qfirst+1)/2)+adjustment_length,)


counter <- 1
for (i in R_chr_order){
  temp <- alignments[alignments$chrR == i,]
  offset <- offset_list[[counter]]
  mtext(text=temp$chrR[1], side=3, at=((max(temp$Rend)-min(temp$Rstart)+1)/2)+offset,)
  counter <- counter + 1
  }



