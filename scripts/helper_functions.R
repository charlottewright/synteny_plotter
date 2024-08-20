### functions ###

read_buscos <- function(file_name, prefix){
  chr_label <- paste0('chr', prefix)
  df <- read.csv(file_name, sep='\t', comment.char = '#', header = FALSE)[,c(0:6)]
  colnames(df) <- c('busco', 'status', chr_label, paste0(prefix, 'start'), paste0(prefix, 'end'), paste0(prefix, 'strand'))
  # these following lines are assuming '|' is only found in the rows with a bug in them
  if ( any(grepl("[|]", df[, chr_label])) ){
    # TODO: add warning here
    df[, chr_label] <- trim_strings(trim_strings(df[, chr_label], ":", 2), "[|]", 2)
  }
  if ( any(grepl("[:]", df[, chr_label])) ){
    # TODO: add warning here
    df[, chr_label] <- trim_strings(df[, chr_label], ":", 1)
  }
  df <- df[df$status == "Complete",]
  df <- subset(df, select=-c(status))
  return(df)
}

####### generate query function #####

get_most_frequent_ref_hit <- function(chr, alignments){
    hit_table <- alignments %>% filter(chrQ == chr) %>% select(chrR) %>% table
    names(hit_table)[which.max(hit_table)]
}

test_invert <- function(chr_row, alignments){
    chr_aln <- alignments %>% filter(chrQ == as.character(chr_row['chrQ'])) %>%  filter(chrR == as.character(chr_row['chrR']))
    if (nrow(chr_aln) <2){
      return (FALSE)
    }
    cortest <- cor.test(unlist(chr_aln[, 'Qstart']), unlist(chr_aln[, 'Rstart']))
    if (cortest$estimate < 0){
        return (TRUE)
    }
    return (FALSE)    
}

generate_auto_query_order <- function(chromosomes_R, busco_R, busco_Q){

    alignments <- merge(busco_Q, busco_R, by='busco') # merging Ref and Query tables into one by BUSCO IDs
    print(alignments)
    chromosomal_correspondences <- sapply(unique(alignments[ ,'chrQ']), get_most_frequent_ref_hit, alignments) # gets the chromosome with the highest number of hits
    chromosomal_correspondences <- data.frame(chrQ = names(chromosomal_correspondences), chrR = as.character(chromosomal_correspondences)) # making a table out of it

    row.names(chromosomes_R) <- chromosomes_R$chr # naming rows for easier extraction
    chromosomal_correspondences$R_order <- chromosomes_R[chromosomal_correspondences$chrR, 'order'] # adding reference order to the table
    chromosomal_correspondences <- chromosomal_correspondences %>% arrange(R_order) # orderining by the final order, here we would add some additional ordering parameters
  #  print(chromosomal_correspondences)
    chromosomal_correspondences$order <- 1:nrow(chromosomal_correspondences) # to give each a unique order number in case there is two chr map to the same ref chr
    chromosomal_correspondences$invert <- apply(chromosomal_correspondences, 1, test_invert, alignments)
    chromosomal_correspondences <- chromosomal_correspondences %>% select(c(chrQ, order, invert)) 
    colnames(chromosomal_correspondences)[1] <- 'chr' 

    return(chromosomal_correspondences)
}


#######

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

trim_strings <- function(values, delimiter, index){
  sapply(strsplit(values, delimiter), function(x){ x[index] })
  
}
perform_inverts <- function(df, Q_chromosomes){
  Q_chr <- Q_chromosomes$chr
  reorientated_alignments <- data.frame(matrix(ncol = length(df), nrow = 0))
  colnames(reorientated_alignments) <- colnames(df)
  for (q in Q_chr){
    df_subset <- df[df$chrQ == q,]
    if (Q_chromosomes[Q_chromosomes$chr == q,]$invert == "TRUE"){
      Q_end_value <- max(df_subset$Qend)
      df_subset$Qstart <- (df_subset$Qstart - Q_end_value)*-1 # should technically be length of chr not R_end value
      df_subset$Qend <- (df_subset$Qend - Q_end_value)*-1 # should technically be length of chr not R_end value
      df_subset$Qstrand[df_subset$Qstrand == '-'] <- '--'
      df_subset$Qstrand[df_subset$Qstrand == '+'] <- '-'
      df_subset$Qstrand[df_subset$Qstrand == '--'] <- '+'
      df_subset$Qstart_temp <- df_subset$Qstart
      df_subset$Qstart <- df_subset$Qend
      df_subset$Qend <- df_subset$Qstart_temp
      df_subset <- subset(df_subset, select = -c(Qstart_temp))
      reorientated_alignments <- rbind(reorientated_alignments, df_subset)
    }
    else{ # no flip!
      reorientated_alignments <- rbind(reorientated_alignments, df_subset)
    }
    
  }
  return(reorientated_alignments)
}


## check order is matching in the two chr i.e. that one chr isn't flipped relative to the other
auto_invert_chr <- function(df, ref_chr){
  df <- df[df$chrR == ref_chr,]
  df <- df %>% arrange(Qstart) # sort
  num_rows <- nrow(df) /20  # 5% of rows
  R_start_value <- mean(head(df, n=num_rows)$Rstart)
  R_end_value <- mean(tail(df, n=num_rows)$Rend)
  Q_start_value <- mean(head(df, n=num_rows)$Qstart)
  Q_end_value <- mean(tail(df, n=num_rows)$Qend)
  R_diff <- sign(R_end_value - R_start_value) # positive means (+) direction 
  Q_diff <- sign(Q_end_value - Q_start_value) # negative means (-) direction
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

offset_chr <- function(df, q_or_r, chr_offset, chr_order){
  counter = 0
  offset = 0
  offset_list <- list()
  offset_alignments <- data.frame(matrix(ncol = 10, nrow = 0))
  chr_prefix <- paste0('chr', q_or_r)
  #chr_order <- unique(df[[chr_prefix]]) # no longer needed as specify order in function
  chrStart <- paste0(q_or_r, 'start')
  chrEnd <- paste0(q_or_r, 'end')
  for (i in chr_order){
    chr_df <- df[df[[chr_prefix]] == i,]
    max_end <- max(chr_df[[chrStart]])
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

plot_one_ref_chr <- function(df, adjustment_length_R, adjustment_length_Q, y_offset, busco_2_colour, alpha){ 
  df$Qstart <- df$Qstart + adjustment_length_Q
  df$Qend <- df$Qend + adjustment_length_Q
  df$Rstart <- df$Rstart + adjustment_length_R
  df$Rend <- df$Rend + adjustment_length_R
  Qstarts <-df$Qstart
  Qends <- df$Qend
  Rstarts <- df$Rstart # was alignments
  Rends <- df$Rend # was alignments
  #cols =  c("blue", "red")
#  cols = c(col1, col2)
 # border = ifelse(sign(Rends - Rstarts) == sign(Qends - Qstarts), cols[1], cols[2]) # if R&Q are same sign, use blue, else use red
 # col = border
  for (i in 1:nrow(df)){ # curved lines- sigmoid connector = x1,y1,x2,y2
    busco <-df$busco[i]
    col1 <- busco_2_colour$chr_colour[busco_2_colour$busco == busco]
    #col2 <- t_col(col1, alpha) # this is occasionally erroring - I don't know why?
    col2 <- col1 # for now lets keep both strands same colour rather than reverse strand being half as transparent (see above line)
    cols = c(col1, col2)
    border = ifelse(sign(Rends - Rstarts) == sign(Qends - Qstarts), cols[1], cols[2]) # if R&Q are same sign, use blue, else use red
    col = border
    lines.to.poly(sigmoid.connector(Qstarts[i], 1-gap-y_offset, Rstarts[i], 0+gap-y_offset, vertical=T),
                  sigmoid.connector(Qends[i], 1-gap-y_offset, Rends[i], 0+gap-y_offset, vertical=T),
                  col = col[i], border=ifelse(show_outline==FALSE, NA, border[i]), lwd=lwd)
  }
}

# stolem from "the internet"; DOES NOT TAKE VECTORS
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

## Get RGB values for named color
rgb.val <- col2rgb(color)

## Make new color using input color as base and alpha set by transparency
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)

## Save the color
return(t.col)
}
