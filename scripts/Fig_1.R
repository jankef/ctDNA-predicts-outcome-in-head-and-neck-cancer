#   Script name:  Fig_1.R
#   
#   Author:       Florian Janke
#   Last update:  17th July 2024
#
#-------------------------------

cpa_and_cnv <- function(data, cytoband, input){
  
  #----
  # Function
  #----
  rounding <- function(x){
    if(x < 0.001) y <- sprintf("%.2e", x) else y <- round(x, 3)
    return(y)
  }
  boxplot <- function(tmp, type){
    
    # Adjust colnames
    colnames(tmp) <- c("condition", "value")
    
    # Extract cut-off
    cut <- max(tmp[tmp$condition == "control", ]$value)
    
    # Add box plot order
    if(type == "timepoint") name <- names(table(tmp$condition))[c(2, 1, 3, 4)] else name <- names(table(tmp$condition))
    if(type == "group") name <- rev(name)
    numb <- c(1:length(table(tmp$condition)))
    names(numb) <- name
    tmp$box_order <- NA
    for(i in 1:nrow(tmp)) tmp[i, ]$box_order <- as.character(which(names(numb) == tmp[i, ]$condition))
    
    # Add point order
    tmp <- tmp[order(tmp$value), ]
    tmp <- tmp[order(tmp$box_order), ]
    tmp$point_order <- NA
    for(i in 1:length(numb)) tmp[tmp$box_order == as.character(i), ]$point_order <- seq((-0.3 +i), (0.3 +i), length.out = nrow(tmp[tmp$box_order == as.character(i), ]))
    
    # Draw boxplot
    base_size <- 14
    plot <- ggplot() +
      stat_boxplot(data = tmp, aes(x = box_order, y = value), geom = "errorbar", width = 0.15, lwd = 0.3) +
      geom_boxplot(data = tmp, aes(x = box_order, y = value), outlier.shape = NA, color = "black") +
      geom_point(data = tmp, aes(x = point_order, y = value), size = 2.5, shape = 21, fill = "black") +
      scale_y_continuous(name = "CPA score fold-change", trans = "log10", limits = c(0.4, 10)) +
      annotation_logticks(sides = "l") +
      scale_x_discrete(labels = c("Control", "Case")) +
      geom_hline(yintercept = cut, linetype = "dashed", alpha = 0.5) +
      theme(
        text = element_text(family = "Arial"),
        axis.line =         element_line(colour = "black"),
        axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
        axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
        axis.title.x =      element_blank(),
        axis.title.y =      element_text(size = base_size * 1, angle = 90, vjust = 0.5),
        axis.ticks.length = unit(0.3, "lines"),
        panel.background =  element_rect(fill = "transparent", colour = NA), 
        panel.border =      element_rect(fill = NA, colour = "transparent", linewidth = 0.75), 
        panel.grid.major.y =  element_line(colour = "grey90"),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.y =  element_blank(),
        panel.grid.minor.x =  element_blank(),
        strip.background =  element_rect(color = "black", linewidth = 0.75, fill = "#9ec5c8"),
        strip.text.x =      element_text(size = base_size * 1, angle = 0, vjust = 1),
        plot.background =   element_rect(colour = NA, fill = "transparent"),
        plot.title =        element_text(size = base_size * 1),
        legend.position = "none",
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.55, 'cm'),
        legend.text = element_text(size=10)
      )
    
    # Summarize statistics
    comb <- expand.grid(numb, numb)
    comb <- comb[comb$Var1 != comb$Var2, ]
    comb <- subset(comb, Var1 <= Var2)
    for(i in 1:nrow(comb)){
      
      stat.tmp <- data.frame(Group_01 = names(comb$Var1[i]), Group_02 = names(comb$Var2[i]),
                             Positive_01 = nrow(tmp[tmp$box_order == comb$Var1[i] & tmp$positive == 1, ]),
                             Total_01 = nrow(tmp[tmp$box_order == comb$Var1[i], ]),
                             Positive_02 = nrow(tmp[tmp$box_order == comb$Var2[i] & tmp$positive == 1, ]),
                             Total_02 = nrow(tmp[tmp$box_order == comb$Var2[i], ]),
                             p = rounding(wilcox.test(tmp[tmp$box_order == comb$Var2[i], ]$value, tmp[tmp$box_order == comb$Var1[i], ]$value, alternative = "two.sided", paired = FALSE)$p.value))
      nrow(tmp[tmp$box_order == comb$Var2[i] & tmp$positive == 1, ])
      nrow(tmp[tmp$box_order == comb$Var2[i], ])
      
      if(i == 1) stat <- stat.tmp
      if(i != 1) stat <- rbind(stat, stat.tmp)
      
    }
    
    list <- list()
    list[["plot"]] <- plot
    list[["statistics"]] <- stat
    
    # Return <stat>
    return(list)
    
  }
  cnv_localization <- function(input, pID, cytoband){
    
    # Initiate per-patient-loop
    for(i in 1:length(pID)){
      
      # List samples of current patient
      sampID <- list.files(input, full.names = FALSE, pattern = pID[i])
      
      # Initiate per-sample-loop
      for(n in 1:length(sampID)){
        
        # Open '*_bins.bed' file
        tmp <- data.frame(fread(file.path(input, sampID[n])))[, c(1:4, 7)]
        
        # Replace empty cells with 'neutral'
        tmp$V7 <- ifelse(tmp$V7 == "", "neutral", tmp$V7)
        
        # Re-name columns
        colnames(tmp) <- c("id", "chr", "start", "end", sampID[n])
        
        # Combine copy number state from all samples of one patient
        if(n == 1) state <- tmp
        if(n != 1) state <- dplyr::left_join(state, tmp[, c(1, 5)], by = "id")
        
      }
      
      # Combine copy number states of longitudinal samples
      for(n in 1:nrow(state)){
        
        # Get copy number state frequency
        tmp <- data.frame(table(t(state[n, c(5:ncol(state))])[, 1]))
        
        # Make call
        call <- ifelse(nrow(tmp) == 1, as.character(tmp$Var1),
                       ifelse("gain" %in% tmp$Var1 & "loss" %in% tmp$Var1, "neutral",
                              ifelse(tmp[tmp$Var1 == "neutral", ]$Freq == (ncol(state) -4), "neutral",
                                     ifelse("gain" %in% tmp$Var1, "gain", "loss"))))
        
        # Add call to <state>
        if(n == 1){
          
          state.tmp <- state
          state.tmp$state <- NA
          
        }
        state.tmp[n, ]$state <- call
        if(n == nrow(state)){
          
          state <- state.tmp
          rm(state.tmp)
          state <- state[, c(1:4, ncol(state))]
          colnames(state)[5] <- pID[i]
          
        }
        
      }
      
      # Append <state>
      if(i == 1) cnv <- state
      if(i != 1) cnv <- dplyr::left_join(cnv, state[, c(1, ncol(state))], by = "id")
      
    }
    
    # Exclude sex chromosomes
    cnv <- cnv[!cnv$chr %in% c("X", "Y"), ]
    
    # Summarize 'gain' and 'loss' occurrences
    for(i in 1:nrow(cnv)){
      
      # Count copy number state occurrences
      tmp <- data.frame(table(t(cnv[i, c(5:ncol(cnv))])[, 1]))
      
      # Summarize
      tmp <- data.frame(id = cnv$id[i], chr = cnv$chr[i], start = cnv$start[i], end = cnv$end[i],
                        gain = ifelse("gain" %in% tmp$Var1, tmp[tmp$Var1 == "gain", ]$Freq, 0),
                        loss = ifelse("loss" %in% tmp$Var1, tmp[tmp$Var1 == "loss", ]$Freq, 0),
                        neutral = ifelse("neutral" %in% tmp$Var1, tmp[tmp$Var1 == "neutral", ]$Freq, 0))
      
      # Append
      if(i == 1) state <- tmp
      if(i != 1) state <- rbind(state, tmp)
      
    }
    
    # Order by chromosome and start position
    state <- state[order(as.numeric(state$start)), ]
    state <- state[order(as.numeric(state$chr)), ]
    
    # Introduce running bin number
    state$bin <- c(1:nrow(state))
    
    # Load cytobands
    cyto <- data.frame(fread(cytoband))[,c(1:4)]
    colnames(cyto) <- c("chr", "start", "end", "arm")
    cyto$arm <- str_sub(cyto$arm, 1, 1)
    cyto <- cyto[!cyto$chr %in% c("chrX", "chrY"), ]
    cyto$chr <- as.character(gsub("chr", "", cyto$chr))
    cyto <- cyto[cyto$arm == "q", ]
    cyto <- cyto %>% dplyr::group_by(chr) %>% dplyr::summarize(end = min(start))
    
    # Get closest bin to cytoband
    for(i in 1:22){
      
      # Extract chromosome-of-interest
      tmp <- state[state$chr == i, ]
      
      # Find nearst bin
      tmp$distance <- abs(tmp$end -cyto[cyto$chr == i, ]$end)
      bin <- tmp[tmp$distance == min(tmp$distance), ]$bin[1]
      
      # Add to cytoband
      if(i == 1) cyto$cyto <- NA
      cyto[cyto$chr == i, ]$cyto <- bin
      
    }
    
    # Create chromosome annotation
    chr <- state %>% dplyr::group_by(chr) %>% dplyr::summarize(count = n())
    chr$chr <- factor(chr$chr, levels = c(1:22))
    chr <- chr[order(chr$chr), ]
    chr$counts_half <- chr$count / 2
    chr$cumsum <- cumsum(chr$count)
    chr$tick <- chr$cumsum -chr$count + chr$counts_half
    chr <- merge(chr, cyto[, c("chr", "cyto")], by = "chr")
    
    # Get relative number of gains and losses
    state$gain_rel <- state$gain /rowSums(state[, c(5:7)])
    state$loss_rel <- (state$loss /rowSums(state[, c(5:7)])) *-1
    
    # Adjust state to include gains and losses per region
    gain <- state[, c("id", "bin", "gain_rel")]
    colnames(gain)[3] <- "freq"
    gain$state <- "gain"
    loss <- state[, c("id", "bin", "loss_rel")]
    colnames(loss)[3] <- "freq"
    loss$state <- "loss"
    state <- rbind(gain, loss)
    rm(gain, loss)
    
    # Plotting parameters
    colors <- c("gain" = "#8a0000", "loss" = "#006300")
    base_size <- 14
    x.limit <- nrow(state) /2 +1
    
    # Plot CNV frequency
    plot <- ggplot() +
      geom_rect(data = chr, aes(xmin = cyto, xmax = cumsum, ymin = -Inf, ymax = Inf), fill = "gray90") +
      geom_bar(data = state, aes(x = bin, y = freq, color = state, group = state), stat = "identity") +
      scale_color_manual(values = colors) +
      geom_vline(xintercept = c(0, chr$cumsum), linetype = "dashed", alpha = 0.5, linewidth = 0.4) +
      scale_x_continuous(name = "Chromosomes", limits = c(0, x.limit), breaks = chr[!chr$chr %in% c(15,17,19,21, 13, 20, 11), ]$tick, labels = chr[!chr$chr %in% c(15,17,19,21, 13, 20, 11), ]$chr) +
      geom_hline(yintercept = 0, color = "grey20") +
      scale_y_continuous(limits = c(-1, 1), name = "CNV frequency (%)", labels = c("1.0", "0.5", "0.0", "0.5", "1.0")) +
      theme(
        text = element_text(family = "Arial"),
        axis.line =         element_blank(),
        axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 0.5, angle = 0),
        axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
        axis.title.x =      element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 0.5),
        axis.title.y =      element_text(size = base_size * 0.9, angle = 90, vjust = 2),
        axis.ticks.x =      element_blank(), 
        axis.ticks.length.x = unit(0.3, "lines"),
        axis.ticks.length.y = unit(0.2, "lines"),
        panel.background =  element_rect(fill = "transparent", colour = NA),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "grey20"),
        panel.border =      element_blank(), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank(),
        strip.background =  element_rect(fill = "grey80", colour = "grey50"), 
        strip.text.x =      element_text(size = base_size * 0.8),
        strip.text.y =      element_text(size = base_size * 0.8, angle = -90),
        plot.background =   element_rect(colour = NA, fill = "transparent"),
        plot.title =        element_text(size = base_size * 1.1),
        legend.position = "none"
      )
    
    # Return <plot>
    return(plot)
    
  }
  
  #----
  # Main script
  #----
  
  # Draw boxplot
  box_plot <- boxplot(tmp = data[, c("group", "CPA")], type = "group")
  
  # Plot recurren CNVs
  cnv_plot <- cnv_localization(input = input,
                               pID = as.character(data.frame(table(data[data$positive == "positive", ]$Patient_ID))[[1]]),
                               cytoband = cytoband)
  
  # Summarize
  list <- list()
  list[["figA"]] <- box_plot
  list[["figB"]] <- cnv_plot
  
  # Return <list>
  return(list)
  
}




