#   Script name:  Supplementary_Fig_4.R
#   
#   Author:       Florian Janke
#   Last update:  16th July 2024
#
#-------------------------------

gtv_comp <- function(table, info, timepoint){
  
  #----
  # Functions
  #----
  rounding <- function(x){
    if(x < 0.001) y <- sprintf("%.2e", x) else y <- round(x, 3)
    return(y)
  }
  boxplot <- function(tmp, type, labels, cut){
    
    # Adjust colnames
    colnames(tmp) <- c("condition", "value")
    
    # Add box plot order
    if(type == "timepoint") name <- names(table(tmp$condition))[c(2, 1, 3, 4)] else name <- names(table(tmp$condition))
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
    if(type == "timepoint") color <- "black" else color <- c(lighten("#23505f", 0.75), "#e3ce97")
    
    base_size <- 14
    plot <- ggplot() +
      stat_boxplot(data = tmp, aes(x = box_order, y = value), geom = "errorbar", width = 0.15, lwd = 0.3) +
      {if(type == "timepoint") geom_boxplot(data = tmp, aes(x = box_order, y = value), outlier.shape = NA, color = "black")} +
      {if(type != "timepoint") geom_boxplot(data = tmp, aes(x = box_order, y = value), outlier.shape = NA, color = "black")} +
      geom_point(data = tmp, aes(x = point_order, y = value), size = 2.5, shape = 21, fill = "black") +
      scale_y_continuous(name = "CPA score", trans = "log10", limits = c(0.4, 10)) +
      annotation_logticks(sides = "l") +
      scale_x_discrete(labels = labels) +
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
  
  #----
  # Main script
  #----
  
  # Extract detectability threshold
  cut <- max(table[table$group == "control", ]$CPA)
  
  # Get the cohorts median GTV
  gtv <- round(median(info$GTV_ccm), 0)
  
  # Add GTV to <table>
  table <- merge(table, info[, c("Patient_ID", "GTV_ccm")], by = "Patient_ID")
  
  # Separate groups by cohort median GTV
  table$type <- ifelse(table$GTV_ccm >=gtv, paste0("≥", gtv), paste0("<", gtv))
  
  # Factorize <type>
  table$type <- factor(table$type, levels = c(paste0("<", gtv), paste0("≥", gtv)))
  
  # Draw box plot for all samples
  figA <- boxplot(tmp = table[, c("type", "CPA")], cut = cut, type = "none", labels = c("<19ccm", ">19ccm"))
  
  # Draw box plot per timepoint
  figB <- list()
  for(i in 1:length(timepoint)){
    
    # Subset <table>
    tmp <- table[grep(pattern = timepoint[i], x = table$timepoint), ]
    
    # Draw box plot
    tmp <- boxplot(tmp = tmp[, c("type", "CPA")], cut = cut, type = "none", labels = c("<19ccm", ">19ccm"))
    
    # Summarize
    figB[["plot"]][[timepoint[i]]] <- tmp$plot
    if(i == 1) figB[["statistics"]] <- cbind(data.frame(Timepoint = timepoint[i], tmp$statistics))
    if(i != 1) figB[["statistics"]] <- rbind(figB[["statistics"]], cbind(data.frame(Timepoint = timepoint[i], tmp$statistics)))
    
  }
  
  # Plot correlation
  corr <- dplyr::left_join(table[table$timepoint == "Baseline", c("Patient_ID", "CPA")], info[, c("Patient_ID", "GTV_ccm", "N", "M", "Histology")], by = "Patient_ID")
  
  # Exclude outliers
  corr_outlier <- corr[!corr$Patient_ID %in% c("P004", "P014"), ]
  
  # Add TNM labels
  corr$M <- ifelse(is.na(corr$M), "0", ifelse(corr$M == "1", "M1", "0"))
  corr$N <- ifelse(corr$N >0, "N≥1", "0")
  corr$distant <- ifelse(corr$N != "0", corr$N, ifelse(corr$M != "0", corr$M, "0"))
  
  # Plot xy-plot
  base_size <- 14
  plot <- ggplot() +
    geom_point(data = corr, aes(x = GTV_ccm, y = CPA, fill = distant), size = 3.5, shape = 21, color = "black") +
    geom_point(data = corr, aes(x = GTV_ccm, y = CPA, shape = Histology), size = 1, color = "white") +
    geom_smooth(data = corr, aes(x = GTV_ccm, y = CPA), method = "lm", se = FALSE, col = "black") +
    geom_smooth(data = corr_outlier, aes(x = GTV_ccm, y = CPA), method = "lm", se = FALSE, col = "gray75")+
    scale_shape_manual(values = c(3, 4)) +
    scale_fill_manual(values = c("black", "#23505f", lighten("#23505f", 0.75))) +
    scale_y_continuous(limits = c(0.4, 3.5), breaks = seq(0, 3.5, 1), name = "CPA score") +
    scale_x_continuous(limits = c(0, 150), name = "GTV (ccm)") +
    theme(
      text = element_text(family = "Arial"),
      axis.line =         element_line(colour = "black"),
      axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
      axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
      axis.title.x =      element_text(size = base_size * 1, angle = 0, vjust = 0, hjust = 0.5),
      axis.title.y =      element_text(size = base_size * 1, angle = 90, vjust = 1),
      axis.ticks.length = unit(0.3, "lines"),
      panel.background =  element_rect(fill = "transparent", colour = NA), 
      panel.border =      element_rect(fill = NA, colour = "transparent", linewidth = 0.75), 
      panel.grid.major.y =  element_blank(),
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
  
  # Get correlation statistics
  type <- c("total", "outliers")
  for(i in 1:length(type)){
    
    # Combine <table> and <info>
    corr <- dplyr::left_join(table[table$timepoint == "Baseline", c("Patient_ID", "CPA")], info[, c("Patient_ID", "GTV_ccm", "N", "M", "Histology")], by = "Patient_ID")
    
    # Exclude outliers
    if(type[i] == "outliers") corr <- corr[!corr$Patient_ID %in% c("P004", "P014"), ]
    
    # Get statistics
    tmp <- data.frame(Type = type[i],
                      Pearson_coeff = sprintf("%.3f", round(cor.test(corr$GTV_ccm, corr$CPA, method = "pearson")$estimate, 3)),
                      Pearson_p = sprintf("%.3f", round(cor.test(corr$GTV_ccm, corr$CPA, method = "pearson")$p.value, 3)),
                      Spearman_coeff = sprintf("%.3f", round(cor.test(corr$GTV_ccm, corr$CPA, method = "spearman")$estimate, 3)),
                      Spearman_p = sprintf("%.3f", round(cor.test(corr$GTV_ccm, corr$CPA, method = "spearman")$p.value, 3)))
    
    # Append
    if(i == 1) stat <- tmp
    if(i != 1) stat <- rbind(stat, tmp)
    
  }
  
  # Summarize correlation
  figC <- list()
  figC[["plot"]] <- plot
  figC[["statistics"]] <- stat
  
  # Summarize everything
  output <- list()
  output[["FigA"]] <- figA
  output[["FigB"]] <- figB
  output[["FigC"]] <- figC
  
  # Return <output>
  return(output)
  
}

