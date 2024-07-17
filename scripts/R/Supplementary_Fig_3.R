#   Script name:  Supplementary_Fig_3.R
#   
#   Author:       Florian Janke
#   Last update:  16th July 2024
#
#-------------------------------

modality_diff <- function(table, info, timepoint){
  
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
      scale_y_continuous(name = "CPA score fold-change", trans = "log10", limits = c(0.1, 10)) +
      annotation_logticks(sides = "l") +
      scale_x_discrete(labels = c("VMAT", "CIRT")) +
      #labs(title = name) +
      #geom_hline(yintercept = cpa$thresholds$`100%`, linetype = "dashed", alpha = 0.5) +
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
  
  # Initiate list
  list <- list()
  
  # Loop through <timepoint>
  for(i in 1:length(timepoint)){
    
    # Subset <table> by <timepoint>
    tmp <- table[table$timepoint %in% timepoint[[i]], ]
    
    # Exclude patients with only one timepoint
    tmp <- tmp[tmp$Patient_ID %in% names(table(tmp$Patient_ID)[table(tmp$Patient_ID) == 2]), ]
    
    # Calculate difference
    tmp <- tmp %>% pivot_wider(names_from = timepoint, values_from = CPA, id_cols = Patient_ID)
    
    # Add irradiation modality information
    tmp <- merge(info[, c("Patient_ID", "Modality")], tmp, by = "Patient_ID")
    
    # Calculate relative difference
    tmp$diff <- tmp[[4]] /tmp[[3]]
    
    # Set order
    tmp$Modality <- factor(tmp$Modality, levels = c("VMAT", "CIRT"))
    tmp <- tmp[order(tmp$Modality), ]
    
    # Draw boxplot
    stat.tmp <- boxplot(tmp = tmp[, c("Modality", "diff")], type = "none")
    
    # Append
    list[["plot"]][[paste(timepoint[[i]], collapse = "|")]] <- stat.tmp$plot
    if(i == 1) list[["statistics"]] <- cbind(data.frame(Timepoint = paste(timepoint[[i]], collapse = "|")), stat.tmp$statistics)
    if(i != 1 ) list[["statistics"]] <- rbind(list[["statistics"]], cbind(data.frame(Timepoint = paste(timepoint[[i]], collapse = "|")), stat.tmp$statistics))
    
  }
  
  # Return <list>
  return(list)
  
}