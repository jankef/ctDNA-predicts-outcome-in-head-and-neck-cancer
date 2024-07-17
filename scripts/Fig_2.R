#   Script name:  Fig_2.R
#   
#   Author:       Florian Janke
#   Last update:  17th July 2024
#
#-------------------------------

oncoprint <- function(table, volume){
  
  # Subset and order by GTV
  volume <- volume[, c("Patient_ID", "GTV_ccm")]
  volume <- volume[order(volume$GTV_ccm, decreasing = FALSE), ]
  volume$Patient_ID <- factor(volume$Patient_ID, levels = volume$Patient_ID)
  
  # Create bar plot of GTVs
  base_size <- 14
  volume_plot <- ggplot(data = volume, aes(x = Patient_ID, y = GTV_ccm), fill = "gray30") +
    geom_bar(stat = "identity", color = "black") +
    theme(
      text = element_text(family = "Arial"),
      axis.line =         element_line(colour = "black"),
      axis.text.x =       element_blank(),
      axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
      axis.title.x =      element_blank(),
      axis.title.y =      element_blank(),
      axis.ticks.length.x = unit(0.0, "lines"),
      axis.ticks.length.y = unit(0.3, "lines"),
      panel.background =  element_rect(fill = "transparent", colour = NA), 
      panel.grid.major.y =  element_blank(),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.y =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      strip.background =  element_rect(color = "black", linewidth = 0.75, fill = "#9ec5c8"),
      strip.text.x =      element_text(size = base_size * 1, angle = 0, vjust = 1),
      plot.background =   element_rect(colour = NA, fill = "transparent"),
      plot.title =        element_text(size = base_size * 1.5),
      legend.position = "none",
      legend.key = element_rect(colour = "transparent", fill = "transparent"),
      legend.title = element_text(size = 12),
      legend.key.size = unit(0.55, 'cm'),
      legend.text = element_text(size=10)
    )
  
  
  # Extract naive CPA-scores of 'cases'
  naive <- ctCPA$naive[ctCPA$naive$group == "case", ]
  
  # Extract positive samples per analysis
  comp <- names(table)
  for(i in 1:length(comp)){
    
    # Extract positive samples
    tmp <- ctCPA[[comp[i]]]
    tmp <- tmp[tmp$positive == "positive", ]
    
    # Extract time point of interest
    #if(comp[i] == "naive") timepoint <- "Baseline"
    if(comp[i] == "Baseline") timepoint <- c("Baseline", "5-Fx")
    if(comp[i] == "Baseline|5-Fx") timepoint <- c("5-Fx", "10-Fx")
    if(comp[i] == "Baseline|5-Fx|10-Fx") timepoint <- c("RT-end", "MRI_01", "MRI_02", "MRI_03")
    tmp <- tmp[tmp$timepoint %in% timepoint, ]
    
    # Subset columns
    tmp <- tmp[, c("Patient_ID", "timepoint", "positive")]
    
    # Append
    if(i == 1) informed <- tmp
    if(i != 1) informed <- rbind(informed, tmp)
    
  }
  
  # Ordering
  naive$timepoint <- factor(naive$timepoint, levels = rev(c("Baseline", "5-Fx", "10-Fx", "RT-end", "MRI_01", "MRI_02", "MRI_03", "control")))
  naive$Patient_ID <- factor(naive$Patient_ID, levels = volume$Patient_ID)
  
  # Plot oncoprint
  onco_plot <- ggplot(data = naive, aes(x = Patient_ID, y = timepoint)) +
    geom_tile(fill="gray90", colour="white", size = 1.1) +
    geom_tile(data = naive[naive$positive == "positive", ], aes(x = Patient_ID, y = timepoint), inherit.aes=FALSE, width = 0.7, height = 0.9, fill = "#93c3ad") +
    geom_tile(data = informed, aes(x = Patient_ID, y = timepoint), inherit.aes=FALSE, width = 0.7, height = 0.3, fill = "#23505f") +
    theme(
      text = element_text(family = "Arial"),
      axis.line =         element_blank(),
      axis.text.x =       element_text(size = base_size * 1, lineheight = 0.9, hjust = 1, vjust = 0.5, angle = 90),
      axis.text.y =       element_blank(),
      axis.title.x =      element_blank(),
      axis.title.y =      element_blank(),
      axis.ticks.length = unit(0.0, "lines"),
      panel.background =  element_rect(fill = "transparent", colour = NA), 
      panel.border =      element_blank(),
      panel.grid.major =  element_blank(),
      panel.grid.minor =  element_blank(),
      strip.background =  element_blank(), 
      strip.text.x =      element_blank(),
      strip.text.y =      element_blank(),
      strip.text =        element_blank(),
      plot.background =   element_rect(colour = NA, fill = "transparent"),
      plot.title =        element_text(size = base_size * 0.95),
      legend.position = "none"
    )
  
  # Summarize
  list <- list()
  list[["volume"]] <- volume_plot
  list[["oncoprint"]] <- onco_plot
  
  # Return <list>
  return(list)
  
}
