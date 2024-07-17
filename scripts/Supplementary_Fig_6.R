#   Script name:  Supplementary_Fig_6.R
#   
#   Author:       Florian Janke
#   Last update:  17th July 2024
#
#-------------------------------

kinetics <- function(patID, ctCPA){
  
  #----
  # Functions
  #----
  kinetic_plot <- function(ctCPA, patient, type){
    
    # Extract (ct)CPA scores
    if(type == "informed") tmp <- ctCPA$`Baseline|5-Fx|10-Fx`
    if(type == "naive") tmp <- ctCPA$naive
    
    # Retrieve detectability threshold
    if(type == "informed") cut <- max(tmp[tmp$informed_by == patient & tmp$group == "control", ]$CPA)
    if(type == "naive") cut <- max(tmp[tmp$group == "control", ]$CPA)
    
    # Subset patient specified by <patID>
    tmp <- tmp[tmp$Patient_ID == patient, ]
    
    # Exit function if there are no ctCPA scores for the given patient
    if(length(tmp$time) == 0) return()
    
    # Transform time to months
    tmp$time <- tmp$time /30.5
    
    # Set colors for plotting
    colors <- c("white", "black")
    if(all(tmp$positive == "positive")) colors <- "black"
    if(all(tmp$positive == "negative")) colors <- "white"
    if(info[info$Patient_ID == patient, ]$Modality == "CIRT") rt_color <- "#e3ce97" else rt_color <- "#93c3ad"
    
    # Adjust very low values for better visibility
    if(any(tmp$CPA < 0.2)) tmp[tmp$CPA < 0.2, ]$CPA <- 0.2
    
    # Plot kinetic
    base_size <- 14
    plot <- ggplot() +
      geom_rect(aes(xmin = 0, xmax = tmp[tmp$timepoint == "RT-end", ]$time, ymin = 0.2, ymax = Inf), fill = rt_color) +
      geom_line(data = tmp, aes(x = time, y = CPA), linetype = "dashed") +
      geom_point(data = tmp, aes(x = time, y = CPA, fill = positive), size = 3, shape = 21) +
      scale_fill_manual(values = colors) +
      scale_x_continuous(limits = c(-0.1, 12), breaks = seq(0, 12, 3), name = "Months since re-RT start") +
      scale_y_continuous(trans = "log10", limits = c(0.2, 13), name = ifelse(type == "informed", "ctCPA score", "CPA score"), breaks = c(0.2, 0.5, 1, 2, 5, 10)) +
      annotation_logticks(sides = "l") +
      geom_hline(yintercept = cut, linetype = "dashed", alpha = 0.5) +
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
    
    # Return <plot>
    return(plot)
    
  }
  
  #----
  # Main script
  #----
  
  # Get kinetics
  list <- list()
  for(i in 1:length(patID)){
    
    # Set <type>
    if(patID[i] %in% c("P016", "P001", "P003", "P006", "P007", "P008", "P018")) type <- "informed" else type <- "naive"
    
    # Draw kinetic
    plot <- kinetic_plot(ctCPA = ctCPA, patient = patID[i], type = type)
    
    # Append
    list[[patID[i]]] <- plot
    
  }
  
  # Return <list>
  return(list)
  
}
