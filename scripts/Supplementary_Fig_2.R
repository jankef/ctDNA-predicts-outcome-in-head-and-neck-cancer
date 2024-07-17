#   Script name:  Swimmer_plot.R
#   
#   Author:       Florian Janke
#   Last update:  16th July 2024
#
#
#   Description:  Function used to draw the swimmer plot illustrated in Supplementary Fig. 2.
#
#
#-------------------------------

swimmer <- function(input, data){
  
  # Load info
  info <- readxl::read_excel(input, sheet = "Characteristics")
  
  # Load therapy info
  ctx <- readxl::read_excel(input, sheet = "Chemotherapy")
  ctx$Start <- as.numeric(ctx$Start)
  ctx$End <- as.numeric(ctx$End)
  io <- readxl::read_excel(input, sheet = "Immunotherapy")
  io$Start <- as.numeric(io$Start)
  io$End <- as.numeric(io$End)
  combi <- readxl::read_excel(input, sheet = "Combination")
  combi$Start <- as.numeric(combi$Start)
  combi$End <- as.numeric(combi$End)
  surgery <- readxl::read_excel(input, sheet = "Surgery")
  
  # Extract 2nd progression
  second_pd <- info[!is.na(info$PFS2), c("Patient_ID", "PFS2")]
  
  # Extract time point of death
  death <- info[info$OS_event == 1 & info$OS <12, c("Patient_ID", "OS")]
  colnames(death)[2] <- "time_start"
  death$time_end <- death$time_start +0.09
  
  # Extract relevant information from <info>
  info <- data.frame(Patient_ID = info$Patient_ID,
                     reRT_start = 0,
                     reRT_end = as.numeric(as.Date(info$RT_end, format = "%m/%d/%Y") -as.Date(info$RT_start, format = "%m/%d/%Y")) /30.5,
                     type = info$Modality,
                     PD_start = info$PFS,
                     PD_end = info$PFS +0.03,
                     PD_info = info$PD_info)
  
  # Order <info>
  info <- info[order(info$PD_start, decreasing = FALSE), ]
  info <- info[order(info$type, decreasing = TRUE), ]
  info$Patient_ID <- factor(info$Patient_ID, levels = info$Patient_ID)
  

  
  # Extract sampling time points
  swimmer_data <- data[str_sub(data$Sample_ID, 1, 4) %in% info$Patient_ID, ]
  
  # Plot
  base_size <- 14
  plot <- ggplot() +
    geom_segment(data = info, aes(x = reRT_start, xend = reRT_end, y = Patient_ID, yend = Patient_ID, color = type), linewidth = 4) +
    geom_segment(data = io, aes(x = Start, xend = End, y = Patient_ID, yend = Patient_ID), color = "#AEDFF2", linewidth = 4) +
    geom_segment(data = combi, aes(x = Start, xend = End, y = Patient_ID, yend = Patient_ID), color = "#F2D7D7", linewidth = 4) +
    geom_segment(data = ctx, aes(x = Start, xend = End, y = Patient_ID, yend = Patient_ID), color = "gray90", linewidth = 4) +
    geom_segment(data = info, aes(x = 0, xend = PD_end, y = Patient_ID, yend = Patient_ID), linewidth = 0.25, color = "black") +
    geom_segment(data = second_pd, aes(x = PFS2, xend = PFS2+0.03, y = Patient_ID, yend = Patient_ID), linewidth = 4) +
    geom_segment(data = death, aes(x = time_start, xend = time_end, y = Patient_ID, yend = Patient_ID), linewidth = 4, color = "red") +
    geom_segment(data = info, aes(x = PD_start, xend = PD_end, y = Patient_ID, yend = Patient_ID), linewidth = 4) +
    scale_color_manual(values = c("#e3ce97", "#93c3ad")) +
    scale_x_continuous(limits = c(-0.1, 12), breaks = seq(0, 12, 3)) +
    scale_y_discrete(labels = info$Patient_ID) +
    geom_point(data = swimmer_data, aes(x = time /30.5, y = Patient_ID), size = 2, shape = 21, fill = "white") +
    geom_point(data = surgery, aes(x = Start, y = Patient_ID), shape = 4, size = 2) +
    theme(
      text = element_text(family = "Arial"),
      axis.line =         element_line(colour = "black"),
      axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
      axis.text.y =       element_text(size = base_size * 0.6, lineheight = 0.9, hjust = 1),
      axis.title.x =      element_blank(),
      axis.title.y =      element_blank(),
      axis.ticks.length = unit(0.3, "lines"),
      panel.background =  element_rect(fill = "transparent", colour = NA), 
      panel.border =      element_rect(fill = NA, colour = "transparent", linewidth = 0.75), 
      panel.grid.major.y =  element_line(colour = "grey90"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.y =  element_line(colour = "grey90"),
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
  
  
  # Return <plot>
  return(plot)
  
}