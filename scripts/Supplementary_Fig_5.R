#   Script name:  Supplementary_Fig_5.R
#   
#   Author:       Florian Janke
#   Last update:  17th July 2024
#
#-------------------------------

ctCPA_changes <- function(timepoint, ctCPA){
  
  #----
  # Functions
  #----
  rounding <- function(x){
    if(x < 0.001) y <- sprintf("%.2e", x) else y <- round(x, 3)
    return(y)
  }
  pairwise <- function(ctCPA, timepoint){
    
    # Extract ctCPA scores
    table <- ctCPA$`Baseline|5-Fx|10-Fx`[, c(1:8)]
    
    # Extract relevant time points
    timepoint <- unlist(str_split(string = timepoint, pattern = "\\|"))
    for(i in 1:length(timepoint)){
      
      # Select relevant time points
      if(timepoint[i] %in% c("RT-end", "Baseline", "5-Fx", "10-Fx")){
        
        tmp <- table[table$group == "case", ]
        tmp <- tmp[tmp$timepoint == timepoint[i], ]
        
      }
      if(timepoint[i] == "PD"){
        
        # Subset <info> to include PFS and Patient_ID
        pd <- info[, c("Patient_ID", "PFS")]
        
        # Merge PFS information with sampling time points
        pd <- merge(table[, c("Patient_ID", "Sample_ID", "time")], pd, by = "Patient_ID")
        
        # Get time points clostest to PD (but not after PD)
        pd$delta <- (pd$time /30.5) - pd$PFS
        pd$delta <- abs(pd$delta)
        pd.tmp <- pd %>% dplyr::group_by(Patient_ID) %>% dplyr::summarize(delta = min(delta))
        pd <- dplyr::left_join(pd.tmp, pd, by = c("Patient_ID", "delta"))
        
        # Subset <tmp> by closest time point to PD
        tmp <- dplyr::left_join(pd[, c("Patient_ID", "time")], table, by = c("Patient_ID", "time"))
        
        # Re-name 'timepoint' entries
        tmp$timepoint <- "PD"
        
        # Adjust to 2nd progression for patient 'CARE_006'
        tmp[tmp$Patient_ID == "P006", ]$CPA <- table[table$Sample_ID == "P006_05-P1", ]$CPA
        
      }
      
      # Append
      if(i == 1) tab <- tmp
      if(i != 1){
        
        diff <- merge(tab[, c("Patient_ID", "CPA")], tmp[, c("Patient_ID", "CPA")], by = "Patient_ID")
        diff <- data.frame(Patient_ID = diff$Patient_ID, change = round(diff$CPA.y /diff$CPA.x -1, 3))
        tab <- rbind(tab, tmp)
        
      }
      
    }
    
    # Calculate statistics
    stat <- data.frame(Group_A = timepoint[1],
                       Group_B = timepoint[2],
                       MedianCPA_A = paste0(round(median(tab[tab$timepoint == timepoint[1], ]$CPA), 3), " (", round(min(tab[tab$timepoint == timepoint[1], ]$CPA), 3), " - ", round(max(tab[tab$timepoint == timepoint[1], ]$CPA), 3), ")"),
                       MedianCPA_B = paste0(round(median(tab[tab$timepoint == timepoint[2], ]$CPA), 3), " (", round(min(tab[tab$timepoint == timepoint[2], ]$CPA), 3), " - ", round(max(tab[tab$timepoint == timepoint[2], ]$CPA), 3), ")"),
                       MedianChange = paste0(round(median(diff$change), 3), " (", round(min(diff$change), 3), " - ", round(max(diff$change), 3), ")"),
                       Decreases = paste0(nrow(diff[diff$change <= 0, ]), "/", nrow(diff)),
                       Wilcox_paired_p = rounding(wilcox.test(tab[tab$timepoint == timepoint[1], ]$CPA, tab[tab$timepoint == timepoint[2], ]$CPA, paired = TRUE)$p.value))
    
    # Factorize
    if("PD" %in% timepoint) tab$timepoint <- factor(tab$timepoint, levels = c("RT-end", "PD")) else tab$timepoint <- factor(tab$timepoint, levels = c("Baseline", "RT-end"))
    
    # Draw plot
    base_size <- 14
    plot <- ggplot() +
      stat_boxplot(data = tab, aes(x = timepoint, y = CPA), geom = "errorbar", width = 0.15, lwd = 0.3) +
      geom_boxplot(data = tab, aes(x = timepoint, y = CPA), outlier.shape = NA, color = "black") +
      geom_line(data = tab, aes(x = timepoint, y = CPA, group = Patient_ID), color = "gray90") +
      geom_point(data = tab, aes(x = timepoint, y = CPA), fill = "black", size = 2, shape = 21) +
      scale_y_continuous(trans = "log10", limits = c(0.2, 13), name = "ctCPA score", breaks = c(0.2, 0.5, 1, 2, 5, 10)) +
      scale_x_discrete(labels = ifelse(timepoint == "RT-end", "re-RT end", timepoint)) +
      annotation_logticks(sides = "l") +
      theme(
        text = element_text(family = "Arial"),
        axis.line =         element_line(colour = "black"),
        axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
        axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
        axis.title.x =      element_blank(),
        axis.title.y =      element_text(size = base_size * 1, angle = 90, vjust = 2),
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
    
    # Summarize
    list <- list()
    list[["plot"]][[paste0(timepoint[1], "|", timepoint[2])]] <- plot
    list[["statistics"]] <- stat
    
    # Return <list>
    return(list)
    
  }
  
  #----
  # Main script
  #----
  list <- list()
  for(i in 1:length(timepoint)){
    
    # Get plot and statistics
    tmp <- pairwise(ctCPA = ctCPA, timepoint = timepoint[i])
    
    # Append
    list[["plot"]][[timepoint[i]]] <- tmp$plot
    if(i == 1) list[["statistics"]] <- tmp$statistics
    if(i != 1) list[["statistics"]] <- rbind(list[["statistics"]], tmp$statistics)
    
  }
  
  # Return <list>
  return(list)
  
}
