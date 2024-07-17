#   Script name:  Fig_3.R
#   
#   Author:       Florian Janke
#   Last update:  17th July 2024
#
#-------------------------------

survival_summary <- function(table, ctCPA, timepoint){
  
  #----
  # Functions
  #----
  survival <- function(table, comp, metric, type){
    
    if(type == "univariate"){
      
      # Subset <info>
      if(metric == "OS") tmp <- table[, c(metric, comp, paste0(metric, "_event"))]
      if(metric == "PFS"){
        
        tmp <- table[, c(metric, comp)]
        tmp$event <- 1
        
      }
      
      # Exclude NAs
      tmp <- tmp[!is.na(tmp[2]), ]
      
      # Re-name columns
      colnames(tmp) <- c("time", "condition", "event")
      
      # Adjust 'condition' column to separate conditions in 1 and 0 
      if(comp == "Age") {tmp$condition <- ifelse(tmp$condition >= 59, "old", "young"); tmp$condition <- factor(tmp$condition, levels = c("young", "old"))}
      if(comp == "GTV_ccm"){
        
        tmp$condition <- ifelse(tmp$condition >= quantile(info$GTV_ccm, 0.55), "large", "small")
        tmp$condition <- factor(tmp$condition, levels = c("small", "large"))
        
      }
      if(comp == "Sex") tmp$condition <- factor(tmp$condition, levels = c("female", "male"))
      if(comp == "Histology") tmp$condition <- factor(tmp$condition, levels = c("SCC", "AC"))
      if(comp == "Modality") tmp$condition <- factor(tmp$condition, levels = c("CIRT", "VMAT"))
      if(comp == "Smoking_status") tmp$condition <- factor(tmp$condition, levels = c("none-smoker", "current-former"))
      if(comp == "PD_info") tmp$condition <- factor(tmp$condition, levels = c("local", "distant"))
      if(!comp %in% c("Smoking_status", "Modality", "Histology", "Sex", "GTV_ccm", "Age")) tmp$condition <- as.factor(tmp$condition)
      
    }
    
    if(type == "multivariable"){
      
      if(metric == "OS") tmp <- info[, c(metric, unlist(strsplit(comp, "|", fixed = TRUE)), paste0(type, "_event"))]
      if(metric == "PFS") tmp <- info[, c(metric, unlist(strsplit(comp, "|", fixed = TRUE)))]; tmp$event <- 1
      
      # Re-name columns
      colnames(tmp)[c(1, (ncol(tmp)-1):ncol(tmp))] <- c("time", "condition", "event")
      
      # Factorize variables
      if("Age" %in% unlist(strsplit(comp, "|", fixed = TRUE))) {tmp$Age <- ifelse(tmp$Age >= 59, "old", "young"); tmp$Age <- factor(tmp$Age, levels = c("young", "old"))}
      if("GTV_ccm" %in% unlist(strsplit(comp, "|", fixed = TRUE))){
        
        tmp$GTV_ccm <- ifelse(tmp$GTV_ccm >= quantile(info$GTV_ccm, 0.5), "large", "small")
        tmp$GTV_ccm <- factor(tmp$GTV_ccm, levels = c("small", "large"))
        
      }
      if("Sex" %in% unlist(strsplit(comp, "|", fixed = TRUE))) tmp$Sex <- factor(tmp$Sex, levels = c("female", "male"))
      if("Histology" %in% unlist(strsplit(comp, "|", fixed = TRUE))) tmp$Histology <- factor(tmp$Histology, levels = c("SCC", "AC"))
      if("reRT_type" %in% unlist(strsplit(comp, "|", fixed = TRUE))) tmp$reRT_type <- factor(tmp$reRT_type, levels = c("CIRT", "VMAT"))
      if("smoking_status" %in% unlist(strsplit(comp, "|", fixed = TRUE))) tmp$smoking_status <- factor(tmp$smoking_status, levels = c("no", "yes"))
      if("metastasis" %in% unlist(strsplit(comp, "|", fixed = TRUE))) tmp$metastasis <- factor(tmp$metastasis, levels = c("local", "distant"))
      tmp$condition <- as.factor(tmp$condition)
      
    }
    
    # Exclude NAs
    tmp <- tmp[!is.na(tmp$condition), ]
    
    # Summarize statistics
    if(type == "univariate"){
      
      # Get log-rank p-value and group name
      log.rank <- survdiff(Surv(time, event) ~ condition, data = tmp)
      group <- paste0(gsub("condition=", "", names(log.rank$n[1])), "_V_", gsub("condition=", "", names(log.rank$n[2])))
      if(group == "0_V_1") group <- "noMRD_V_MRD"
      log.rank <- rounding(log.rank$pvalue)
      
      # Get median survival and range
      surv.median.01 <- paste0(round(median(tmp[tmp$condition == levels(tmp$condition)[2], ]$time), 2),
                               " (",
                               round(min(tmp[tmp$condition == levels(tmp$condition)[2], ]$time), 2),
                               " - ",
                               round(max(tmp[tmp$condition == levels(tmp$condition)[2], ]$time), 2),
                               ")")
      
      surv.median.02 <- paste0(round(median(tmp[tmp$condition == levels(tmp$condition)[1], ]$time), 2),
                               " (",
                               round(min(tmp[tmp$condition == levels(tmp$condition)[1], ]$time), 2),
                               " - ",
                               round(max(tmp[tmp$condition == levels(tmp$condition)[1], ]$time), 2),
                               ")")
      
      # Get hazard ratio
      cox <- coxph(Surv(time, event) ~ condition, data = tmp)
      cox <- summary(cox)
      
      # Summarize
      stat <- data.frame(survival_metric = metric,
                         analysis_type = type,
                         variable = comp,
                         group = group,
                         log_rank = log.rank,
                         HR = rounding(cox$conf.int[1]),
                         HR_low = rounding(cox$conf.int[3]),
                         HR_high = rounding(cox$conf.int[4]),
                         median_survival_01 = surv.median.01,
                         median_survival_02 = surv.median.02, 
                         n_01 = nrow(tmp[tmp$condition == levels(tmp$condition)[2], ]),
                         n_02 = nrow(tmp[tmp$condition == levels(tmp$condition)[1], ]))
      
    }
    
    if(type == "multivariable"){
      
      # Get hazard ratio
      cox <- coxph(Surv(time, event) ~ condition + GTV_ccm + metastasis + reRT_type + Age, data = tmp)
      cox <- summary(cox)
      
      # Extract statistics
      stat <- data.frame(cox$conf.int)[, c(1, 3:4)]
      stat <- cbind(stat, data.frame(cox$coefficients)[5])
      stat$Parameter <- rownames(stat)
      stat <- stat[, c(5, 1:4)]
      colnames(stat) <- c("Parameter", "HR", "HR_low", "HR_high", "P")
      
    }
    
    # Plot Kaplan-Meier curve
    if(type == "univariate"){
      
      # Parameters
      x.limit <- ifelse(metric == "OS", 36, 12)
      x.steps <- ifelse(metric == "OS", 6, 2)
      censoring <- ifelse(all(tmp$event == 1), FALSE, TRUE)
      y.label <- ifelse(metric == "PFS", "Progression-free survival (%)", "Overall survival (%)")
      
      # Plotting
      s.fit <- survfit(Surv(time, event) ~ condition, data = tmp)
      base_size <- 14
      plot <- ggsurv(s = s.fit, plot.cens = censoring) +
        scale_color_manual(values = c("black", "gray65")) +
        scale_y_continuous(labels = seq(0, 100, 25), limits = c(0, 1), name = y.label) +
        scale_x_continuous(limits = c(0, x.limit), breaks = seq(0, x.limit, x.steps), name = "Months since re-RT start") +
        theme(
          text = element_text(family = "Arial"),
          axis.line =         element_line(colour = "black"),
          axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 0.5, angle = 0, hjust = 0.5),
          axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
          axis.title.x =      element_text(size = base_size * 0.8, angle = 0, vjust = 0, hjust = 0.5),
          axis.title.y =      element_text(size = base_size * 0.8, angle = 90, vjust = 2),
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
          plot.title =        element_text(size = base_size * 1.5),
          legend.position = "none",
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.title = element_text(size = 12),
          legend.key.size = unit(0.55, 'cm'),
          legend.text = element_text(size=10)
        )
      
      
      # Get individuals at risk
      for(i in seq(0, x.limit, x.steps)){
        
        # Number at risk
        if(i >max(tmp[tmp$condition == levels(tmp$condition)[1], ]$time)) no.risk.A <- 0 else no.risk.A <- nrow(tmp[tmp$condition == levels(tmp$condition)[1], ][tmp[tmp$condition == levels(tmp$condition)[1], ]$time >i, ])
        if(i >max(tmp[tmp$condition == levels(tmp$condition)[2], ]$time)) no.risk.B <- 0 else no.risk.B <- nrow(tmp[tmp$condition == levels(tmp$condition)[2], ][tmp[tmp$condition == levels(tmp$condition)[2], ]$time >i, ])
        
        # Summarize
        risk.tmp <- data.frame(Time = i,
                               Group_A = no.risk.A,
                               Group_B = no.risk.B)
        
        # Append
        if(i == 0) no.risk <- risk.tmp
        if(i != 0) no.risk <- rbind(no.risk, risk.tmp)
        
      }
      
      
    }
    
    # Summarize
    list <- list()
    list[["statistics"]] <- stat
    list[["plot"]] <- plot
    list[["no-risk"]] <- no.risk
    
    
    # Return <list>
    return(list)
    
  }
  rounding <- function(x){
    if(x < 0.001) y <- sprintf("%.2e", x) else y <- round(x, 3)
    return(y)
  }
  
  
  #----
  # Main script
  #----
  
  # Run survival analysis
  list <- list()
  for(i in 1:length(timepoint)){
    
    tmp <- survival(table = info,
                    comp = timepoint[i],
                    metric = "PFS",
                    type = "univariate")
    
    # Append
    if(i == 1) list[["statistics"]][[timepoint[i]]] <- tmp$statistics
    if(i != 1) list[["statistics"]][[timepoint[i]]] <- rbind(list[["statistics"]][[timepoint[i]]], tmp$statistics)
    list[["KM_plot"]][[timepoint[i]]] <- tmp$plot
    list[["Number_at_risk"]][[timepoint[i]]] <- tmp$`no-risk`
    
  }
  
  # Return <list>
  return(list)
  
}






