#   Script name:  ctCPA_score.R
#   
#   Author:       Florian Janke
#   Last update:  16th July 2024
#
#-------------------------------

ctDNA_informed_CPA <- function(input, timepoints, data){
  
  #----
  # Functions
  #----
  intersect_custom <- function(table){
    
    table$intersect <- NA
    
    COND01 <- between(table$i.start, table$start, table$i.start)
    COND02 <- table$i.end >= table$end
    COND03 <- table$i.end < table$end
    
    table$intersect <- ifelse(COND01,
                              ifelse(COND02, table$end - table$i.start, table$i.end - table$i.start),
                              ifelse(COND02, table$end - table$start,
                                     ifelse(COND03, table$i.end - table$start, table$intersect)))
    
    table$intersect <- table$intersect +1
    return(table)
    
    
  }
  duplicate_removal <- function(table){
    
    # Identify duplicated segments
    duplicates <- data.frame(table(table$id))
    
    # Subset for duplicated segments
    duplicates <- as.character(duplicates[duplicates$Freq > 1, ][[1]])
    
    # Returns table if no duplicates are present
    if(rlang::is_empty(duplicates)) return(table)
    
    # Store non-duplicates for merging
    non_duplicates <- table[!table$id %in% duplicates, ]
    
    # Reduce duplicated segments to the one with the largest overlap
    for(n in 1:length(duplicates)){
      
      tmp <- table[table$id == duplicates[n], ]
      tmp <- tmp[order(tmp$intersect, decreasing = TRUE), ]
      tmp <- tmp[1, ]
      
      if(n == 1) tmp_sum <- tmp
      if(n != 1) tmp_sum <- rbind(tmp_sum, tmp)
      if(n == length(duplicates)) table <- rbind(non_duplicates, tmp_sum)
      
    }
    
    # Return <table>
    return(table)
    
  }
  ctDNA_informed_CPA_calc <- function(input, timepoints = "Baseline|5-Fx|10-Fx", data){
    
    # Defines which longituidnal samples should be included for ctCPA score calculation
    if(timepoints == "Baseline") tp <- 1
    if(timepoints == "Baseline|5-Fx") tp <- c(1:2)
    if(timepoints == "Baseline|5-Fx|10-Fx") tp <- c(1:3)
    
    # Extract patients with at least one positive sample before the end of RT
    samp <- data[data$timepoint %in% unlist(strsplit(timepoints, "|", fixed = TRUE)) & data$positive == "positive", ]
    samp <- samp[!duplicated(samp$Patient_ID), ]$Patient_ID
    
    # Calculate ctDNA-informed CPA-scores
    for(i in 1:length(samp)){
      
      # Load copy number states of time points 'Baseline', '5-Fx' and '10-Fx'
      states <- list.files(file.path(input), pattern = samp[i], full.names = TRUE)
      states <- states[as.numeric(str_sub(basename(states), 6, 7)) %in% tp]
      
      for(n in 1:length(states)){
        
        # Load file
        tmp <- data.frame(fread(states[n]))
        
        # Subset relevant columns
        tmp <- tmp[, c(1:4, 7)]
        
        # Add column names
        colnames(tmp) <- c("id", "chr", "start", "end", basename(states[n]))
        
        # Append
        if(n == 1) state <- tmp
        if(n != 1) state <- merge(state, tmp[, c(1, 5)], by = "id", all = TRUE)
        if(n == length(states)) rm(states, tmp, n)
        
      }
      
      # Exclude sex chromosomes
      state <- state[!state$chr %in% c("X", "Y"), ]
      
      # Find consensus copy number state; i.e., exclude conflicting regions (gain & loss), neutral regions and empty regions
      if(timepoints == "Baseline"){
        
        state$neutral <- ifelse(state[[5]] == "neutral", 1, 0)
        state$loss <- ifelse(state[[5]] == "loss", 1, 0)
        state$gain <- ifelse(state[[5]] == "gain", 1, 0)
        
      } else {
        
        state$neutral <- rowSums(state[, c(5:ncol(state))] == "neutral", na.rm = TRUE)
        state$loss <- rowSums(state[, c(5:ncol(state))] == "loss", na.rm = TRUE)
        state$gain <- rowSums(state[, c(5:ncol(state))] == "gain", na.rm = TRUE)
        
      }
      state <- state[state$neutral != length(tp), ]
      state <- state[rowSums(state[, c((ncol(state) -2):ncol(state))]) != 0, ]
      state$neutral <- NULL
      state$loss <- ifelse(state$loss == 0, 0, -1)
      state$gain <- ifelse(state$gain == 0, 0, 1)
      state$state <- rowSums(state[, c((ncol(state) -1):ncol(state))])
      state <- state[state$state != 0, ]
      state <- state[, c(2:4, ncol(state))]
      state <- data.table(state)
      
      # Calculate ctDNA-informed CPA-scores
      files <- list.files(input, full.names = TRUE)
      files <- files[grep(paste(c(samp[i], "C0"), collapse = "|"), files)]
      for(n in 1:length(files)){
        
        # Load copy number state file
        tmp <- fread(files[n])
        
        # Re-format
        tmp <- tmp[, c(1:4, 6)]
        colnames(tmp) <- c("id", "chr", "start", "end", "zscore")
        
        # Set keys for merging
        setkey(tmp, chr, start, end)
        setkey(state, chr, start, end)
        
        # Merge
        tmp <- foverlaps(state, tmp)
        
        # Calculate overlap (in bp)
        tmp <- intersect_custom(tmp)
        
        # Exclude regions w/o z-score in plasma file
        tmp <- tmp[!is.na(tmp$zscore), ]
        
        # Reduce duplicates to the one with the largest overlap
        tmp <- duplicate_removal(table = tmp)
        
        # Multiple z-score by directionality
        tmp$zscore_direction <- tmp$zscore *tmp$state
        
        # Sum directionality-corrected z-scores
        score <- sum(tmp$zscore_direction) / sum(tmp$intersect) *10^8
        
        # Summarize
        tmp <- data.frame(Sample_ID = basename(files[n]),
                          Patient_ID = basename(samp[i]),
                          group = ifelse(str_sub(basename(files[n]), 1, 1) == "C", "control", "case"),
                          score = score)
        
        # Append
        if(n == 1) informed_scores <- tmp
        if(n != 1) informed_scores <- rbind(informed_scores, tmp)
        
      }
      
      # Define positivity
      informed_scores$positive <- ifelse(informed_scores$score >max(informed_scores[informed_scores$group == "control", ]$score), "positive", "negative")
      
      # Adjust Sample IDs
      informed_scores$Sample_ID <- gsub("_bins.bed", "", informed_scores$Sample_ID)
      
      # Append
      if(i == 1) scores <- informed_scores
      if(i != 1) scores <- rbind(scores, informed_scores)
      if(i == length(samp)) rm(informed_scores, files, n, i, input, tmp, score, samp, state)
      
    }
    
    # Return <scores>
    return(scores)
    
  }
  
  #----
  # Main script
  #----
  for(i in 1:length(timepoints)){
    
    # Initiate list for data storage
    if(i == 1) scores <- list()
    
    # Calculate CNV scores
    tmp <- ctDNA_informed_CPA_calc(input = input, data = data, timepoints = timepoints[i])
    
    # Add time point information
    tmp <- merge(tmp, data[, c("Sample_ID", "timepoint", "date", "time")], by = "Sample_ID", all.x = TRUE)
    
    # Add 'Patient_ID'
    colnames(tmp)[2] <- "informed_by"
    tmp$Patient_ID <- str_sub(tmp$Sample_ID, 1, 4)
    
    # Re-order
    colnames(tmp)[4] <- "CPA"
    tmp <- tmp[, c(colnames(data), "informed_by")]
    
    # Store
    scores[[timepoints[i]]] <- tmp
    
    # Add <data> to <scores>
    if(i == length(timepoints)) scores[["naive"]] <- data
    
  }
  
  # Return <scores>
  return(scores)
  
}




