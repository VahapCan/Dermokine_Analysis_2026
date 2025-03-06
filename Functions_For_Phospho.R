
reformat_fasta_df <- function(fasta_df) {
  reformatted_lines <- c()
  current_sequence <- ""
  
  for (i in 1:nrow(fasta_df)) {
    line <- as.character(fasta_df[i, 1])
    if (!is.na(line)) {
      if (startsWith(line, ">")) {
        if (nchar(current_sequence) > 0) {
          reformatted_lines <- c(reformatted_lines, current_sequence)
          current_sequence <- ""
        }
        reformatted_lines <- c(reformatted_lines, line)
      } else {
        current_sequence <- paste0(current_sequence, line)
      }
    }
  }} 


fasta_to_named_list <- function(fasta_lines) {
  named_list <- list()
  current_header <- ""
  
  for (line in fasta_lines) {
    if (startsWith(line, ">")) {
      uniprot_id <- sub("^.*\\|(.+?)\\|.*$", "\\1", line)
      gene_id <- sub("^.*GN=(\\S+).*?$", "\\1", line)
      current_header <- paste(uniprot_id, gene_id, sep = ";")
    } else {
      named_list[[current_header]] <- line
    }
  }
  
  return(named_list)
}

extract_flanking_seq <- function(named_FASTA, rownames_ppe) {
  updated_rownames <- c()
  
  for (rowname in rownames_ppe) {
    split_name <- strsplit(rowname, ";")[[1]]
    
    # Remove the first and third terms if there are more than 4 terms
    if (length(split_name) > 4) {
      split_name <- split_name[-c(1, 3)]
    }
    # Remove any terms after the third semicolon
    split_name <- split_name[1:4]
    
    rowname <- remove_terms_after_third_semicolon(rowname)
    gene_id <- split_name[2]
    position <- as.integer(gsub("\\D", "", split_name[3]))
    matches <- grep(paste0(".*;", gene_id, "$"), names(named_FASTA), value = TRUE)
    
    if (length(matches) > 0) {
      protein_seq <- named_FASTA[matches[1]]
      
      # Correct the position indexing
      start <- max(1, position - 7)
      stop <- min(nchar(protein_seq), position + 7)
      
      flanking_seq <- substr(protein_seq, start, stop)
      new_name <- paste(split_name[1], split_name[2], split_name[3], flanking_seq, sep = ";")
      updated_rownames <- c(updated_rownames, new_name)
    } else {
      updated_rownames <- c(updated_rownames, rowname)
    }
  }
  
  return(updated_rownames)
}

remove_terms_after_third_semicolon <- function(rowname) {
  split_name <- strsplit(rowname, ";")[[1]]
  new_name <- paste(split_name[1:4], collapse = ";")
  return(new_name)
}


generate_volcano_plot <- function(df) {
  # Extract the name of the data frame
  df_name <- deparse(substitute(df))
  
  # Define the title based on the name of the data frame
  title <- switch(df_name,
                  "results_df_CS_Rescue_sorted" = "DMKN\U03B1\U03B2-/- Rescue vs DMKN\U03B1\U03B2-/- Endogenous",
                  "results_df_CS_Endogenous_sorted" = "WT vs DMKN\U03B1\U03B2-/-",
                  "results_df_DTU_Rescue_sorted" = "DMKN\U03B2\U03B3-/- Rescue vs DMKN\U03B2\U03B3-/- Endogenous",
                  "results_df_DTU_Endogenous_sorted" = "WT vs DMKN\U03B2\U03B3-/-",
                  "results_df_WT_Endogenous_sorted" = "WT Rescue vs WT endogenous",
                  "Unknown data frame")
  
  log2FC_threshold = log2(1.5)
  point_color = "black"
  hl.col.down = "blue"
  hl.col.up = "red"
  # Generate the volcano plot
  volcano_plot <- ggplot(df, aes(x = diff_coeff, y = -log10(phos_pvals_significant), label = phos_ID)) +
    geom_vline(xintercept = log2FC_threshold, linetype = 3) + geom_vline(xintercept = -log2FC_threshold, linetype = 3) +
    geom_point(aes(color = ifelse(diff_coeff >= log2FC_threshold, hl.col.up, ifelse(diff_coeff <= -log2FC_threshold, hl.col.down, point_color)))) +
    scale_alpha_identity(guide=FALSE) +
    scale_color_identity(guide = FALSE) +
    ylim(0, 3) + xlim(-3, 3) +
    labs(x = "Difference in Coefficients (Log2 fold change)", y = "-Log10(limma moderated p value)") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'transparent'),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          text = element_text(size = 12) ,
          legend.position = "none") +
    geom_label_repel(aes(label = phos_ID, fill = ifelse(startsWith(phos_ID, "CTNND1"), "orange", "black")), 
                     nudge_y = 0.25, nudge_x = .15, segment.curvature = -1e-20,
                     arrow = arrow(length = unit(0.015, "npc")),
                     box.padding   = 0.35, max.overlaps = Inf,
                     point.padding = 0.5,
                     segment.color = 'grey50', size  = 2)+
    ggtitle(paste("Difference in Log2 Fold Change of \nPhosphorylated Peptide and Protein \nAbundance", title))
  
  return(volcano_plot)
}


volcano_plot_ggplot_multicol_overlay <- function(fit, coefs = 1L, style = "p-value", 
                                                 highlight = 0L, names = rownames(fit$coef), 
                                                 hl.col.up = "red", hl.col.down = "blue", 
                                                 point_color = "darkgray", xlab = "Log2 Fold Change",
                                                 ylab = NULL, pch = 16, cex = 0.5, 
                                                 p_value_threshold = 0.05, log2FC_threshold = 1, ...) {
  if (!is(fit, "MArrayLM")) 
    stop("fit must be an MArrayLM")
  
  style <- match.arg(tolower(style), c("p-value", "b-statistic"))
  
  data_list <- list()
  
  for (coef in coefs) {
    x <- as.matrix(fit$coef)[, coef]
    
    if (style == "p-value") {
      if (is.null(fit$p.value)) 
        stop("No p-values found in linear model fit object. Please run eBayes.")
      y <- as.matrix(fit$p.value)[, coef]
      y <- -log10(y)
      if (is.null(ylab)) 
        ylab = "-log10(limma moderated p value)"
    }
    else {
      if (is.null(fit$lods)) 
        stop("No B-statistics found in linear model fit object. Please run eBayes.")
      y <- as.matrix(fit$lods)[, coef]
      if (is.null(ylab)) 
        ylab = "Log Odds of Differential Expression"
    }
    
    if (is.null(names)) {
      names <- 1:length(x)
    }
    
    coef_data <- data.frame(x = x, y = y, names = names, coef = factor(coef))
    data_list[[as.character(coef)]] <- coef_data
  }
  
  # ... (rest of your function) ...
  
  data <- do.call(rbind, data_list)
  
  p <- ggplot(data, aes(x, y)) + labs(x = xlab, y = ylab)
  
  # Manually specify each comparison
  coef1_data <- data[data$coef == coefs[1],]
  coef2_data <- data[data$coef == coefs[2],]
  
  # Determine labels based on the coefficients provided
  labels <- c()
  for (coef in coefs) {
    coef_name <- colnames(fit$coef)[coef]
    labels <- c(labels, paste0(coef_name, " up"))
  }
  labels <- c(labels, "Not significant")
  for (coef in coefs) {
    coef_name <- colnames(fit$coef)[coef]
    labels <- c(labels, paste0(coef_name, " down"))
  }
  
  
  
  p <- p + geom_point(data = coef1_data, aes(color = ifelse(x >= log2FC_threshold, "blue", 
                                                            ifelse(x <= -log2FC_threshold, "darkred",
                                                                   point_color)), alpha = ifelse(abs(x) >= log2FC_threshold & y >=-log10(p_value_threshold), 1, 0.01)), 
                      shape = pch, size = cex)
  p <- p + geom_point(data = coef2_data, aes(color = ifelse(x >= log2FC_threshold, "cyan", 
                                                            ifelse(x <= -log2FC_threshold, "pink", 
                                                                   point_color)), alpha = ifelse(abs(x) >= log2FC_threshold & y >= -log10(p_value_threshold), 1, 0.01)), 
                      shape = 18, size = cex)
  # Getting 2 colors from each palette
  red_spectrum <- brewer.pal(9, "YlOrRd")[c(8, 5)]  # Choosing darker and somewhat lighter reds
  blue_spectrum <- brewer.pal(9, "GnBu")[c(8, 5)]  # Choosing darker and somewhat lighter blues
  
  
  # ... (rest of your function) ...
  color_order <- c(red_spectrum[1], red_spectrum[2], "darkgrey", blue_spectrum[1], blue_spectrum[2])
  shape_order  <- c("pch", 18, "pch", "pch", 18)
  p = p + scale_alpha_identity(guide="none") + scale_color_manual(values = color_order, labels =labels,  name = "Comparisons")  + scale_shape_manual(values = shape_order, labels = labels) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent'),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),legend.text=element_text(size=10),
          text = element_text(size = 12),legend.position="bottom") + 
    guides(colour = guide_legend(ncol = 2,title.position="top", title.hjust = 0.5)) + ggtitle("Comparison of phosphorylation \nsite abundances") 
  nudge_directions <- c(0.3, -0.3)  # replace these with your preferred values
  
  
  print(p)
  invisible()
  return(p)
}

volcano_plot_ggplot_multicol_overlay_CS <- function(fit, coefs = 1L, style = "p-value", 
                                                    highlight = 0L, names = rownames(fit$coef), 
                                                    hl.col.up = "red", hl.col.down = "blue", 
                                                    point_color = "darkgray", xlab = "Log2 Fold Change",
                                                    ylab = NULL, pch = 16, cex = 0.5, 
                                                    p_value_threshold = 0.05, log2FC_threshold = 1, ...) {
  if (!is(fit, "MArrayLM")) 
    stop("fit must be an MArrayLM")
  
  style <- match.arg(tolower(style), c("p-value", "b-statistic"))
  
  data_list <- list()
  
  for (coef in coefs) {
    x <- as.matrix(fit$coef)[, coef]
    
    if (style == "p-value") {
      if (is.null(fit$p.value)) 
        stop("No p-values found in linear model fit object. Please run eBayes.")
      y <- as.matrix(fit$p.value)[, coef]
      y <- -log10(y)
      if (is.null(ylab)) 
        ylab = "-log10(limma moderated p value)"
    }
    else {
      if (is.null(fit$lods)) 
        stop("No B-statistics found in linear model fit object. Please run eBayes.")
      y <- as.matrix(fit$lods)[, coef]
      if (is.null(ylab)) 
        ylab = "Log Odds of Differential Expression"
    }
    
    if (is.null(names)) {
      names <- 1:length(x)
    }
    
    coef_data <- data.frame(x = x, y = y, names = names, coef = factor(coef))
    data_list[[as.character(coef)]] <- coef_data
  }
  
  # ... (rest of your function) ...
  
  data <- do.call(rbind, data_list)
  
  p <- ggplot(data, aes(x, y)) + labs(x = xlab, y = ylab)
  
  # Manually specify each comparison
  coef1_data <- data[data$coef == coefs[1],]
  coef2_data <- data[data$coef == coefs[2],]
  
  # Determine labels based on the coefficients provided
  labels <- c()
  for (coef in coefs) {
    coef_name <- colnames(fit$coef)[coef]
    labels <- c(labels, paste0(coef_name, " up"))
  }
  labels <- c(labels, "Not significant")
  for (coef in coefs) {
    coef_name <- colnames(fit$coef)[coef]
    labels <- c(labels, paste0(coef_name, " down"))
  }
  
  
  
  p <- p + geom_point(data = coef1_data, aes(color = ifelse(x >= log2FC_threshold, "blue", 
                                                            ifelse(x <= -log2FC_threshold, "darkred",
                                                                   point_color)), alpha = ifelse(abs(x) >= log2FC_threshold & y >=-log10(p_value_threshold), 1, 0.01)), 
                      shape = pch, size = cex)
  p <- p + geom_point(data = coef2_data, aes(color = ifelse(x >= log2FC_threshold, "cyan", 
                                                            ifelse(x <= -log2FC_threshold, "pink", 
                                                                   point_color)), alpha = ifelse(abs(x) >= log2FC_threshold & y >= -log10(p_value_threshold), 1, 0.01)), 
                      shape = 18, size = cex)
  # Getting 2 colors from each palette
  red_spectrum <- brewer.pal(9, "YlOrRd")[c(8, 5)]  # Choosing darker and somewhat lighter reds
  blue_spectrum <- brewer.pal(9, "GnBu")[c(8, 5)]  # Choosing darker and somewhat lighter blues
  
  
  # ... (rest of your function) ...
  color_order <- c(red_spectrum[1],"darkgrey", blue_spectrum[1])
  shape_order  <- c("pch", "pch", "pch")
  p = p + scale_alpha_identity(guide="none") + scale_color_manual(values = color_order, labels =labels,  name = "Comparisons")  + scale_shape_manual(values = shape_order, labels = labels) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent'),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),legend.text=element_text(size=10),
          text = element_text(size = 12),legend.position="bottom") + 
    guides(colour = guide_legend(ncol = 2,title.position="top", title.hjust = 0.5)) + ggtitle("Comparison of endogenous phosphorylation site \nabundances") 
  nudge_directions <- c(0.3, -0.3)  # replace these with your preferred values
  
  
  print(p)
  invisible()
  return(p)
}

volcano_plot_ggplot_multicol_overlay_DTU <- function(fit, coefs = 1L, style = "p-value", 
                                                     highlight = 0L, names = rownames(fit$coef), 
                                                     hl.col.up = "red", hl.col.down = "blue", 
                                                     point_color = "darkgray", xlab = "Log2 Fold Change",
                                                     ylab = NULL, pch = 16, cex = 0.5, 
                                                     p_value_threshold = 0.05, log2FC_threshold = 1, ...) {
  if (!is(fit, "MArrayLM")) 
    stop("fit must be an MArrayLM")
  
  style <- match.arg(tolower(style), c("p-value", "b-statistic"))
  
  data_list <- list()
  
  for (coef in coefs) {
    x <- as.matrix(fit$coef)[, coef]
    
    if (style == "p-value") {
      if (is.null(fit$p.value)) 
        stop("No p-values found in linear model fit object. Please run eBayes.")
      y <- as.matrix(fit$p.value)[, coef]
      y <- -log10(y)
      if (is.null(ylab)) 
        ylab = "-log10(limma moderated p value)"
    }
    else {
      if (is.null(fit$lods)) 
        stop("No B-statistics found in linear model fit object. Please run eBayes.")
      y <- as.matrix(fit$lods)[, coef]
      if (is.null(ylab)) 
        ylab = "Log Odds of Differential Expression"
    }
    
    if (is.null(names)) {
      names <- 1:length(x)
    }
    
    coef_data <- data.frame(x = x, y = y, names = names, coef = factor(coef))
    data_list[[as.character(coef)]] <- coef_data
  }
  
  # ... (rest of your function) ...
  
  data <- do.call(rbind, data_list)
  
  p <- ggplot(data, aes(x, y)) + labs(x = xlab, y = ylab)
  
  # Manually specify each comparison
  coef1_data <- data[data$coef == coefs[1],]
  coef2_data <- data[data$coef == coefs[2],]
  
  # Determine labels based on the coefficients provided
  labels <- c()
  for (coef in coefs) {
    coef_name <- colnames(fit$coef)[coef]
    labels <- c(labels, paste0(coef_name, " up"))
  }
  labels <- c(labels, "Not significant")
  for (coef in coefs) {
    coef_name <- colnames(fit$coef)[coef]
    labels <- c(labels, paste0(coef_name, " down"))
  }
  
  
  
  p <- p + geom_point(data = coef1_data, aes(color = ifelse(x >= log2FC_threshold, "blue", 
                                                            ifelse(x <= -log2FC_threshold, "darkred",
                                                                   point_color)), alpha = ifelse(abs(x) >= log2FC_threshold & y >=-log10(p_value_threshold), 1, 0.01)), 
                      shape = pch, size = cex)
  p <- p + geom_point(data = coef2_data, aes(color = ifelse(x >= log2FC_threshold, "cyan", 
                                                            ifelse(x <= -log2FC_threshold, "pink", 
                                                                   point_color)), alpha = ifelse(abs(x) >= log2FC_threshold & y >= -log10(p_value_threshold), 1, 0.01)), 
                      shape = 18, size = cex)
  # Getting 2 colors from each palette
  red_spectrum <- brewer.pal(9, "YlOrRd")[c(8, 5)]  # Choosing darker and somewhat lighter reds
  blue_spectrum <- brewer.pal(9, "GnBu")[c(8, 5)]  # Choosing darker and somewhat lighter blues
  
  
  # ... (rest of your function) ...
  color_order <- c(red_spectrum[2], "darkgrey", blue_spectrum[2])
  shape_order  <- c(18, "pch", 18)
  p = p + scale_alpha_identity(guide="none") + scale_color_manual(values = color_order, labels =labels,  name = "Comparisons")  + scale_shape_manual(values = shape_order, labels = labels) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent'),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),legend.text=element_text(size=10),
          text = element_text(size = 12),legend.position="bottom") + 
    guides(colour = guide_legend(ncol = 2,title.position="top", title.hjust = 0.5)) + ggtitle("Comparison of phosphorylation \nsite abundances") 
  nudge_directions <- c(0.3, -0.3)  # replace these with your preferred values
  
  
  print(p)
  invisible()
  return(p)
}
volcano_plot_ggplot_multicol <- function(fit, coefs = 1L, style = "p-value", 
                                         highlight = 0L, names = rownames(fit$coef), 
                                         hl.col.up = "red", hl.col.down = "blue", 
                                         point_color = "darkgray", xlab = "Log2 Fold Change",
                                         ylab = NULL, pch = 16, cex = 0.5, facet = "grid", 
                                         p_value_threshold = 0.05, log2FC_threshold = 1, ...) {
  if (!is(fit, "MArrayLM")) 
    stop("fit must be an MArrayLM")
  
  style <- match.arg(tolower(style), c("p-value", "b-statistic"))
  
  data_list <- list()
  
  for (coef in coefs) {
    x <- as.matrix(fit$coef)[, coef]
    
    if (style == "p-value") {
      if (is.null(fit$p.value)) 
        stop("No p-values found in linear model fit object. Please run eBayes.")
      y <- as.matrix(fit$p.value)[, coef]
      y <- -log10(y)
      if (is.null(ylab)) 
        ylab = "-log10(limma moderated p value)"
    }
    else {
      if (is.null(fit$lods)) 
        stop("No B-statistics found in linear model fit object. Please run eBayes.")
      y <- as.matrix(fit$lods)[, coef]
      if (is.null(ylab)) 
        ylab = "Log Odds of Differential Expression"
    }
    
    if (is.null(names)) {
      names <- 1:length(x)
    }
    
    coef_data <- data.frame(x = x, y = y, names = names, coef = factor(coef))
    data_list[[as.character(coef)]] <- coef_data
  }
  
  data <- do.call(rbind, data_list)
  
  p <- ggplot(data, aes(x, y)) + 
    geom_point(aes(color = ifelse(x >= log2FC_threshold, hl.col.up, ifelse(x <= -log2FC_threshold, hl.col.down, point_color)), alpha = ifelse(abs(x) >= log2FC_threshold & y >= -log10(p_value_threshold), 1, 0.01)),
               shape = pch, size = cex) +
    scale_color_identity(guide = FALSE) + scale_alpha_identity(guide=FALSE) +
    #theme_classic() + 
    labs(x = xlab, y = ylab) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent'),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          text = element_text(size = 12))
  
  # ... (rest of the function remains the same) ...
  
  if (highlight > 0) {
    for (coef in coefs) {
      coef_data <- data[data$coef == coef,]
      o <- order(coef_data$y, decreasing = TRUE)
      i <- o[1:highlight]
      p <- p + ggrepel::geom_text_repel(data = coef_data[i,], aes(label = names, 
                                                                  color = ifelse(x >= log2FC_threshold, hl.col.up, ifelse(x <= -log2FC_threshold, hl.col.down, point_color)),
                                                                  alpha = ifelse(abs(x) >= log2FC_threshold & y >= -log10(p_value_threshold), 1, 0.5)),
                                        size = 2, nudge_y = 0.3)
    }
  }
  panel_labels <- vector(mode = "character", length = length(coefs))
  names(panel_labels) <- as.character(coefs)
  
  for (i in seq_along(coefs)) {
    panel_labels[as.character(coefs[i])] <- colnames(fit$coef)[coefs[i]]
  }
  
  
  p <- p + facet_wrap(. ~ coef, ncol =2,labeller = as_labeller(panel_labels))
  
  
  print(p)
  invisible()
  return(p)
}

volcano_plot_ggplot_multicol_phosSites <- function(fit, coefs = 1L, style = "p-value", 
                                                   highlight = 20, names = rownames(fit$coef), 
                                                   hl.col.up = "red", hl.col.down = "blue", 
                                                   point_color = "black", xlab = "Log2 (WT vs KO) ratio",
                                                   ylab = NULL, pch = 16, cex = 0.5, facet = "grid", 
                                                   p_value_threshold = 0.05, log2FC_threshold = 1,title, ...) {
  if (!is(fit, "MArrayLM")) 
    stop("fit must be an MArrayLM")
  
  style <- match.arg(tolower(style), c("p-value", "b-statistic"))
  
  data_list <- list()
  
  for (coef in coefs) {
    x <- as.matrix(fit$coef)[, coef]
    
    if (style == "p-value") {
      if (is.null(fit$p.value)) 
        stop("No p-values found in linear model fit object. Please run eBayes.")
      y <- as.matrix(fit$p.value)[, coef]
      y <- -log10(y)
      if (is.null(ylab)) 
        ylab = "-log10(limma moderated p value)"
    }
    else {
      if (is.null(fit$lods)) 
        stop("No B-statistics found in linear model fit object. Please run eBayes.")
      y <- as.matrix(fit$lods)[, coef]
      if (is.null(ylab)) 
        ylab = "Log Odds of Differential Expression"
    }
    
    if (is.null(names)) {
      names <- 1:length(x)
    }
    
    coef_data <- data.frame(x = x, y = y, names = names, coef = factor(coef))
    data_list[[as.character(coef)]] <- coef_data
  }
  
  data <- do.call(rbind, data_list)
  
  p <- ggplot(data, aes(x, y,label=phos_Site       )) + geom_vline(xintercept=log2(1.5),linetype=3) +
    geom_point(aes(color = ifelse(x >= log2FC_threshold, hl.col.up, ifelse(x <= -log2FC_threshold, hl.col.down, point_color)), alpha = ifelse(abs(x) >= log2FC_threshold & y >= -log10(p_value_threshold), 1, 1)),
               shape = pch, size = cex) +
    scale_color_identity(guide = FALSE) + scale_alpha_identity(guide=FALSE) + ylim(0,5) + xlim(-1.0,2)+
    #theme_classic() + 
    labs(x = xlab, y = ylab) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent'),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          text = element_text(size = 12))
  
  # ... (rest of the function remains the same) ...
  
  p <- p + geom_label_repel(aes(label = phos_Site), nudge_y = 0.25,nudge_x = .15,segment.curvature = -1e-20,
                            arrow = arrow(length = unit(0.015, "npc")),
                            box.padding   = 0.35, max.overlaps = Inf,
                            point.padding = 0.5,
                            segment.color = 'grey50', size  = 2) + ggtitle(title)
  
  #print(p)
  invisible()
  return(p)
}

combine_GOplots_with_shared_legend <- function(plot1, plot2) {
  # Remove the legend from both plots
  plot1_no_legend <- plot1 + theme(legend.position = "none")
  plot2_no_legend <- plot2 + theme(legend.position = "none")
  
  # Combine the plots using ggarrange without legends
  combined_plots <- ggarrange(plot1_no_legend, plot2_no_legend, ncol = 2)
  
  # Extract the legend from the first plot
  plot1_legend <- get_legend(plot1)
  plot2_legend <- get_legend(plot2)
  
  # Get the number of legend elements in both plots
  n_legend_elements_plot1 <- max(plot1$data$n)
  n_legend_elements_plot2 <- max(plot2$data$n)
  # Choose the legend with more elements
  selected_legend <- if (n_legend_elements_plot1 >= n_legend_elements_plot2) plot1_legend else plot2_legend
  
  # Combine the plots and the shared legend using the cowplot package
  final_plot <- cowplot::plot_grid(combined_plots, selected_legend, ncol = 1, rel_heights = c(1, 0.4))
  
  return(final_plot)
}

combine_plots_with_shared_legend <- function(plot1, plot2) {
  # Remove the legend from both plots
  plot1_no_legend <- plot1 + theme(legend.position = "none")
  plot2_no_legend <- plot2 + theme(legend.position = "none")
  
  # Combine the plots using ggarrange without legends
  combined_plots <- ggarrange(plot1_no_legend, plot2_no_legend, ncol = 2)
  
  # Extract the legend from the first plot
  plot1_legend <- get_legend(plot1)
  plot2_legend <- get_legend(plot2)
  
  # Get the number of legend elements in both plots
  n_legend_elements_plot1 <- max(plot1$data$n)
  n_legend_elements_plot2 <- max(plot2$data$n)
  # Choose the legend with more elements
  selected_legend <- if (n_legend_elements_plot1 >= n_legend_elements_plot2) plot1_legend else plot2_legend
  
  # Combine the plots and the shared legend using the cowplot package
  final_plot <- cowplot::plot_grid(combined_plots, selected_legend, ncol = 1, rel_heights = c(0.9, 0.1))
  
  return(final_plot)
}
process_signalomes_results <- function(Signalomes_results_DTU_Up) {
  
  # Convert the character vector into a data frame
  signalomes_results_DTU_Up_df <- data.frame(gene_and_phos_site = names(Signalomes_results_DTU_Up$proteinModules),
                                             module_number = as.numeric(Signalomes_results_DTU_Up$proteinModules),
                                             stringsAsFactors = FALSE)
  
  # Create a list to store the new data frames
  module_dfs <- list()
  
  # Iterate through the unique module numbers and create a data frame for each module
  for (module_number in unique(signalomes_results_DTU_Up_df$module_number)) {
    module_subset <- signalomes_results_DTU_Up_df[signalomes_results_DTU_Up_df$module_number == module_number, "gene_and_phos_site"]
    module_dfs[[paste0("module_", module_number)]] <- data.frame(gene_and_phos_site = module_subset)
  }
  
  convert_to_entrez <- function(gene_and_phos_site) {
    # Extract gene names
    gene_names <- gsub(" .*", "", gene_and_phos_site)
    
    # Convert gene names to Entrez IDs
    entrez_ids <- mapIds(org.Hs.eg.db,
                         keys = gene_names,
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")
    
    # Remove NA values
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    
    return(entrez_ids)
  }
  
  # Create a new list to store the Entrez IDs as character vectors
  module_entrez_ids_chr <- list()
  
  # Iterate through each module and convert gene names to Entrez IDs
  for (module_name in names(module_dfs)) {
    module_entrez_ids_chr[[module_name]] <- as.character(convert_to_entrez(module_dfs[[module_name]]$gene_and_phos_site))
  }
  
  return(module_entrez_ids_chr)
}
#Remove irregular sequences
short_fourth_term <- function(row_name) {
  terms <- unlist(strsplit(row_name, ";"))
  if (length(terms) >= 4) {
    return (nchar(terms[4]) < 4)
  }
  return (FALSE)
}
remove_nonunique_rownames_matrix <- function(mat) {
  # Convert the matrix to a data frame
  df <- as.data.frame(mat)
  
  # Create a new column 'rownames' with the current row names
  df$rownames <- rownames(mat)
  
  # Count the occurrences of each row name
  rowname_counts <- table(df$rownames)
  
  # Find the row names that occur only once
  unique_rownames <- names(rowname_counts[rowname_counts == 1])
  
  # Keep only the rows with unique row names
  df_unique <- df[df$rownames %in% unique_rownames, ]
  
  # Reset the row names
  rownames(df_unique) <- df_unique$rownames
  
  # Remove the 'rownames' column
  df_unique$rownames <- NULL
  
  # Convert the data frame back to a matrix
  mat_unique <- as.matrix(df_unique)
  
  return(mat_unique)
}

kinaseSubstrateHeatmap_modified = function (phosScoringMatrices, top = 3, printPlot = NULL, filePath = "./kinaseSubstrateHeatmap.pdf", 
                                            width = 10, height = 10) 
{
  utils::data("KinaseFamily", envir = environment())
  sites <- c()
  for (i in seq_len(ncol(phosScoringMatrices$combinedScoreMatrix))) {
    sites <- union(sites, names(sort(phosScoringMatrices$combinedScoreMatrix[, 
                                                                             i], decreasing = TRUE)[seq_len(top)]))
  }
  o <- intersect(colnames(phosScoringMatrices$combinedScoreMatrix), 
                 rownames(KinaseFamily))
  annotation_col = data.frame(group = KinaseFamily[o, "kinase_group"], 
                              family = KinaseFamily[o, "kinase_family"])
  rownames(annotation_col) <- o
  if (is.null(printPlot) == TRUE) {
    
    pheatmap(phosScoringMatrices$combinedScoreMatrix[sites, 
    ], annotation_col = annotation_col, cluster_rows = TRUE, 
    cluster_cols = TRUE, fontsize = 7)
    
  }
  else {
    pdf(file = filePath, width = width, height = height)
    
    pheatmap(phosScoringMatrices$combinedScoreMatrix[sites, 
    ], annotation_col = annotation_col, cluster_rows = TRUE, 
    cluster_cols = TRUE, fontsize = 7)
    
    dev.off()
  }
}
siteAnnotate_modified <- function(site, phosScoringMatrices, predMatrix) {
  od <- order(predMatrix[site, ], decreasing = FALSE)
  kinases <- colnames(predMatrix)[od]
  
  par(mfrow = c(1, 4), mar = c(1, 5, 5, 1))
  
  barplot(predMatrix[site, kinases], las = 1, xlab = "Prediction score", 
          col = "red3", main = "Prediction score", xlim = c(0, 1), horiz = TRUE, xaxt= "n")
  axis(3)
  barplot(phosScoringMatrices$combinedScoreMatrix[site, kinases], 
          las = 1, main = "Combined score", col = "orange2", xlim = c(0, 1), horiz = TRUE,xaxt= "n")
  axis(3)
  barplot(phosScoringMatrices$motifScoreMatrix[site, kinases], 
          las = 1, main = "Motif score", col = "green4", xlim = c(0, 1), horiz = TRUE,xaxt= "n")
  axis(3)
  barplot(phosScoringMatrices$profileScoreMatrix[site, kinases], 
          las = 1, main = "Profile score", col = "lightblue3", 
          xlim = c(0, 1), horiz = TRUE,xaxt= "n")
  axis(3)
  mtext(paste("Site =", site),                   # Add main title
        side = 3,
        line = - 2,
        outer = TRUE)
}

siteAnnotate_modified_ggplot <- function(site, phosScoringMatrices, predMatrix) {
  od <- order(predMatrix[site, ], decreasing = FALSE)
  kinases <- colnames(predMatrix)[od]
  
  plot_data <- list(predMatrix[site, kinases], phosScoringMatrices$combinedScoreMatrix[site, kinases],
                    phosScoringMatrices$motifScoreMatrix[site, kinases], phosScoringMatrices$profileScoreMatrix[site, kinases])
  plot_types <- c("Prediction score", "Combined score", "Motif score", "Profile score")
  plot_colors <- c("red3", "orange2", "green4", "lightblue3")
  
  plots <- list()
  
  for (i in 1:4) {
    df <- data.frame(Kinase = factor(kinases, levels = kinases), Value = plot_data[[i]])
    
    p <- ggplot(df, aes(x = Kinase, y = Value, fill = factor(plot_types[i]))) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = plot_colors[i], guide = "none") +
      theme_classic() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(size = 7),
            axis.title.x = element_text(size = 0),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x.top = element_text(angle = 45, hjust = 1)) +
      labs(title = plot_types[i], x = NULL) +
      scale_x_discrete(position = "top") +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0)), limits = c(0, 1),breaks = c(0,1))
    
    plots[[i]] <- p
  }
  title1=text_grob(paste(site), size = 9, face = "bold")   #### this worked for me
  
  grid.arrange(grobs = plots, ncol = 4, top = title1)
}


plotSignalomeMap_modified <- function(signalomes, color,module_colors, max_percent = 100,size =10) {
  df <- stack(signalomes$kinaseSubstrates)
  modules <- signalomes$proteinModule
  names(modules) <- unlist(lapply(strsplit(as.character(names(signalomes$proteinModules)), ";"), "[[", 1))
  df$cluster <- modules[df$values]
  df_balloon <- df
  df_balloon <- na.omit(df_balloon) %>% dplyr::count(.data$cluster, .data$ind)
  df_balloon$ind <- as.factor(df_balloon$ind)
  df_balloon$cluster <- as.factor(df_balloon$cluster)
  df_balloon <- tidyr::spread(df_balloon, .data$ind, .data$n)[, -1]
  df_balloon[is.na(df_balloon)] <- 0
  df_balloon <- do.call(rbind, lapply(seq(nrow(df_balloon)), function(x) {
    res <- unlist(lapply(df_balloon[x, ], function(y) y/sum(df_balloon[x, ]) * 100))
  }))
  df_balloon <- reshape2::melt(as.matrix(df_balloon))
  colnames(df_balloon) <- c("cluster", "ind", "n")
  df_balloon <- df_balloon[df_balloon$n > 0, ]
  
  # Normalize the n column based on the max_percent value
  max_percentage <- max(df_balloon$n)
  normalized_n <- (df_balloon$n / max_percentage) * max_percent
  df_balloon$normalized_n <- normalized_n
  
  # Function to find the last balloon for each group
  last_balloon_y <- function(group_data) {
    max(group_data$y[group_data$x == max(group_data$x)])
  }
  
  # Add a new column to the dataframe containing the last balloon's y value for each group
  df_balloon$last_y <- sapply(unique(df_balloon$cluster), function(cluster) {
    last_balloon_y(df_balloon[df_balloon$cluster == cluster, ])
  })[df_balloon$cluster]
  
  
  # Function to find the first balloon for each group
  first_balloon_y <- function(group_data) {
    min(group_data$y[group_data$x == min(group_data$x)])
  }
  
  # Add a new column to the dataframe containing the first balloon's y value for each group
  df_balloon$first_y <- sapply(unique(df_balloon$cluster), function(cluster) {
    first_balloon_y(df_balloon[df_balloon$cluster == cluster, ])
  })[df_balloon$cluster]  
  
  g <- ggplot2::ggplot(df_balloon, aes(y = .data$ind, x = .data$cluster)) + 
    
    
    #geom_text(aes(label = .data$cluster, x = .data$cluster, y = 0), angle = 90, hjust = 1, size = 3) +  # Show x-axis labels using geom_text
    
    theme_classic() + 
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_text(angle = 0),
          plot.margin = margin(, 0.5, , , "cm"),aspect.ratio = 1,
          axis.line = element_blank(), 
          axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 14),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank()) +
    guides(col = FALSE)  # Remove the legend for the kinases
  # Add rectangles and text labels for module names on the x-axis
  
  
  g <- g + geom_line(aes(group = .data$cluster, col = .data$ind), size = 0.5, linetype = "solid", color = "grey") +  # Connect the balloons within the horizontal line
    geom_segment(data = df_balloon, aes(x = .data$cluster, xend = .data$cluster, y = .data$first_y, yend = .data$last_y, col = .data$ind), size = 0.5, linetype = "solid", show.legend = FALSE, color = "grey") 
  
  rectangle_height <- 30
  # Add rectangles for module names on the x-axis
  for (i in unique(df_balloon$cluster)) {
    g <- g + annotate("rect", xmin = as.numeric(i) - 0.45, xmax = as.numeric(i) + 0.45, ymin = -1.5, ymax = 0, alpha = 0.7, fill = module_colors[i]) + annotate("text", x = i, y = -0.025 * rectangle_height, label = i, size = 4, fontface = "bold")
  }
  g= g + scale_color_manual(values = color) +
    scale_size_continuous(name = "% of phosphosites regulated by a kinase", range = c(1, size),
                          guide = guide_legend(title.position = "top", title.hjust = 0.5),
                          limits = c(0, max_percent))
  g = g + geom_point(aes(col = .data$ind, size = .data$normalized_n))
  
  #Set the limits using the y_axis_limits argument
  #g <- g + scale_y_discrete(limits = y_axis_limits)
  # Use geom_segment() to draw lines from the x-axis to the last balloon for each group
  #g <- g + 
  g
}

plot_kinase_activity_NoSerum <- function(input_data, group,title) {
  if (group == "DTU") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("WT_NoSerum|DTU", colnames(input_data))]
    label <- "DMKN \U03B2\U03B3-/-"
  } else if (group == "CS") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("WT_NoSerum|CS", colnames(input_data))]
    label <- "DMKN \U03B1\U03B2-/-"
  } else {
    stop("Invalid group. Please specify either 'DTU' or 'CS'")
  }
  
  long_data <- melt(filtered_data, variable.name = "Condition", value.name = "Score")
  
  long_data$Genotype <- ifelse(str_detect(long_data$Var2, group), label, "WT")
  long_data$Treatment <- ifelse(str_detect(long_data$Var2, "Treated"), "Treated", "Control")
  long_data$Genotype_Treatment <- paste(long_data$Genotype, long_data$Treatment, sep = " ")
  
  # Define the list of comparisons based on the input group
  if (group == "DTU") {
    my_comparisons <- list(c("WT Control", "WT Treated"),
                           c("DMKN \U03B2\U03B3-/- Control", "DMKN \U03B2\U03B3-/- Treated"),
                           c("WT Control","DMKN \U03B2\U03B3-/- Control"),
                           c("WT Treated","DMKN \U03B2\U03B3-/- Treated"),
                           c("WT Control", "DMKN \U03B2\U03B3-/- Treated")
    )
    
    long_data$Genotype_Treatment <- factor(long_data$Genotype_Treatment,
                                           levels = c("WT Control", "WT Treated",
                                                      "DMKN \U03B2\U03B3-/- Control", 
                                                      "DMKN \U03B2\U03B3-/- Treated"))
  } else if (group == "CS") {
    my_comparisons <- list(c("WT Control", "WT Treated"),
                           c("DMKN \U03B1\U03B2-/- Control", "DMKN \U03B1\U03B2-/- Treated"),
                           c("WT Control","DMKN \U03B1\U03B2-/- Control"),
                           c("WT Treated","DMKN \U03B1\U03B2-/- Treated"),
                           c("WT Control", "DMKN \U03B1\U03B2-/- Treated")
    )
    
    long_data$Genotype_Treatment <- factor(long_data$Genotype_Treatment,
                                           levels = c("WT Control", "WT Treated", 
                                                      "DMKN \U03B1\U03B2-/- Control", 
                                                      "DMKN \U03B1\U03B2-/- Treated"))
  } else {
    stop("Invalid group specified. Please use either 'DTU' or 'CS'.")
  }
  
  plot <- ggplot(long_data, aes(x = Genotype_Treatment, y = Score, fill = Genotype)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.7) + 
    #geom_errorbar(aes(ymin = Score - sd(Score), ymax = Score + sd(Score)), position = position_dodge(0.8), width = 0.3) +
    facet_wrap(~Var1, scales = "free_x", ncol = 5) +
    labs(x = "Comparison", y = "Kinase Activity Score") +
    theme_bw()+ ylim(-3,4.5)+ theme(
      legend.position = "bottom", 
      axis.text.x = element_text(size = 3),
      axis.text.y = element_text(size = 7), 
      axis.title = element_text(size=9, face = "bold"),
      strip.text = element_text(size = 8),
      plot.margin = margin(, 0.5, , , "cm"),
      strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")))+
    stat_compare_means(comparisons = my_comparisons, tip.length = 0.01, vjust = 0.55,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("****", "***", "**", "*", "ns")), 
                       label = "p.signif", method = "t.test", label.y = max(long_data$Score) * 1.0,
                       step.increase = 0.08, hide.ns	= TRUE) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
  plot = plot + ggtitle(title)
  return(plot)
}
plot_kinase_activity_NoSerum_all_combos <- function(input_data, group) {
  if (group == "DTU") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("", colnames(input_data))]
    labelA <- "DMKN \U03B2\U03B3-/-"
    labelB <- "DMKN \U03B1\U03B2-/-"
  } else if (group == "CS") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("", colnames(input_data))]
    labelA <- "DMKN \U03B2\U03B3-/-"
    labelB <- "DMKN \U03B1\U03B2-/-"
  } else {
    stop("Invalid group. Please specify either 'DTU' or 'CS'")
  }
  
  long_data <- melt(filtered_data, variable.name = "Condition", value.name = "Score")
  
  long_data$Genotype <- ifelse(str_detect(long_data$Var2, "DTU"), "DMKN ????-/-",
                               ifelse(str_detect(long_data$Var2, "CS"), "DMKN \U03B1\U03B2-/-", "WT"))
  long_data$Treatment <- ifelse(str_detect(long_data$Var2, "Treated"), "Treated", "Control")
  long_data$Genotype_Treatment <- paste(long_data$Genotype, long_data$Treatment, sep = " ")
  
  # Define the list of comparisons based on the input group
  if (group == "DTU") {
    my_comparisons <- list(c("WT Control", "WT Treated"),
                           c("DMKN \U03B2\U03B3-/- Control", "DMKN \U03B2\U03B3-/- Treated"),
                           c("DMKN \U03B1\U03B2-/- Control", "DMKN \U03B1\U03B2-/- Treated"),
                           c("DMKN \U03B1\U03B2-/- Control", "DMKN \U03B2\U03B3-/- Control"),
                           c("DMKN \U03B1\U03B2-/- Treated", "DMKN \U03B2\U03B3-/- Treated"),
                           c("WT Control","DMKN \U03B2\U03B3-/- Control"),
                           c("WT Treated","DMKN \U03B1\U03B2-/- Treated"),
                           c("WT Treated","DMKN \U03B2\U03B3-/- Treated"),
                           c("WT Control", "DMKN \U03B2\U03B3-/- Treated"),
                           c("WT Control","DMKN \U03B1\U03B2-/- Treated")
    )
    
    long_data$Genotype_Treatment <- factor(long_data$Genotype_Treatment,
                                           levels = c("WT Control", "WT Treated",
                                                      "DMKN \U03B2\U03B3-/- Control", 
                                                      "DMKN \U03B2\U03B3-/- Treated",
                                                      "DMKN \U03B1\U03B2-/- Control", 
                                                      "DMKN \U03B1\U03B2-/- Treated"))
  } else if (group == "CS") {
    my_comparisons <- list(c("WT Control", "WT Treated"),
                           c("DMKN \U03B2\U03B3-/- Control", "DMKN \U03B2\U03B3-/- Treated"),
                           c("DMKN \U03B1\U03B2-/- Control", "DMKN \U03B1\U03B2-/- Treated"),
                           c("DMKN \U03B1\U03B2-/- Control", "DMKN \U03B2\U03B3-/- Control"),
                           c("DMKN \U03B1\U03B2-/- Treated", "DMKN \U03B2\U03B3-/- Treated"),
                           c("WT Control","DMKN \U03B2\U03B3-/- Control"),
                           c("WT Treated","DMKN \U03B1\U03B2-/- Treated"),
                           c("WT Treated","DMKN \U03B2\U03B3-/- Treated"),
                           c("WT Control", "DMKN \U03B2\U03B3-/- Treated"),
                           c("WT Control","DMKN \U03B1\U03B2-/- Treated")
    )
    
    long_data$Genotype_Treatment <- factor(long_data$Genotype_Treatment,
                                           levels = c("WT Control", "WT Treated",
                                                      "DMKN \U03B2\U03B3-/- Control", 
                                                      "DMKN \U03B2\U03B3-/- Treated",
                                                      "DMKN \U03B1\U03B2-/- Control", 
                                                      "DMKN \U03B1\U03B2-/- Treated"))
  } else {
    stop("Invalid group specified. Please use either 'DTU' or 'CS'.")
  }
  
  plot <- ggplot(long_data, aes(x = Genotype_Treatment, y = Score, fill = Genotype)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.7) + 
    #geom_errorbar(aes(ymin = Score - sd(Score), ymax = Score + sd(Score)), position = position_dodge(0.8), width = 0.3) +
    facet_wrap(~Var1, scales = "free_x", ncol = 5) +
    labs(x = "Comparison", y = "Kinase Activity Score") +
    theme_bw()+ ylim(-3,4.5)+ theme(
      legend.position = "bottom", 
      axis.text.x = element_text(size = 3),
      axis.text.y = element_text(size = 7), 
      axis.title = element_text(size=9, face = "bold"),
      strip.text = element_text(size = 8),
      plot.margin = margin(, 0.5, , , "cm"),
      strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")))+
    stat_compare_means(comparisons = my_comparisons, tip.length = 0.01, vjust = 0.55,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("****", "***", "**", "*", "ns")), 
                       label = "p.signif", method = "t.test", label.y = max(long_data$Score) * 1.0,
                       step.increase = 0.08, hide.ns	= TRUE) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
  
  return(plot)
}
plot_kinase_activity_NoSerum_Genotypes <- function(input_data, group, title) {
  if (group == "DTU") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("_control|_Control", colnames(input_data))]
    label1 <- "DMKN \U03B2\U03B3-/-"
  } else if (group == "CS") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("_control|_Control", colnames(input_data))]
    label2 <- "DMKN \U03B1\U03B2-/-"
  } else {
    stop("Invalid group. Please specify either 'DTU' or 'CS'")
  }
  
  long_data <- melt(filtered_data, variable.name = "Condition", value.name = "Score")
  
  long_data$Genotype <- ifelse(str_detect(long_data$Var2, "DTU"), "DMKN ????-/-",
                               ifelse(str_detect(long_data$Var2, "CS"), "DMKN \U03B1\U03B2-/-", "WT"))
  long_data$Treatment <- ifelse(str_detect(long_data$Var2, "Treated"), "Treated", "Control")
  long_data$Genotype_Treatment <- paste(long_data$Genotype, long_data$Treatment, sep = " ")
  
  # Define the list of comparisons based on the input group
  if (group == "DTU") {
    my_comparisons <- list(c("WT", "DMKN \U03B2\U03B3-/-"),
                           c("WT","DMKN \U03B1\U03B2-/-"),
                           c("DMKN \U03B2\U03B3-/-", "DMKN \U03B2\U03B3-/-")
    )
    
    long_data$Genotype_Treatment <- factor(long_data$Genotype,
                                           levels = c("WT", "DMKN \U03B2\U03B3-/-",
                                                      "DMKN \U03B1\U03B2-/-"))
  } else if (group == "CS") {
    my_comparisons <- list(c("WT", "DMKN \U03B2\U03B3-/-"),
                           c("WT","DMKN \U03B1\U03B2-/-"),
                           c("DMKN \U03B2\U03B3-/-", "DMKN \U03B2\U03B3-/-")
                           
    )
    
    long_data$Genotype_Treatment <- factor(long_data$Genotype_Treatment,
                                           levels = c("WT", "DMKN \U03B2\U03B3-/-",
                                                      "DMKN \U03B1\U03B2-/-"))
  } else {
    stop("Invalid group specified. Please use either 'DTU' or 'CS'.")
  }
  
  plot <- ggplot(long_data, aes(x = Genotype, y = Score, fill = Genotype)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.7) + 
    #geom_errorbar(aes(ymin = Score - sd(Score), ymax = Score + sd(Score)), position = position_dodge(0.8), width = 0.3) +
    facet_wrap(~Var1, scales = "free_x", ncol = 5) +
    labs(x = "Comparison", y = "Kinase Activity Score") +
    theme_bw()+ ylim(-3,6)+ theme(
      legend.position = "bottom", 
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7), 
      axis.title = element_text(size=9, face = "bold"),
      strip.text = element_text(size = 8),
      plot.margin = margin(, 0.5, , , "cm"),
      strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")))+
    stat_compare_means(comparisons = my_comparisons, tip.length = 0.01, vjust = 0.55,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 2), 
                                          symbols = c("****", "***", "**", "*", "ns")), 
                       label = "p.signif", method = "t.test", label.y = max(long_data$Score) * 1.0,
                       step.increase = 0.25, hide.ns	= TRUE) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
  plot = plot + ggtitle(title)
  return(plot)
}
plot_kinase_activity_NoSerum_Genotypes_forFIT2 <- function(input_data, group,title) {
  if (group == "DTU") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("_control|_Control", colnames(input_data))]
    label1 <- "DMKN \U03B2\U03B3-/-"
  } else if (group == "CS") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("_control|_Control", colnames(input_data))]
    label2 <- "DMKN \U03B1\U03B2-/-"
  } else {
    stop("Invalid group. Please specify either 'DTU' or 'CS'")
  }
  
  long_data <- melt(filtered_data, variable.name = "Condition", value.name = "Score")
  
  long_data$Genotype <- ifelse(str_detect(long_data$Var2, "DTU"), "DMKN ????-/-",
                               ifelse(str_detect(long_data$Var2, "CS"), "DMKN \U03B1\U03B2-/-", "WT"))
  long_data$Treatment <- ifelse(str_detect(long_data$Var2, "Treated"), "Treated", "Control")
  long_data$Genotype_Treatment <- paste(long_data$Genotype, long_data$Treatment, sep = " ")
  
  # Define the list of comparisons based on the input group
  if (group == "DTU") {
    my_comparisons <- list(c("WT", "DMKN \U03B2\U03B3-/-"),
                           c("WT","DMKN \U03B1\U03B2-/-"),
                           c("DMKN \U03B2\U03B3-/-", "DMKN \U03B2\U03B3-/-")
    )
    
    long_data$Genotype_Treatment <- factor(long_data$Genotype,
                                           levels = c("WT", "DMKN \U03B2\U03B3-/-",
                                                      "DMKN \U03B1\U03B2-/-"))
  } else if (group == "CS") {
    my_comparisons <- list(c("WT", "DMKN \U03B2\U03B3-/-"),
                           c("WT","DMKN \U03B1\U03B2-/-"),
                           c("DMKN \U03B2\U03B3-/-", "DMKN \U03B2\U03B3-/-")
                           
    )
    
    long_data$Genotype_Treatment <- factor(long_data$Genotype_Treatment,
                                           levels = c("WT", "DMKN \U03B2\U03B3-/-",
                                                      "DMKN \U03B1\U03B2-/-"))
  } else {
    stop("Invalid group specified. Please use either 'DTU' or 'CS'.")
  }
  
  plot <- ggplot(long_data, aes(x = Genotype, y = Score, fill = Genotype)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.7) + 
    #geom_errorbar(aes(ymin = Score - sd(Score), ymax = Score + sd(Score)), position = position_dodge(0.8), width = 0.3) +
    facet_wrap(~Var1, scales = "free_x", ncol = 5) +
    labs(x = "Comparison", y = "Kinase Activity Score") +
    theme_bw()+ ylim(-3,6)+ theme(
      legend.position = "bottom", 
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7), 
      axis.title = element_text(size=9, face = "bold"),
      strip.text = element_text(size = 8),
      plot.margin = margin(, 0.5, , , "cm"),
      strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")))+
    stat_compare_means(comparisons = my_comparisons, tip.length = 0.01, vjust = 0.55,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 2), 
                                          symbols = c("****", "***", "**", "*", "ns")), 
                       label = "p.signif", method = "t.test", label.y = max(long_data$Score) * 1.0,
                       step.increase = 0.25, hide.ns	= TRUE) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
  plot = plot + ggtitle(title)
  return(plot)
}
plot_kinase_activity_NoSerum_Genotypes_forFIT3 <- function(input_data, group,title) {
  if (group == "DTU") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("WT_NoSerum_Control|DTU_NoSerum_Control", colnames(input_data))]
    label <- "DMKN \U03B2\U03B3-/-"
  } else if (group == "CS") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("WT_NoSerum_Control|CS_NoSerum_Control", colnames(input_data))]
    label <- "DMKN \U03B1\U03B2-/-"
  } else {
    stop("Invalid group. Please specify either 'DTU' or 'CS'")
  }
  
  long_data <- melt(filtered_data, variable.name = "Condition", value.name = "Score")
  
  long_data$Genotype <- ifelse(str_detect(long_data$Var2, "DTU"), "DMKN ????-/-",
                               ifelse(str_detect(long_data$Var2, "CS"), "DMKN ????-/-", "WT"))
  
  long_data$Genotype <- factor(long_data$Genotype,
                               levels = c("WT", "DMKN ????-/-", "DMKN ????-/-"))
  
  long_data$Treatment <- ifelse(str_detect(long_data$Var2, "Treated"), "Treated", "Control")
  long_data$Genotype_Treatment <- paste(long_data$Genotype, long_data$Treatment, sep = " ")
  
  # Define the list of comparisons based on the input group
  if (group == "DTU") {
    my_comparisons <- list(c("DMKN ????-/-", "WT"))
  } else if (group == "CS") {
    my_comparisons <- list(c("DMKN ????-/-","WT"))
  }
  else {
    stop("Invalid group specified. Please use either 'DTU' or 'CS'.")
  }
  
  plot <- ggplot(long_data, aes(x = Genotype, y = Score, fill = Genotype)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.7) + 
    #geom_errorbar(aes(ymin = Score - sd(Score), ymax = Score + sd(Score)), position = position_dodge(0.8), width = 0.3) +
    facet_wrap(~Var1, scales = "free_x", ncol = 5) +
    labs(x = "Comparison", y = "Kinase Activity Score") +
    theme_bw()+ ylim(-3,5)+ theme(
      legend.position = "bottom", 
      axis.text.x=element_blank(),
      axis.text.y = element_text(size = 7), 
      axis.title = element_text(size=9, face = "bold"),
      strip.text = element_text(size = 8),
      plot.margin = margin(, 0.5, , , "cm"),
      legend.key.size = unit(1, 'cm'),
      legend.title = element_text(size=14), 
      legend.text = element_text(size=14),
      strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")))+
    stat_compare_means(comparisons = my_comparisons, tip.length = 0.01, vjust = 0.55,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 2), 
                                          symbols = c("****", "***", "**", "*", "ns")), 
                       label = "p.signif", method = "t.test", label.y = max(long_data$Score) * 1.0,
                       step.increase = 0.25, hide.ns	= TRUE) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
  
  if (group == "DTU") {
    plot <- plot + scale_fill_manual(values = c("DMKN ????-/-" = "#1b9e77", "WT" = "#d95f02"))
  } else if (group == "CS") {
    plot <- plot + scale_fill_manual(values = c("DMKN ????-/-" = "#7570b3", "WT" = "#d95f02"))
  }
  plot = plot + ggtitle(title)
  return(plot)
}

plot_kinase_activity_Serum_Genotypes_forFIT4 <- function(input_data, group,title) {
  if (group == "DTU") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("WT_Serum_Control|DTU_Serum_Control", colnames(input_data))]
    label <- "DMKN \U03B2\U03B3-/-"
  } else if (group == "CS") {
    #Change to WT_NoSerum_Control if I dont want to show treated
    filtered_data <- input_data[, grep("WT_Serum_Serum_Control|CS_Serum_Control", colnames(input_data))]
    label <- "DMKN \U03B1\U03B2-/-"
  } else {
    stop("Invalid group. Please specify either 'DTU' or 'CS'")
  }
  
  long_data <- melt(filtered_data, variable.name = "Condition", value.name = "Score")
  
  long_data$Genotype <- ifelse(str_detect(long_data$Var2, "DTU"), "DMKN ????-/-",
                               ifelse(str_detect(long_data$Var2, "CS"), "DMKN ????-/-", "WT"))
  
  long_data$Genotype <- factor(long_data$Genotype,
                               levels = c("WT", "DMKN ????-/-", "DMKN ????-/-"))
  
  long_data$Treatment <- ifelse(str_detect(long_data$Var2, "Treated"), "Treated", "Control")
  long_data$Genotype_Treatment <- paste(long_data$Genotype, long_data$Treatment, sep = " ")
  
  # Define the list of comparisons based on the input group
  if (group == "DTU") {
    my_comparisons <- list(c("DMKN ????-/-", "WT"))
  } else if (group == "CS") {
    my_comparisons <- list(c("DMKN ????-/-","WT"))
  }
  else {
    stop("Invalid group specified. Please use either 'DTU' or 'CS'.")
  }
  
  plot <- ggplot(long_data, aes(x = Genotype, y = Score, fill = Genotype)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.7) + 
    #geom_errorbar(aes(ymin = Score - sd(Score), ymax = Score + sd(Score)), position = position_dodge(0.8), width = 0.3) +
    facet_wrap(~Var1, scales = "free_x") +
    labs(x = "Comparison", y = "Kinase Activity Score") +
    theme_bw()+ ylim(-3,5)+ theme(
      legend.position = "bottom", 
      axis.text.x=element_blank(),
      axis.text.y = element_text(size = 7), 
      axis.title = element_text(size=9, face = "bold"),
      strip.text = element_text(size = 8),
      plot.margin = margin(, 0.5, , , "cm"),
      legend.key.size = unit(1, 'cm'),
      legend.title = element_text(size=14), 
      legend.text = element_text(size=14),
      strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")))+
    stat_compare_means(comparisons = my_comparisons, tip.length = 0.01, vjust = 0.55,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 2), 
                                          symbols = c("****", "***", "**", "*", "ns")), 
                       label = "p.signif", method = "t.test", label.y = max(long_data$Score) * 1.0,
                       step.increase = 0.25, hide.ns	= TRUE) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
  
  if (group == "DTU") {
    plot <- plot + scale_fill_manual(values = c("DMKN ????-/-" = "#1b9e77", "WT" = "#d95f02"))
  } else if (group == "CS") {
    plot <- plot + scale_fill_manual(values = c("DMKN ????-/-" = "#7570b3", "WT" = "#d95f02"))
  }
  plot = plot + ggtitle(title)
  return(plot)
}


remove_terms_after_third_semicolon <- function(rowname) {
  split_name <- strsplit(rowname, ";")[[1]]
  new_name <- paste(split_name[1:4], collapse = ";")
  return(new_name)
}
extract_flanking_seq <- function(named_FASTA, rownames_ppe) {
  updated_rownames <- c()
  
  for (rowname in rownames_ppe) {
    split_name <- strsplit(rowname, ";")[[1]]
    
    # Remove the first and third terms if there are more than 4 terms
    if (length(split_name) > 4) {
      split_name <- split_name[-c(1, 3)]
    }
    # Remove any terms after the third semicolon
    split_name <- split_name[1:4]
    
    rowname <- remove_terms_after_third_semicolon(rowname)
    gene_id <- split_name[2]
    position <- as.integer(gsub("\\D", "", split_name[3]))
    matches <- grep(paste0(".*;", gene_id, "$"), names(named_FASTA), value = TRUE)
    
    if (length(matches) > 0) {
      protein_seq <- named_FASTA[matches[1]]
      
      # Correct the position indexing
      start <- max(1, position - 7)
      stop <- min(nchar(protein_seq), position + 7)
      
      flanking_seq <- substr(protein_seq, start, stop)
      new_name <- paste(split_name[1], split_name[2], split_name[3], flanking_seq, sep = ";")
      updated_rownames <- c(updated_rownames, new_name)
    } else {
      updated_rownames <- c(updated_rownames, rowname)
    }
  }
  
  return(updated_rownames)
}



fasta_to_named_list <- function(fasta_lines) {
  named_list <- list()
  current_header <- ""
  
  for (line in fasta_lines) {
    if (startsWith(line, ">")) {
      uniprot_id <- sub("^.*\\|(.+?)\\|.*$", "\\1", line)
      gene_id <- sub("^.*GN=(\\S+).*?$", "\\1", line)
      current_header <- paste(uniprot_id, gene_id, sep = ";")
    } else {
      named_list[[current_header]] <- line
    }
  }
  
  return(named_list)
}


reformat_fasta_df <- function(fasta_df) {
  reformatted_lines <- c()
  current_sequence <- ""
  
  for (i in 1:nrow(fasta_df)) {
    line <- as.character(fasta_df[i, 1])
    if (!is.na(line)) {
      if (startsWith(line, ">")) {
        if (nchar(current_sequence) > 0) {
          reformatted_lines <- c(reformatted_lines, current_sequence)
          current_sequence <- ""
        }
        reformatted_lines <- c(reformatted_lines, line)
      } else {
        current_sequence <- paste0(current_sequence, line)
      }
    }
  }
  
  if (nchar(current_sequence) > 0) {
    reformatted_lines <- c(reformatted_lines, current_sequence)
  }
  
  return(reformatted_lines)
}
rename_columns <- function(df) { new_colnames <- c("WT_NoSerum_Control", "WT_NoSerum_Treated", "WT_Serum_Control", "WT_Serum_Treated", "KO_DTU_NoSerum_Control", "KO_DTU_NoSerum_Treated", "KO_DTU_Serum_Control", "KO_DTU_Serum_Treated", "KO_CS_NoSerum_Control", "KO_CS_NoSerum_Treated", "KO_CS_Serum_Control", "KO_CS_Serum_Treated") 
for (i in 1:length(new_colnames)) 
{ 
  start_col <- (i - 1) * 3 + 1 
  end_col <- start_col + 2 
  colnames(df)[start_col:end_col] <- paste(new_colnames[i], rep(1:3, 1), sep = "_Rep") } 
return(df) }


rename_columns_2 <- function(df) {
  new_colnames <- c("WT_NoSerum_Control", "WT_NoSerum_Treated",
                    "KO_DTU_NoSerum_Control", "KO_DTU_NoSerum_Treated",
                    "KO_CS_NoSerum_Control", "KO_CS_NoSerum_Treated")
  for (i in 1:length(new_colnames)) {
    start_col <- (i - 1) * 3 + 1
    end_col <- start_col + 2
    colnames(df)[start_col:end_col] <- paste(new_colnames[i], rep(1:3, 1), sep = "_Rep")
  }
  return(df)
}

calculate_log2fc_combined <- function(df) {
  log2fc_df <- data.frame(matrix(ncol = 14, nrow = nrow(df)))
  colnames(log2fc_df) <- c("WT_NoSerum", "WT_Serum", "KO_DTU_NoSerum", "KO_DTU_Serum", "KO_CS_NoSerum", "KO_CS_Serum",
                           "KO_DTU_NoSerum_treated_WT_NoSerum_control", "KO_DTU_Serum_treated_WT_Serum_control",
                           "KO_CS_NoSerum_treated_WT_NoSerum_control", "KO_CS_Serum_treated_WT_Serum_control",
                           "KO_DTU_NoSerum_control_WT_NoSerum_control", "KO_DTU_Serum_control_WT_Serum_control",
                           "KO_CS_NoSerum_control_WT_NoSerum_control", "KO_CS_Serum_control_WT_Serum_control")
  
  
  # Define control and treated columns for each group
  WT_NoSerum_control <- c(1, 2, 3)
  WT_NoSerum_treated <- c(4, 5, 6)
  
  WT_Serum_control <- c(7, 8, 9)
  WT_Serum_treated <- c(10, 11, 12)
  
  KO_DTU_NoSerum_control <- c(13, 14, 15)
  KO_DTU_NoSerum_treated <- c(16, 17, 18)
  
  KO_DTU_Serum_control <- c(19, 20, 21)
  KO_DTU_Serum_treated <- c(22, 23, 24)
  
  KO_CS_NoSerum_control <- c(25, 26, 27)
  KO_CS_NoSerum_treated <- c(28, 29, 30)
  
  KO_CS_Serum_control <- c(31, 32, 33)
  KO_CS_Serum_treated <- c(34, 35, 36)
  
  group_list <- list(WT_NoSerum = list(control = WT_NoSerum_control, treated = WT_NoSerum_treated),
                     WT_Serum = list(control = WT_Serum_control, treated = WT_Serum_treated),
                     KO_DTU_NoSerum = list(control = KO_DTU_NoSerum_control, treated = KO_DTU_NoSerum_treated),
                     KO_DTU_Serum = list(control = KO_DTU_Serum_control, treated = KO_DTU_Serum_treated),
                     KO_CS_NoSerum = list(control = KO_CS_NoSerum_control, treated = KO_CS_NoSerum_treated),
                     KO_CS_Serum = list(control = KO_CS_Serum_control, treated = KO_CS_Serum_treated),
                     KO_DTU_NoSerum_treated_WT_NoSerum_control = list(control = WT_NoSerum_control, treated = KO_DTU_NoSerum_treated),
                     KO_DTU_Serum_treated_WT_Serum_control = list(control = WT_Serum_control, treated = KO_DTU_Serum_treated),
                     KO_CS_NoSerum_treated_WT_NoSerum_control = list(control = WT_NoSerum_control, treated = KO_CS_NoSerum_treated),
                     KO_CS_Serum_treated_WT_Serum_control = list(control = WT_Serum_control, treated = KO_CS_Serum_treated),
                     KO_DTU_NoSerum_control_WT_NoSerum_control = list(control = WT_NoSerum_control, treated = KO_DTU_NoSerum_control),
                     KO_DTU_Serum_control_WT_Serum_control = list(control = WT_Serum_control, treated = KO_DTU_Serum_control),
                     KO_CS_NoSerum_control_WT_NoSerum_control = list(control = WT_NoSerum_control, treated = KO_CS_NoSerum_control),
                     KO_CS_Serum_control_WT_Serum_control = list(control = WT_Serum_control, treated = KO_CS_Serum_control))
  
  for (i in seq(1, length(group_list))) {
    treated_cols <- group_list[[names(group_list)[i]]][["treated"]]
    control_cols <- group_list[[names(group_list)[i]]][["control"]]
    # Filter out rows with NA values in either control or treated groups
    no_na_rows <- !rowSums(is.na(df[, c(control_cols, treated_cols)]))
    filtered_df <- df[no_na_rows, ]
    
    log2fc_df[no_na_rows, i] <- log2(rowMeans(filtered_df[, treated_cols]) / rowMeans(filtered_df[, control_cols]))
  }
  return(log2fc_df)
}

calculate_log2fc_combined_2 <- function(df) {
  log2fc_df <- data.frame(matrix(ncol = 6, nrow = nrow(df)))
  colnames(log2fc_df) <- c("WT_NoSerum", "KO_DTU_NoSerum", "KO_CS_NoSerum",
                           "KO_DTU_NoSerum_treated_WT_NoSerum_control",
                           "KO_CS_NoSerum_treated_WT_NoSerum_control",
                           "KO_DTU_NoSerum_control_WT_NoSerum_control",
                           "KO_CS_NoSerum_control_WT_NoSerum_control")
  
  # Define control and treated columns for each group
  WT_NoSerum_control <- c(1, 2, 3)
  WT_NoSerum_treated <- c(4, 5, 6)
  
  KO_DTU_NoSerum_control <- c(7, 8, 9)
  KO_DTU_NoSerum_treated <- c(10, 11, 12)
  
  KO_CS_NoSerum_control <- c(13, 14, 15)
  KO_CS_NoSerum_treated <- c(16, 17, 18)
  
  group_list <- list(WT_NoSerum = list(control = WT_NoSerum_control, treated = WT_NoSerum_treated),
                     KO_DTU_NoSerum = list(control = KO_DTU_NoSerum_control, treated = KO_DTU_NoSerum_treated),
                     KO_CS_NoSerum = list(control = KO_CS_NoSerum_control, treated = KO_CS_NoSerum_treated),
                     KO_DTU_NoSerum_treated_WT_NoSerum_control = list(control = WT_NoSerum_control, treated = KO_DTU_NoSerum_treated),
                     KO_CS_NoSerum_treated_WT_NoSerum_control = list(control = WT_NoSerum_control, treated = KO_CS_NoSerum_treated),
                     KO_DTU_NoSerum_control_WT_NoSerum_control = list(control = WT_NoSerum_control, treated = KO_DTU_NoSerum_control),
                     KO_CS_NoSerum_control_WT_NoSerum_control = list(control = WT_NoSerum_control, treated = KO_CS_NoSerum_control))
  
  for (i in seq(1, length(group_list))) {
    treated_cols <- group_list[[names(group_list)[i]]][["treated"]]
    control_cols <- group_list[[names(group_list)[i]]][["control"]]
    # Filter out rows with NA values in either control or treated groups
    no_na_rows <- !rowSums(is.na(df[, c(control_cols, treated_cols)]))
    filtered_df <- df[no_na_rows, ]
    
    log2fc_df[no_na_rows, i] <- log2(rowMeans(filtered_df[, treated_cols]) / rowMeans(filtered_df[, control_cols]))
  }
  return(log2fc_df)
}
create_histogram_subplots <- function(df) {
  plot_list <- list()
  for (i in 1:ncol(df)) {
    p <- ggplot(df, aes_string(x = colnames(df)[i])) +
      geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
      theme_minimal() + annotate("text", x = max(3), y =max(10), 
                                 label = print(nrow(na.omit(df)))) +
      labs(x = colnames(df)[i], y = "Frequency", title = colnames(df)[i]) +
      theme(plot.title = element_text(hjust = 0.5))
    plot_list[[i]] <- p
  }
  grid.arrange(grobs = plot_list, ncol = 2)
}

display_min_max_values <- function(df) {
  for (i in 1:ncol(df)) {
    col_name <- colnames(df)[i]
    min_value <- min(df[[col_name]], na.rm = TRUE)
    max_value <- max(df[[col_name]], na.rm = TRUE)
    cat(sprintf("Column: %s\nMin value: %f\nMax value: %f\n\n", col_name, min_value, max_value))
  }
}

#Updated create_volcano
create_volcano_plot <- function(data, calculate_log2fc) {
  # Calculate mean intensities for each group
  data <- data %>%
    mutate(WT_NoSerum_Mean = rowMeans(cbind(WT_NoSerum_Control_Rep1, 
                                            WT_NoSerum_Control_Rep2, 
                                            WT_NoSerum_Control_Rep3,
                                            WT_NoSerum_Treated_Rep1, 
                                            WT_NoSerum_Treated_Rep2, 
                                            WT_NoSerum_Treated_Rep3)),
           WT_Serum_Mean = rowMeans(cbind(WT_Serum_Control_Rep1, 
                                          WT_Serum_Control_Rep2, 
                                          WT_Serum_Control_Rep3,
                                          WT_Serum_Treated_Rep1, 
                                          WT_Serum_Treated_Rep2, 
                                          WT_Serum_Treated_Rep3)),
           KO_DTU_NoSerum_Mean = rowMeans(cbind(KO_DTU_NoSerum_Control_Rep1, 
                                                KO_DTU_NoSerum_Control_Rep2,
                                                KO_DTU_NoSerum_Control_Rep3,
                                                KO_DTU_NoSerum_Treated_Rep1, 
                                                KO_DTU_NoSerum_Treated_Rep2, 
                                                KO_DTU_NoSerum_Treated_Rep3)),
           KO_DTU_Serum_Mean = rowMeans(cbind(KO_DTU_Serum_Control_Rep1, 
                                              KO_DTU_Serum_Control_Rep2, 
                                              KO_DTU_Serum_Control_Rep3,
                                              KO_DTU_Serum_Treated_Rep1, 
                                              KO_DTU_Serum_Treated_Rep2, 
                                              KO_DTU_Serum_Treated_Rep3)),
           KO_CS_NoSerum_Mean = rowMeans(cbind(KO_CS_NoSerum_Control_Rep1, 
                                               KO_CS_NoSerum_Control_Rep2,
                                               KO_CS_NoSerum_Control_Rep3,
                                               KO_CS_NoSerum_Treated_Rep1, 
                                               KO_CS_NoSerum_Treated_Rep2, 
                                               KO_CS_NoSerum_Treated_Rep3)),
           KO_CS_Serum_Mean = rowMeans(cbind(KO_CS_Serum_Control_Rep1, 
                                             KO_CS_Serum_Control_Rep2, 
                                             KO_CS_Serum_Control_Rep3,
                                             KO_CS_Serum_Treated_Rep1, 
                                             KO_CS_Serum_Treated_Rep2, 
                                             KO_CS_Serum_Treated_Rep3)),
           KO_DTU_NoSerum_treated_WT_NoSerum_control_Mean = rowMeans(cbind(KO_DTU_NoSerum_Treated_Rep1, 
                                                                           KO_DTU_NoSerum_Treated_Rep2, 
                                                                           KO_DTU_NoSerum_Treated_Rep3,
                                                                           WT_NoSerum_Control_Rep1, 
                                                                           WT_NoSerum_Control_Rep2, 
                                                                           WT_NoSerum_Control_Rep3)),
           KO_DTU_Serum_treated_WT_Serum_control_Mean = rowMeans(cbind(KO_DTU_Serum_Treated_Rep1, 
                                                                       KO_DTU_Serum_Treated_Rep2, 
                                                                       KO_DTU_Serum_Treated_Rep3,
                                                                       WT_Serum_Control_Rep1, 
                                                                       WT_Serum_Control_Rep2, 
                                                                       WT_Serum_Control_Rep3)),
           KO_CS_NoSerum_treated_WT_NoSerum_control_Mean = rowMeans(cbind(KO_CS_NoSerum_Treated_Rep1, 
                                                                          KO_CS_NoSerum_Treated_Rep2, 
                                                                          KO_CS_NoSerum_Treated_Rep3,
                                                                          WT_NoSerum_Control_Rep1, 
                                                                          WT_NoSerum_Control_Rep2, 
                                                                          WT_NoSerum_Control_Rep3)),
           KO_CS_Serum_treated_WT_Serum_control_Mean = rowMeans(cbind(KO_CS_Serum_Treated_Rep1, 
                                                                      KO_CS_Serum_Treated_Rep2, 
                                                                      KO_CS_Serum_Treated_Rep3,
                                                                      WT_Serum_Control_Rep1, 
                                                                      WT_Serum_Control_Rep2, 
                                                                      WT_Serum_Control_Rep3)),
           KO_DTU_NoSerum_control_WT_NoSerum_control_Mean = rowMeans(cbind(KO_DTU_NoSerum_Control_Rep1, 
                                                                           KO_DTU_NoSerum_Control_Rep2, 
                                                                           KO_DTU_NoSerum_Control_Rep3,
                                                                           WT_NoSerum_Control_Rep1, 
                                                                           WT_NoSerum_Control_Rep2, 
                                                                           WT_NoSerum_Control_Rep3)),
           KO_DTU_Serum_control_WT_Serum_control_Mean = rowMeans(cbind(KO_DTU_Serum_Control_Rep1, 
                                                                       KO_DTU_Serum_Control_Rep2, 
                                                                       KO_DTU_Serum_Control_Rep3,
                                                                       WT_Serum_Control_Rep1, 
                                                                       WT_Serum_Control_Rep2, 
                                                                       WT_Serum_Control_Rep3)),
           KO_CS_NoSerum_control_WT_NoSerum_control_Mean = rowMeans(cbind(KO_CS_NoSerum_Control_Rep1, 
                                                                          KO_CS_NoSerum_Control_Rep2, 
                                                                          KO_CS_NoSerum_Control_Rep3,
                                                                          WT_NoSerum_Control_Rep1, 
                                                                          WT_NoSerum_Control_Rep2, 
                                                                          WT_NoSerum_Control_Rep3)),
           KO_CS_Serum_control_WT_Serum_control_Mean= rowMeans(cbind(KO_CS_Serum_Control_Rep1, 
                                                                     KO_CS_Serum_Control_Rep2, 
                                                                     KO_CS_Serum_Control_Rep3,
                                                                     WT_Serum_Control_Rep1, 
                                                                     WT_Serum_Control_Rep2, 
                                                                     WT_Serum_Control_Rep3)))
  
  
  # Calculate log2 fold changes within each group
  data <- data %>%
    mutate(WT_NoSerum_Mean = log10(WT_NoSerum_Mean),
           WT_Serum_Mean = log10(WT_Serum_Mean),
           KO_DTU_NoSerum_Mean = log10(KO_DTU_NoSerum_Mean),
           KO_DTU_Serum_Mean = log10(KO_DTU_Serum_Mean),
           KO_CS_NoSerum_Mean = log10(KO_CS_NoSerum_Mean),
           KO_CS_Serum_Mean = log10(KO_CS_Serum_Mean),
           KO_DTU_NoSerum_treated_WT_NoSerum_control_Mean = log10(KO_DTU_NoSerum_treated_WT_NoSerum_control_Mean),
           KO_DTU_Serum_treated_WT_Serum_control_Mean = log10(KO_DTU_Serum_treated_WT_Serum_control_Mean),
           KO_CS_NoSerum_treated_WT_NoSerum_control_Mean = log10(KO_CS_NoSerum_treated_WT_NoSerum_control_Mean),
           KO_CS_Serum_treated_WT_Serum_control_Mean = log10(KO_CS_Serum_treated_WT_Serum_control_Mean),
           KO_DTU_NoSerum_control_WT_NoSerum_control_Mean = log10(KO_DTU_NoSerum_control_WT_NoSerum_control_Mean),
           KO_DTU_Serum_control_WT_Serum_control_Mean = log10(KO_DTU_Serum_control_WT_Serum_control_Mean),
           KO_CS_NoSerum_control_WT_NoSerum_control_Mean = log10(KO_CS_NoSerum_control_WT_NoSerum_control_Mean),
           KO_CS_Serum_control_WT_Serum_control_Mean = log10(KO_CS_Serum_control_WT_Serum_control_Mean),
           
           WT_NoSerum_FC = calculate_log2fc$WT_NoSerum,
           WT_Serum_FC = calculate_log2fc$WT_Serum,
           KO_DTU_NoSerum_FC = calculate_log2fc$KO_DTU_NoSerum,
           KO_DTU_Serum_FC = calculate_log2fc$KO_DTU_Serum,
           KO_CS_NoSerum_FC = calculate_log2fc$KO_CS_NoSerum,
           KO_CS_Serum_FC = calculate_log2fc$KO_CS_Serum,
           KO_DTU_NoSerum_treated_WT_NoSerum_control_FC = calculate_log2fc$KO_DTU_NoSerum_treated_WT_NoSerum_control,
           KO_DTU_Serum_treated_WT_Serum_control_FC = calculate_log2fc$KO_DTU_Serum_treated_WT_Serum_control,
           KO_CS_NoSerum_treated_WT_NoSerum_control_FC = calculate_log2fc$KO_CS_NoSerum_treated_WT_NoSerum_control,
           KO_CS_Serum_treated_WT_Serum_control_FC = calculate_log2fc$KO_CS_Serum_treated_WT_Serum_control,
           KO_DTU_NoSerum_control_WT_NoSerum_control_FC = calculate_log2fc$KO_DTU_NoSerum_control_WT_NoSerum_control,
           KO_DTU_Serum_control_WT_Serum_control_FC = calculate_log2fc$KO_DTU_Serum_control_WT_Serum_control,
           KO_CS_NoSerum_control_WT_NoSerum_control_FC = calculate_log2fc$KO_CS_NoSerum_control_WT_NoSerum_control,
           KO_CS_Serum_control_WT_Serum_control_FC = calculate_log2fc$KO_CS_Serum_control_WT_Serum_control)
  
  
  
  create_group_plot <- function(data, group_mean_col, group_fc_col, title) {
    
    group_fc_values <- data[[group_fc_col]]
    not_significant_count <- sum(group_fc_values >= -0.5849 & group_fc_values <= 0.5849, na.rm = TRUE)
    downregulated_count  <- sum(group_fc_values < -0.5849, na.rm = TRUE)
    upregulated_count <- sum(group_fc_values > 0.5849, na.rm = TRUE)
    
    ggplot(data, aes_string(x = group_fc_col, y = group_mean_col)) +
      geom_point(aes_string(color = paste("ifelse(", group_fc_col, " > 0.5849, 'red', ifelse(", group_fc_col, " < -0.5849, 'blue', 'black'))")), alpha = 0.3, shape =16, size = 1) +
      scale_color_identity() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA), 
            axis.title=element_text(size=8)) +
      labs(x = title, y = "Log10 XIC", size =1.0) +
      geom_point(aes(x = -6.75, y = log10(10^(10.0)), color = 'black'), size = 1.0, shape = 16) +
      geom_point(aes(x = -6.75, y = log10(10^(9.5)), color = 'blue'), size = 1.0, shape = 16) +
      geom_point(aes(x = -6.75, y = log10(10^(9.00)), color = 'red'), size = 1.0, shape = 16) +
      annotate("text", x = -6.5, y = log10(10^(10.0)),
               label = paste("Not significant: n=", not_significant_count),
               hjust = 0, size = 1.75) +
      annotate("text", x = -6.5, y = log10(10^(9.5)),
               label = paste("Sig. down-regulated: n=", downregulated_count),
               hjust = 0, size = 1.75) +
      annotate("text", x = -6.5, y = log10(10^(9.00)),
               label = paste("Sig. up-regulated: n=", upregulated_count),
               hjust = 0, size = 1.75) +
      annotate("text", x = 4.0, y = log10(10^(9.8)),
               label = paste("N = 3"),
               hjust = 0, size = 2.5) +
      coord_cartesian(xlim = c(-6.5, 5), ylim = c(log10(10^(4.00)), log10(10^(10.00))))
  }
  
  # Create subplots for each group
  plot1 <- create_group_plot(data, "WT_NoSerum_Mean", "WT_NoSerum_FC", 
                             expression(paste("Log2 ratio (WT"," (Treated)/WT ","(Control), Starved)")))
  plot2 <- create_group_plot(data, "KO_DTU_NoSerum_Mean", "KO_DTU_NoSerum_FC", expression(paste("Log2 ratio (DMKN ",beta,gamma, atop(scriptscriptstyle("-/-")), " (Treated) vs DMKN ",beta,gamma, atop(scriptscriptstyle("-/-")), " (Control), Starved)")))
  plot3 <- create_group_plot(data, "KO_CS_NoSerum_Mean", "KO_CS_NoSerum_FC", expression(paste("Log2 ratio (DMKN ",alpha,beta, atop(scriptscriptstyle("-/-")), " (Treated) vs DMKN ",alpha,beta, atop(scriptscriptstyle("-/-")), " (Control), Starved)")))
  plot4 <- create_group_plot(data, "WT_Serum_Mean", "WT_Serum_FC", expression(paste("Log2 ratio (DMKN ",beta,gamma, atop(scriptscriptstyle("-/-")), " (Treated) vs DMKN ",beta,gamma, atop(scriptscriptstyle("-/-")), " (Control), Serum)")))
  plot5 <- create_group_plot(data, "KO_DTU_Serum_Mean", "KO_DTU_Serum_FC", expression(paste("Log2 ratio (DMKN ",beta,gamma, atop(scriptscriptstyle("-/-")), " (Treated) vs DMKN ",beta,gamma, atop(scriptscriptstyle("-/-")), " (Control), Serum)")))
  plot6 <- create_group_plot(data, "KO_CS_Serum_Mean", "KO_CS_Serum_FC", expression(paste("Log2 ratio (DMKN ",alpha,beta, atop(scriptscriptstyle("-/-")), " (Treated) vs DMKN ",alpha,beta, atop(scriptscriptstyle("-/-")), " (Control), Serum)")))
  # Create subplots for each new group
  plot7 <- create_group_plot(data, "KO_DTU_NoSerum_treated_WT_NoSerum_control_Mean", "KO_DTU_NoSerum_treated_WT_NoSerum_control_FC",
                             expression(paste("Log2 ratio (DMKN ",beta,gamma, atop(scriptscriptstyle("-/-")), " (Treated) vs WT (Control), Starved)")))
  plot8 <- create_group_plot(data, "KO_DTU_Serum_treated_WT_Serum_control_Mean", "KO_DTU_Serum_treated_WT_Serum_control_FC", 
                             expression(paste("Log2 ratio (DMKN ",beta,gamma, atop(scriptscriptstyle("-/-")), " (Treated) vs WT (Control), Serum)")))
  plot9 <- create_group_plot(data, "KO_CS_NoSerum_treated_WT_NoSerum_control_Mean", "KO_CS_NoSerum_treated_WT_NoSerum_control_FC",
                             expression(paste("Log2 ratio (DMKN ",alpha,beta, atop(scriptscriptstyle("-/-")), " (Treated) vs WT (Control), Starved)")))
  plot10 <- create_group_plot(data, "KO_CS_Serum_treated_WT_Serum_control_Mean", "KO_CS_Serum_treated_WT_Serum_control_FC", 
                              expression(paste("Log2 ratio (DMKN ",alpha,beta, atop(scriptscriptstyle("-/-")), " (Treated) vs WT (Control), Serum)")))
  plot11 <- create_group_plot(data, "KO_DTU_NoSerum_control_WT_NoSerum_control_Mean", "KO_DTU_NoSerum_control_WT_NoSerum_control_FC",
                              expression(paste("Log2 ratio (DMKN ",beta,gamma, atop(scriptscriptstyle("-/-")), " (Control) vs WT (Control), Starved)")))
  plot12 <- create_group_plot(data, "KO_DTU_Serum_control_WT_Serum_control_Mean", "KO_DTU_Serum_control_WT_Serum_control_FC", 
                              expression(paste("Log2 ratio (DMKN ",beta,gamma, atop(scriptscriptstyle("-/-")), " (Control) vs WT (Control), Serum)")))
  plot13 <- create_group_plot(data, "KO_CS_NoSerum_control_WT_NoSerum_control_Mean", "KO_CS_NoSerum_control_WT_NoSerum_control_FC",
                              expression(paste("Log2 ratio (DMKN ",alpha,beta, atop(scriptscriptstyle("-/-")), " (Control) vs WT (Control), Starved)")))
  plot14 <- create_group_plot(data, "KO_CS_Serum_control_WT_Serum_control_Mean", "KO_CS_Serum_control_WT_Serum_control_FC", 
                              expression(paste("Log2 ratio (DMKN ",alpha,beta, atop(scriptscriptstyle("-/-")), " (Control) vs WT (Control), Serum)")))
  
  # Combine subplots into a grid
  library(gridExtra)
  #grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, plot13, plot14, ncol = 3)
  #  options(repr.plot.width = 14, repr.plot.height =16)
  NoSerum_plot = grid.arrange(plot11, plot13, plot7, plot9, plot2, plot3, ncol = 2)
  Serum_plot = grid.arrange(plot12, plot14, plot8, plot10, plot5, plot6, ncol = 2)
  NoSerum_plot ; Serum_plot ; ggsave(filename=paste("NoSerum_plot", ".pdf"), plot = NoSerum_plot, device = "pdf") ; ggsave(filename=paste("Serum_plot", ".pdf"), plot = Serum_plot, device = "pdf")
}

fuzzPlot <- function (Tc, clustObj, mfrow = c(1, 1), cols, min.mem = 0, new.window = FALSE, llwd = 3, x_labels) {
  clusterindex <- clustObj$cluster
  memship <- clustObj$membership
  memship[memship < min.mem] <- -1
  
  colorindex <- integer(dim(Tc)[[1]])
  if (missing(cols)) {
    cols <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", "#FF0030", "#FF0018")
  }
  colorseq <- seq(0, 1, length = length(cols))
  for (j in 1:max(clusterindex)) {
    tmp <- Tc[clusterindex == j, ]
    tmpmem <- memship[clusterindex == j, j]
    if (((j - 1) %% (mfrow[1] * mfrow[2])) == 0) {
      if (new.window) 
        dev.new()
      par(mfrow = mfrow)
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      plot(x = NA, xlim = c(1, dim(Tc)[[2]]), ylim = c(ymin, ymax), xaxt = "n", xlab = "", ylab = "Standardized Profile", main = paste("Cluster", j, "; size=", nrow(tmp)))
      axis(1, at = 1:ncol(Tc), labels = x_labels)
      
      mean_tmp <- apply(tmp, 2, mean)
      
      # Fit a linear model
      fit <- lm(mean_tmp ~ log(1:ncol(Tc)))
      
      
      # Add the linear trend line to the plot
      lines(predict(fit), col = 'black')
    }
    else {
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      
      plot(x = NA, xlim = c(1, dim(Tc)[[2]]), ylim = c(ymin, ymax), xaxt = "n", xlab = "", ylab = "Standardized Profile", main = paste("Cluster", j, "; size=", nrow(tmp)))
      axis(1, at = 1:ncol(Tc), labels = x_labels, las = 2, cex.axis = 1, tck = -0.01)
      # Add the angled labels with custom size
      
      # Calculate the overall trend line
      y_all <- c(tmp)
      x_all <- c(rep(1:ncol(Tc), nrow(tmp)))
      loess_fit <- loess(y_all ~ x_all)
      # Add the overall trend line to the plot
      lines(1:ncol(Tc), predict(loess_fit, data.frame(x_all = 1:ncol(Tc))), col = "black", lwd = 2, lty = 2)
      
    }
    if (!(sum(clusterindex == j) == 0)) {
      for (jj in 1:(length(colorseq) - 1)) {
        tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= colorseq[jj + 1])
        if (sum(tmpcol) > 0) {
          tmpind <- which(tmpcol)
          for (k in 1:length(tmpind)) {
            lines(tmp[tmpind[k], ], col = cols[jj], lwd = llwd)
          }
        }
      }
    }
  }
}

custom_clustOptimal <- function  (clueObj, rep = 5, user.maxK = NULL, visualize = TRUE, x_labels,...) {
  bstPvalue <- 1
  bst.clustObj <- c()
  bst.evaluation <- c()
  for (i in 1:rep) {
    clustObj <- c()
    if (clueObj$clustAlg == "cmeans") {
      if (is.null(user.maxK)) {
        clustObj <- cmeans(clueObj$Tc, centers = clueObj$maxK, 
                           iter.max = 50, m = 1.25)
      }
      else {
        clustObj <- cmeans(clueObj$Tc, centers = user.maxK, 
                           iter.max = 50, m = 1.25)
      }
    }
    else {
      if (is.null(user.maxK)) {
        clustObj <- kmeans(clueObj$Tc, centers = clueObj$maxK, 
                           iter.max = 50)
        clustObj$membership <- matrix(1, nrow = nrow(clueObj$Tc), 
                                      ncol = nrow(clustObj$centers))
        rownames(clustObj$membership) <- names(clustObj$cluster)
      }
      else {
        clustObj <- kmeans(clueObj$Tc, centers = user.maxK, 
                           iter.max = 50)
        clustObj$membership <- matrix(1, nrow = nrow(clueObj$Tc), 
                                      ncol = nrow(clustObj$centers))
        rownames(clustObj$membership) <- names(clustObj$cluster)
      }
    }
    evaluation <- clustEnrichment(clustObj, clueObj$annotation, 
                                  clueObj$effectiveSize, clueObj$pvalueCutoff)
    currentPvalue <- evaluation$fisher.pvalue
    if (currentPvalue < bstPvalue) {
      bstPvalue <- currentPvalue
      bst.clustObj <- clustObj
      bst.enrichList <- evaluation$enrich.list
    }
  }
  if (visualize) {
    fuzzPlot(clueObj$Tc, clustObj = bst.clustObj, x_labels = x_labels,...)
  }
  results <- list()
  results$clustObj <- bst.clustObj
  results$enrichList <- bst.enrichList
  return(results)
}

custom_runClue = function (Tc, annotation, rep = 5, kRange = 2:10, clustAlg = "cmeans", 
                           effectiveSize = c(5, 100), pvalueCutoff = 0.05, alpha = 0.5) 
{
  means <- apply(Tc, 1, mean)
  stds <- apply(Tc, 1, sd)
  tmp <- sweep(Tc, 1, means, FUN = "-")
  Tc <- sweep(tmp, 1, stds, FUN = "/")
  annotation.intersect <- lapply(annotation, intersect, rownames(Tc))
  annotation.filtered <- annotation.intersect[lapply(annotation.intersect, 
                                                     length) > 0]
  # ... (same code as before, up to the repeat.list)
  
  repeat.list <- mclapply(1:rep, function(rp) {
    cat("repeat", rp, "\n")
    enrichment <- c()
    for (k in kRange) {
      clustered <- c()
      
      # Perform na.omit on Tc before clustering
      Tc_no_na <- na.omit(Tc)
      
      if (clustAlg == "cmeans") {
        clustered <- cmeans(Tc_no_na, centers = k, iter.max = 50, m = 1.25)
      } else if (clustAlg == "kmeans") {
        clustered <- kmeans(Tc_no_na, centers = k, iter.max = 50)
      } else {
        print("Unknown clustering algorithm specified. Using cmeans clustering instead")
        clustered <- cmeans(Tc_no_na, centers = k, iter.max = 50, m = 1.25)
      }
      evaluate <- clustEnrichment(clustered, annotation.filtered, effectiveSize, pvalueCutoff)
      fisher.pvalue <- evaluate$fisher.pvalue
      escore <- -log10(fisher.pvalue) - alpha * nrow(clustered$centers)
      enrichment <- c(enrichment, escore)
    }
    enrichment
  })
  
  # ... (same code as before, after the repeat.list)
  x <- do.call(rbind, repeat.list)
  x.normalize <- (x - min(x))/(max(x) - min(x))
  rownames(x.normalize) <- paste("repeat", 1:rep, sep = "")
  colnames(x.normalize) <- paste("k", kRange, sep = "=")
  maxK <- which.max(apply(x.normalize, 2, median)) + (kRange[1] - 
                                                        1)
  result <- list()
  result$Tc <- Tc
  result$annotation <- annotation.filtered
  result$clustAlg <- clustAlg
  result$effectiveSize <- effectiveSize
  result$pvalueCutoff <- pvalueCutoff
  result$evlMat <- x.normalize
  result$maxK <- maxK
  return(result)
}


custom_kinaseSubstratePred = function (phosScoringMatrices, ensembleSize = 10, top = 50, cs = 0.8, 
                                       inclusion = 20, iter = 5, verbose = TRUE) 
{
  substrate.list = substrateList(phosScoringMatrices, top, 
                                 cs, inclusion)
  if (verbose) 
    message("Predicting kinases for phosphosites:")
  featureMat <- phosScoringMatrices$combinedScoreMatrix
  predMatrix <- matrix(0, nrow = nrow(featureMat), ncol = length(substrate.list))
  colnames(predMatrix) <- names(substrate.list)
  rownames(predMatrix) <- rownames(featureMat)
  tmp.list = lapply(seq(length(substrate.list)), function(i) {
    positive.train <- featureMat[substrate.list[[i]], ]
    positive.cls <- rep(1, length(substrate.list[[i]]))
    negative.pool <- featureMat[!(rownames(featureMat) %in% 
                                    substrate.list[[i]]), ]
    if (verbose) 
      message(paste(i, ".", sep = ""))
    tmp_col = predMatrix[, i]
    for (e in seq_len(ensembleSize)) {
      negativeSize <- length(substrate.list[[i]])
      idx <- sample(seq_len(nrow(negative.pool)), size = negativeSize, 
                    replace = TRUE)
      negative.samples <- rownames(negative.pool)[idx]
      negative.train <- featureMat[negative.samples, ]
      negative.cls <- rep(2, length(negative.samples))
      train.mat <- rbind(positive.train, negative.train)
      cls <- as.factor(c(positive.cls, negative.cls))
      names(cls) <- rownames(train.mat)
      pred <- multiAdaSampling(train.mat, test.mat = featureMat, 
                               label = cls, kernelType = "radial", iter = iter)
      tmp_col <- tmp_col[names(pred[, 1])] + pred[, 1]
    }
    tmp_col
  })
  predMatrix = matrix(unlist(tmp.list), ncol = ncol(predMatrix))
  colnames(predMatrix) = names(substrate.list)
  rownames(predMatrix) = rownames(featureMat)
  predMatrix <- predMatrix/ensembleSize
  if (verbose) 
    message("done")
  return(predMatrix)
}

multiAdaSampling <- function(train.mat, test.mat,
                             label, kernelType, iter = 5) {
  
  X <- train.mat
  Y <- label
  
  model <- c()
  prob.mat <- c()
  
  for (i in seq_len(iter)) {
    tmp <- X
    rownames(tmp) <- NULL
    model <- e1071::svm(tmp, factor(Y),
                        kernel = kernelType, probability = TRUE)
    prob.mat <- attr(predict(model, train.mat,
                             decision.values = FALSE, probability = TRUE),
                     "probabilities")
    
    X <- c()
    Y <- c()
    for (j in seq_len(ncol(prob.mat))) {
      voteClass <- prob.mat[label == colnames(prob.mat)[j], ]
      idx <- c()
      idx <- sample(seq_len(nrow(voteClass)),
                    size = nrow(voteClass), replace = TRUE,
                    prob = voteClass[, j])
      X <- rbind(X, train.mat[rownames(voteClass)[idx],])
      Y <- c(Y, label[rownames(voteClass)[idx]])
    }
  }
  
  pred <- attr(predict(model, newdata = test.mat,
                       probability = TRUE), "prob")
  return(pred)
}

kinaseSubstrateHeatmap_modified <- function (phosScoringMatrices, top = 3, printPlot = NULL, filePath = "./kinaseSubstrateHeatmap.pdf", 
                                             width = 10, height = 10) 
{
  utils::data("KinaseFamily", envir = environment())
  sites <- c()
  for (i in seq_len(ncol(phosScoringMatrices$combinedScoreMatrix))) {
    sites <- union(sites, names(sort(phosScoringMatrices$combinedScoreMatrix[, 
                                                                             i], decreasing = TRUE)[seq_len(top)]))
  }
  o <- intersect(colnames(phosScoringMatrices$combinedScoreMatrix), 
                 rownames(KinaseFamily))
  annotation_col = data.frame(group = KinaseFamily[o, "kinase_group"], 
                              family = KinaseFamily[o, "kinase_family"])
  rownames(annotation_col) <- o
  if (is.null(printPlot) == TRUE) {
    
    pheatmap::pheatmap(phosScoringMatrices$combinedScoreMatrix[sites, 
    ], annotation_col = annotation_col, cluster_rows = TRUE, 
    cluster_cols = TRUE, fontsize = 6, show_row_dendrogram = FALSE, show_col_dendrogram = FALSE)
    
  }
  else {
    pdf(file = filePath, width = width, height = height)
    
    pheatmap::pheatmap(phosScoringMatrices$combinedScoreMatrix[sites, 
    ], annotation_col = annotation_col, cluster_rows = TRUE, 
    cluster_cols = TRUE, fontsize = 6, show_row_dendrogram = FALSE, show_col_dendrogram = FALSE)
    
    dev.off()
  }
}
kinaseSubstrateHeatmap_modified <- function (phosScoringMatrices, top = 3, printPlot = NULL, filePath = "./kinaseSubstrateHeatmap.pdf", 
                                             width = 10, height = 10) 
{
  utils::data("KinaseFamily", envir = environment())
  sites <- c()
  for (i in seq_len(ncol(phosScoringMatrices$combinedScoreMatrix))) {
    sites <- union(sites, names(sort(phosScoringMatrices$combinedScoreMatrix[, 
                                                                             i], decreasing = TRUE)[seq_len(top)]))
  }
  if (is.null(printPlot) == TRUE) {
    pheatmap::pheatmap(phosScoringMatrices$combinedScoreMatrix[sites, 
    ], cluster_rows = TRUE, 
    cluster_cols = TRUE, fontsize = 6, show_row_dendrogram = FALSE, show_col_dendrogram = FALSE, show_rownames = FALSE, legend = FALSE)
  }
  else {
    pdf(file = filePath, width = width, height = height)
    pheatmap::pheatmap(phosScoringMatrices$combinedScoreMatrix[sites, 
    ], cluster_rows = TRUE, 
    cluster_cols = TRUE, fontsize = 6, show_row_dendrogram = FALSE, show_col_dendrogram = FALSE, show_rownames = FALSE, legend = FALSE)
    dev.off()
  }
}
siteAnnotate_modified <- function(site, phosScoringMatrices, predMatrix) {
  od <- order(predMatrix[site, ], decreasing = FALSE)
  kinases <- colnames(predMatrix)[od]
  
  par(mfrow = c(1, 4), mar = c(1, 5, 5, 1))
  
  barplot(predMatrix[site, kinases], las = 1, xlab = "Prediction score", 
          col = "red3", main = "Prediction score", xlim = c(0, 1), horiz = TRUE, xaxt= "n")
  axis(3)
  barplot(phosScoringMatrices$combinedScoreMatrix[site, kinases], 
          las = 1, main = "Combined score", col = "orange2", xlim = c(0, 1), horiz = TRUE,xaxt= "n")
  axis(3)
  barplot(phosScoringMatrices$motifScoreMatrix[site, kinases], 
          las = 1, main = "Motif score", col = "green4", xlim = c(0, 1), horiz = TRUE,xaxt= "n")
  axis(3)
  barplot(phosScoringMatrices$profileScoreMatrix[site, kinases], 
          las = 1, main = "Profile score", col = "lightblue3", 
          xlim = c(0, 1), horiz = TRUE,xaxt= "n")
  axis(3)
  mtext(paste("Site =", site),                   # Add main title
        side = 3,
        line = - 2,
        outer = TRUE)
}

siteAnnotate_modified_ggplot <- function(site, phosScoringMatrices, predMatrix) {
  od <- order(predMatrix[site, ], decreasing = FALSE)
  kinases <- colnames(predMatrix)[od]
  
  plot_data <- list(predMatrix[site, kinases], phosScoringMatrices$combinedScoreMatrix[site, kinases],
                    phosScoringMatrices$motifScoreMatrix[site, kinases], phosScoringMatrices$profileScoreMatrix[site, kinases])
  plot_types <- c("Prediction score", "Combined score", "Motif score", "Profile score")
  plot_colors <- c("red3", "orange2", "green4", "lightblue3")
  
  plots <- list()
  
  for (i in 1:4) {
    df <- data.frame(Kinase = factor(kinases, levels = kinases), Value = plot_data[[i]])
    
    p <- ggplot(df, aes(x = Kinase, y = Value, fill = factor(plot_types[i]))) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = plot_colors[i], guide = "none") +
      theme_classic() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(size = 4),
            axis.title.x = element_text(size = 0),
            plot.title = element_text(hjust = 0.5, size = 7),
            axis.text.x.top = element_text(angle = 45, hjust = 1)) +
      labs(title = plot_types[i], x = NULL) +
      scale_x_discrete(position = "top") +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0)), limits = c(0, 1),breaks = c(0,1))
    
    plots[[i]] <- p
  }
  title1=text_grob(paste(site), size = 9, face = "bold")   #### this worked for me
  
  grid.arrange(grobs = plots, ncol = 4, top = title1)
}


plotSignalomeMap_modified <- function(signalomes, color) {
  df <- stack(signalomes$kinaseSubstrates)
  modules <- signalomes$proteinModule
  names(modules) <- unlist(lapply(strsplit(as.character(names(signalomes$proteinModules)), ";"), "[[", 1))
  df$cluster <- modules[df$values]
  df_balloon <- df
  df_balloon <- na.omit(df_balloon) %>% dplyr::count(.data$cluster, .data$ind)
  df_balloon$ind <- as.factor(df_balloon$ind)
  df_balloon$cluster <- as.factor(df_balloon$cluster)
  df_balloon <- tidyr::spread(df_balloon, .data$ind, .data$n)[, -1]
  df_balloon[is.na(df_balloon)] <- 0
  df_balloon <- do.call(rbind, lapply(seq(nrow(df_balloon)), function(x) {
    res <- unlist(lapply(df_balloon[x, ], function(y) y/sum(df_balloon[x, ]) * 100))
  }))
  df_balloon <- reshape2::melt(as.matrix(df_balloon))
  colnames(df_balloon) <- c("cluster", "ind", "n")
  df_balloon <- df_balloon[df_balloon$n >= 0, ]
  
  # Function to find the last balloon for each group
  last_balloon_y <- function(group_data) {
    max(group_data$y[group_data$x == max(group_data$x)])
  }
  
  # Add a new column to the dataframe containing the last balloon's y value for each group
  df_balloon$last_y <- sapply(unique(df_balloon$cluster), function(cluster) {
    last_balloon_y(df_balloon[df_balloon$cluster == cluster, ])
  })[df_balloon$cluster]
  
  
  # Function to find the first balloon for each group
  first_balloon_y <- function(group_data) {
    min(group_data$y[group_data$x == min(group_data$x)])
  }
  
  # Add a new column to the dataframe containing the first balloon's y value for each group
  df_balloon$first_y <- sapply(unique(df_balloon$cluster), function(cluster) {
    first_balloon_y(df_balloon[df_balloon$cluster == cluster, ])
  })[df_balloon$cluster]  
  
  g <- ggplot2::ggplot(df_balloon, aes(y = .data$ind, x = .data$cluster)) + 
    
    
    #geom_text(aes(label = .data$cluster, x = .data$cluster, y = 0), angle = 90, hjust = 1, size = 3) +  # Show x-axis labels using geom_text
    scale_color_manual(values = color) + 
    scale_size_continuous(name = "% of phosphosites regulated by a kinase", range = c(2, 17)) + 
    theme_classic() + 
    theme(legend.position = "right", legend.direction = "vertical", legend.title = element_text(angle = 90),
          axis.line = element_blank(), 
          axis.title = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank()) +
    guides(col = FALSE)  # Remove the legend for the kinases
  # Add rectangles and text labels for module names on the x-axis
  
  
  g <- g + geom_line(aes(group = .data$cluster, col = .data$ind), size = 0.5, linetype = "solid", color = "grey") +  # Connect the balloons within the horizontal line
    geom_segment(data = df_balloon, aes(x = .data$cluster, xend = .data$cluster, y = .data$first_y, yend = .data$last_y, col = .data$ind), size = 0.5, linetype = "solid", show.legend = FALSE, color = "grey") 
  
  # Add rectangles for module names on the x-axis
  for (i in unique(df_balloon$cluster)) {
    g <- g + annotate("rect", xmin = as.numeric(i) - 0.45, xmax = as.numeric(i) + 0.45, ymin = -0.05 * max(df_balloon$n), ymax = 0, alpha = 0.7, fill = color[i])
  }
  
  g = g +geom_point(aes(col = .data$ind, size = .data$n)) 
  
  
  # Use geom_segment() to draw lines from the x-axis to the last balloon for each group
  #g <- g + 
  g
}

