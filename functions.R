# load("required_data.RData")

#-------------------Analysis for unprocessed data-------------------
# Parameters
## method (drop down menu): t.test, wilcox.test, mod.t.test
## ctrl: blank or ?
## exp: blank or ?
## unmapped_FA: blank or ?
## exo_lipid: blank or ?
## species (drop down menu): human, mouse, rat
## progress: progress object
## update_progress: progress bar update function

FA_substructure_analysis <- function(exp_raw, method, ctrl, exp,
                                     unmapped_FA = c("w9-18:2;0", "w3-20:4;0"),
                                     exo_lipid = "w3-22:6;0", species = "rat",
                                     progress = NULL, update_progress = NULL) {
  exp_data <- build_char_table(exp_raw, network_node = network_node)[[1]]

  char_data <- build_char_table(exp_raw, network_node = network_node)[[2]]

  no_sub_t <- unprocessed_data_test(
    exp_data = exp_data,
    char_table = char_data,
    method = method,
    significant = "adj_p_value",
    ctrl_group = ctrl, exp_group = exp
  )

  #-------------------FA substructure analysis-------------------

  # FA biosynthetic network data transformation

  FA_network_new <- build_FA_net(
    FA_network = FA_network,
    unprocessed_data_result = no_sub_t
  )

  # Decompose lipids into FA substructures
  # 18:2 and 20:4 are majorly omega-6 FAs, so we only kept omega-6 forms of them

  FA_substructure <- FA_sub_transform(
    FA_network = FA_network_new,
    unprocessed_data_result = no_sub_t,
    unmapped_FA = unmapped_FA
  )

  # Extract FA substructures using fold changes

  FA_sub_stop <- FA_sub_extract(
    char_table = char_data,
    FA_substructure = FA_substructure,
    unprocessed_data_result = no_sub_t,
    exact_FA = "no", exo_lipid = exo_lipid
  )

  # Transform FA exp into substructure exp

  FA_sub_exp <- lipid_sub_matrix(
    exp_data = exp_data, sub_data = FA_sub_stop,
    sub_type = "FA"
  )


  # Differential expression analysis for FA substructures

  FA_sub_exp_t <- t_test(
    data = FA_sub_exp[[3]], ctrl = ctrl, exp = exp,
    method = method, significant = "adj_p_value"
  )

  print("Substructure transformation complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Substructure transformation complete.")
  }

  # Essential pathway analysis for FA substructures

  set.seed(1)

  path_score_FA <- path_scoring(
    network = FA_network_new, sub_t = FA_sub_exp_t,
    calibrate = T, data_type = "FA"
  )
  path_score_FA_sel <-
    path_score_FA %>%
    filter(Significant == "yes") %>%
    mutate(Type = ifelse(Type == "Active", "Increase", "Decrease"))

  if (nrow(path_score_FA_sel) != 0) {
    path_data <- rbind(
      path_score_FA_sel %>%
        filter(Type == "Decrease") %>%
        arrange(cal_score) %>% .[!duplicated(.$rep_sub_path), ] %>% .[1:5, ],
      path_score_FA_sel %>%
        filter(Type == "Increase") %>%
        .[!duplicated(.$rep_sub_path), ] %>% .[1:5, ]
    ) %>%
      filter(Significant == "yes") %>%
      dplyr::select(-score)

    add_suffix <- function(strings) {
      counts <- table(strings)
      duplicated_indices <- which(duplicated(strings))
      count <- 2
      for (index in duplicated_indices) {
        strings[index] <- paste(strings[index], paste0("(", count, ")"), sep = "")
        count <- count + 1
      }

      return(strings)
    }
    path_data$path <- str_split(path_data$path, " --> ") %>% map_chr(~ str_c(.x[[1]], " --> ", last(.x)))
    path_data$path <- add_suffix(path_data$path)

    path_data_fig <- path_data %>%
      mutate(path = factor(.$path, levels = .$path)) %>%
      ggplot(aes(x = reorder(path, cal_score), y = cal_score, fill = Type)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0) +
      geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "gray") +
      coord_flip() +
      theme_bw() +
      scale_fill_manual(values = c("Increase" = "red", "Decrease" = "blue")) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "top"
      ) +
      labs(
        x = "", y = "Pathway score",
        title = "Top 5 representative pathways"
      )
  } else {
    path_data_fig <- NA
    path_score_FA_sel <- NA
  }

  path_score_FA_sel <- path_score_FA[, c("path", "from", "to", "cal_score", "Significant", "Type", "rep_sub_path")]

  colnames(path_score_FA_sel) <- c("Pathway", "From", "To", "Score", "Significant", "Type", "Representative pathway")

  print("Pathway analysis complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Pathway analysis complete.")
  }

  # Essential edges (reactions) analysis for FA substructures

  reaction_score_FA <- reaction_scoring(
    network = FA_network_new,
    sub_exp = FA_sub_exp[[3]],
    sub_t = FA_sub_exp_t,
    ctrl = ctrl, exp = exp,
    Species = species
  )

  reaction_score_FA_sel <-
    reaction_score_FA[, c(
      "edge_name", "p_value", "mlog10p",
      "perturbation_score", "Mode", "genes"
    )] %>%
    filter(p_value < 0.05)

  if (nrow(reaction_score_FA_sel) != 0) {
    reaction_data <- rbind(
      reaction_score_FA %>% filter(perturbation_score > 0) %>% .[1:5, ],
      reaction_score_FA %>% filter(perturbation_score < 0) %>%
        arrange(perturbation_score) %>% .[1:5, ]
    ) %>%
      filter(p_value < 0.05)

    reaction_data <- reaction_data %>%
      mutate(
        node1 = str_split(.$edge_name, " --> ") %>% map_chr(.f = function(x) {
          x[1]
        }),
        node2 = str_split(.$edge_name, " --> ") %>% map_chr(.f = function(x) {
          x[2]
        })
      ) %>%
      mutate(node1_color = ifelse(node1_log2FC > 0, paste0("<i style='color:#FF0000'>", node1, " --> ", "</i>"),
        paste0("<i style='color:#0000FF'>", node1, " --> ", "</i>")
      )) %>%
      mutate(node2_color = ifelse(node2_log2FC > 0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
        paste0("<i style='color:#0000FF'>", node2, "</i>")
      )) %>%
      mutate(edge_color = paste0(node1_color, node2_color))

    reaction_data_fig <- reaction_data %>%
      ggplot(aes(
        x = perturbation_score, y = reorder(edge_name, perturbation_score),
        fill = Mode
      )) +
      geom_bar(stat = "identity", size = 0.8) +
      scale_y_discrete(
        labels = rev(reaction_data$edge_color)
      ) +
      geom_vline(xintercept = 0) +
      theme_bw() +
      theme(
        legend.position = "top",
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_fill_manual(values = c("Increase" = "red", "Decrease" = "blue")) +
      labs(
        y = "", fill = "Reaction", x = "Perturbation score",
        title = "Top 5 significant reations"
      )
  } else {
    reaction_score_FA_sel <- NA
    reaction_data_fig <- NA
  }

  reaction_score_FA_sel <- reaction_score_FA %>% mutate(Significant = ifelse(p_value < 0.05, "yes", "no"))
  reaction_score_FA_sel <- reaction_score_FA_sel[, c("edge_name", "p_value", "perturbation_score", "Significant", "Mode", "genes")]

  colnames(reaction_score_FA_sel) <- c("Reaction", "P-value", "Perturbation score", "Significant", "Type", "Gene")

  print("Reaction analysis complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Reaction analysis complete.")
  }

  # FA biosynthetic network construction
  FA_network_data <- draw_network(
    network_data = FA_network_new,
    DE_data = FA_sub_exp_t,
    if_species = F, significant = "adj_p_value",
    path_scoring_result = path_score_FA,
    reaction_scoring_result = reaction_score_FA,
    top_n = 5, path_type = "both"
  )


  network <- visNetwork(FA_network_data[[1]], FA_network_data[[2]]) %>%
    visIgraphLayout(
      layout = "layout_with_sugiyama", type = "square",
      physics = F, smooth = TRUE, randomSeed = 5
    )

  print("Network analysis complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Network analysis complete.")
  }

  #---------------------new code--------------------------------

  sub_result <- FA_sub_exp_t[, c(
    "lipid", "mean_ctrl", "mean_exp", "FC",
    "log2FC", "p_value", "adj_p_value", "sig"
  )] %>%
    `colnames<-`(c(
      "Substructure", "Mean(ctrl)", "Mean(exp)", "FC", "Log2(FC)",
      "P-value", "Adjusted p-value", "Significance"
    )) %>%
    arrange(`Adjusted p-value`, desc(Substructure))

  network_node <- FA_network_data[[1]] %>% left_join(sub_result, by = c("id" = "Substructure"))
  network_edge <- FA_network_data[[2]] %>%
    mutate(Reaction = str_c(from, " --> ", to)) %>%
    left_join(
      reaction_score_FA[, c(
        "edge_name", "p_value",
        "perturbation_score", "Mode", "genes"
      )] %>%
        `colnames<-`(c("Reaction", "P-value", "Perturbation score", "Type", "Gene")),
      by = "Reaction"
    ) %>%
    arrange(desc(label))

  DE_lipid_data <- FA_sub_exp_t %>%
    mutate(log2FC = ifelse(is.infinite(log2FC), 10 * sign(log2FC), log2FC)) %>%
    mutate(Significance = ifelse(log2FC > 0, "Increase", "Decrease")) %>%
    mutate(Significance = ifelse(sig == "yes", Significance, "No change")) %>%
    ggplot(aes(x = log2FC, y = mlog10padj, col = Significance)) +
    geom_jitter(width = 0.3) +
    # scale_x_continuous(limits = c(-10,10), labels = c('-Inf','-5','0','5', 'Inf'))+
    scale_color_manual(values = c("Increase" = "red", "Decrease" = "blue", "No change" = "gray")) +
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
    theme_classic() +
    theme(legend.position = "top") +
    labs(y = "-Log10 (Adjusted p-value)", x = "Log2 (Fold change)")

  return(list(
    path_score_FA_sel, path_data_fig, reaction_score_FA_sel,
    reaction_data_fig, network, sub_result, network_node, network_edge, DE_lipid_data
  ))
}

#-------------------Analysis for unprocessed data-------------------
# Parameters
## method (drop down menu): t.test, wilcox.test, mod.t.test
## ctrl: blank or ?
## exp: blank or ?
## exo_lipid: blank or ?
## species (drop down menu): human, mouse, rat
## progress: progress object
## update_progress: progress bar update function

lipid_species_substructure_analysis <- function(exp_raw, method, ctrl, exp,
                                                non_missing_pct = 0.3,
                                                exo_lipid = NULL, species = "rat",
                                                progress = NULL, update_progress = NULL) {
  exp_data <- build_char_table(exp_raw, network_node = network_node)[[1]]

  char_data <- build_char_table(exp_raw, network_node = network_node)[[2]]

  no_sub_t <- unprocessed_data_test(
    exp_data = exp_data,
    char_table = char_data,
    method = method,
    significant = "adj_p_value",
    ctrl_group = ctrl, exp_group = exp
  )

  #-------------------Lipid species substructure analysis-------------------


  # Decompose lipids into species substructures

  species_substructure <- species_sub_transform(
    char = char_data,
    lipid_substructure = lipid_substructure,
    network_node = network_node
  )

  # Extract species substructures using fold changes

  species_sub_stop <- species_sub_extract(
    lipid_substructure = species_substructure,
    unprocessed_data_result = no_sub_t,
    type = "species", pct_limit = non_missing_pct,
    exo_lipid = exo_lipid
  )

  # Transform lipid exp into substructure exp

  species_sub_exp <- lipid_sub_matrix(
    exp_data = exp_data,
    sub_data = species_sub_stop,
    sub_type = "species"
  )


  # Differential expression analysis for substructures

  species_sub_exp_t <- t_test(
    data = species_sub_exp[[3]], ctrl = ctrl, exp = exp,
    method = method, significant = "adj_p_value"
  )

  print("Substructure transformation complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Substructure transformation complete.")
  }

  # Species biosynthetic network data transformation

  species_network <- build_species_net(species_substructure = species_substructure)


  # Essential pathway analysis for species substructures

  set.seed(1)
  path_score <- path_scoring(
    network = species_network,
    sub_t = species_sub_exp_t,
    calibrate = T, data_type = "Species"
  )


  path_score_sel <- path_score %>%
    filter(Significant == "yes") %>%
    mutate(Type = ifelse(Type == "Active", "Increase", "Decrease"))

  if (nrow(path_score_sel) == 0) {
    path_data_fig <- NA
    path_score_sel <- NA
  } else {
    path_data <- rbind(
      path_score_sel %>%
        filter(Type == "Decrease") %>%
        arrange(cal_score) %>% .[!duplicated(.$rep_sub_path), ] %>% .[1:5, ],
      path_score_sel %>%
        filter(Type == "Increase") %>%
        .[!duplicated(.$rep_sub_path), ] %>% .[1:5, ]
    ) %>%
      filter(Significant == "yes") %>%
      dplyr::select(-score)

    add_suffix <- function(strings) {
      counts <- table(strings)
      duplicated_indices <- which(duplicated(strings))
      count <- 2
      for (index in duplicated_indices) {
        strings[index] <- paste(strings[index], paste0("(", count, ")"), sep = "")
        count <- count + 1
      }

      return(strings)
    }
    path_data$path <- str_split(path_data$path, " --> ") %>% map_chr(~ str_c(.x[[1]], " --> ", last(.x)))
    path_data$path <- add_suffix(path_data$path)

    path_data_fig <- path_data %>%
      mutate(path = factor(.$path, levels = .$path)) %>%
      ggplot(aes(x = reorder(path, cal_score), y = cal_score, fill = Type)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0) +
      geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "gray") +
      coord_flip() +
      theme_bw() +
      scale_fill_manual(values = c("Increase" = "red", "Decrease" = "blue")) +
      theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5)
      ) +
      labs(
        x = "", y = "Pathway score",
        title = "Top 5 representative pathways"
      )
  }
  path_score_sel <- path_score[, c("path", "from", "to", "cal_score", "Significant", "Type", "rep_sub_path")]

  colnames(path_score_sel) <- c("Pathway", "From", "To", "Score", "Significant", "Type", "Representative pathway")




  print("Pathway analysis complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Pathway analysis complete.")
  }


  # Essential edges (reactions) analysis for species substructures

  species_net_w_rev <- add_rev_reaction(
    network_edge = network_edge,
    species_net = species_network
  )

  reaction_score <- reaction_scoring(
    network = species_net_w_rev,
    sub_exp = species_sub_exp[[3]],
    sub_t = species_sub_exp_t,
    ctrl = ctrl, exp = exp,
    Species = species
  )


  reaction_score_sel <-
    reaction_score[, c(
      "edge_name", "p_value", "mlog10p",
      "perturbation_score", "Mode", "genes"
    )] %>%
    filter(p_value < 0.05)


  if (nrow(reaction_score_sel) == 0) {
    reaction_data_fig <- NA
    reaction_score_sel <- NA
  } else {
    reaction_data <- rbind(
      reaction_score %>% filter(perturbation_score > 0) %>% .[1:5, ],
      reaction_score %>% filter(perturbation_score < 0) %>%
        arrange(perturbation_score) %>% .[1:5, ]
    ) %>%
      filter(p_value < 0.05)


    reaction_data <- reaction_data %>%
      mutate(
        node1 = str_split(.$edge_name, " --> ") %>% map_chr(.f = function(x) {
          x[1]
        }),
        node2 = str_split(.$edge_name, " --> ") %>% map_chr(.f = function(x) {
          x[2]
        })
      ) %>%
      mutate(node1_color = ifelse(node1_log2FC > 0, paste0("<i style='color:#FF0000'>", node1, " --> ", "</i>"),
        paste0("<i style='color:#0000FF'>", node1, " --> ", "</i>")
      )) %>%
      mutate(node2_color = ifelse(node2_log2FC > 0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
        paste0("<i style='color:#0000FF'>", node2, "</i>")
      )) %>%
      mutate(edge_color = paste0(node1_color, node2_color))



    reaction_data_fig <- reaction_data %>%
      ggplot(aes(
        x = perturbation_score, y = reorder(edge_name, perturbation_score),
        fill = Mode
      )) +
      geom_bar(stat = "identity", size = 0.8) +
      scale_y_discrete(
        labels = rev(reaction_data$edge_color)
      ) +
      geom_vline(xintercept = 0) +
      theme_bw() +
      theme(
        legend.position = "top",
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_fill_manual(values = c("Increase" = "red", "Decrease" = "blue")) +
      labs(
        y = "", fill = "Reaction", x = "Perturbation score",
        title = "Top 5 significant reactions"
      )
  }
  reaction_score_sel <- reaction_score %>% mutate(Significant = ifelse(p_value < 0.05, "yes", "no"))
  reaction_score_sel <- reaction_score_sel[, c("edge_name", "p_value", "perturbation_score", "Significant", "Mode", "genes")]

  colnames(reaction_score_sel) <- c("Reaction", "P-value", "Perturbation score", "Significant", "Type", "Gene")

  print("Reaction analysis complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Reaction analysis complete.")
  }

  # Lipid species biosynthetic network construction

  species_network_data <- draw_network(
    network_data = species_net_w_rev,
    DE_data = species_sub_exp_t,
    if_species = T, significant = "adj_p_value",
    path_scoring_result = path_score,
    reaction_scoring_result = reaction_score,
    top_n = 5, path_type = "both"
  )


  network <- visNetwork(species_network_data[[1]], species_network_data[[2]])
  print("Network analysis complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Network analysis complete.")
  }

  #---------------------new code--------------------------------

  sub_result <- species_sub_exp_t[, c(
    "lipid", "mean_ctrl", "mean_exp", "FC",
    "log2FC", "p_value", "adj_p_value", "sig"
  )] %>%
    `colnames<-`(c(
      "Substructure", "Mean(ctrl)", "Mean(exp)", "FC", "Log2(FC)",
      "P-value", "Adjusted p-value", "Significance"
    )) %>%
    arrange(`Adjusted p-value`, desc(Substructure))

  network_node <- species_network_data[[1]] %>% left_join(sub_result, by = c("id" = "Substructure"))
  network_edge <- species_network_data[[2]] %>%
    mutate(Reaction = str_c(from, " --> ", to)) %>%
    left_join(
      reaction_score[, c(
        "edge_name", "p_value",
        "perturbation_score", "Mode", "genes"
      )] %>%
        `colnames<-`(c("Reaction", "P-value", "Perturbation score", "Type", "Gene")),
      by = "Reaction"
    ) %>%
    arrange(desc(label))

  DE_lipid_data <- species_sub_exp_t %>%
    mutate(log2FC = ifelse(is.infinite(log2FC), 10 * sign(log2FC), log2FC)) %>%
    mutate(Significance = ifelse(log2FC > 0, "Increase", "Decrease")) %>%
    mutate(Significance = ifelse(sig == "yes", Significance, "No change")) %>%
    ggplot(aes(x = log2FC, y = mlog10padj, col = Significance)) +
    geom_jitter(width = 0.3) +
    # scale_x_continuous(limits = c(-10,10), labels = c('-Inf','-5','0','5', 'Inf'))+
    scale_color_manual(values = c("Increase" = "red", "Decrease" = "blue", "No change" = "gray")) +
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
    theme_classic() +
    theme(legend.position = "top") +
    labs(y = "-Log10 (Adjusted p-value)", x = "Log2 (Fold change)")

  return(list(
    path_score_sel, path_data_fig, reaction_score_sel,
    reaction_data_fig, network, sub_result, network_node, network_edge,
    DE_lipid_data
  ))
}

#-------------------Analysis for unprocessed data-------------------
# Parameters
## method (drop down menu): t.test, wilcox.test, mod.t.test
## ctrl: blank or ?
## exp: blank or ?
## exo_lipid: blank or ?
## species (drop down menu): human, mouse, rat
## progress: progress object
## update_progress: progress bar update function

lipid_class_substructure_analysis <- function(exp_raw, method, ctrl, exp,
                                              exo_lipid = NULL, species = "rat",
                                              progress = NULL, update_progress = NULL) {
  exp_data <- build_char_table(exp_raw, network_node = network_node)[[1]]

  char_data <- build_char_table(exp_raw, network_node = network_node)[[2]]

  no_sub_t <- unprocessed_data_test(
    exp_data = exp_data,
    char_table = char_data,
    method = method,
    significant = "adj_p_value",
    ctrl_group = ctrl, exp_group = exp
  )

  #-------------------Lipid class substructure analysis-------------------

  # Extract class substructures using fold changes

  class_sub_stop <- species_sub_extract(
    lipid_substructure = lipid_substructure,
    unprocessed_data_result = no_sub_t,
    type = "class", pct_limit = 0.01,
    exo_lipid = exo_lipid
  )

  # Transform lipid exp into substructures exp

  class_exp <- no_sub_t[[1]] %>%
    filter(type == "class") %>%
    dplyr::select(-type)


  class_sub_exp <- lipid_sub_matrix(
    exp_data = class_exp,
    sub_data = class_sub_stop,
    sub_type = "Class"
  )



  # Differential expression analysis for substructures

  class_sub_exp_t <- t_test(
    data = class_sub_exp[[3]], ctrl = ctrl, exp = exp,
    method = method, significant = "adj_p_value"
  )

  print("Substructure transformation complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Substructure transformation complete.")
  }

  # Class biosynthetic network data transformation

  class_network <- network_edge[c("S1", "P1")] %>%
    filter(S1 %in% class_sub_exp_t$lipid, P1 %in% class_sub_exp_t$lipid)

  # Essential pathway analysis for class substructures

  set.seed(1)
  path_score <- path_scoring(
    network = class_network,
    sub_t = class_sub_exp_t,
    calibrate = T, data_type = "Class"
  )


  path_score_sel <- path_score %>%
    filter(Significant == "yes") %>%
    mutate(Type = ifelse(Type == "Active", "Increase", "Decrease"))

  if (nrow(path_score_sel) == 0) {
    path_data_fig <- NA
    path_score_sel <- NA
  } else {
    path_data <- rbind(
      path_score_sel %>%
        filter(Type == "Decrease") %>%
        arrange(cal_score) %>% .[!duplicated(.$rep_sub_path), ] %>% .[1:5, ],
      path_score_sel %>%
        filter(Type == "Increase") %>%
        .[!duplicated(.$rep_sub_path), ] %>% .[1:5, ]
    ) %>%
      filter(Significant == "yes") %>%
      dplyr::select(-score)


    add_suffix <- function(strings) {
      counts <- table(strings)
      duplicated_indices <- which(duplicated(strings))
      count <- 2
      for (index in duplicated_indices) {
        strings[index] <- paste(strings[index], paste0("(", count, ")"), sep = "")
        count <- count + 1
      }

      return(strings)
    }
    path_data$path <- str_split(path_data$path, " --> ") %>% map_chr(~ str_c(.x[[1]], " --> ", last(.x)))
    path_data$path <- add_suffix(path_data$path)

    path_data_fig <- path_data %>%
      mutate(path = factor(.$path, levels = .$path)) %>%
      ggplot(aes(x = reorder(path, cal_score), y = cal_score, fill = Type)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0) +
      geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "gray") +
      coord_flip() +
      theme_bw() +
      scale_fill_manual(values = c("Increase" = "red", "Decrease" = "blue")) +
      theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5)
      ) +
      labs(
        x = "", y = "Pathway score",
        title = "Top 5 representative pathways"
      )
  }

  path_score_sel <- path_score[, c("path", "from", "to", "cal_score", "Significant", "Type", "rep_sub_path")]

  colnames(path_score_sel) <- c("Pathway", "From", "To", "Score", "Significant", "Type", "Representative pathway")


  print("Pathway analysis complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Pathway analysis complete.")
  }


  # Essential edges (reactions) analysis for class substructures

  reaction_score <- reaction_scoring(
    network = class_network,
    sub_exp = class_sub_exp[[3]],
    sub_t = class_sub_exp_t,
    ctrl = ctrl, exp = exp,
    Species = species
  )

  reaction_score_sel <-
    reaction_score[, c(
      "edge_name", "p_value", "mlog10p",
      "perturbation_score", "Mode", "genes"
    )] %>%
    filter(p_value < 0.05)

  if (nrow(reaction_score_sel) == 0) {
    reaction_data_fig <- NA
    reaction_score_sel <- NA
  } else {
    reaction_data <- rbind(
      reaction_score %>% filter(perturbation_score > 0) %>% .[1:5, ],
      reaction_score %>% filter(perturbation_score < 0) %>%
        arrange(perturbation_score) %>% .[1:5, ]
    ) %>%
      filter(p_value < 0.05)


    reaction_data <- reaction_data %>%
      mutate(
        node1 = str_split(.$edge_name, " --> ") %>% map_chr(.f = function(x) {
          x[1]
        }),
        node2 = str_split(.$edge_name, " --> ") %>% map_chr(.f = function(x) {
          x[2]
        })
      ) %>%
      mutate(node1_color = ifelse(node1_log2FC > 0, paste0("<i style='color:#FF0000'>", node1, " --> ", "</i>"),
        paste0("<i style='color:#0000FF'>", node1, " --> ", "</i>")
      )) %>%
      mutate(node2_color = ifelse(node2_log2FC > 0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
        paste0("<i style='color:#0000FF'>", node2, "</i>")
      )) %>%
      mutate(edge_color = paste0(node1_color, node2_color))


    reaction_data_fig <- reaction_data %>%
      ggplot(aes(
        x = perturbation_score, y = reorder(edge_name, perturbation_score),
        fill = Mode
      )) +
      geom_bar(stat = "identity", size = 0.8) +
      scale_y_discrete(
        labels = rev(reaction_data$edge_color)
      ) +
      geom_vline(xintercept = 0) +
      theme_bw() +
      theme(
        legend.position = "top",
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_fill_manual(values = c("Increase" = "red", "Decrease" = "blue")) +
      labs(
        y = "", fill = "Reaction", x = "Perturbation score",
        title = "Top 5 significant reactions"
      )
  }
  reaction_score_sel <- reaction_score %>% mutate(Significant = ifelse(p_value < 0.05, "yes", "no"))
  reaction_score_sel <- reaction_score_sel[, c("edge_name", "p_value", "perturbation_score", "Significant", "Mode", "genes")]

  colnames(reaction_score_sel) <- c("Reaction", "P-value", "Perturbation score", "Significant", "Type", "Gene")

  print("Reaction analysis complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Reaction analysis complete.")
  }

  # Lipid class biosynthetic network construction

  class_network_data <- draw_network(
    network_data = class_network,
    DE_data = class_sub_exp_t,
    if_species = F, significant = "adj_p_value",
    path_scoring_result = path_score,
    reaction_scoring_result = reaction_score,
    top_n = 5, path_type = "both"
  )

  network <- visNetwork(class_network_data[[1]], class_network_data[[2]])
  print("Network analysis complete.")
  if (is.function(update_progress)) {
    update_progress(progress = progress, detail = "Network analysis complete.")
  }

  #---------------------new code--------------------------------

  sub_result <- class_sub_exp_t[, c(
    "lipid", "mean_ctrl", "mean_exp", "FC",
    "log2FC", "p_value", "adj_p_value", "sig"
  )] %>%
    `colnames<-`(c(
      "Substructure", "Mean(ctrl)", "Mean(exp)", "FC", "Log2(FC)",
      "P-value", "Adjusted p-value", "Significance"
    )) %>%
    arrange(`Adjusted p-value`, desc(Substructure))

  network_node <- class_network_data[[1]] %>% left_join(sub_result, by = c("id" = "Substructure"))
  network_edge <- class_network_data[[2]] %>%
    mutate(Reaction = str_c(from, " --> ", to)) %>%
    left_join(
      reaction_score[, c(
        "edge_name", "p_value",
        "perturbation_score", "Mode", "genes"
      )] %>%
        `colnames<-`(c("Reaction", "P-value", "Perturbation score", "Type", "Gene")),
      by = "Reaction"
    ) %>%
    arrange(desc(label))

  DE_lipid_data <- class_sub_exp_t %>%
    mutate(log2FC = ifelse(is.infinite(log2FC), 10 * sign(log2FC), log2FC)) %>%
    mutate(Significance = ifelse(log2FC > 0, "Increase", "Decrease")) %>%
    mutate(Significance = ifelse(sig == "yes", Significance, "No change")) %>%
    ggplot(aes(x = log2FC, y = mlog10padj, col = Significance)) +
    geom_jitter(width = 0.3) +
    # scale_x_continuous(limits = c(-10,10), labels = c('-Inf','-5','0','5', 'Inf'))+
    scale_color_manual(values = c("Increase" = "red", "Decrease" = "blue", "No change" = "gray")) +
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
    theme_classic() +
    theme(legend.position = "top") +
    labs(y = "-Log10 (Adjusted p-value)", x = "Log2 (Fold change)")

  return(list(
    path_score_sel, path_data_fig, reaction_score_sel,
    reaction_data_fig, network, sub_result, network_node, network_edge,
    DE_lipid_data
  ))
}

#-------------------Check for correct data format-------------------
# Parameters
## exp_data: check format of this file

check_data_format <- function(exp_data) {
  library(tidyverse)
  warning_message <- character()
  # check variable class
  ckeck_feature <- map_chr(exp_data, ~ class(.x))[1] != "character"

  ckeck_value <- !map_chr(exp_data, ~ class(.x))[-1] %in% c("numeric", "integer")

  ckeck_feature_name <- colnames(exp_data)[1] != "feature"

  if (ckeck_feature || sum(ckeck_value) != 0 || ckeck_feature_name) {
    warning_message <- "! Please ensure that the first column's name is 'feature' and it is a string variable, while the remaining columns are numeric variables."
  }

  lipid_class <- str_replace(str_extract(exp_data$feature, "(.+?(_))|(.+)"), "_", "")



  colnames(exp_data)[1] <- "feature"
  exp_data <- exp_data[, c(T, !ckeck_value)]
  exp_data <- exp_data %>% filter(!feature %in% c(c("G3P", "DHAP")))

  if (sum(unique(lipid_class) %in% network_node$Abbreviation) < 5 || nrow(exp_data) < 30) {
    warning_message <- c(warning_message, "! The dataset must consist of a minimum of 5 lipid classes and 30 lipid species.")
  }

  FA_format <- str_split(exp_data$feature, "_")
  FA_format1 <- FA_format %>% map_lgl(~ length(.x) < 2)
  FA_format2 <- FA_format %>% map_lgl(~ sum(!str_detect(.x[-1], "[0-9]+:[0-9]+;[0-9]")) != 0)
  FA_format <- FA_format1 | FA_format2

  if (sum(FA_format) != 0) {
    warning_text <- " with wrong lipid format."

    warning_message <- c(warning_message, str_c("! ", str_c(exp_data$feature[FA_format],
      collapse = ", "
    ), warning_text))
  }

  exp_data <- exp_data[!FA_format, ]


  lipid_class_not_support <- unique(lipid_class)[!unique(lipid_class) %in%
    network_node$Abbreviation]

  if (length(lipid_class_not_support) != 0) {
    warning_text <- " are not supported by iLipidome"
    warning_message <- c(warning_message, str_c("! ", str_c(lipid_class_not_support,
      collapse = ", "
    ), warning_text))
  }
  exp_data <- exp_data[lipid_class %in% network_node$Abbreviation, ]
  lipid_class <- lipid_class[lipid_class %in% network_node$Abbreviation]


  FA_num_ref <- network_node$FA[match(lipid_class, network_node$Abbreviation)]

  FA_num <- map_int(str_split(exp_data$feature, "_"), ~ length(.x))

  FA_num_ckeck <- (FA_num_ref == FA_num - 1) | (FA_num - 1 == 1)


  if (sum(FA_num_ckeck) != nrow(exp_data)) {
    warning_text <- " with wrong FA number."

    warning_message <- c(warning_message, str_c("! ", str_c(exp_data$feature[!FA_num_ckeck],
      collapse = ", "
    ), warning_text))
  }
  if (length(warning_message) != 0) {
    warning_message <- str_c(
      str_c(warning_message, collapse = "\n"),
      "\n! Please refer to section 2.2.1 for instructions on how to format your data"
    )
  }

  return(warning_message)
}
