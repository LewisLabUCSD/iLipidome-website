load('required_data.RData')

#-------------------Analysis for unprocessed data-------------------
# Parameters
## method (drop down menu): t.test, wilcox.test, mod.t.test
## ctrl: blank or ?
## exp: blank or ?
## unmapped_FA: blank or ?
## exo_lipid: blank or ?
## species (drop down menu): human, mouse, rat

FA_substructure_analysis <- function(exp_raw, method, ctrl, exp,
                                     unmapped_FA = c('w9-18:2;0','w3-20:4;0'),
                                     exo_lipid='w3-22:6;0',species='rat'){
  
  exp_data <- build_char_table(exp_raw, network_node = network_node)[[1]]
  
  char_data <- build_char_table(exp_raw, network_node = network_node)[[2]]
  
  no_sub_t <- unprocessed_data_test(exp_data = exp_data,
                                    char_table = char_data,
                                    method = method,
                                    significant='adj_p_value',
                                    ctrl_group = ctrl, exp_group = exp)
  
  #-------------------FA substructure analysis-------------------
  
  #FA biosynthetic network data transformation
  
  FA_network_new <- build_FA_net(FA_network = FA_network,
                                 unprocessed_data_result = no_sub_t)

  #Decompose lipids into FA substructures
  #18:2 and 20:4 are majorly omega-6 FAs, so we only kept omega-6 forms of them
  
  FA_substructure <- FA_sub_transform(FA_network = FA_network_new,
                                      unprocessed_data_result = no_sub_t,
                                      unmapped_FA = unmapped_FA)
  
  #Extract FA substructures using fold changes
  
  FA_sub_stop <- FA_sub_extract(char_table = char_data,
                                FA_substructure = FA_substructure,
                                unprocessed_data_result = no_sub_t,
                                exact_FA='no', exo_lipid=exo_lipid)
  
  #Transform FA exp into substructure exp
  
  FA_sub_exp <- lipid_sub_matrix(exp_data = exp_data, sub_data = FA_sub_stop,
                                 sub_type = 'FA')
  
  
  #Differential expression analysis for FA substructures
  
  FA_sub_exp_t <- t_test(data = FA_sub_exp[[3]], ctrl = ctrl, exp = exp,
                         method = method, significant = 'adj_p_value')
  
  print('Substructure transformation complete.')
  #Essential pathway analysis for FA substructures
  
  set.seed(1)
  
  path_score_FA <- path_scoring(network = FA_network_new, sub_t = FA_sub_exp_t, 
                                calibrate = T, data_type = 'FA')
  path_score_FA_sel <- 
    path_score_FA %>% filter(Significant=='yes')
  
  
  
  path_data <- rbind(path_score_FA_sel %>% 
                       filter(Type=='Suppressed') %>% 
                       arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                     path_score_FA_sel %>% 
                       filter(Type=='Active') %>% 
                       .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
    filter(Significant=='yes')
  
  
  
  
  path_data_fig <- path_data %>% 
    mutate(path=factor(.$path, levels = .$path)) %>% 
    ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=Type))+
    geom_bar(stat='identity')+
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
    coord_flip()+
    theme_bw()+
    #scale_y_continuous(limits = c(-7,7))+
    scale_fill_manual(values = bluered(100)[c(1,100)])+
    theme(legend.position='none',
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size=8))+
    labs(x='', y='Path score',
         title='Top 5 representative pathways')
  
  print('Pathway analysis complete.')
  
  #Essential edges (reactions) analysis for FA substructures
  
  reaction_score_FA <- reaction_scoring(network = FA_network_new, 
                                        sub_exp = FA_sub_exp[[3]],
                                        sub_t = FA_sub_exp_t, 
                                        ctrl = ctrl, exp = exp, 
                                        Species = species)
  
  reaction_score_FA_sel <- 
    reaction_score_FA[,c("edge_name", "p_value", "mlog10p",
                         'perturbation_score' ,'edge_type', 'FA_change', 'genes')] %>% 
    filter(p_value<0.05)
  

  reaction_data <- rbind(reaction_score_FA %>% filter(perturbation_score>0) %>% .[1:5,],
                         reaction_score_FA %>% filter(perturbation_score<0) %>% 
                           arrange(perturbation_score) %>%.[1:5,]) %>% 
    filter(p_value<0.05)
  

  reaction_data <- reaction_data %>% 
    mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
    mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
           node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
    mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                              paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
    mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                              paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
    mutate(edge_color=paste0(node1_color,node2_color))
  

  
  reaction_data_fig <- reaction_data %>%
    mutate(Mode=factor(Mode, levels = c("Increase","Decrease"))) %>% 
    mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
    ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
               fill=Mode, color=Edge_direction))+
    geom_bar(stat='identity', size=0.8)+
    scale_y_discrete(
      labels=rev(reaction_data$edge_color)
    ) +
    geom_vline(xintercept = 0)+
    theme_bw()+
    theme(legend.position = 'right', 
          axis.text.y = element_markdown(),
          plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(values = rev(pal_lancet()(2)))+
    scale_color_manual(values = c('gold','white'))+
    labs(y='', fill='Reaction', x='Perturbation score', 
         title='Top 5 significant edges',color='Edge type')+
    guides(fill=guide_legend(order=1),
           color=guide_legend(order=2))
  
  
  print('Reaction analysis complete.')
  
  #FA biosynthetic network construction
  FA_network_data <- draw_network(network_data = FA_network_new,
                                  DE_data = FA_sub_exp_t,
                                  if_species = F, significant = 'adj_p_value',
                                  path_scoring_result = path_score_FA,
                                  reaction_scoring_result = reaction_score_FA,
                                  top_n = 5, path_type = 'both')
  
  
  network <- visNetwork(FA_network_data[[1]],FA_network_data[[2]]) %>% 
    visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                    physics = F, smooth = TRUE, randomSeed =5) 
  
  print('Network analysis complete.')
  
  return(list(path_score_FA_sel, path_data_fig, reaction_score_FA_sel,
              reaction_data_fig, network))
}

#-------------------Analysis for unprocessed data-------------------
# Parameters
## method (drop down menu): t.test, wilcox.test, mod.t.test
## ctrl: blank or ?
## exp: blank or ?
## exo_lipid: blank or ?
## species (drop down menu): human, mouse, rat

lipid_species_substructure_analysis <- function(exp_raw, method, ctrl, exp,
                                              non_missing_pct=0.3,
                                              exo_lipid=NULL,species='rat'){
  
  exp_data <- build_char_table(exp_raw, network_node = network_node)[[1]]
  
  char_data <- build_char_table(exp_raw, network_node = network_node)[[2]]
  
  no_sub_t <- unprocessed_data_test(exp_data = exp_data,
                                    char_table = char_data,
                                    method = method,
                                    significant='adj_p_value',
                                    ctrl_group = ctrl, exp_group = exp)
  
  #-------------------Lipid species substructure analysis-------------------
  

  #Decompose lipids into species substructures
  
  species_substructure <- species_sub_transform(char = char_data,
                                                lipid_substructure = lipid_substructure,
                                                network_node = network_node)
  
  #Extract species substructures using fold changes
  
  species_sub_stop <- species_sub_extract(lipid_substructure = species_substructure,
                                          unprocessed_data_result =  no_sub_t,
                                          type = 'species', pct_limit = non_missing_pct,
                                          exo_lipid=exo_lipid)
  
  #Transform lipid exp into substructure exp
  
  species_sub_exp <- lipid_sub_matrix(exp_data = exp_data, 
                                      sub_data = species_sub_stop,
                                      sub_type = 'species')
  
  
  #Differential expression analysis for substructures
  
  species_sub_exp_t <- t_test(data = species_sub_exp[[3]], ctrl = ctrl, exp = exp,
                              method = method, significant = 'adj_p_value')
  
  print('Substructure transformation complete.')
  
  #Species biosynthetic network data transformation
  
  species_network <- build_species_net(species_substructure = species_substructure)
  
  
  #Essential pathway analysis for species substructures
  
  set.seed(1)
  path_score <-  path_scoring(network = species_network,
                                      sub_t = species_sub_exp_t,
                                      calibrate = T, data_type = 'Species')
  
  
  path_score_sel <- path_score %>% filter(Significant=='yes')
  
  if(nrow(path_score_sel)==0){
    path_data_fig <- NULL
  }
  else{
    path_data <- rbind(path_score_sel %>% 
                         filter(Type=='Suppressed') %>% 
                         arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                       path_score_sel %>% 
                         filter(Type=='Active') %>% 
                         .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
      filter(Significant=='yes')
    
    
    path_data_fig <- path_data %>% 
      mutate(path=factor(.$path, levels = .$path)) %>% 
      ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=Type))+
      geom_bar(stat='identity')+
      geom_hline(yintercept = 0)+
      geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
      coord_flip()+
      theme_bw()+
      #scale_y_continuous(limits = c(-7,7))+
      scale_fill_manual(values = bluered(100)[c(1,100)])+
      theme(legend.position='none',
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size=8))+
      labs(x='', y='Path score',
           title='Top 5 representative pathways')
  }
  
  
  
  print('Pathway analysis complete.')
  
  
  #Essential edges (reactions) analysis for species substructures
  
  species_net_w_rev <- add_rev_rection(network_edge = network_edge,
                                       species_net = species_network)
  
  reaction_score <- reaction_scoring(network = species_net_w_rev,
                                             sub_exp = species_sub_exp[[3]],
                                             sub_t = species_sub_exp_t,
                                             ctrl=ctrl, exp=exp,
                                             Species = species)
  
  
  reaction_score_sel <- 
    reaction_score[,c("edge_name", "p_value", "mlog10p",
                      'perturbation_score' ,'edge_type', 'genes')] %>% 
    filter(p_value<0.05)
  
  
  if(nrow(reaction_score_sel)==0){
    reaction_data_fig <- NULL
  }
  else{
    reaction_data <- rbind(reaction_score %>% filter(perturbation_score>0) %>% .[1:5,],
                           reaction_score %>% filter(perturbation_score<0) %>% 
                             arrange(perturbation_score) %>%.[1:5,]) %>% 
      filter(p_value<0.05)
    
    
    reaction_data <- reaction_data %>% 
      mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
      mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
             node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
      mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                                paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
      mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                                paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
      mutate(edge_color=paste0(node1_color,node2_color))
    
    
    
    reaction_data_fig <- reaction_data %>%
      mutate(Mode=factor(Mode, levels = c("Increase","Decrease"))) %>% 
      mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
      ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
                 fill=Mode, color=Edge_direction))+
      geom_bar(stat='identity', size=0.8)+
      scale_y_discrete(
        labels=rev(reaction_data$edge_color)
      ) +
      geom_vline(xintercept = 0)+
      theme_bw()+
      theme(legend.position = 'right', 
            axis.text.y = element_markdown(),
            plot.title = element_text(hjust = 0.5))+
      scale_fill_manual(values = rev(pal_lancet()(2)))+
      scale_color_manual(values = c('gold','white'))+
      labs(y='', fill='Reaction', x='Perturbation score', 
           title='Top 5 significant edges',color='Edge type')+
      guides(fill=guide_legend(order=1),
             color=guide_legend(order=2))
  }
  print('Reaction analysis complete.')
  
  
  #Lipid species biosynthetic network construction
  
  species_network_data <- draw_network(network_data = species_net_w_rev,
                                       DE_data = species_sub_exp_t,
                                       if_species = T,significant = 'adj_p_value',
                                       path_scoring_result = path_score,
                                       reaction_scoring_result = reaction_score,
                                       top_n = 5, path_type = 'both')
  
  
  network <- visNetwork(species_network_data[[1]], species_network_data[[2]])
  print('Network analysis complete.')
  
  return(list(path_score_sel, path_data_fig, reaction_score_sel,
              reaction_data_fig,network))  
  
}

#-------------------Analysis for unprocessed data-------------------
# Parameters
## method (drop down menu): t.test, wilcox.test, mod.t.test
## ctrl: blank or ?
## exp: blank or ?
## exo_lipid: blank or ?
## species (drop down menu): human, mouse, rat

lipid_class_substructure_analysis <- function(exp_raw, method, ctrl, exp,
                                     exo_lipid=NULL,species='rat'){
  
  exp_data <- build_char_table(exp_raw, network_node = network_node)[[1]]
  
  char_data <- build_char_table(exp_raw, network_node = network_node)[[2]]
  
  no_sub_t <- unprocessed_data_test(exp_data = exp_data,
                                    char_table = char_data,
                                    method = method,
                                    significant='adj_p_value',
                                    ctrl_group = ctrl, exp_group = exp)

  #-------------------Lipid class substructure analysis-------------------
  
  #Extract class substructures using fold changes
  
  class_sub_stop <- species_sub_extract(lipid_substructure =lipid_substructure,
                                        unprocessed_data_result = no_sub_t,
                                        type = 'class', pct_limit = 0.01,
                                        exo_lipid=exo_lipid)
  
  #Transform lipid exp into substructures exp
  
  class_exp <- no_sub_t[[1]] %>% filter(type=='class') %>% 
    dplyr::select(-type)
  
  
  class_sub_exp <- lipid_sub_matrix(exp_data = class_exp, 
                                    sub_data = class_sub_stop,
                                    sub_type = 'Class')
  
  
  
  #Differential expression analysis for substructures
  
  class_sub_exp_t <- t_test(data = class_sub_exp[[3]], ctrl = ctrl, exp = exp,
                            method = 't.test', significant = 'adj_p_value')
  
  
  #Class biosynthetic network data transformation
  
  class_network <- network_edge[c('S1','P1')] %>% 
    filter(S1 %in% class_sub_exp_t$lipid, P1 %in% class_sub_exp_t$lipid)
  
  #Essential pathway analysis for class substructures
  
  set.seed(1)
  path_score <-  path_scoring(network = class_network,
                                    sub_t = class_sub_exp_t,
                                    calibrate = T, data_type = 'Class')
  
  
  path_score_sel <- path_score %>% filter(Significant=='yes')
  
  if(nrow(path_score_sel)==0){
    path_data_fig <- NULL
  }
  else{
    path_data <- rbind(path_score_sel %>% 
                         filter(Type=='Suppressed') %>% 
                         arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                       path_score_sel %>% 
                         filter(Type=='Active') %>% 
                         .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
      filter(Significant=='yes')
    
    
    
    path_data_fig <- path_data %>% 
      mutate(path=factor(.$path, levels = .$path)) %>% 
      ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=Type))+
      geom_bar(stat='identity')+
      geom_hline(yintercept = 0)+
      geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
      coord_flip()+
      theme_bw()+
      #scale_y_continuous(limits = c(-7,7))+
      scale_fill_manual(values = bluered(100)[c(1,100)])+
      theme(legend.position='none',
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size=8))+
      labs(x='', y='Path score',
           title='Top 5 representative pathways')
  }
  
  
  
  print('Pathway analysis complete.')
  
  
  #Essential edges (reactions) analysis for class substructures
  
  reaction_score <- reaction_scoring(network = class_network,
                                           sub_exp = class_sub_exp[[3]],
                                           sub_t = class_sub_exp_t,
                                           ctrl=ctrl, exp=exp,
                                           Species = species)
  
  reaction_score_sel <- 
    reaction_score[,c("edge_name", "p_value", "mlog10p",
                         'perturbation_score' ,'edge_type', 'genes')] %>% 
    filter(p_value<0.05)
  
  
  if(nrow(reaction_score_sel)==0){
    reaction_data_fig <- NULL
  }
  else{
    reaction_data <- rbind(reaction_score %>% filter(perturbation_score>0) %>% .[1:5,],
                           reaction_score %>% filter(perturbation_score<0) %>% 
                             arrange(perturbation_score) %>%.[1:5,]) %>% 
      filter(p_value<0.05)
    
    
    reaction_data <- reaction_data %>% 
      mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
      mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
             node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
      mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                                paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
      mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                                paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
      mutate(edge_color=paste0(node1_color,node2_color))
    
    
    
    reaction_data_fig <- reaction_data %>%
      mutate(Mode=factor(Mode, levels = c("Increase","Decrease"))) %>% 
      mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
      ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
                 fill=Mode, color=Edge_direction))+
      geom_bar(stat='identity', size=0.8)+
      scale_y_discrete(
        labels=rev(reaction_data$edge_color)
      ) +
      geom_vline(xintercept = 0)+
      theme_bw()+
      theme(legend.position = 'right', 
            axis.text.y = element_markdown(),
            plot.title = element_text(hjust = 0.5))+
      scale_fill_manual(values = rev(pal_lancet()(2)))+
      scale_color_manual(values = c('gold','white'))+
      labs(y='', fill='Reaction', x='Perturbation score', 
           title='Top 5 significant edges',color='Edge type')+
      guides(fill=guide_legend(order=1),
             color=guide_legend(order=2))
  }
  print('Reaction analysis complete.')
  
  
  #Lipid class biosynthetic network construction
  
  class_network_data <- draw_network(network_data = class_network,
                                     DE_data = class_sub_exp_t,
                                     if_species = F,significant = 'adj_p_value',
                                     path_scoring_result = path_score,
                                     reaction_scoring_result = reaction_score,
                                     top_n = 5, path_type = 'both')
  
  network <- visNetwork(class_network_data[[1]], class_network_data[[2]])
  print('Network analysis complete.')
  
  return(list(path_score_sel, path_data_fig, reaction_score_sel,
              reaction_data_fig,network))
}

#-------------------Check for correct data format-------------------
# Parameters
## exp_data: check format of this file

check_data_format <- function(exp_data){
  library(tidyverse)
  warning_message <- character()
  #check variable class
  ckeck_feature <- map_chr(exp_data, ~class(.x))[1]!='character'
  
  ckeck_value <- !map_chr(exp_data, ~class(.x))[-1] %in% c('numeric','integer')
  
  ckeck_feature_name <- colnames(exp_data)[1]!='feature'
  
  if(ckeck_feature || sum(ckeck_value)!=0 || ckeck_feature_name){
    warning_message <- "! Please ensure that the first column's name is 'feature' and it is a string variable, while the remaining columns are numeric variables."
  }
  
  lipid_class <- str_replace(str_extract(exp_data$feature, '(.+?(_))|(.+)'), '_','')
  

  
  colnames(exp_data)[1] <- 'feature'
  exp_data <- exp_data[,c(T,!ckeck_value)]
  exp_data <- exp_data %>% filter(!feature %in% c(c('G3P', 'DHAP')))
  
  if(sum(unique(lipid_class) %in% network_node$Abbreviation)<5 || nrow(exp_data)<30){
    warning_message <- c(warning_message, '! The dataset must consist of a minimum of 5 lipid classes and 30 lipid species.')
  }
  
  FA_format <- str_split(exp_data$feature, '_')
  FA_format1 <- FA_format %>% map_lgl(~length(.x)<2)
  FA_format2 <- FA_format %>% map_lgl(~sum(!str_detect(.x[-1],'[0-9]+:[0-9]+;[0-9]'))!=0)
  FA_format <- FA_format1 | FA_format2
  
  if(sum(FA_format)!=0){
    warning_text <- ' with wrong lipid format.'
    
    warning_message <- c(warning_message, str_c('! ',str_c(exp_data$feature[FA_format],
                                                      collapse = ', '), warning_text))
  }
  
  exp_data <- exp_data[!FA_format,]


  lipid_class_not_support <- unique(lipid_class)[!unique(lipid_class) %in% 
                                                   network_node$Abbreviation]
  
  if(length(lipid_class_not_support)!=0){
    warning_text <- ' are not supported by iLipidome'
    warning_message <- c(warning_message, str_c('! ',str_c(lipid_class_not_support,
                                                      collapse = ', '),warning_text))
  }
  exp_data <- exp_data[lipid_class %in%network_node$Abbreviation,]
  lipid_class <- lipid_class[lipid_class %in%network_node$Abbreviation]
  

  FA_num_ref <- network_node$FA[match(lipid_class, network_node$Abbreviation)]
  
  FA_num <- map_int(str_split(exp_data$feature, '_'), ~length(.x))
  
  FA_num_ckeck <- (FA_num_ref==FA_num-1)|(FA_num-1==1)

  
  if(sum(FA_num_ckeck)!=nrow(exp_data)){
    warning_text <- ' with wrong FA number.'
    
    warning_message <- c(warning_message, str_c('! ',str_c(exp_data$feature[!FA_num_ckeck],
                                                      collapse = ', '), warning_text))
  }
  if(length(warning_message)!=0){
    warning_message <- str_c(str_c(warning_message, collapse = '\n'),
          '\n! Please check the help page for instructions on how to format the data')
  }
  
  return(warning_message)
}