#-------------------Source function and required data-------------------

file <- '/Users/linwj/waxlos987@gmail.com - Google Drive/我的雲端硬碟/ucsd/Research/manuscript/submission/main text/revise2/code/'

# Source function
suppressWarnings(suppressPackageStartupMessages(source(file.path(file,'Required_function/required_function.R'))))

# Load required data
load(file.path(file,'Required_data/required_data.RData'))

#-------------------Data upload and process-------------------


exp_raw <- read.csv(file.path(file, 'Documentation/exp.csv'))


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



#-------------------Lipid class substructure analysis-------------------



lipid_class_substructure_result <- lipid_class_substructure_analysis(exp_raw, method='t.test',
                                                   ctrl=1:7, exp=8:13,
                                                   exo_lipid=NULL, species='rat')



