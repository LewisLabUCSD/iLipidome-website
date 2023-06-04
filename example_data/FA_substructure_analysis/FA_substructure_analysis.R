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
  
  FA_sub_stop <- FA_sub_extract(char_table = char_sel,
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
              reaction_data_fig,network))
}



FA_substructure_result <- FA_substructure_analysis(exp_raw, method='mod.t.test',
                                                   ctrl=1:7, exp=8:13,
                                                   unmapped_FA = c('w9-18:2;0','w3-20:4;0'),
                                                   exo_lipid='w3-22:6;0', species='rat')



