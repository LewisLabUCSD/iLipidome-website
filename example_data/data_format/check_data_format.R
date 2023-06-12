#-------------------Source function and required data-------------------

# file <- '/Users/linwj/Desktop/iLipidome_web/Code'

# Load required data
load('../..required_data.RData')

#-------------------Data upload and process-------------------


exp_raw <- read.csv('exp.csv')
exp_wrong <- exp_raw
exp_wrong$feature[1:8] <- c('CL_12:0', 'CL_12:0;0_13:1;3', 'CL=12:0','CL0_12:0;0', 
                            'CL__12:0;0','G3P','/', '_GPVD_13L_13:0;2')
colnames(exp_wrong)[1] <- 'sad'

exp_wrong$Ctrl2 <- 'a'

check_data_format <- function(exp_data){
  library(tidyverse)
  warning_message <- character()
  #check variable class
  ckeck_feature <- map_chr(exp_data, ~class(.x))[1]!='character'
  
  ckeck_value <- !map_chr(exp_data, ~class(.x))[-1] %in% c('numeric','integer')
  
  ckeck_feature_name <- colnames(exp_data)!='feature'
  
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

#correct
check_data_format(exp_raw)

#show error message and stop
check_data_format(exp_wrong) %>% cat()


