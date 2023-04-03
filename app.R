library(shiny)

library(tidyverse)
library(dplyr)
library(igraph)
library(visNetwork)
library(xlsx)
library(data.table)
library(MKmisc)
library(gplots)
library(gtools)

ui <- fluidPage(
    # app title
    navbarPage(
        "iLipidome",
        tabPanel("Try it out!", 
            fluidRow(
                column(12,
                    fileInput(
                            "file1",
                            "Choose files", 
                            multiple = FALSE, 
                            accept = c("text/csv", "text/comma-separated-values", ".csv")
                        ),
                    checkboxGroupInput(inputId = "things to run",
                                        label = "Select Analyses to run:",
                                        choices = c("option 1", "option 2", "option 3"),
                        ),
                    actionButton("runcode", "run code")
                )
            ),
            fluidRow(
                #plotOutput("test_plot", click = "plot_click")
                tableOutput("t1")
            )
        ),
        tabPanel("About",
            value = 4,
            fluidRow(
                column(6, 
                    p("about the project")
                ),
                column(6, 
                    p("authors")
                )
            ),
            fluidRow(
                h3("paper:"),
                br(),
                downloadButton("downloadData", "Download")
            ),
            fluidRow(
                tags$iframe(
                    style = "height:400px; width:100%; scrolling=yes",
                    src = "C:/Users/Evan/Code/iLipidome-website/iLipidome-paper.pdf"
                    #edit path
                )
            )
        )
    ),

    
)

server <- function(input, output) {


    build_char_table <- function(raw_data, network_node){
  browser()
  class <- raw_data$feature %>% str_extract('[A-Za-z]+( O-)*')
  class_included <- class %in% network_node$Abbreviation
  
  
  totallength <- raw_data$feature %>% str_extract_all('\\d+:') %>% 
    map_int(~str_sub(.x, end = -2) %>% as.integer() %>% sum)
  
  totaldb <- raw_data$feature %>% str_extract_all('\\d+;') %>% 
    map_int(~str_sub(.x, end = -2) %>% as.integer() %>% sum)
  
  totaloh <- raw_data$feature %>% str_extract_all(';\\d+') %>% 
    map_int(~str_sub(.x, start = 2) %>% as.integer() %>% sum)
  
  FA_sum <- str_c(totallength, ':', totaldb, ';', totaloh)
  
  FA_split <- raw_data$feature %>% str_extract_all('\\d+:\\d+;\\d+') %>% 
    map_chr(~str_c(.x, collapse = '_'))

    FA_num <- vector("numeric", length(FA_sum))
    #browser()
  #FA_num <- FA_sum

  char_table <- data.frame(feature=raw_data$feature, class=class, totallength=totallength,
                           totaldb=totaldb, totaloh=totaloh, FA_sum=FA_sum, FA_split=FA_split, FA_num=FA_num) %>% 
    .[class_included,] %>% left_join(network_node[c('Abbreviation','FA')], by=c('class'='Abbreviation'))
  
  #colnames(char_table)[8] <- 'FA_num'
  FA_exact <- map_int(str_split(char_table$FA_split, '_'), ~length(.x))==char_table$FA_num
  char_table$FA_split[!FA_exact] <- ''
  
  each_FA <- str_extract_all(char_table$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})
  
  char_table <- char_table %>% mutate(each_FA=each_FA)
  
  raw_data <- raw_data[class_included,]
  
  return(list(raw_data, char_table))
}

    # observeEvent(input$runcode,
    #     {
    #         output$test_plot <- renderPlot({
    #             plot(c(1, 2, 3), c(2, 3, 4))
    #         })
    #     }
    # )
    # currently crashes when no file is selected
    # fix pls
    observeEvent(input$runcode,
        {
            exp <- read.csv(input$file1$datapath,
                header = TRUE,
                sep = ",",
                quote = input$quote
            )
            output$table <- renderTable(
                head(exp)
            )
            load('~/Code/iLipidome/Documentation/required_data.RData')
            
            # do prelim data processing
            # req(input$file1)

            

            built <- build_char_table(raw_data = exp, network_node = network_node)
            
            exp_sel <- built[[1]]
            char_sel <- built[[2]]

            # no_sub_t <- unprocessed_data_test(
            #     exp_data = exp_sel,
            #     char_table = char_sel,
            #     method = 't.test',
            #     significant='adj_p_value',
            #     ctrl_group = 1:7, 
            #     exp_group = 8:13
            # )

            output$t1 = renderTable(
                char_sel
                # exp
            )
        }
    )

    output$downloadData <- downloadHandler(
        #change file names
        filename = "iLipidome-paper.pdf",
        content = function(fileDownload) {
            file.copy("C:/Users/Evan/Code/iLipidome-website/iLipidome-paper.pdf", fileDownload)
        }
    )
}

shinyApp(ui, server)