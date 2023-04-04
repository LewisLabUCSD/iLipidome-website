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

library(iLipidome)
#library(test1)

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
                #tableOutput("t1")
                visNetworkOutput("visNet")
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
            #browser()

            if (!is.null(input$file1$datapath)) {
              exp <- read.csv(input$file1$datapath,
                              header = TRUE,
                              sep = ",",
                              quote = "\""
              )
              output$table <- renderTable(
                head(exp)
              )
              load('~/iLipidome/data/required_data.RData')

              # do prelim data processing
              # req(input$file1)



              built <- build_char_table(raw_data = exp, network_node = network_node)

              exp_sel <- built[[1]]
              char_sel <- built[[2]]

              no_sub_t <- unprocessed_data_test(
                exp_data = exp_sel,
                char_table = char_sel,
                method = 't.test',
                significant='adj_p_value',
                ctrl_group = 1:7,
                exp_group = 8:13
              )

              # no_sub_t[[1]] %>% head()

              # no_sub_t[[2]] %>% head()

              #browser()

              FA_network_new <- build_FA_net(FA_network = FA_network,
                                             unprocessed_data_result = no_sub_t)

              FA_substructure <- FA_sub_transform(FA_network = FA_network_new,
                                                  unprocessed_data_result = no_sub_t,
                                                  unmapped_FA = c('w9-18:2;0','w3-20:4;0'))

              FA_sub_stop <- FA_sub_extract(char_table = char_sel,
                                            FA_substructure = FA_substructure,
                                            unprocessed_data_result = no_sub_t,
                                            exact_FA='no', exo_lipid='w3-22:6;0')

              FA_sub_exp <- lipid_sub_matrix(exp_data = exp_sel, sub_data = FA_sub_stop,
                                             sub_type = 'FA')

              FA_sub_exp_t <- t_test(data = FA_sub_exp[[3]], ctrl = 1:7, exp = 8:13,
                                     method = 't.test', significant = 'adj_p_value')

              set.seed(1)
              path_score_FA <- path_scoring(network = FA_network_new, sub_t = FA_sub_exp_t,
                                            calibrate = T, data_type = 'FA')

              reaction_score_FA <- reaction_scoring(network = FA_network_new,
                                                    sub_exp = FA_sub_exp[[3]],
                                                    sub_t = FA_sub_exp_t,
                                                    ctrl = 1:7, exp = 8:13,
                                                    Species = 'rat')

              FA_network_data <- draw_network(network_data = FA_network_new,
                                              DE_data = FA_sub_exp_t,
                                              if_species = F, significant = 'adj_p_value',
                                              path_scoring_result = path_score_FA,
                                              reaction_scoring_result = reaction_score_FA,
                                              top_n = 5, path_type = 'both')

              output$visNet <- renderVisNetwork({
                visNetwork(FA_network_data[[1]],FA_network_data[[2]]) %>%
                  visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                                  physics = F, smooth = TRUE, randomSeed =5)
              })
            }

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
