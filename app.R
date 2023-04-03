library(shiny)

# library(tidyverse)
# library(dplyr)
# library(igraph)
# library(visNetwork)
# library(xlsx)
# library(data.table)
# library(MKmisc)
# library(gplots)
# library(gtools)

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
                plotOutput("test_plot", click = "plot_click")
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

    observeEvent(input$runcode,
        {
            output$test_plot <- renderPlot({
                plot(c(1, 2, 3), c(2, 3, 4))
            })
        }
    )
    # currently crashes when no file is selected
    # fix pls
    # observeEvent(input$runcode,
    #     {
    #         exp <- read.csv(input$file1$datapath,
    #             header = TRUE,
    #             sep = input$sep,
    #             quote = input$quote
    #         )
    #         output$table <- renderTable(
    #             head(exp)
    #         )
    #         load('~/Code/iLipidome/Documentation/required_data.RData')
            
    #         # do prelim data processing
    #         # req(input$file1)

            

    #         # built <- build_char_table(raw_data = exp, network_node = network_node)
            
    #         # exp_sel <- built[[1]]
    #         # char_sel <- built[[2]]

    #         # no_sub_t <- unprocessed_data_test(
    #         #     exp_data = exp_sel,
    #         #     char_table = char_sel,
    #         #     method = 't.test',
    #         #     significant='adj_p_value',
    #         #     ctrl_group = 1:7, 
    #         #     exp_group = 8:13
    #         # )

    #         output$table = renderTable(
    #             head(exp)
    #             # exp
    #         )
    #     }
    # )

    output$downloadData <- downloadHandler(
        #change file names
        filename = "iLipidome-paper.pdf",
        content = function(fileDownload) {
            file.copy("C:/Users/Evan/Code/iLipidome-website/iLipidome-paper.pdf", fileDownload)
        }
    )
}

shinyApp(ui, server)