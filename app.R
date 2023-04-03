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
        tabPanel("About"),
        tabPanel("Try it out!")
    ),

    sidebarLayout(
        sidebarPanel(
            fileInput(
                "file1",
                "Choose files", 
                multiple = FALSE, 
                accept = c("text/csv", "text/comma-separated-values", ".csv")
            ),
            tags$hr(),
            checkboxGroupInput(
                "things to run", 
                "things to run", 
                choices = c(
                    "option 1",
                    "option 2",
                    "option 3"
                ),
            ),
            actionButton(
                "run",
                "run code"
            )
            
        ),
        mainPanel(
            tableOutput("table")
        )
    ),
)

server <- function(input, output) {
    

    observeEvent(input$run,
        {
            load('~/Code/iLipidome/Documentation/required_data.RData')
            
            # do prelim data processing
            req(input$file1)

            exp <- read.csv(input$file1$datapath,
                header = TRUE,
                sep = input$sep,
                quote = input$quote
            )

            # built <- build_char_table(raw_data = exp, network_node = network_node)
            
            # exp_sel <- built[[1]]
            # char_sel <- built[[2]]

            # no_sub_t <- unprocessed_data_test(
            #     exp_data = exp_sel,
            #     char_table = char_sel,
            #     method = 't.test',
            #     significant='adj_p_value',
            #     ctrl_group = 1:7, 
            #     exp_group = 8:13
            # )

            output$table = renderTable(
                head(exp)
                # exp
            )
        }
    )
}

shinyApp(ui, server)