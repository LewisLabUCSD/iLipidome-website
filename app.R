library(shiny)
# library(shinydisconnect)

library(tidyverse)
library(dplyr)
library(igraph)
library(visNetwork)
# library(xlsx)
library(data.table)
# library(MKmisc)
library(gplots)
library(gtools)
library(ggsci)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(ggtext)

library(ComplexHeatmap)
library(fpc)
library(cowplot)

library(ggplot2)
# library(ggvenn)

# install.packages('iLipidome_0.1.0.tar.gz', repos=NULL, type='source')
# library(iLipidome)

source("functions.R")
source("required_function.R")

load('required_data.RData')

ui <- fluidPage(
    # app title
    navbarPage(
        "iLipidome",
        tabPanel("Tutorial",

        ),
        tabPanel("Lipid Substructure Analysis",
            tabsetPanel(
                tabPanel("Fatty Acid Analysis",
                    fluidRow(
                        column(4, 
                            radioButtons("FAData", "Data Source",
                                c(
                                    "Example dataset (<dataset name>)" = "FAExample",
                                    "Upload your own data" = "FACustom"
                                )
                            ),
                            tabsetPanel(id = "FAFileIn", type = "hidden",
                                tabPanel("FAExample"
                                ),
                                tabPanel("FACustom",
                                    fileInput("FAfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c("text/csv", "text/comma-separated-values", ".csv")
                                    )
                                ),
                            ),
                        ),
                        column(8, style = "background-color:#defae0; padding: 5px",
                            h3("Parameter Selection"),
                            radioButtons("FAMethod", 
                                label = "method:",
                                choices = c(
                                    "t.test" = "t.test",
                                    "wilcox.test" = "wilcox.test",
                                    "mod.t.test" = "mod.t.test"
                                ),
                                selected = "mod.t.test",
                                inline = TRUE
                            ),
                            sliderInput("FActrl", 
                                label = "ctrl:", 
                                min = 1,
                                max = 20,
                                value = c(1, 7)
                            ),
                            sliderInput("FAexp", 
                                label = "exp:", 
                                min = 1,
                                max = 20,
                                value = c(8, 13)
                            ),
                            textInput("FAunmappedFA", 
                                label = "unmapped FA:", 
                                "w9-18:2;0, w3-20:4;0"
                            ),
                            textInput("FAexolipid", 
                                label = "exo lipid:", 
                                "w3-22:6;0"
                            ),
                            radioButtons("FAspecies", 
                                label = "species:", 
                                choices = c(
                                    "human" = "human",
                                    "mouse" = "mouse",
                                    "rat" = "rat"
                                ),
                                selected = "rat",
                                inline = TRUE
                            ),
                        )
                    ),
                    fluidRow(
                        actionButton("FARun", "run code", padding = "8px")
                    ),
                    fluidRow( # for visualizations
                        tableOutput("FAPathScore"),
                        plotOutput("FAPathScorePlot"),
                        tableOutput("FAReactionScore"),
                        plotOutput("FAReactionScorePlot"),
                        visNetworkOutput("FANetwork")
                    )
                ),
                tabPanel("Lipid Species Analysis",
                    fluidRow(
                        column(4, 
                            radioButtons("LSData", "Data Source",
                                c(
                                    "Example dataset (<dataset name>)" = "LSExample",
                                    "Upload your own data" = "LSCustom"
                                )
                            ),
                            tabsetPanel(
                                id = "LSFileIn", type = "hidden",
                                tabPanel("LSExample"
                                ),
                                tabPanel("LSCustom",
                                    fileInput("LSfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c("text/csv", "text/comma-separated-values", ".csv")
                                    )
                                ),
                            )
                        ),
                        column(8, style = "background-color:#defae0;",
                            selectInput("LSParams", "Parameter Selection:", 
                                c(
                                    # put params for lipid species here
                                )
                            ),
                        )
                    ),
                    fluidRow(
                        actionButton("LSRun", "run code")
                    ),
                    fluidRow( # for visualizations

                    )
                ),
                tabPanel("Lipid Class Analysis",
                    fluidRow(
                        column(4, 
                            radioButtons("LCData", "Data Source",
                                c(
                                    "Example dataset (<dataset name>)" = "LCExample",
                                    "Upload your own data" = "LCCustom"
                                )
                            ),
                            tabsetPanel(
                                id = "LCFileIn",
                                type = "hidden",
                                tabPanel("LCExample"
                                ),
                                tabPanel("LCCustom",
                                    fileInput(
                                        "LCfile",
                                        "Choose file",
                                        multiple = FALSE,
                                        accept = c("text/csv", "text/comma-separated-values", ".csv")
                                    )
                                ),
                            )
                        ),
                        column(8, style = "background-color:#defae0;",
                            selectInput("LCParams", "Parameter Selection:", 
                                c(
                                    # put params for lipid class here
                                )
                            ),
                        )
                    ),
                    fluidRow(
                        actionButton("LCRun", "run code")
                    ),
                    fluidRow( # for visualizations

                    )
                )
            )
        ),
        navbarMenu("Other",
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
                        src = "iLipidome-paper.pdf"
                        #edit path
                    )
                )
            ),
            tabPanel("FAQ",

            ),
            tabPanel("Report an Issue",

            ),
        ),
    ),
    # disconnectMessage()
)

server <- function(input, output, session) {
    # observeEvent(input$disconnect, {
    #     session$close()
    # })

    observeEvent(input$FAData, {
        updateTabsetPanel(inputId = "FAFileIn", selected = input$FAData)
    })
    observeEvent(input$LSData, {
        updateTabsetPanel(inputId = "LSFileIn", selected = input$LSData)
    })
    observeEvent(input$LCData, {
        updateTabsetPanel(inputId = "LCFileIn", selected = input$LCData)
    })

    observeEvent(input$FARun, {
        if (input$FAData == "FAExample") {
            exp_raw <- read.csv("example_data/FA_substructure_analysis/exp.csv",
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )

            print(input$FAunmappedFA)
            print(str_trim(strsplit(input$FAunmappedFA, ",")[[1]]))
            print(input$FAexolipid)
            
            FA_substructure_result <- FA_substructure_analysis(exp_raw, method=input$FAMethod,
                                                   ctrl=input$FActrl[1]:input$FActrl[2], exp=input$FAexp[1]:input$FAexp[2],
                                                   unmapped_FA = str_trim(strsplit(input$FAunmappedFA, ",")[[1]]),
                                                   exo_lipid=input$FAexolipid, species=input$FAspecies)
        }
        if (input$FAData == "FACustom") {
            exp_raw <- read.csv(input$FAfile$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )

            FA_substructure_result <- FA_substructure_analysis(exp_raw, method=input$FAMethod,
                                        ctrl=input$FActrl[1]:input$FActrl[2], exp=input$FAexp[1]:input$FAexp[2],
                                        unmapped_FA = str_trim(strsplit(input$FAunmappedFA, ",")[[1]]),
                                        exo_lipid=input$FAexolipid, species=input$FAspecies)
        }

        output$FAPathScore <- renderTable(head(FA_substructure_result[[1]]))
        output$FAPathScorePlot <- renderPlot(plot(FA_substructure_result[[2]]))
        output$FAReactionScore <- renderTable(head(FA_substructure_result[[3]]))
        output$FAReactionScorePlot <- renderPlot(plot(FA_substructure_result[[4]]))
        output$FANetwork <- renderVisNetwork(FA_substructure_result[[5]])
    })
    observeEvent(input$LSRun, {
        
    })
    observeEvent(input$LCRun, {
        
    })

    output$downloadData <- downloadHandler(
        #change file paths
        filename = "iLipidome-paper.pdf",
        content = function(fileDownload) {
            file.copy("iLipidome-paper.pdf", fileDownload)
        }
    )
}

shinyApp(ui, server)
