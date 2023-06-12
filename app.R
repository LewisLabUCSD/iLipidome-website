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

library(DT)
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
                    sidebarLayout(
                        sidebarPanel(
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
                            h3("Parameter Selection"),
                            radioButtons("FAMethod", 
                                label = "method:",
                                choices = c(
                                    "t.test" = "t.test",
                                    "wilcox.test" = "wilcox.test"
                                    # "mod.t.test" = "mod.t.test"
                                ),
                                selected = "t.test",
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
                            selectInput("FAunmappedFA", 
                                label = "unmapped FA:", 
                                choices = c(
                                    "w9-18:2;0, w3-20:4;0"
                                )
                            ),
                            selectInput("FAexolipid", 
                                label = "exo lipid:", 
                                choices = c(
                                    "w3-22:6;0"
                                )
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
                            checkboxInput("FADownload",
                                "Download Tables and Plots"
                            ),
                            actionButton("FARun", "run code", padding = "8px")
                        ),
                        mainPanel(
                            tabsetPanel( # tabPanels for visualizations
                                tabPanel("Input Data",
                                    span(textOutput("FA_error"), style="color:red"),
                                    DT::dataTableOutput("FAInData")
                                ),
                                tabPanel("Path Score",
                                    DT::dataTableOutput("FAPathScoreDT"),
                                    plotOutput("FAPathScorePlot"),
                                ),
                                tabPanel("Reaction Score",
                                    DT::dataTableOutput("FAReactionScoreDT"),
                                    plotOutput("FAReactionScorePlot"),
                                ),
                                tabPanel("Network Graph",
                                    visNetworkOutput("FANetworkGraph", height = "700px")
                                )
                            )
                        )
                    )
                ),
                tabPanel("Lipid Species Analysis",
                    sidebarLayout(
                        sidebarPanel(
                            radioButtons("LSData", "Data Source",
                                c(
                                    "Example dataset (<dataset name>)" = "LSExample",
                                    "Upload your own data" = "LSCustom"
                                )
                            ),
                            tabsetPanel(id = "LSFileIn", type = "hidden",
                                tabPanel("LSExample"
                                ),
                                tabPanel("LSCustom",
                                    fileInput("LSfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c("text/csv", "text/comma-separated-values", ".csv")
                                    )
                                ),
                            ),
                            h3("Parameter Selection"),
                            radioButtons("LSMethod", 
                                label = "method:",
                                choices = c(
                                    "t.test" = "t.test",
                                    "wilcox.test" = "wilcox.test"
                                ),
                                selected = "t.test",
                                inline = TRUE
                            ),
                            sliderInput("LSctrl", 
                                label = "ctrl:", 
                                min = 1,
                                max = 20,
                                value = c(1, 7)
                            ),
                            sliderInput("LSexp", 
                                label = "exp:", 
                                min = 1,
                                max = 20,
                                value = c(8, 13)
                            ),
                            numericInput("LSnonMissingPCT", #ASK ABOUT THIS
                                label = "non missing pct", 
                                value = 0.3,
                                min = 0,
                                max = 1
                            ),
                            selectInput("LSexolipid", 
                                label = "exo lipid:", 
                                choices = c(
                                    NULL
                                )
                            ),
                            radioButtons("LSspecies", 
                                label = "species:", 
                                choices = c(
                                    "human" = "human",
                                    "mouse" = "mouse",
                                    "rat" = "rat"
                                ),
                                selected = "rat",
                                inline = TRUE
                            ),
                            checkboxInput("LSDownload",
                                "Download Tables and Plots"
                            ),
                            actionButton("LSRun", "run code", padding = "8px")
                        ),
                        mainPanel(
                            tabsetPanel( # tabPanels for visualizations
                                tabPanel("Input Data",
                                    span(textOutput("LS_error"), style="color:red"),
                                    DT::dataTableOutput("LSInData")
                                ),
                                tabPanel("Path Score",
                                    DT::dataTableOutput("LSPathScoreDT"),
                                    plotOutput("LSPathScorePlot"),
                                ),
                                tabPanel("Reaction Score",
                                    DT::dataTableOutput("LSReactionScoreDT"),
                                    plotOutput("LSReactionScorePlot"),
                                ),
                                tabPanel("Network Graph",
                                    visNetworkOutput("LSNetworkGraph", height = "700px")
                                )
                            )
                        )
                    )
                ),
                tabPanel("Lipid Class Analysis",
                    sidebarLayout(
                        sidebarPanel(
                            radioButtons("LCData", "Data Source",
                                c(
                                    "Example dataset (<dataset name>)" = "LCExample",
                                    "Upload your own data" = "LCCustom"
                                )
                            ),
                            tabsetPanel(id = "LCFileIn", type = "hidden",
                                tabPanel("LCExample"
                                ),
                                tabPanel("LCCustom",
                                    fileInput("LCfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c("text/csv", "text/comma-separated-values", ".csv")
                                    )
                                ),
                            ),
                            h3("Parameter Selection"),
                            radioButtons("LCMethod", 
                                label = "method:",
                                choices = c(
                                    "t.test" = "t.test",
                                    "wilcox.test" = "wilcox.test"
                                ),
                                selected = "t.test",
                                inline = TRUE
                            ),
                            sliderInput("LCctrl", 
                                label = "ctrl:", 
                                min = 1,
                                max = 20,
                                value = c(1, 7)
                            ),
                            sliderInput("LCexp", 
                                label = "exp:", 
                                min = 1,
                                max = 20,
                                value = c(8, 13)
                            ),
                            selectInput("LCexolipid", 
                                label = "exo lipid:", 
                                choices = c(
                                    NULL
                                )
                            ),
                            radioButtons("LCspecies", 
                                label = "species:", 
                                choices = c(
                                    "human" = "human",
                                    "mouse" = "mouse",
                                    "rat" = "rat"
                                ),
                                selected = "rat",
                                inline = TRUE
                            ),
                            checkboxInput("LCDownload",
                                "Download Tables and Plots"
                            ),
                            actionButton("LCRun", "run code", padding = "8px")
                        ),
                        mainPanel(
                            tabsetPanel( # tabPanels for visualizations
                                tabPanel("Input Data",
                                    span(textOutput("LC_error"), style="color:red"),
                                    DT::dataTableOutput("LCInData")
                                ),
                                tabPanel("Path Score",
                                    DT::dataTableOutput("LCPathScoreDT"),
                                    plotOutput("LCPathScorePlot"),
                                ),
                                tabPanel("Reaction Score",
                                    DT::dataTableOutput("LCReactionScoreDT"),
                                    plotOutput("LCReactionScorePlot"),
                                ),
                                tabPanel("Network Graph",
                                    visNetworkOutput("LCNetworkGraph", height = "700px")
                                )
                            )
                        )
                    )
                )
            ),
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
                downloadButton("downloadPaper", "Download")
            ),
            fluidRow(
                tags$iframe(
                    style = "height:400px; width:100%; scrolling=yes",
                    src = "iLipidome-paper.pdf"
                    #edit path
                )
            )
        ),
        tabPanel("Download Example Datasets",
            downloadButton("downloadExampleData", "Download Example Dataset"),
        ),
        tabPanel("FAQ",

        ),
        tabPanel("Contact Us",
            p("If you would like to reach out to use, please send an email to <email address> with the
                topic of the email in the subject line. "
            )
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
        # reset plots to default state
        output$FAInData <- NULL
        output$FAPathScoreDT <- NULL
        output$FAPathScorePlot <- NULL
        output$FAReactionScoreDT <- NULL
        output$FAReactionScorePlot <- NULL
        output$FANetworkGraph <- NULL

        if (input$FAData == "FAExample") {
            FA_exp_raw <- read.csv("example_data/FA_substructure_analysis/exp.csv",
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
        }
        if (input$FAData == "FACustom") {
            FA_exp_raw <- read.csv(input$FAfile$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
        }

        FA_format <- check_data_format(FA_exp_raw)
        print(FA_format)

        if (length(FA_format) == 0) {
            FA_substructure_result <- FA_substructure_analysis(FA_exp_raw, method=input$FAMethod,
                                        ctrl=input$FActrl[1]:input$FActrl[2], exp=input$FAexp[1]:input$FAexp[2],
                                        unmapped_FA = str_trim(strsplit(input$FAunmappedFA, ",")[[1]]),
                                        exo_lipid=input$FAexolipid, species=input$FAspecies)
            
            output$FAInData <- DT::renderDataTable({
                DT::datatable(FA_exp_raw, options = list(orderClasses = TRUE))
            })
            output$FAPathScoreDT <- DT::renderDataTable({
                DT::datatable(FA_substructure_result[[1]], options = list(orderClasses = TRUE))
            })
            output$FAPathScorePlot <- renderPlot(plot(FA_substructure_result[[2]]))

            output$FAReactionScoreDT <- DT::renderDataTable({
                DT::datatable(FA_substructure_result[[3]], options = list(orderClasses = TRUE))
            })
            output$FAReactionScorePlot <- renderPlot(plot(FA_substructure_result[[4]]))

            output$FANetworkGraph <- renderVisNetwork(FA_substructure_result[[5]])

            # if download box is selected, download tables and plots
            
        }
        else {
            #figure out how to show FA_format error to user
            output$FA_error <- renderText(FA_format)
        }
    })
    observeEvent(input$LSRun, {
        # reset plots to default state
        output$LSInData <- NULL
        output$LSPathScoreDT <- NULL
        output$LSPathScorePlot <- NULL
        output$LSReactionScoreDT <- NULL
        output$LSReactionScorePlot <- NULL
        output$LSNetworkGraph <- NULL

        if (input$LSData == "LSExample") {
            LS_exp_raw <- read.csv("example_data/lipid_species_substructure_analysis/exp.csv",
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
        }
        if (input$LSData == "LSCustom") {
            LS_exp_raw <- read.csv(input$LSfile$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
        }

        LS_format <- check_data_format(LS_exp_raw)
        print(LS_format)

        if (length(LS_format) == 0) {
            LS_substructure_result <- lipid_species_substructure_analysis(LS_exp_raw, method=input$LSMethod,
                                                                        ctrl=input$LSctrl[1]:input$LSctrl[2], exp=input$LSexp[1]:input$LSexp[2],
                                                                        non_missing_pct = input$LSnonMissingPCT,
                                                                        exo_lipid=input$LSexolipid, species=input$LSspecies) 
                                                                        #FIGURE OUT WHAT NON_MISSING_PCT IS AND HOW TO PASS NULL

            output$LSInData <- DT::renderDataTable({
                DT::datatable(LS_exp_raw, options = list(orderClasses = TRUE))
            })
            output$LSPathScoreDT <- DT::renderDataTable({
                DT::datatable(LS_substructure_result[[1]], options = list(orderClasses = TRUE))
            })
            output$LSPathScorePlot <- renderPlot(plot(LS_substructure_result[[2]]))

            output$LSReactionScoreDT <- DT::renderDataTable({
                DT::datatable(LS_substructure_result[[3]], options = list(orderClasses = TRUE))
            })
            output$LSReactionScorePlot <- renderPlot(plot(LS_substructure_result[[4]]))

            output$LSNetworkGraph <- renderVisNetwork(LS_substructure_result[[5]])
        }
        else {
            output$LS_error <- renderText(LS_format)
        }
    })
    observeEvent(input$LCRun, {
        # reset plots to default state
        output$LCInData <- NULL
        output$LCPathScoreDT <- NULL
        output$LCPathScorePlot <- NULL
        output$LCReactionScoreDT <- NULL
        output$LCReactionScorePlot <- NULL
        output$LCNetworkGraph <- NULL
        
        if (input$LCData == "LCExample") {
            LC_exp_raw <- read.csv("example_data/lipid_class_substructure_analysis/exp.csv",
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
        }
        if (input$LCData == "LCCustom") {
            LC_exp_raw <- read.csv(input$LCfile$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
        }

        LC_format <- check_data_format(LC_exp_raw)
        print(LC_format)

        if (length(LC_format) == 0) {
            LC_substructure_result <- lipid_class_substructure_analysis(LC_exp_raw, method=input$LCMethod,
                                                   ctrl=input$LCctrl[1]:input$LCctrl[2], exp=input$LCexp[1]:input$LCexp[2],
                                                   exo_lipid=NULL, species=input$LCspecies)

            output$LCInData <- DT::renderDataTable({
                DT::datatable(LC_exp_raw, options = list(orderClasses = TRUE))
            })
            output$LCPathScoreDT <- DT::renderDataTable({
                DT::datatable(LC_substructure_result[[1]], options = list(orderClasses = TRUE))
            })
            output$LCPathScorePlot <- renderPlot(plot(LC_substructure_result[[2]]))

            output$LCReactionScoreDT <- DT::renderDataTable({
                DT::datatable(LC_substructure_result[[3]], options = list(orderClasses = TRUE))
            })
            output$LCReactionScorePlot <- renderPlot(plot(LC_substructure_result[[4]]))

            output$LCNetworkGraph <- renderVisNetwork(LC_substructure_result[[5]])
        }
        else {
            output$LC_error <- renderText(LC_format)
        }
    })

    output$downloadPaper <- downloadHandler(
        #change file paths
        filename = "iLipidome-paper.pdf",
        content = function(fileDownload) {
            file.copy("iLipidome-paper.pdf", fileDownload)
        }
    )

    output$downloadExampleData <- downloadHandler(
        filename = "example_data/FA_substructure_analysis/exp.csv",
        content = function(fileDownload) {
            file.copy("example_data/FA_substructure_analysis/exp.csv", fileDownload)
        }
    )
}

shinyApp(ui, server)
