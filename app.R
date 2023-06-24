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
# library(ggsci)
# library(ggpubr)
# library(gridExtra)
# library(ggrepel)
library(ggtext)

# library(ComplexHeatmap)
# library(fpc)
# library(cowplot)

# library(ggplot2)

# library(DT)
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
        # imageOutput("logo"),
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
                                tabPanel("FAExample", 
                                    downloadButton("FAdownload", "Download FA Example Dataset"),
                                ),
                                tabPanel("FACustom",
                                    fileInput("FAfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c("text/csv", "text/comma-separated-values", ".csv")
                                    )
                                ),
                            ),
                            tabsetPanel(id = "FAparams", type = "hidden",
                                tabPanel("FAExample"
                                ),
                                tabPanel("FACustom",
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
                                        multiple = TRUE, 
                                        choices = c(
                                            "w3-22:6;0",
                                            "16:0;0",
                                            "18:0;0",
                                            "20:0;0",
                                            "22:0;0",
                                            "24:0;0",
                                            "26:0;0",
                                            "16:0;0",
                                            "w7-16:1;0",
                                            "16:0;0",
                                            "18:0;0",
                                            "w9-18:1;0",
                                            "w9-18:1;0",
                                            "w9-20:1;0",
                                            "w9-18:2;0",
                                            "w9-20:2;0",
                                            "w9-20:1;0",
                                            "w9-22:1;0",
                                            "w9-24:1;0",
                                            "w6-18:2;0",
                                            "w6-18:3;0",
                                            "w6-20:3;0",
                                            "w6-20:4;0",
                                            "w6-22:4;0",
                                            "w6-24:4;0",
                                            "w6-26:4;0",
                                            "w6-24:4;0",
                                            "w6-24:5;0",
                                            "w6-24:5;0",
                                            "w6-22:5;0",
                                            "w6-22:4;0",
                                            "w6-26:5;0",
                                            "w3-18:3;0",
                                            "w3-18:4;0",
                                            "w3-20:4;0",
                                            "w3-20:5;0",
                                            "w3-22:5;0",
                                            "w3-24:5;0",
                                            "w3-26:5;0",
                                            "w3-24:5;0",
                                            "w3-24:6;0",
                                            "w3-24:6;0",
                                            "w3-26:6;0",
                                            "w3-22:6;0",
                                            "w3-22:5;0",
                                            "4:0;0",
                                            "6:0;0",
                                            "8:0;0",
                                            "10:0;0",
                                            "12:0;0",
                                            "14:0;0",
                                            "16:0;0",
                                            "18:0;0",
                                            "20:0;0",
                                            "22:0;0",
                                            "24:0;0",
                                            "26:0;0",
                                            "28:0;0",
                                            "w7-16:1;0",
                                            "w7-18:1;0",
                                            "w10-16:1;0",
                                            "w9-18:1;0",
                                            "w9-18:2;0",
                                            "w9-20:1;0",
                                            "w9-20:2;0",
                                            "w9-20:2;0",
                                            "w9-20:3;0",
                                            "w9-22:1;0",
                                            "w9-24:1;0",
                                            "w9-26:1;0",
                                            "w6-18:3;0",
                                            "w6-20:3;0",
                                            "w6-20:4;0",
                                            "w6-22:4;0",
                                            "w6-24:4;0",
                                            "w6-26:4;0",
                                            "w6-28:4;0",
                                            "w6-24:5;0",
                                            "w6-22:5;0",
                                            "w6-26:5;0",
                                            "w6-24:5;0",
                                            "w6-22:5;0",
                                            "w6-28:5;0",
                                            "w3-18:4;0",
                                            "w3-20:4;0",
                                            "w3-20:5;0",
                                            "w3-22:5;0",
                                            "w3-24:5;0",
                                            "w3-26:5;0",
                                            "w3-28:5;0",
                                            "w3-24:6;0",
                                            "w3-22:6;0",
                                            "w3-26:6;0",
                                            "w3-28:6;0",
                                            "w3-24:6;0",
                                            "w3-22:6;0"
                                        )
                                    ),
                                    textInput("FAexolipid", 
                                        label = "exo lipid:", 
                                        value = "w9-18:2;0, w3-20:4;0"
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
                                ),
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
                                    span(textOutput("FA_nosig_path"), style="color:red"),
                                    DT::dataTableOutput("FAPathScoreDT"),
                                    plotOutput("FAPathScorePlot", width = "70%", height = "600"),
                                ),
                                tabPanel("Reaction Score",
                                    span(textOutput("FA_nosig_reaction"), style="color:red"),
                                    DT::dataTableOutput("FAReactionScoreDT"),
                                    plotOutput("FAReactionScorePlot", width = "70%", height = "600"),
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
                                tabPanel("LSExample",
                                    downloadButton("LSdownload", "Download LS Example Dataset"),
                                ),
                                tabPanel("LSCustom",
                                    fileInput("LSfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c("text/csv", "text/comma-separated-values", ".csv")
                                    )
                                ),
                            ),
                            tabsetPanel(id = "LSparams", type = "hidden",
                                tabPanel("LSExample"
                                ),
                                tabPanel("LSCustom",
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
                                    textInput("LSexolipid", 
                                        label = "exo lipid:", 
                                        value = NULL
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
                                ),
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
                                    span(textOutput("LS_nosig_path"), style="color:red"),
                                    DT::dataTableOutput("LSPathScoreDT"),
                                    plotOutput("LSPathScorePlot", width = "70%", height = "600"),
                                ),
                                tabPanel("Reaction Score",
                                    span(textOutput("LS_nosig_reaction"), style="color:red"),
                                    DT::dataTableOutput("LSReactionScoreDT"),
                                    plotOutput("LSReactionScorePlot", width = "70%", height = "600"),
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
                                tabPanel("LCExample",
                                    downloadButton("LCdownload", "Download LC Example Dataset"),
                                ),
                                tabPanel("LCCustom",
                                    fileInput("LCfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c("text/csv", "text/comma-separated-values", ".csv")
                                    )
                                ),
                            ),
                            tabsetPanel(id = "LCparams", type = "hidden",
                                tabPanel("LCExample"
                                ),
                                tabPanel("LCCustom",
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
                                    textInput("LCexolipid", 
                                        label = "exo lipid:", 
                                        value = NULL
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
                                ),
                            ),
                            checkboxInput("LCDownload",
                                "Download Tables and Plots"
                            ),
                            actionButton("LCRun", "run code", padding = "8px")
                        ),
                        mainPanel(
                            tabsetPanel( # tabPanels for visualizations
                                tabPanel("Input Data",
                                    span(textOutput("LC_error"), style = "color:red"),
                                    DT::dataTableOutput("LCInData")
                                ),
                                tabPanel("Path Score",
                                    span(textOutput("LC_nosig_path"), style = "color:red"),
                                    DT::dataTableOutput("LCPathScoreDT"),
                                    plotOutput("LCPathScorePlot", width = "70%", height = "600"),
                                ),
                                tabPanel("Reaction Score",
                                    span(textOutput("LC_nosig_reaction"), style = "color:red"),
                                    DT::dataTableOutput("LCReactionScoreDT"),
                                    plotOutput("LCReactionScorePlot", width = "70%", height = "600"),
                                ),
                                tabPanel("Network Graph",
                                    tags$div(
                                        visNetworkOutput("LCNetworkGraph", height = "700px"),
                                        style = "color:red"
                                    )
                                    
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
        # tabPanel("Download Example Datasets",
            
        # ),
        tabPanel("FAQ",

        ),
        tabPanel("Contact Us",
            p("If you would like to reach out to us, please send an email to <u104001424 at cmu dot edu dot tw> with the
                topic of the email in the subject line. "
            ),
            br(),
            p("  "),
            tags$div("If there is an issue with the website you would like to report, please create an issue ", 
                tags$a(href = "https://github.com/LewisLabUCSD/iLipidome-website/issues", "here.")
            )
        ),
    ),
    # disconnectMessage()
)

server <- function(input, output, session) {
    # observeEvent(input$disconnect, {
    #     session$close()
    # })

    # output$logo <- renderImage({
    #     list(
    #         src = file.path("www/logo.png"),
    #         contentType = "image/png",
    #         width = 80,
    #         height = 30
    #     )
    # }, deleteFile = FALSE)

    observeEvent(input$FAData, {
        updateTabsetPanel(inputId = "FAFileIn", selected = input$FAData)
        updateTabsetPanel(inputId = "FAparams", selected = input$FAData)
    })
    observeEvent(input$LSData, {
        updateTabsetPanel(inputId = "LSFileIn", selected = input$LSData)
        updateTabsetPanel(inputId = "LSparams", selected = input$LSData)
    })
    observeEvent(input$LCData, {
        updateTabsetPanel(inputId = "LCFileIn", selected = input$LCData)
        updateTabsetPanel(inputId = "LCparams", selected = input$LCData)
    })

    observeEvent(input$FARun, {
        # reset plots to default state
        output$FAInData <- NULL
        output$FAPathScoreDT <- NULL
        output$FAPathScorePlot <- NULL
        output$FAReactionScoreDT <- NULL
        output$FAReactionScorePlot <- NULL
        output$FANetworkGraph <- NULL

        FA_substructure_result <- NULL

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

        if (length(FA_format) == 0) {
            if (input$FAData == "FAExample") {
                FA_substructure_result <- FA_substructure_analysis(FA_exp_raw, method = 't.test',
                                                   ctrl = 1:7, exp = 8:13,
                                                   unmapped_FA = c('w9-18:2;0','w3-20:4;0'),
                                                   exo_lipid = 'w3-22:6;0', species = 'rat')
            }
            if (input$FAData == "FACustom") {
                FA_substructure_result <- FA_substructure_analysis(FA_exp_raw, method = input$FAMethod,
                                        ctrl = input$FActrl[1]:input$FActrl[2], exp = input$FAexp[1]:input$FAexp[2],
                                        unmapped_FA = input$FAunmappedFA,
                                        exo_lipid = str_trim(strsplit(input$FAexolipid, ",")[[1]]), species = input$FAspecies)
            }

            output$FAInData <- DT::renderDataTable({
                DT::datatable(format(FA_exp_raw, digits = 1, justify = "none"), options = list(orderClasses = TRUE))
            })

            if (FA_substructure_result[[1]] != "NA" && FA_substructure_result[[2]] != "NA") {
                output$FAPathScoreDT <- DT::renderDataTable({
                    DT::datatable(format(FA_substructure_result[[1]], digits = 1, justify = "none"), options = list(orderClasses = TRUE))
                })
                output$FAPathScorePlot <- renderPlot(plot(FA_substructure_result[[2]]))
            }
            else {
                output$FA_nosig_path <- renderText("No Significant Pathways Found")
            }
            
            if (FA_substructure_result[[3]] != "NA" && FA_substructure_result[[4]] != "NA") {
                output$FAReactionScoreDT <- DT::renderDataTable({
                    DT::datatable(format(FA_substructure_result[[3]], digits = 1, justify = "none"), options = list(orderClasses = TRUE))
                })
                output$FAReactionScorePlot <- renderPlot(plot(FA_substructure_result[[4]]))
            }
            else {
                output$FA_nosig_reaction <- renderText("No Significant Reactions Found")
            }

            output$FANetworkGraph <- renderVisNetwork(FA_substructure_result[[5]]) # ask what the output for the visnetwork would be if nothing is outputted
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

        LS_substructure_result <- NULL

        if (input$LSData == "LSExample") {
            LS_exp_raw <- read.csv("example_data/lipid_species_substructure_analysis/exp2.csv",
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

        if (length(LS_format) == 0) {
            if (input$LSData == "LSExample") {
                LS_substructure_result <- lipid_species_substructure_analysis(LS_exp_raw, method = 't.test',
                                                                         ctrl = 1:7, exp = 8:13,
                                                                         non_missing_pct = 0.3,
                                                                         exo_lipid = NULL, species = 'rat')
            }
            if (input$LSData == "LSCustom") {
                LS_substructure_result <- lipid_species_substructure_analysis(LS_exp_raw, method = input$LSMethod,
                                                                        ctrl = input$LSctrl[1]:input$LSctrl[2], exp = input$LSexp[1]:input$LSexp[2],
                                                                        non_missing_pct = input$LSnonMissingPCT,
                                                                        exo_lipid = input$LSexolipid, species = input$LSspecies) 
                                                                        #FIGURE OUT WHAT NON_MISSING_PCT IS AND HOW TO PASS NULL
            }

            output$LSInData <- DT::renderDataTable({
                DT::datatable(format(LS_exp_raw, digits = 1, justify = "none"), options = list(orderClasses = TRUE))
            })

            if (LS_substructure_result[[1]] != "NA" && LS_substructure_result[[2]] != "NA") {
                output$LSPathScoreDT <- DT::renderDataTable({
                    DT::datatable(format(LS_substructure_result[[1]], digits = 1, justify = "none"), options = list(orderClasses = TRUE))
                })
                output$LSPathScorePlot <- renderPlot(plot(LS_substructure_result[[2]]))
            }
            else {
                output$LS_nosig_path <- renderText("No Significant Pathways Found")
            }

            if (LS_substructure_result[[3]] != "NA" && LS_substructure_result[[4]] != "NA") {
                output$LSReactionScoreDT <- DT::renderDataTable({
                    DT::datatable(format(LS_substructure_result[[3]], digits = 1, justify = "none"), options = list(orderClasses = TRUE))
                })
                output$LSReactionScorePlot <- renderPlot(plot(LS_substructure_result[[4]]))
            }
            else {
                output$LS_nosig_reaction <- renderText("No Significant Reactions Found")
            }

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

        LC_substructure_result <- NULL
        
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

        if (length(LC_format) == 0) {
            if (input$LCData == "LCExample") {
                LC_substructure_result <- lipid_class_substructure_analysis(LC_exp_raw, method = 't.test',
                                                   ctrl = 1:7, exp = 8:13,
                                                   exo_lipid = NULL, species = 'rat')
            }
            if (input$LCData == "LCCustom") {
                LC_substructure_result <- lipid_class_substructure_analysis(LC_exp_raw, method = input$LCMethod,
                                                   ctrl = input$LCctrl[1]:input$LCctrl[2], exp = input$LCexp[1]:input$LCexp[2],
                                                   exo_lipid = NULL, species = input$LCspecies)
            }

            output$LCInData <- DT::renderDataTable({
                DT::datatable(format(LC_exp_raw, digits = 1, justify = "none"), options = list(orderClasses = TRUE))
            })

            if (length(LC_substructure_result[[1]]) > 1 && length(LC_substructure_result[[2]]) > 1) {
                output$LCPathScoreDT <- DT::renderDataTable({
                    DT::datatable(format(LC_substructure_result[[1]], digits = 1, justify = "none"), options = list(orderClasses = TRUE))
                })
                output$LCPathScorePlot <- renderPlot(plot(LC_substructure_result[[2]]))
            }
            else {
                output$LC_nosig_path <- renderText("No Significant Pathways Found")
            }
            
            if (length(LC_substructure_result[[3]]) > 1 && length(LC_substructure_result[[4]]) > 1) {
                output$LCReactionScoreDT <- DT::renderDataTable({
                    DT::datatable(format(LC_substructure_result[[3]], digits = 1, justify = "none"), options = list(orderClasses = TRUE))
                })
                output$LCReactionScorePlot <- renderPlot(plot(LC_substructure_result[[4]]))
            }
            else {
                output$LC_nosig_reaction <- renderText("No Significant Reactions Found")
            }

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

    output$FAdownload <- downloadHandler(
        filename = "example_data/FA_substructure_analysis/exp.csv",
        content = function(fileDownload) {
            file.copy("example_data/FA_substructure_analysis/exp.csv", fileDownload)
        }
    )

    output$LSdownload <- downloadHandler(
        filename = "example_data/lipid_species_substructure_analysis/exp.csv",
        content = function(fileDownload) {
            file.copy("example_data/lipid_species_substructure_analysis/exp2.csv", fileDownload)
        }
    )

    output$LCdownload <- downloadHandler(
        filename = "example_data/lipid_class_substructure_analysis/exp.csv",
        content = function(fileDownload) {
            file.copy("example_data/lipid_class_substructure_analysis/exp.csv", fileDownload)
        }
    )
}

shinyApp(ui, server)
