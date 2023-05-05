library(shiny)
library(shinydisconnect)

library(tidyverse)
library(dplyr)
library(igraph)
library(visNetwork)
library(xlsx)
library(data.table)
library(MKmisc)
library(gplots)
library(gtools)
library(ggsci)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(ggtext)

# install.packages('~/iLipidome_0.1.0.tar.gz', repos=NULL, type='source')

library(iLipidome)

source("functions.R")

load('data/required_data.RData')

file_inputs <- tabsetPanel(
    id = "fileInputs",
    type = "hidden",
    tabPanel("default",
        # default file input
        fileInput(
                "file1",
                "Choose file",
                multiple = FALSE,
                accept = c("text/csv", "text/comma-separated-values", ".csv")
        )
    ),
    tabPanel("CVD", 
    ),
    tabPanel("DHA", 
    ),
    tabPanel("KO", 
    ),
    tabPanel("LPCAT", 
    ),
    tabPanel("test",
    ),
)

plots <- tabsetPanel(
    id = "plots",
    type = "hidden",
    tabPanel("default",
        visNetworkOutput("visNet1"),
        visNetworkOutput("visNet2"),
        visNetworkOutput("visNet3"),
    ),
    tabPanel("CVD", 
    ),
    tabPanel("DHA", 
    ),
    tabPanel("KO", 
    ),
    tabPanel("LPCAT", 
        #c(SF3a, F5c, F5d, F5b, SF5b, SF3b, F5e, SF3c, SF3d)
        plotOutput("SF3a"),
        plotOutput("F5c"),
        plotOutput("F5d"),
        visNetworkOutput("F5b"),
        visNetworkOutput("SF5b"),
        plotOutput("SF3b"),
        plotOutput("F5e"),
        plotOutput("SF3c"),
        plotOutput("SF3d"),
    ),
    tabPanel("test",
    ),
)

ui <- fluidPage(
    # app title
    navbarPage(
        "iLipidome",
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
        ),
        tabPanel("default",
            fluidRow(
                column(4,
                    radioButtons("datasetselector", "Select dataset:",
                        c(
                            "default" = "default",
                            "DHA" = "DHA",
                            "CVD" = "CVD",
                            "Lipid KO" = "KO",
                            "LPCAT" = "LPCAT",
                            "test" = "test"
                        )
                    ),
                    file_inputs,
                    radioButtons("deffuncselector", "Select functions:",
                        c("func1" = "f1",
                        "func2" = "f2",
                        "func3" = "f3",
                        "func4" = "f4")
                    ),
                    checkboxInput("downloadPDF", "Download Plots as PDF", FALSE),
                    actionButton("runcode", "run code"),
                    textOutput("loading")
                ),
                column(8,
                    plots
                ),
            )
        ),
    ),
    disconnectMessage()
)

server <- function(input, output, session) {
    observeEvent(input$disconnect, {
        session$close()
    })

    observeEvent(input$datasetselector, {
        updateTabsetPanel(inputId = "fileInputs", selected = input$datasetselector)
        updateTabsetPanel(inputId = "plots", selected = input$datasetselector)
    }) 

    observeEvent(input$runcode,
        {
            if (input$datasetselector == "default") {
                if (is.null(input$file1$datapath)) {
                    output$loading <- renderText("please select a file") 
                    # in the future, could change this to an r shiny error, instead of this jank custom one
                }
                if (!is.null(input$file1$datapath)) {
                    withProgress(message = "doing computation", value = 0, {
                        exp <- read.csv(input$file1$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
                        incProgress(1/2, detail = paste("processing"))

                        def1 <- default1(exp)
                        output$visNet1 <- renderVisNetwork({
                            visNetwork(def1[[1]], def1[[2]]) %>%
                                visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                                                physics = F, smooth = TRUE, randomSeed =5)
                        })

                        def2 <- default2(exp)
                        output$visNet2 <- renderVisNetwork({
                            visNetwork(def2[[1]], def2[[2]])
                        })

                        def3 <- default3(exp)
                        output$visNet3 <- renderVisNetwork({
                            visNetwork(def3[[1]], def3[[2]])
                        })

                        # add pdf logic here
                        # if (input$dow)

                        incProgress(1/2, detail = paste("done"))
                    })
                    
                    output$loading <- renderText("")
                }
            }
            if (input$datasetselector == "DHA") {
                withProgress(message = "doing computation", value = 0, {
                    exp <- read.csv("~/Code/iLipidome-website/DHA/lipidome_data/exp_DHA_raw.csv",
                        header = TRUE,
                        sep = ",",
                        quote = "\""
                    )
                    char <- read.csv("~/Code/iLipidome-website/DHA/lipidome_data/char_DHA.csv",
                        header = TRUE,
                        sep = ",",
                        quote = "\""
                    )
                    
                    DHAres <- processDHA(exp, char)
                    output$loadingDHA <- renderText("done")
                })
            }

            if (input$datasetselector == "CVD") {
                withProgress(message = "doing computation", value = 0, {
                    exp <- read.csv("~/Code/iLipidome-website/CVD/lipidome_data/exp_CVD_raw.csv",
                        header = TRUE,
                        sep = ",",
                        quote = "\"")
                    char <- read.csv("~/Code/iLipidome-website/CVD/lipidome_data/char_CVD.csv",
                        header = TRUE,
                        sep = ",",
                        quote = "\""
                    )

                    incProgress(1/5, detail = paste("loaded data"))

                    CVDres <- processCVD(exp, char)

                    #c(F6b, SF4c, SF4e, SF5c, SF4b, SF4d, SF4a, F6c, F6d1, F6d2, F6d3, F6d4, SF4g, SF4f))

                    output$loadingCVD <- renderText("done")
                })
            }
            if (input$datasetselector == "KO") {
                withProgress(message = "doing computation", value = 0, {

                    raw1 <- read.csv("~/Code/iLipidome-website/Lipid_gene_KO/lipidome_data/raw_data_deletion_t.csv",
                        header = TRUE,
                        sep = ",",
                        quote = "\"")
                    raw2 <- read.csv("~/Code/iLipidome-website/Lipid_gene_KO/lipidome_data/raw_data_trap_t.csv",
                        header = TRUE,
                        sep = ",",
                        quote = "\""
                    )
                    
                    KOres <- processKO(raw1, raw2)
                    output$loadingKO <- renderText("done")
                })
            }
            if (input$datasetselector == "LPCAT") {
                withProgress(message = "doing computation", value = 0, {
                    exp <- read.csv("~/Code/iLipidome-website/LPCAT1/lipidome_data/exp_LPCAT1_raw.csv",
                        header = TRUE,
                        sep = ",",
                        quote = "\""
                    )
                    char <- read.csv("~/Code/iLipidome-website/LPCAT1/lipidome_data/char_LPCAT1_raw.csv",
                        header = TRUE,
                        sep = ",",
                        quote = "\""
                    )
                    ogname <- read.csv("~/Code/iLipidome-website/LPCAT1/lipidome_data/lipid_original_name.csv",
                        header = TRUE,
                        sep = ",",
                        quote = "\""
                    )

                    LPCATres <- processLPCAT(exp, char, ogname)

                    #c(SF3a, F5c, F5d, F5b, SF5b, SF3b, F5e, SF3c, SF3d)

                    output$SF3a <- renderPlot(
                        plot(LPCATres[[1]])
                    )
                    output$F5c <- renderPlot(
                        plot(LPCATres[[2]])
                    )
                    output$F5d <- renderPlot(
                        plot(LPCATres[[3]])
                    )
                    # asdf <- LPCATres[[4]]
                    # output$F5b <- renderVisNetwork({
                    #         visNetwork(asdf[[1]], asdf[[2]])
                    #     })
                    output$SF3b <- renderPlot(
                        plot(LPCATres[[6]])
                    )
                    output$F5e <- renderPlot(
                        plot(LPCATres[[7]])
                    )
                    output$SF3c <- renderPlot(
                        plot(LPCATres[[8]])
                    )
                    output$SF3d <- renderPlot(
                        plot(LPCATres[[9]])
                    )

                    output$loadingLPCAT <- renderText("done")
                })
            }
            if (input$datasetselector == "test") {
            }
        }
    )

    output$downloadData <- downloadHandler(
        #change file paths
        filename = "iLipidome-paper.pdf",
        content = function(fileDownload) {
            file.copy("C:/Users/Evan/Code/iLipidome-website/iLipidome-paper.pdf", fileDownload)
        }
    )
}

shinyApp(ui, server)
