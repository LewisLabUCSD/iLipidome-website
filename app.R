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

#install.packages('~/iLipidome_0.1.0.tar.gz', repos=NULL, type='source')

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
        # CVD file input
        fileInput(
            "fileCVDexp",
            "Choose EXP file",
            multiple = FALSE,
            accept = c("text/csv", "text/comma-separated-values", ".csv")
        ),
        fileInput(
            "fileCVDchar",
            "Choose CHAR file",
            multiple = FALSE,
            accept = c("text/csv", "text/comma-separated-values", ".csv")
        ),
    ),
    tabPanel("DHA", 
        # DHA file input
        fileInput(
            "fileDHAexp",
            "Choose EXP file",
            multiple = FALSE,
            accept = c("text/csv", "text/comma-separated-values", ".csv")
        ),
        fileInput(
            "fileDHAchar",
            "Choose CHAR file",
            multiple = FALSE,
            accept = c("text/csv", "text/comma-separated-values", ".csv")
        ),
    ),
    tabPanel("KO", 
        # KO file input 
        fileInput(
            "fileKOraw1",
            "Choose file RAW DATA 1",
            multiple = FALSE,
            accept = c("text/csv", "text/comma-separated-values", ".csv")
        ),
        fileInput(
            "fileKOraw2",
            "Choose file RAW DATA 2",
            multiple = FALSE,
            accept = c("text/csv", "text/comma-separated-values", ".csv")
        ),
    ),
    tabPanel("LPCAT", 
        # LPCAT file input
        fileInput(
            "fileLPCATchar",
            "Choose CHAR file",
            multiple = FALSE,
            accept = c("text/csv", "text/comma-separated-values", ".csv")
        ),
        fileInput(
            "fileLPCATexp",
            "Choose EXP file",
            multiple = FALSE,
            accept = c("text/csv", "text/comma-separated-values", ".csv")
        ),
        fileInput(
            "fileLPCATogname",
            "Choose ORIGINAL NAME file",
            multiple = FALSE,
            accept = c("text/csv", "text/comma-separated-values", ".csv")
        ),
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
                column(2,
                    radioButtons("datasetselector", "Select dataset:",
                        c(
                            "default" = "default",
                            "DHA" = "DHA",
                            "CVD" = "CVD",
                            "Lipid KO" = "KO",
                            "LPCAT" = "LPCAT"
                        )
                    ),
                    radioButtons("deffuncselector", "Select functions:",
                        c("func1" = "f1",
                        "func2" = "f2",
                        "func3" = "f3",
                        "func4" = "f4")
                    ),
                    file_inputs,
                    
                    actionButton("runcode", "run code"),
                    textOutput("loading")
                ),
                column(5,
                    visNetworkOutput("visNet1", height = "600px")
                ),
                column(5,
                    visNetworkOutput("visNet2", height = "600px"),
                    visNetworkOutput("visNet3", height = "600px")
                )
            )
        ),

        # commented out while i work on default with everything 

        # tabPanel("CVD",
        #     fluidRow(
        #         column(2,
        #             fileInput(
        #                 "fileCVDexp",
        #                 "Choose EXP file",
        #                 multiple = FALSE,
        #                 accept = c("text/csv", "text/comma-separated-values", ".csv")
        #             ),
        #             fileInput(
        #                 "fileCVDchar",
        #                 "Choose CHAR file",
        #                 multiple = FALSE,
        #                 accept = c("text/csv", "text/comma-separated-values", ".csv")
        #             ),
        #             actionButton("runcodeCVD", "run code"),
        #             textOutput("loadingCVD")
        #         ),
        #         column(5,
        #             visNetworkOutput("visNetCVD", height = "600px")
        #         ),
        #     )
        # ),
        # tabPanel("DHA",
        #     fluidRow(
        #         column(2,
        #             fileInput(
        #                 "fileDHAexp",
        #                 "Choose EXP file",
        #                 multiple = FALSE,
        #                 accept = c("text/csv", "text/comma-separated-values", ".csv")
        #             ),
        #             fileInput(
        #                 "fileDHAchar",
        #                 "Choose CHAR file",
        #                 multiple = FALSE,
        #                 accept = c("text/csv", "text/comma-separated-values", ".csv")
        #             ),
        #             actionButton("runcodeDHA", "run code"),
        #             textOutput("loadingDHA")
        #         ),
        #         column(5,
        #             visNetworkOutput("visNetDHA", height = "600px")
        #         ),
        #     )
        # ),
        # tabPanel("Lipid gene KO",
        #     fluidRow(
        #         column(2,
        #             fileInput(
        #                 "fileKOraw1",
        #                 "Choose file RAW DATA 1",
        #                 multiple = FALSE,
        #                 accept = c("text/csv", "text/comma-separated-values", ".csv")
        #             ),
        #             fileInput(
        #                 "fileKOraw2",
        #                 "Choose file RAW DATA 2",
        #                 multiple = FALSE,
        #                 accept = c("text/csv", "text/comma-separated-values", ".csv")
        #             ),
        #             actionButton("runcodeKO", "run code"),
        #             textOutput("loadingKO")
        #         ),
        #         column(5,
        #                 visNetworkOutput("visNetKO", height = "600px")
        #         ),
        #     )
        # ),
        # tabPanel("LPCAT1",
        #     fluidRow(
        #         column(2,
        #             fileInput(
        #                 "fileLPCATchar",
        #                 "Choose CHAR file",
        #                 multiple = FALSE,
        #                 accept = c("text/csv", "text/comma-separated-values", ".csv")
        #             ),
        #             fileInput(
        #                 "fileLPCATexp",
        #                 "Choose EXP file",
        #                 multiple = FALSE,
        #                 accept = c("text/csv", "text/comma-separated-values", ".csv")
        #             ),
        #             fileInput(
        #                 "fileLPCATogname",
        #                 "Choose ORIGINAL NAME file",
        #                 multiple = FALSE,
        #                 accept = c("text/csv", "text/comma-separated-values", ".csv")
        #             ),
        #             actionButton("runcodeLPCAT", "run code"),
        #             textOutput("loadingLPCAT")
        #         ),
        #         column(5,
        #                 visNetworkOutput("visNetLPCAT", height = "600px")
        #         ),
        #     )
        # )
    ),
    disconnectMessage()
    # actionButton("disconnect", "Disconnect the app")
)

server <- function(input, output, session) {
    observeEvent(input$disconnect, {
        session$close()
    })

    observeEvent(input$datasetselector, {
        updateTabsetPanel(inputId = "fileInputs", selected = input$datasetselector)
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

                        incProgress(1/2, detail = paste("done"))
                    })
                    
                    output$loading <- renderText("")
                }
            }
            if (input$datasetselector == "DHA") {
                if (is.null(input$fileDHAexp$datapath) || is.null(input$fileDHAchar$datapath)) {
                    output$loadingDHA <- renderText("please select a file")
                }
                if (!is.null(input$fileDHAexp$datapath) && !is.null(input$fileDHAchar$datapath)) {
                    withProgress(message = "doing computation", value = 0, {

                        exp <- read.csv(input$fileDHAexp$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
                        char <- read.csv(input$fileDHAchar$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
                        
                        DHAres <- processDHA(exp, char)

                        output$loadingDHA <- renderText("done")

                    })
                }
            }
            if (input$datasetselector == "CVD") {
                if (is.null(input$fileCVDexp$datapath) || is.null(input$fileCVDchar$datapath)) {
                    output$loadingCVD <- renderText("please select a file")
                }
                if (!is.null(input$fileCVDexp$datapath) && !is.null(input$fileCVDchar$datapath)) {
                    withProgress(message = "doing computation", value = 0, {
                        exp <- read.csv(input$fileCVDexp$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\"")
                        char <- read.csv(input$fileCVDchar$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )

                        incProgress(1/5, detail = paste("loaded data"))

                        CVDres <- processCVD(exp, char)

                        output$loadingCVD <- renderText("done")
                    })
                }
            }
            if (input$datasetselector == "KO") {
                if (is.null(input$fileKOraw1$datapath) || is.null(input$fileKOraw2$datapath)) {
                    output$loadingKO <- renderText("please select a file")
                }
                if (!is.null(input$fileKOraw1$datapath) && !is.null(input$fileKOraw2$datapath)) {
                    withProgress(message = "doing computation", value = 0, {

                        raw1 <- read.csv(input$fileKOraw1$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\"")
                        raw2 <- read.csv(input$fileKOraw2$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
                        
                        KOres <- processKO(raw1, raw2)

                        output$loadingKO <- renderText("done")
                    })
                }
            }
            if (input$datasetselector == "LPCAT") {
                if (is.null(input$fileLPCATexp$datapath) || is.null(input$fileLPCATchar$datapath) || is.null(input$fileLPCATogname$datapath)) {
                    output$loadingLPCAT <- renderText("please select a file")
                }
                if (!is.null(input$fileLPCATexp$datapath) && !is.null(input$fileLPCATchar$datapath) && !is.null(input$fileLPCATogname$datapath)) {
                    withProgress(message = "doing computation", value = 0, {
                        exp <- read.csv(input$fileLPCATexp$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\"")
                        char <- read.csv(input$fileLPCATchar$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
                        ogname <- read.csv(input$fileLPCATogname$datapath,
                            header = TRUE,
                            sep = ",",
                            quote = "\""
                        )
                        
                        LPCATres <- processLPCAT(exp, char, ogname)

                        output$loadingLPCAT <- renderText("done")
                    })
                }
            }
        }
    )

    # observeEvent(input$runcode,
    #     {
    #         if (is.null(input$file1$datapath)) {
    #             output$loading <- renderText("please select a file")
    #         }
    #         if (!is.null(input$file1$datapath)) {
    #             # do a switch case depending on which functions and demos are selected

    #             withProgress(message = "doing computation", value = 0, {

    #                 exp <- read.csv(input$file1$datapath,
    #                     header = TRUE,
    #                     sep = ",",
    #                     quote = "\""
    #                 )

    #                 incProgress(1/2, detail = paste("processing"))

    #                 def1 <- default1(exp)

    #                 output$visNet1 <- renderVisNetwork({
    #                     visNetwork(def1[[1]], def1[[2]]) %>%
    #                         visIgraphLayout(layout = "layout_with_sugiyama", type='square',
    #                                         physics = F, smooth = TRUE, randomSeed =5)
    #                 })

    #                 def2 <- default2(exp)

    #                 output$visNet2 <- renderVisNetwork({
    #                     visNetwork(def2[[1]], def2[[2]])
    #                 })

    #                 def3 <- default3(exp)

    #                 output$visNet3 <- renderVisNetwork({
    #                     visNetwork(def3[[1]], def3[[2]])
    #                 })

    #                 incProgress(1/2, detail = paste("done"))
    #             })
                
    #             output$loading <- renderText("")
    #         }

    #     }
    # )

    # observeEvent(input$runcodeCVD,
    #     {
    #         #CVD processing
    #         if (is.null(input$fileCVDexp$datapath) || is.null(input$fileCVDchar$datapath)) {
    #             output$loadingCVD <- renderText("please select a file")
    #         }
    #         if (!is.null(input$fileCVDexp$datapath) && !is.null(input$fileCVDchar$datapath)) {
    #             withProgress(message = "doing computation", value = 0, {
    #                 exp <- read.csv(input$fileCVDexp$datapath,
    #                     header = TRUE,
    #                     sep = ",",
    #                     quote = "\"")
    #                 char <- read.csv(input$fileCVDchar$datapath,
    #                     header = TRUE,
    #                     sep = ",",
    #                     quote = "\""
    #                 )

    #                 # load('data/required_data.RData')

    #                 incProgress(1/5, detail = paste("loaded data"))

    #                 CVDres <- processCVD(exp, char)

    #                 output$loadingCVD <- renderText("done")
    #             })
    #         }
    #     }
    # )
    # observeEvent(input$runcodeDHA,
    #     {
    #         if (is.null(input$fileDHAexp$datapath) || is.null(input$fileDHAchar$datapath)) {
    #             output$loadingDHA <- renderText("please select a file")
    #         }
    #         if (!is.null(input$fileDHAexp$datapath) && !is.null(input$fileDHAchar$datapath)) {
    #             withProgress(message = "doing computation", value = 0, {

    #                 exp <- read.csv(input$fileDHAexp$datapath,
    #                     header = TRUE,
    #                     sep = ",",
    #                     quote = "\""
    #                 )
    #                 char <- read.csv(input$fileDHAchar$datapath,
    #                     header = TRUE,
    #                     sep = ",",
    #                     quote = "\""
    #                 )

    #                 # load('data/required_data.RData')
                    
    #                 DHAres <- processDHA(exp, char)

    #                 output$loadingDHA <- renderText("done")

    #             })
    #         }
    #     }
    # )
    # observeEvent(input$runcodeKO,
    #     {
    #         if (is.null(input$fileKOraw1$datapath) || is.null(input$fileKOraw2$datapath)) {
    #             output$loadingKO <- renderText("please select a file")
    #         }
    #         if (!is.null(input$fileKOraw1$datapath) && !is.null(input$fileKOraw2$datapath)) {
    #             withProgress(message = "doing computation", value = 0, {

    #                 raw1 <- read.csv(input$fileKOraw1$datapath,
    #                     header = TRUE,
    #                     sep = ",",
    #                     quote = "\"")
    #                 raw2 <- read.csv(input$fileKOraw2$datapath,
    #                     header = TRUE,
    #                     sep = ",",
    #                     quote = "\""
    #                 )

    #                 # load('data/required_data.RData')
                    
    #                 KOres <- processKO(raw1, raw2)

    #                 output$loadingKO <- renderText("done")
    #             })
    #         }
    #     }
    # )
    # observeEvent(input$runcodeLPCAT,
    #     {
    #        #LPCAT processing
    #         if (is.null(input$fileLPCATexp$datapath) || is.null(input$fileLPCATchar$datapath) || is.null(input$fileLPCATogname$datapath)) {
    #             output$loadingLPCAT <- renderText("please select a file")
    #         }
    #         if (!is.null(input$fileLPCATexp$datapath) && !is.null(input$fileLPCATchar$datapath) && !is.null(input$fileLPCATogname$datapath)) {
    #             withProgress(message = "doing computation", value = 0, {
    #                 exp <- read.csv(input$fileLPCATexp$datapath,
    #                     header = TRUE,
    #                     sep = ",",
    #                     quote = "\"")
    #                 char <- read.csv(input$fileLPCATchar$datapath,
    #                     header = TRUE,
    #                     sep = ",",
    #                     quote = "\""
    #                 )
    #                 ogname <- read.csv(input$fileLPCATogname$datapath,
    #                     header = TRUE,
    #                     sep = ",",
    #                     quote = "\""
    #                 )

    #                 # load('data/required_data.RData')
                    
    #                 LPCATres <- processLPCAT(exp, char, ogname)

    #                 output$loadingLPCAT <- renderText("done")
    #             })
    #         }
    #     }
    # )

    output$downloadData <- downloadHandler(
        #change file paths
        filename = "iLipidome-paper.pdf",
        content = function(fileDownload) {
            file.copy("C:/Users/Evan/Code/iLipidome-website/iLipidome-paper.pdf", fileDownload)
        }
    )
}

shinyApp(ui, server)
