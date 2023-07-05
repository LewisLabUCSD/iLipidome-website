library(shiny)
library(shinyhelper)
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

install.packages('iLipidome_0.1.0.tar.gz', repos=NULL, type='source')
library(iLipidome)

source("functions.R")

# source("required_function.R")
# load('required_data.RData')

data_info <- "<ol><li>Lipid dataset can be uploaded by users or using example datasets. Data needs to be uploaded in CSV or TSV format. The maximum file size is 30MB.
                <li>Once the file is chosen and shown 'Upload complete' then press 'Run analysis'.</ol>"

ctrl_info <- "<ol><li>Enter the column number separated by comma to assign control and experimental groups.
                <li>Note that the first column contains the names of the lipids, so the grouping information should be counted from the second column onwards.
                <li>For instance, in a dataset where the first column contains lipid names, columns 2 to 4 represent the control group, and columns 5 to 7 indicate the experimental group, you should fill '1,2,3' in the Control group and '4,5,6' in the Experimental group.</ol>"

unmapped_info <- "<ol><li>Select the low-expressed fatty acid isomers for exclusion
                <li>Due to the limitations of mass spectrometry, the exact double bond locations for fatty acids are often not provided in most lipidomics data. Consequently, certain fatty acids may be mapped to multiple candidates in the fatty acid network (e.g., FA 20:4 could be omega-3 or omega-6). This parameter aids in the specific removal of low-expressed fatty acid isomers, thereby enabling more precise calculations.</ol>"

FA_exo_info <- "<ol><li>Select the exogenous fatty acids in the study to prevent substructure decomposition.
                <li>If an exogenous treatment is present in the study, it can significantly impact the calculation results. This parameter allows users to exclude the effects of exogenous treatment. </ol>"

lipid_exo_info <- "<ol><li>Enter the exogenous lipids separated by comma to prevent substructure decomposition. For example: 'PC_16:0;0_22:6;0,PC_18:0;0_22:6;0'
                    <li>The lipid names or classes must be present in the uploaded dataset.
                    <li>If an exogenous treatment is present in the study, it can significantly impact the calculation results. This parameter allows users to exclude the effects of exogenous treatment.</ol>"

pct_info <- "<ol><li>Enter a value between 0 and 1 to set the threshold for the percentage of non-missing values in a biosynthetic pathway. Increasing this value will result in fewer biosynthetic pathways being retained. This parameter enables users to regulate the substructure decomposition process, reducing artifacts that may arise from excessive decomposition.
                <li>Usually, values between 0.3 and 0.7 are commonly used for this parameter.</ol>"

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
                                    "Example dataset (Levental KR, et al. Nat Commun. 2020)" = "FAExample",
                                    "Upload your own data" = "FACustom"
                                )
                            ) %>%
                            helper(
                                type = "inline",
                                title = "Data Source:",
                                size = "l",
                                content = data_info
                            ),
                            tabsetPanel(id = "FAFileIn", type = "hidden",
                                tabPanel("FAExample", 
                                    downloadButton("FAdownload", "Download Example"),
                                ),
                                tabPanel("FACustom",
                                    fileInput("FAfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c(".csv", ".tsv")
                                    )
                                ),
                            ),
                            tabsetPanel(id = "FAparams", type = "hidden",
                                tabPanel("FAExample"
                                ),
                                tabPanel("FACustom",
                                    radioButtons("FAMethod", 
                                        label = "Method:",
                                        choices = c(
                                            "t-test" = "t.test",
                                            "Wilcoxon test" = "wilcox.test"
                                            # "mod.t.test" = "mod.t.test"
                                        ),
                                        selected = "t.test",
                                        inline = TRUE
                                    ),
                                    textInput("FActrl", 
                                        label = "Control Group:", 
                                        value = "1:7"
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Assign Group Information:",
                                        size = "l",
                                        content = ctrl_info
                                    ),
                                    textInput("FAexp", 
                                        label = "Experimental Group:",
                                        value = "8:13"
                                    ),
                                    selectInput("FAunmappedFA", 
                                        label = "Remove Low-expressed Fatty Acid Isomers:", 
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
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Remove low-expressed fatty acid isomers:",
                                        size = "l",
                                        content = unmapped_info
                                    ),
                                    textInput("FAexolipid", 
                                        label = "Remove Exogenous Lipid Effect:"
                                        # "w9-18:2;0, w3-20:4;0"
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Remove exogenous lipid effect:",
                                        size = "l",
                                        content = FA_exo_info
                                    ),
                                    radioButtons("FAspecies", 
                                        label = "Species:", 
                                        choices = c(
                                            "Human" = "human",
                                            "Mouse" = "mouse",
                                            "Rat" = "rat"
                                        ),
                                        selected = "rat",
                                        inline = TRUE
                                    ),
                                ),
                            ),
                            checkboxInput("FADownload",
                                "Download Tables and Plots"
                            ),
                            actionButton("FARun", "Run Analysis", padding = "8px")
                        ),
                        mainPanel(
                            tabsetPanel( # tabPanels for visualizations
                                tabPanel("Lipid Expression Data",
                                    span(textOutput("FA_error"), style="color:red"),
                                    DT::dataTableOutput("FAInData")
                                ),
                                tabPanel("Pathway Analysis",
                                    span(textOutput("FA_nosig_path"), style="color:red"),
                                    DT::dataTableOutput("FAPathScoreDT"),
                                    plotOutput("FAPathScorePlot", width = "70%", height = "600"),
                                ),
                                tabPanel("Reaction Analysis",
                                    span(textOutput("FA_nosig_reaction"), style="color:red"),
                                    DT::dataTableOutput("FAReactionScoreDT"),
                                    plotOutput("FAReactionScorePlot", width = "70%", height = "600"),
                                ),
                                tabPanel("Lipid Network",
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
                                    "Example dataset (Levental KR, et al. Nat Commun. 2020)" = "LSExample",
                                    "Upload your own data" = "LSCustom"
                                )
                            ) %>%
                            helper(
                                type = "inline",
                                title = "Data Source:",
                                size = "l",
                                content = data_info
                            ),
                            tabsetPanel(id = "LSFileIn", type = "hidden",
                                tabPanel("LSExample",
                                    downloadButton("LSdownload", "Download Example"),
                                ),
                                tabPanel("LSCustom",
                                    fileInput("LSfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c(".csv", ".tsv")
                                    )
                                ),
                            ),
                            tabsetPanel(id = "LSparams", type = "hidden",
                                tabPanel("LSExample"
                                ),
                                tabPanel("LSCustom",
                                    radioButtons("LSMethod", 
                                        label = "Method:",
                                        choices = c(
                                            "t-test" = "t.test",
                                            "Wilcoxon test" = "wilcox.test"
                                        ),
                                        selected = "t.test",
                                        inline = TRUE
                                    ),
                                    textInput("LSctrl", 
                                        label = "Control Group:", 
                                        value = "1:7"
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Assign Group Information:",
                                        size = "l",
                                        content = ctrl_info
                                    ),
                                    textInput("LSexp", 
                                        label = "Experimental Group:", 
                                        value = "8:13"
                                    ),
                                    numericInput("LSnonMissingPCT",
                                        label = "Proportion of Non-missing Values to Retain a Pathway", 
                                        value = 0.3,
                                        min = 0,
                                        max = 1,
                                        step = 0.1
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Percentage of non-missing values to retain a pathway:",
                                        size = "l",
                                        content = pct_info
                                    ), 
                                    textInput("LSexolipid", 
                                        label = "Remove Exogenous Lipid Effect:", 
                                        value = NULL
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Remove exogenous lipid effect:",
                                        size = "l",
                                        content = lipid_exo_info
                                    ),
                                    radioButtons("LSspecies", 
                                        label = "Species:", 
                                        choices = c(
                                            "Human" = "human",
                                            "Mouse" = "mouse",
                                            "Rat" = "rat"
                                        ),
                                        selected = "rat",
                                        inline = TRUE
                                    ),
                                ),
                            ),
                            checkboxInput("LSDownload",
                                "Download Tables and Plots"
                            ),
                            actionButton("LSRun", "Run Analysis", padding = "8px")
                        ),
                        mainPanel(
                            tabsetPanel( # tabPanels for visualizations
                                tabPanel("Lipid Expression Data",
                                    span(textOutput("LS_error"), style="color:red"),
                                    DT::dataTableOutput("LSInData")
                                ),
                                tabPanel("Pathway Analysis",
                                    span(textOutput("LS_nosig_path"), style="color:red"),
                                    DT::dataTableOutput("LSPathScoreDT"),
                                    plotOutput("LSPathScorePlot", width = "70%", height = "600"),
                                ),
                                tabPanel("Reaction Analysis",
                                    span(textOutput("LS_nosig_reaction"), style="color:red"),
                                    DT::dataTableOutput("LSReactionScoreDT"),
                                    plotOutput("LSReactionScorePlot", width = "70%", height = "600"),
                                ),
                                tabPanel("Lipid Network",
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
                                    "Example dataset (Levental KR, et al. Nat Commun. 2020)" = "LCExample",
                                    "Upload your own data" = "LCCustom"
                                )
                            ) %>%
                            helper(
                                type = "inline",
                                title = "Data Source:",
                                size = "l",
                                content = data_info
                            ),
                            tabsetPanel(id = "LCFileIn", type = "hidden",
                                tabPanel("LCExample",
                                    downloadButton("LCdownload", "Download Example"),
                                ),
                                tabPanel("LCCustom",
                                    fileInput("LCfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c(".csv", ".tsv")
                                    )
                                ),
                            ),
                            tabsetPanel(id = "LCparams", type = "hidden",
                                tabPanel("LCExample"
                                ),
                                tabPanel("LCCustom",
                                    radioButtons("LCMethod", 
                                        label = "Method:",
                                        choices = c(
                                            "t-test" = "t.test",
                                            "Wilcoxon test" = "wilcox.test"
                                        ),
                                        selected = "t.test",
                                        inline = TRUE
                                    ),
                                    textInput("LCctrl", 
                                        label = "Control Group:", 
                                        value = "1:7"
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Assign Group Information:",
                                        size = "l",
                                        content = ctrl_info
                                    ),
                                    textInput("LCexp", 
                                        label = "Experimental Group:", 
                                        value = "8:13"
                                    ),
                                    textInput("LCexolipid", 
                                        label = "Remove Exogenous Lipid Effect:", 
                                        value = NULL
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Remove exogenous lipid effect:",
                                        size = "l",
                                        content = lipid_exo_info
                                    ),
                                    radioButtons("LCspecies", 
                                        label = "Species:", 
                                        choices = c(
                                            "Human" = "human",
                                            "Mouse" = "mouse",
                                            "Rat" = "rat"
                                        ),
                                        selected = "rat",
                                        inline = TRUE
                                    ),
                                ),
                            ),
                            checkboxInput("LCDownload",
                                "Download Tables and Plots"
                            ),
                            actionButton("LCRun", "Run Analysis", padding = "8px")
                        ),
                        mainPanel(
                            tabsetPanel( # tabPanels for visualizations
                                tabPanel("Lipid Expression Data",
                                    span(textOutput("LC_error"), style = "color:red"),
                                    DT::dataTableOutput("LCInData")
                                ),
                                tabPanel("Pathway Analysis",
                                    span(textOutput("LC_nosig_path"), style = "color:red"),
                                    DT::dataTableOutput("LCPathScoreDT"),
                                    plotOutput("LCPathScorePlot", width = "70%", height = "600"),
                                ),
                                tabPanel("Reaction Analysis",
                                    span(textOutput("LC_nosig_reaction"), style = "color:red"),
                                    DT::dataTableOutput("LCReactionScoreDT"),
                                    plotOutput("LCReactionScorePlot", width = "70%", height = "600"),
                                ),
                                tabPanel("Lipid Network",
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

    observe_helpers()

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

    # function to pass to data processing functions to update progress bar
    update_progress <- function(progress = NULL, detail = NULL) {
        progress$inc(amount = 1/6, detail = detail)
    }

    observeEvent(input$FARun, {
        # reset plots and text indicators to default state
        output$FAInData <- NULL
        output$FAPathScoreDT <- NULL
        output$FAPathScorePlot <- NULL
        output$FAReactionScoreDT <- NULL
        output$FAReactionScorePlot <- NULL
        output$FANetworkGraph <- NULL

        FA_substructure_result <- NULL
        output$FA_error <- renderText("")

        # Create a Progress object
        FA_progress <- shiny::Progress$new()
        FA_progress$set(message = "Processing...", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(FA_progress$close())

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
        update_progress(progress = FA_progress, detail = "File uploaded, data processing started")

        FA_format <- check_data_format(FA_exp_raw)

        if (length(FA_format) == 0) {
            if (input$FAData == "FAExample") {
                FA_substructure_result <- FA_substructure_analysis(FA_exp_raw, method = 't.test',
                                                   ctrl = 1:7, exp = 8:13,
                                                   unmapped_FA = c('w9-18:2;0','w3-20:4;0'),
                                                   exo_lipid = 'w3-22:6;0', species = 'rat', 
                                                   progress = FA_progress, update_progress = update_progress)
            }
            if (input$FAData == "FACustom") {
                if (input$FAexolipid == "") { FAexo <- NULL }
                else { FAexo <- str_trim(strsplit(input$FAexolipid, ",")[[1]]) }

                FA_substructure_result <- FA_substructure_analysis(FA_exp_raw, method = input$FAMethod,
                                        ctrl = unlist(lapply(str_trim(unlist(strsplit(input$FActrl, ","))), function(x) eval(parse(text = x)))),
                                        exp = unlist(lapply(str_trim(unlist(strsplit(input$FAexp, ","))), function(x) eval(parse(text = x)))),
                                        unmapped_FA = input$FAunmappedFA,
                                        exo_lipid = FAexo, species = input$FAspecies,
                                        progress = FA_progress, update_progress = update_progress)
            }

            update_progress(progress = FA_progress, detail = "Displaying tables and plots...")

            output$FAInData <- DT::renderDataTable({
                FA_exp_raw %>% 
                    DT::datatable() %>% 
                    DT::formatRound(which(sapply(FA_exp_raw, is.numeric)), digits = 3)
            })

            if (FA_substructure_result[[1]] != "NA" && FA_substructure_result[[2]] != "NA") {
                output$FAPathScoreDT <- DT::renderDataTable({
                    FA_substructure_result[[1]] %>%
                        DT::datatable() %>%
                        DT::formatRound(which(sapply(FA_substructure_result[[1]], is.numeric)), digits = 3)
                })
                output$FAPathScorePlot <- renderPlot(plot(FA_substructure_result[[2]]))
            }
            else { # output of result is "NA"
                output$FA_nosig_path <- renderText("No Significant Pathways Found")
            }
            
            if (FA_substructure_result[[3]] != "NA" && FA_substructure_result[[4]] != "NA") {
                output$FAReactionScoreDT <- DT::renderDataTable({
                    FA_substructure_result[[3]] %>%
                        DT::datatable() %>%
                        DT::formatRound(which(sapply(FA_substructure_result[[3]], is.numeric)), digits = 3)
                })
                output$FAReactionScorePlot <- renderPlot(plot(FA_substructure_result[[4]]))
            }
            else { # output of result is "NA"
                output$FA_nosig_reaction <- renderText("No Significant Reactions Found")
            }

            output$FANetworkGraph <- renderVisNetwork(FA_substructure_result[[5]]) # ask what the output for the visnetwork would be if nothing is outputted
        }
        else { # display errors
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
        output$LS_error <- renderText("")

        # Create a Progress object
        LS_progress <- shiny::Progress$new()
        LS_progress$set(message = "Processing...", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(LS_progress$close())

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
        update_progress(progress = LS_progress, detail = "File uploaded, data processing started")

        LS_format <- check_data_format(LS_exp_raw)

        if (length(LS_format) == 0) {
            if (input$LSData == "LSExample") {
                LS_substructure_result <- lipid_species_substructure_analysis(LS_exp_raw, method = 't.test',
                                                                         ctrl = 1:7, exp = 8:13,
                                                                         non_missing_pct = 0.3,
                                                                         exo_lipid = NULL, species = 'rat',
                                                                         progress = LS_progress, update_progress = update_progress)
            }
            if (input$LSData == "LSCustom") {
                if (input$LSexolipid == "") { LSexo <- NULL }
                else { LSexo <- str_trim(strsplit(input$LSexolipid, ",")[[1]]) }

                LS_substructure_result <- lipid_species_substructure_analysis(LS_exp_raw, method = input$LSMethod,
                                                                        ctrl = unlist(lapply(str_trim(unlist(strsplit(input$LSctrl, ","))), function(x) eval(parse(text = x)))), 
                                                                        exp = unlist(lapply(str_trim(unlist(strsplit(input$LSexp, ","))), function(x) eval(parse(text = x)))),
                                                                        non_missing_pct = input$LSnonMissingPCT,
                                                                        exo_lipid = LSexo, species = input$LSspecies,
                                                                        progress = LS_progress, update_progress = update_progress) 
            }

            update_progress(progress = LS_progress, detail = "Displaying tables and plots...")

            output$LSInData <- DT::renderDataTable({
                LS_exp_raw %>%
                        DT::datatable() %>%
                        DT::formatRound(which(sapply(LS_exp_raw, is.numeric)), digits = 3)
            })

            if (LS_substructure_result[[1]] != "NA" && LS_substructure_result[[2]] != "NA") {
                output$LSPathScoreDT <- DT::renderDataTable({
                    LS_substructure_result[[1]] %>%
                        DT::datatable() %>%
                        DT::formatRound(which(sapply(LS_substructure_result[[1]], is.numeric)), digits = 3)
                })
                output$LSPathScorePlot <- renderPlot(plot(LS_substructure_result[[2]]))
            }
            else { # output of result is "NA"
                output$LS_nosig_path <- renderText("No Significant Pathways Found")
            }

            if (LS_substructure_result[[3]] != "NA" && LS_substructure_result[[4]] != "NA") {
                output$LSReactionScoreDT <- DT::renderDataTable({
                    LS_substructure_result[[3]] %>%
                        DT::datatable() %>%
                        DT::formatRound(which(sapply(LS_substructure_result[[3]], is.numeric)), digits = 3)
                })
                output$LSReactionScorePlot <- renderPlot(plot(LS_substructure_result[[4]]))
            }
            else { # output of result is "NA"
                output$LS_nosig_reaction <- renderText("No Significant Reactions Found")
            }

            output$LSNetworkGraph <- renderVisNetwork(LS_substructure_result[[5]])
        }
        else { # display errors
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
        output$LC_error <- renderText("")

        # Create a Progress object
        LC_progress <- shiny::Progress$new()
        LC_progress$set(message = "Processing...", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(LC_progress$close())
        
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

        update_progress(progress = LC_progress, detail = "File uploaded, data processing started")

        LC_format <- check_data_format(LC_exp_raw)

        if (length(LC_format) == 0) {
            if (input$LCData == "LCExample") {
                LC_substructure_result <- lipid_class_substructure_analysis(LC_exp_raw, method = 't.test',
                                                   ctrl = 1:7, exp = 8:13,
                                                   exo_lipid = NULL, species = 'rat',
                                                   progress = LC_progress, update_progress = update_progress)
            }
            if (input$LCData == "LCCustom") {
                if (input$LCexolipid == "") { LCexo <- NULL }
                else { LCexo <- str_trim(strsplit(input$LCexolipid, ",")[[1]]) }

                LC_substructure_result <- lipid_class_substructure_analysis(LC_exp_raw, method = input$LCMethod,
                                                   ctrl = unlist(lapply(str_trim(unlist(strsplit(input$LCctrl, ","))), function(x) eval(parse(text = x)))), 
                                                   exp = unlist(lapply(str_trim(unlist(strsplit(input$LCexp, ","))), function(x) eval(parse(text = x)))),
                                                   exo_lipid = LCexo, species = input$LCspecies,
                                                   progress = LC_progress, update_progress = update_progress)
            }

            update_progress(progress = LC_progress, detail = "Displaying tables and plots...")

            output$LCInData <- DT::renderDataTable({
                LC_exp_raw %>%
                        DT::datatable() %>%
                        DT::formatRound(which(sapply(LC_exp_raw, is.numeric)), digits = 3)
            })

            if (length(LC_substructure_result[[1]]) > 1 && length(LC_substructure_result[[2]]) > 1) {
                output$LCPathScoreDT <- DT::renderDataTable({
                    LC_substructure_result[[1]] %>%
                        DT::datatable() %>%
                        DT::formatRound(which(sapply(LC_substructure_result[[1]], is.numeric)), digits = 3)
                })
                output$LCPathScorePlot <- renderPlot(plot(LC_substructure_result[[2]]))
            }
            else { # output of result is "NA"
                output$LC_nosig_path <- renderText("No Significant Pathways Found")
            }
            
            if (length(LC_substructure_result[[3]]) > 1 && length(LC_substructure_result[[4]]) > 1) {
                output$LCReactionScoreDT <- DT::renderDataTable({
                    LC_substructure_result[[3]] %>%
                        DT::datatable() %>%
                        DT::formatRound(which(sapply(LC_substructure_result[[3]], is.numeric)), digits = 3)
                })
                output$LCReactionScorePlot <- renderPlot(plot(LC_substructure_result[[4]]))
            }
            else { # output of result is "NA"
                output$LC_nosig_reaction <- renderText("No Significant Reactions Found")
            }

            output$LCNetworkGraph <- renderVisNetwork(LC_substructure_result[[5]])
        }
        else { # display errors
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
