library(shiny)
library(shinyhelper)

library(tidyverse)
library(dplyr)
library(igraph)
library(visNetwork)
library(data.table)
library(gplots)
library(gtools)
library(ggtext)

library(iLipidome)

source("functions.R")

data_info <- "<ol><li>Lipid dataset can be uploaded by users or using example datasets. Data needs to be uploaded in CSV or TSV format. The maximum file size is 30MB.
                <li>Once the file is chosen and shown 'Upload complete' then press 'Run analysis'.</ol>"

ctrl_info <- "<ol><li>Enter individual column numbers (1, 2, 3), or ranges of numbers (1:5), separated by commas to assign control and experimental groups.
                <li>Note that the first column contains the names of the lipids, so the grouping information should be counted from the second column onwards.
                <li>For instance, in a dataset where the first column contains lipid names, columns 2 to 4 represent the control group, and columns 5 to 7 
                indicate the experimental group, you should fill '1,2,3' in the Control group and '4,5,6' in the Experimental group.</ol>"

unmapped_info <- "<ol><li>Select the low-expressed fatty acid isomers for exclusion
                <li>Due to the limitations of mass spectrometry, the exact double bond locations for fatty acids are often not provided in most lipidomics data. 
                Consequently, certain fatty acids may be mapped to multiple candidates in the fatty acid network (e.g., FA 20:4 could be omega-3 or omega-6). This 
                parameter aids in the specific removal of low-expressed fatty acid isomers, thereby enabling more precise calculations.</ol>"

FA_exo_info <- "<ol><li>Select the exogenous fatty acids in the study to prevent substructure decomposition.
                <li>If an exogenous treatment is present in the study, it can significantly impact the calculation results. This parameter allows users to exclude 
                the effects of exogenous treatment. </ol>"

lipid_exo_info <- "<ol><li>Enter the exogenous lipids separated by comma to prevent substructure decomposition. For example: 'PC_16:0;0_22:6;0,PC_18:0;0_22:6;0'
                    <li>The lipid names or classes must be present in the uploaded dataset.
                    <li>If an exogenous treatment is present in the study, it can significantly impact the calculation results. This parameter allows users to exclude 
                    the effects of exogenous treatment.</ol>"

pct_info <- "<ol><li>Enter a value between 0 and 1 to set the threshold for the percentage of non-missing values in a biosynthetic pathway. Increasing this value will 
                result in fewer biosynthetic pathways being retained. This parameter enables users to regulate the substructure decomposition process, reducing artifacts 
                that may arise from excessive decomposition.
                <li>Usually, values between 0.3 and 0.7 are commonly used for this parameter.</ol>"

sig_path_info <- "The figure illustrates the top 5 significant representative pathways within the network, with red and blue indicating increase and decrease, respectively. 
                These pathways are represented by starting and ending lipids. A comprehensive summary of all significant pathways can be found in the accompanying table, 
                providing detailed information."

sig_reaction_info <- "The figure illustrates the top 5 significant reactions within the network, with red and blue indicating increase and decrease, respectively. These 
                    reactions are represented by substrate and product lipids. Red and blue text indicate the fold change of lipids. A comprehensive summary of all significant reactions can 
                    be found in the accompanying table, providing detailed information."

FALC_network_info <- "Top 5 significantly increased/decreased representative pathways and reactions are labeled in the network, with red and blue indicating increase and 
                    decrease, respectively. The line width and color depth indicate the importance of pathways, while the text size represents the importance of reactions. Furthermore, 
                    nodes in the figure are filled based on the log2 (fold change), and their sizes indicate -log10 (adjusted p-value). If a node exhibits significant changes in abundance, 
                    its border will be highlighted in purple."

LS_network_info <- "The network consists of all the significant pathways included in the top 5 increased and decreased representative pathways. Top 5 significantly increased/decreased 
                    representative pathways and reactions are also labeled in the network, with red and blue indicating increase and decrease, respectively. The line width and color 
                    depth indicate the importance of pathways, while the text size represents the importance of reactions. Furthermore, nodes in the figure are filled based on the 
                    log2 (fold change), and their sizes indicate -log10 (adjusted p-value). If a node exhibits significant changes in abundance, its border will be highlighted in purple."

substructure_info <- "In this section, we present the results of the differential expression analysis conducted on the substructure-transformed data."

FA_choices <- c("w3-22:6;0", "16:0;0", "18:0;0", "20:0;0", "22:0;0", "24:0;0", "26:0;0", "16:0;0", "w7-16:1;0", "16:0;0", "18:0;0", "w9-18:1;0", "w9-18:1;0", "w9-20:1;0", "w9-18:2;0", 
                "w9-20:2;0", "w9-20:1;0", "w9-22:1;0", "w9-24:1;0", "w6-18:2;0", "w6-18:3;0", "w6-20:3;0", "w6-22:4;0", "w6-20:4;0", "w6-24:4;0", "w6-26:4;0", "w6-24:4;0", "w6-24:5;0", 
                "w6-24:5;0", "w6-22:5;0", "w6-22:4;0", "w6-26:5;0", "w3-18:3;0", "w3-18:4;0", "w3-20:4;0", "w3-20:5;0", "w3-22:5;0", "w3-24:5;0", "w3-26:5;0", "w3-24:5;0", "w3-24:6;0", 
                "w3-24:6;0", "w3-26:6;0", "w3-22:6;0", "w3-22:5;0", "4:0;0", "6:0;0", "8:0;0", "10:0;0", "12:0;0", "14:0;0", "16:0;0", "18:0;0", "20:0;0", "22:0;0", "24:0;0", "26:0;0", 
                "28:0;0", "w7-16:1;0", "w7-18:1;0", "w10-16:1;0", "w9-18:1;0", "w9-18:2;0", "w9-20:1;0", "w9-20:2;0", "w9-20:2;0", "w9-20:3;0", "w9-22:1;0", "w9-24:1;0", "w9-26:1;0", 
                "w6-18:3;0", "w6-20:3;0", "w6-20:4;0", "w6-22:4;0", "w6-24:4;0", "w6-26:4;0", "w6-28:4;0", "w6-24:5;0", "w6-22:5;0", "w6-26:5;0", "w6-24:5;0", "w6-22:5;0", "w6-28:5;0", 
                "w3-18:4;0", "w3-20:4;0", "w3-20:5;0", "w3-22:5;0", "w3-24:5;0", "w3-26:5;0", "w3-28:5;0", "w3-24:6;0", "w3-22:6;0", "w3-26:6;0", "w3-28:6;0", "w3-24:6;0", "w3-22:6;0")

ui <- fluidPage(
    # image in browser tab
    tags$head(tags$link(rel = "icon", type = "image/png", sizes = "32x32", href = "ilipidome_logo.png")),

    navbarPage(
        # app title
        "iLipidome",
        # imageOutput("logo"),
        tabPanel("Lipid Substructure Analysis",
            tabsetPanel(
                tabPanel("Fatty Acid Analysis",
                    fluidRow(br(),),
                    fluidRow(
                        column(6,
                            radioButtons("FAData", "Data Source",
                                c(
                                    "Example dataset (Levental KR, et al. Nat Commun. 2020)" = "FAExample",
                                    "Upload your own data" = "FACustom"
                                ), 
                                width = "100%"
                            ) %>%
                            helper(
                                type = "inline",
                                title = "Data Source:",
                                size = "l",
                                content = data_info
                            ),
                            tabsetPanel(id = "FAFileIn", type = "hidden",
                                tabPanel("FAExample", 
                                    downloadButton("FAdownload", "Download Example Dataset"),
                                ),
                                tabPanel("FACustom",
                                    fileInput("FAfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c(".csv", ".tsv"), 
                                        width = "100%"
                                    )
                                ),
                            ),
                            br(),
                            uiOutput("FAresDownload"),
                            br(),
                            actionButton("FARun", "Run Analysis", padding = "8px")
                        ),
                        column(6, 
                            tabsetPanel(id = "FAparams", type = "hidden",
                                tabPanel("FAExample"
                                ),
                                tabPanel("FACustom",
                                    radioButtons("FAMethod", 
                                        label = "Method:",
                                        choices = c(
                                            "t-test" = "t.test",
                                            "Wilcoxon test" = "wilcox.test"
                                        ),
                                        selected = "t.test",
                                        inline = TRUE
                                    ),
                                    textInput("FActrl", 
                                        label = "Control Group:", 
                                        value = "1:7", 
                                        width = "100%"
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Assign Group Information:",
                                        size = "l",
                                        content = ctrl_info
                                    ),
                                    textInput("FAexp", 
                                        label = "Experimental Group:",
                                        value = "8:13", 
                                        width = "100%"
                                    ),
                                    selectInput("FAunmappedFA", 
                                        label = "Remove Low-expressed Fatty Acid Isomers:", 
                                        multiple = TRUE, 
                                        choices = FA_choices, 
                                        width = "100%"
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Remove low-expressed fatty acid isomers:",
                                        size = "l",
                                        content = unmapped_info
                                    ),
                                    selectInput("FAexolipid", 
                                        label = "Remove Exogenous Lipid Effect:",
                                        multiple = TRUE, 
                                        choices = FA_choices,
                                        width = "100%"
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
                                        inline = TRUE, 
                                        width = "100%"
                                    ),
                                ),
                            ),
                        ), 
                    ),
                    fluidRow(br(),),
                    fluidRow(
                        tabsetPanel( # tabPanels for visualizations
                            tabPanel("Lipid Expression Data",
                                # span(textOutput("FA_error"), style="color:red"),
                                verbatimTextOutput("FA_error"),
                                DT::dataTableOutput("FAInData")
                            ),
                            tabPanel("Substructure Analysis",
                                h4("Differential Expression Analysis"),
                                p(substructure_info),
                                plotOutput("FAvolcano", width = "700", height = "600"),
                                DT::dataTableOutput("FAsubres")
                            ),
                            tabPanel("Pathway Analysis",
                                h4("Significant pathways in Fatty Acid Network"),
                                p(sig_path_info),
                                span(textOutput("FA_nosig_path"), style="color:red"),
                                plotOutput("FAPathScorePlot", width = "700", height = "600"),
                                DT::dataTableOutput("FAPathScoreDT"),
                            ),
                            tabPanel("Reaction Analysis",
                                h4("Significant reactions in Fatty Acid Network"),
                                p(sig_reaction_info),
                                span(textOutput("FA_nosig_reaction"), style="color:red"),
                                plotOutput("FAReactionScorePlot", width = "700", height = "600"),
                                DT::dataTableOutput("FAReactionScoreDT"),
                            ),
                            tabPanel("Lipid Network",
                                h4("Fatty Acid Network"),
                                p(FALC_network_info),
                                uiOutput("FANodeEdge"),
                                visNetworkOutput("FANetworkGraph", height = "700px")
                            )
                        )
                    ),
                    fluidRow(br(),)
                ),
                tabPanel("Lipid Species Analysis",
                    fluidRow(br(), ),
                    fluidRow(
                        column(6,
                            radioButtons("LSData", "Data Source",
                                c(
                                    "Example dataset (Levental KR, et al. Nat Commun. 2020)" = "LSExample",
                                    "Upload your own data" = "LSCustom"
                                ), 
                                width = "100%"
                            ) %>%
                            helper(
                                type = "inline",
                                title = "Data Source:",
                                size = "l",
                                content = data_info
                            ),
                            tabsetPanel(id = "LSFileIn", type = "hidden",
                                tabPanel("LSExample",
                                    downloadButton("LSdownload", "Download Example Dataset"),
                                ),
                                tabPanel("LSCustom",
                                    fileInput("LSfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c(".csv", ".tsv"), 
                                        width = "100%"
                                    )
                                ),
                            ),
                            br(),
                            uiOutput("LSresDownload"),
                            br(),
                            actionButton("LSRun", "Run Analysis", padding = "8px")
                        ),
                        column(6,
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
                                        inline = TRUE,
                                        width = "100%"
                                    ),
                                    textInput("LSctrl", 
                                        label = "Control Group:", 
                                        value = "1:7",
                                        width = "100%"
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Assign Group Information:",
                                        size = "l",
                                        content = ctrl_info
                                    ),
                                    textInput("LSexp", 
                                        label = "Experimental Group:", 
                                        value = "8:13",
                                        width = "100%"
                                    ),
                                    numericInput("LSnonMissingPCT",
                                        label = "Proportion of Non-missing Values to Retain a Pathway", 
                                        value = 0.3,
                                        min = 0,
                                        max = 1,
                                        step = 0.1,
                                        width = "100%"
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Percentage of non-missing values to retain a pathway:",
                                        size = "l",
                                        content = pct_info
                                    ), 
                                    textInput("LSexolipid", 
                                        label = "Remove Exogenous Lipid Effect:", 
                                        value = NULL,
                                        width = "100%"
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
                                        inline = TRUE,
                                        width = "100%"
                                    ),
                                ),
                            ),
                        ),
                    ),
                    fluidRow(br(),),
                    fluidRow(
                        tabsetPanel( # tabPanels for visualizations
                            tabPanel("Lipid Expression Data",
                                verbatimTextOutput("LS_error"),
                                DT::dataTableOutput("LSInData")
                            ),
                            tabPanel("Substructure Analysis",
                                h4("Differential Expression Analysis"),
                                p(substructure_info),
                                plotOutput("LSvolcano", width = "700", height = "600"),
                                DT::dataTableOutput("LSsubres")
                            ),
                            tabPanel("Pathway Analysis",
                                h4("Significant pathways in Lipid Species Network"),
                                p(sig_path_info),
                                span(textOutput("LS_nosig_path"), style="color:red"),
                                plotOutput("LSPathScorePlot", width = "700", height = "600"),
                                DT::dataTableOutput("LSPathScoreDT"),
                            ),
                            tabPanel("Reaction Analysis",
                                h4("Significant reactions in Lipid Species Network"),
                                p(sig_reaction_info),
                                span(textOutput("LS_nosig_reaction"), style="color:red"),
                                plotOutput("LSReactionScorePlot", width = "700", height = "600"),
                                DT::dataTableOutput("LSReactionScoreDT"),
                            ),
                            tabPanel("Lipid Network",
                                h4("Lipid Species Network"),
                                p(LS_network_info),
                                uiOutput("LSNodeEdge"),
                                visNetworkOutput("LSNetworkGraph", height = "700px")
                            )
                        )
                    ),
                    fluidRow(br(),),
                ),
                tabPanel("Lipid Class Analysis",
                    fluidRow(br(),),
                    fluidRow(
                        column(6,
                            radioButtons("LCData", "Data Source",
                                c(
                                    "Example dataset (Levental KR, et al. Nat Commun. 2020)" = "LCExample",
                                    "Upload your own data" = "LCCustom"
                                ), 
                                width = "100%"
                            ) %>%
                            helper(
                                type = "inline",
                                title = "Data Source:",
                                size = "l",
                                content = data_info
                            ),
                            tabsetPanel(id = "LCFileIn", type = "hidden",
                                tabPanel("LCExample",
                                    downloadButton("LCdownload", "Download Example Dataset"),
                                ),
                                tabPanel("LCCustom",
                                    fileInput("LCfile", "Choose file",
                                        multiple = FALSE,
                                        accept = c(".csv", ".tsv"),
                                        width = "100%"
                                    )
                                ),
                            ),
                            br(),
                            uiOutput("LCresDownload"),
                            br(),
                            actionButton("LCRun", "Run Analysis", padding = "8px")
                        ),
                        column(6, 
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
                                        inline = TRUE,
                                        width = "100%"
                                    ),
                                    textInput("LCctrl", 
                                        label = "Control Group:", 
                                        value = "1:7",
                                        width = "100%"
                                    ) %>% helper(
                                        type = "inline",
                                        title = "Assign Group Information:",
                                        size = "l",
                                        content = ctrl_info
                                    ),
                                    textInput("LCexp", 
                                        label = "Experimental Group:", 
                                        value = "8:13",
                                        width = "100%"
                                    ),
                                    textInput("LCexolipid", 
                                        label = "Remove Exogenous Lipid Effect:", 
                                        value = NULL,
                                        width = "100%"
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
                                        inline = TRUE,
                                        width = "100%"
                                    ),
                                ),
                            ),
                        )
                    ),
                    fluidRow(br(),),
                    fluidRow(
                        tabsetPanel( # tabPanels for visualizations
                            tabPanel("Lipid Expression Data",
                                verbatimTextOutput("LC_error"),
                                DT::dataTableOutput("LCInData")
                            ),
                            tabPanel("Substructure Analysis",
                                h4("Differential Expression Analysis"),
                                p(substructure_info),
                                plotOutput("LCvolcano", width = "700", height = "600"),
                                DT::dataTableOutput("LCsubres")
                            ),
                            tabPanel("Pathway Analysis",
                                h4("Significant pathways in Lipid Class Network"),
                                p(sig_path_info),
                                span(textOutput("LC_nosig_path"), style = "color:red"),
                                plotOutput("LCPathScorePlot", width = "700", height = "600"),
                                DT::dataTableOutput("LCPathScoreDT"),
                            ),
                            tabPanel("Reaction Analysis",
                                h4("Significant reactions in Lipid Class Network"),
                                p(sig_reaction_info),
                                span(textOutput("LC_nosig_reaction"), style = "color:red"),
                                plotOutput("LCReactionScorePlot", width = "700", height = "600"),
                                DT::dataTableOutput("LCReactionScoreDT"),
                            ),
                            tabPanel("Lipid Network",
                                h4("Lipid Class Network"),
                                p(FALC_network_info),
                                uiOutput("LCNodeEdge"),
                                visNetworkOutput("LCNetworkGraph", height = "700px")
                            )
                        )
                    ),
                    fluidRow(br(),),
                )
            ),
        ),
        tabPanel("Tutorial",
            fluidRow(
                column(2),
                column(8,
                    withMathJax(), # load mathjax library to enable latex 

                    tags$div(HTML("<script type='text/x-mathjax-config'>
                                    MathJax.Hub.Config({
                                    tex2jax: {inlineMath: [['$','$']]}
                                    });
                                    </script>
                                    ")),

                    tags$div(
                        tags$img(src = "ilipidome_title_logo.png", width = "400px", height = "100px"),
                        style = "text-align: center;"
                    ),

                    tags$div(
                        h3("Contents"),
                        tags$ul(style = "list-style: none",
                            tags$li(a(href = "#about", "1 What is iLipidome?")),
                            tags$li(a(href = "#data", "2 Data source"),),
                            tags$li(tags$ul(style = "list-style: none",
                                tags$li(a(href = "#tryex", "2.1 Try our example"),),
                                tags$li(a(href = "#uploaddata", "2.2 Upload your data"),),
                                tags$li(tags$ul(style = "list-style: none",
                                    tags$li(a(href = "#prepdata", "2.2.1 How to prepare your dataset?"),),
                                    tags$li(a(href = "#groupinfo", "2.2.2 Assign group information"),),
                                    tags$li(a(href = "#lowFA", "2.2.3 Remove low-expressed fatty acid isomers (only in the fatty acid section)"),),
                                    tags$li(a(href = "#nonmissingpct", "2.2.4 Percentage of non-missing values to retain a pathway (only in the lipid species section):"),),
                                    tags$li(a(href = "#exolipid", "2.2.5 Remove exogenous lipid effect"),),
                                    tags$li(a(href = "#species", "2.2.6 Species"),),
                                ))
                            )),
                            tags$li(a(href = "#results", "3 Analysis Methods and Result Interpretation"),),
                            tags$li(tags$ul(style = "list-style: none",
                                tags$li(a(href = "#pathway", "3.1 Pathway analysis"),),
                                tags$li(tags$ul(style = "list-style: none",
                                    tags$li(a(href = "#pathscore", "3.1.1 Pathway scoring method"),),
                                    tags$li(a(href = "#idpathway", "3.1.2 Identify representative pathways"),)
                                )),
                                tags$li(a(href = "#reactanalysis", "3.2 Reaction analysis"),),
                                tags$li(tags$ul(style = "list-style: none",
                                    tags$li(a(href = "#reactscore", "3.2.1 Reaction scoring method"),)
                                ))
                            )),
                            tags$li(a(href = "#network", "4 Lipid network"),)
                        ),
                    ),

                    h1("About iLipidome"),
                    
                    a(name = "about"),
                    h3("1 What is iLipidome?"),
                    p("iLipidome is an innovative substructure-based method for analyzing lipidomics data using the lipid biosynthetic network, which takes into account the 
                        interdependence and interconnections of measured lipids. It now offers 'Lipid Substructure Analysis' functionality, allowing users to identify significant 
                        altered lipid pathways and link lipidomic changes to their genetic origins or reactions. In addition, iLipidome provides visualizations such as bar plots 
                        and networks, with accompanying statistical values in the table below the figures. The tool is divided into three main sections: “Fatty Acid Analysis”, 
                        “Lipid Species Analysis”, and “Lipid Class Analysis”, enabling comprehensive comparisons of lipid profiles at different levels. Our objective is to empower 
                        researchers with a deeper understanding of the intricate lipidomic alterations observed across various samples."),
                    
                    a(name = "data"),
                    h3("2 Data source"),
                    p("Users have the option to upload their own lipid datasets or utilize example datasets (Levental KR, et al. Nat Commun. 2020). The lipid expression data should 
                        be provided in CSV format. The lipid expression data consists of a feature column, which includes the identification of differential lipids such as PC_16:0;0_18:1;0 
                        or Cer_38:1;2, followed by multiple columns representing the samples and their corresponding values."),
                    
                    a(name = "tryex"),
                    h4("2.1 Try our example"),
                    p("The example dataset from the study 'Lipidomic and biophysical homeostasis of mammalian membranes counteracts dietary lipid perturbations to maintain cellular fitness' 
                        (1) includes the cellular lipidome data from 7 control samples and 6 samples treated with docosahexaenoic acid (DHA). The data was analyzed using MS-based shotgun 
                        lipidomics, revealing significant lipid and fatty acid reprogramming in the presence of DHA."),
                    
                    a(name = "uploaddata"),
                    h4("2.2 Upload your data"),
                    p("Users have the option to upload their own lipid dataset, which should be provided in CSV format. The analysis process offers various parameter options to customize the 
                        analysis, including 'Method', 'Group information', 'Remove low-expressed fatty acid isomers' (available only in the fatty acid section), 'Percentage of non-missing 
                        values to retain a pathway' (available only in the lipid species section), 'Remove exogenous lipid effect', and 'Species'. Users can adjust these parameters according 
                        to their study design to achieve more precise and tailored results."),
                    
                    a(name = "prepdata"),
                    h4("2.2.1 How to prepare your dataset?"),
                    p("To use iLipidome, you need to upload a lipid expression table (CSV format). The table should have lipids as rows and samples as columns. The first column should be named 
                        'feature' and contain the lipid names. The remaining columns should contain the sample names and the corresponding values for each lipid."),
                    tags$img(src = "dataset_format.png", align = "center", width = "100%"),
                    p("Currently, iLipidome supports two-group comparison for analyzing lipidomics data. Ensure that you have at least two samples in each group for calculating statistics. 
                        Depending on the data source, you may need to perform data processing or normalization methods such as missing value imputation or log transformation to improve the 
                        analysis results before uploading."),
                    p("Lipid names in the table can be represented in two formats:"),
                    p("1. When the exact identity of FAs is unknown, the lipids can be represented using the following format: [LipidClassAbbreviation]_[sum of FA chain length] : 
                        [sum of FA double bonds] ; [sum of FA oxygens] For example: 'PC_34:1;0' or 'TAG_52:1;0'"),
                    p("2. When the exact identity of FAs is known, the lipids can be represented using the following format: [LipidClassAbbreviation]_[FA1 chain length] : [FA1 double bonds] 
                        ; [FA1 oxygens]_[FA2 chain length] : [FA2 double bonds] ; [FA2 oxygens]... For example: 'PC_16:0;0_18:1;0' or 'TAG_16:0;0_18:0;0_18:1;0'"),
                    p("You can refer to the 'supported_lipid_class.csv' file for the supported lipid classes, their abbreviations, and the corresponding number of FAs. Note that when using 
                        the exact identity format of FAs, we will verify if the fatty acid numbers match those recorded in the 'supported_lipid_class.csv' file. If they do not match, the 
                        analysis will be interrupted. Also, lipid classes with the same number of FAs (e.g., PC, PE) in the same pathways (e.g., Glycerophospholipid) should have a consistent 
                        lipid naming format. For example, PC_36:0;0 and PE_34:0;0 or PC_18:0;0_18:0;0 and PE_16:0;0_18:0;0. Additionally, dihydrosphingolipids (dh-) specify sphingolipids with 
                        sphingoid bases of 18:0:2 instead of 18:1:2."),
                    # ADD DOWNLOAD BUTTON FOR SUPPORTED LIPID CLASS DOT CSV
                    downloadButton("Supported_lipid_class_download", "Download Supported Lipid Classes"),

                    a(name = "groupinfo"),
                    h4("2.2.2 Assign group information"),
                    p("Users can utilize the 'Control group' and 'Experimental group' parameters to assign group information for conducting a two-group comparison. Each group must include a 
                        minimum of two samples.  When entering the information in the 'Control group' or 'Experimental group' field, users should provide the column numbers separated by commas. 
                        It is important to consider that the first column of the dataset contains the names of the lipids, so the counting for grouping should start from the second column onwards. 
                        For example, if the lipid dataset has the lipid names in the first column, and columns 2 to 4 represent the control group while columns 5 to 7 represent the experimental 
                        group, you would enter '1,2,3' in the 'Control group' field and '4,5,6' in the 'Experimental group' field."),
                    tags$img(src = "ctrl_img.png", align = "center", width = "100%"),

                    a(name = "lowFA"),
                    h4("2.2.3 Remove low-expressed fatty acid isomers (only in the fatty acid section)"),
                    p("Due to limitations in mass spectrometry, precise double bond locations for fatty acids are often not available in lipidomics data. As a result, certain fatty acids may have 
                        multiple candidate mappings in the fatty acid network. However, some fatty acid isomers may be dominant, while others may be negligible. For exmaple, the major isomer of 
                        FA 20:4 is omega-6, not omega-3. Treating all isomers equally in the substructure calculation may not accurately reflect their true abundance. This parameter enables users 
                        to select low-expressed fatty acid isomers to exclude from decomposition into substructures within the fatty acid network, therefore improving the accuracy of calculations."),
                    tags$img(src = "FA_img.png", align = "center", width = "100%"),

                    a(name = "nonmissingpct"),
                    h4("2.2.4 Percentage of non-missing values to retain a pathway (only in the lipid species section)"),
                    p("During the substructure decomposition process in iLipidome, each lipid is decomposed based on its potential biosynthetic pathways. To avoid artifacts, a threshold is set to 
                        control the percentage of missing data within a pathway. You can enter a value between 0 and 1 to determine this threshold. If the proportion of missing substructures exceeds 
                        the threshold for a particular biosynthetic route, the target fatty acid or lipid species will not be decomposed into substructures. Increasing this value will result in fewer 
                        biosynthetic pathways being retained. Adjusting this parameter allows users to regulate the substructure decomposition process and minimize artifacts that may arise from 
                        unlimited decomposition."),

                    a(name = "exolipid"),
                    h4("2.2.5 Remove exogenous lipid effect"),
                    p("If an exogenous lipid treatment is involved in the study, it can significantly influence the results of substructure calculation based on biosynthetic pathways. To address this 
                        issue, iLipidome provides a parameter for users to exclude the effects of the exogenous treatment. In the fatty acid section, users have the option to select the exogenous 
                        fatty acids observed in the study to restrict substructure decomposition. In the lipid species and lipid class sections, users can enter the exogenous lipid names separated 
                        by commas. For example, you can enter 'PC_16:0;0_22:6;0, PC_18:0;0_22:6;0'. Please ensure that the lipid names or classes you enter correspond to those present in the uploaded 
                        dataset."),
                    tags$img(src = "exo_img.png", align = "center", width = "100%"),

                    a(name = "species"),
                    h4("2.2.6 Species"),
                    p("Currently, iLipidome supports the analysis of lipidomics data for three species: human, mouse, and rat."),

                    a(name = "results"),
                    h3("3 Analysis Methods and Result Interpretation"),
                    p("iLipidome offers a range of figures accompanied by corresponding tables to facilitate the interpretation of analysis results."),

                    a(name = "pathway"),
                    h4("3.1 Pathway analysis"),
                    p("In the 'Pathway analysis' section, the figure showcases the top 5 significant representative pathways within the network. Increased pathways are highlighted in red, while 
                        decreased pathways are shown in blue. A pathway is considered significant if its score exceeds $1.96$. The figure represents pathways using starting and ending lipids. 
                        Additionally, a comprehensive summary of all significant pathways can be found in the accompanying table. For a deeper understanding of how we calculate pathway scores, 
                        calibrate pathways, and select representative pathways, detailed information is available below."),
                    tags$img(src = "path_img.png", align = "center", width = "100%"),

                    a(name = "pathscore"),
                    h4("3.1.1 Pathway scoring method"), # USE MATHEMATICAL SYMBOLS FOR THIS SECTION, WITHMATHJAX
                    p("To analyze the pathways within the biosynthetic network, we firstly identify all possible pathways between any two nodes. We then examine the increased and decreased 
                        pathways within this set using a method adapted from previous studies (Nguyen A, et al. Curr Opin Biotechnol. 2017 and Trey Ideker, et al. Bioinformatics 2002)."),
                    p("In this method, we transform the p-value ($P_i$) of each node ($i$) in a pathway into a z-score, which is used to calculate the pathway score. Specifically, we convert each 
                        p-value to a z-score using the formula $Z_i = CDF^{-1}(1 - \\frac{1}{2} * P_i)$, where CDF is the cumulative distribution function. If the fold change between two experimental conditions 
                        is less than 1, we assign a negative sign to $Z_i$. For random data, $Z_i$ follows a standard normal distribution."),
                    p("To calculate the score of a pathway ($a$) with $n$ nodes, we sum all $Z_i$ values within the pathway and divide by $\\sqrt{n}$, yielding $Z_a = \\frac{1}{\\sqrt{n}} \\times \\Sigma Z_i$ for $i \\in a$. Since the variance of a 
                        sum is the sum of the variances for independent random variables, $Z_a$ also follows a standard normal distribution if the $Z_i$ values are independently drawn from a standard 
                        normal distribution. A high $Z_a$ indicates an active pathway, while a low $Z_a$ corresponds to a suppressed pathway. The length $n$ of the pathway is considered in the function 
                        to prevent over-weighting certain pathways."),
                    p("To calibrate $Z_a$ against the background distribution, we employ a Monte Carlo approach. We randomly select $n$ nodes and compute their scores $Z_a'$. These values are used to 
                        calculate the mean $\\mu a'$ and standard deviation $\\sigma a'$ for each $n$. $Z_a$ is then adjusted using these estimates to produce a final pathway score $S_a = (Z_a - Z_a') / \\sigma a'$. This 
                        calibration helps reduce noise and ensures that a randomized dataset has a mean $\\mu$ of $0$ and a standard deviation $\\sigma$ of $1$."),
                    p("To determine if a pathway is significantly increased or decreased, we use a threshold of $1.96$ as the critical value, based on the significance level $\\text{p-value} < 0.05$ from a 
                        two-tailed test. Thus, $S_a > 1.96$ represents a significantly increased pathway, while $S_a < -1.96$ indicates a significantly decreased pathway."),

                    a(name = "idpathway"),
                    h4("3.1.2 Identify representative pathways"),
                    p("To classify the pathways found in the biosynthetic network, we employ two methods based on fatty acids or lipid species and lipid classes. Initially, the pathways are 
                        separated into increased or decreased based on the signs of their scores. We then rank their importance using the absolute value of the pathway scores and perform a 
                        pathway similarity search from top to bottom. The pathway with the highest score is designated as the first representative pathway."),
                    p("In the FA analysis, we calculate the proportion of overlapping edges between each candidate pathway and the set of representative pathways. If the overlapping proportion 
                        of a pathway exceeds 50% with one representative pathway, it is assigned to that representative pathway. If none of the representative pathways have an overlapping 
                        proportion over 50%, the candidate pathway becomes a new representative pathway."),
                    p("On the other hand, since FA information is inherited in lipid species through biosynthetic connections, we can identify the dominant FA composition with the highest 
                        frequency (e.g., 16:0-18:1 or 32:0) within each lipid species pathway. Similar to the FA analysis, pathways with the same dominant FA composition are categorized into 
                        the same representative pathway. If a dominant FA appears for the first time during this process, the pathway is considered a new representative pathway."),

                    a(name = "reactanalysis"),
                    h4("3.2 Reaction analysis"),
                    p("In the 'Reaction analysis' section, the figure showcases the top 5 significant reactions within the network, where red and blue colors indicate an increase and decrease, 
                        respectively. A reaction is deemed significant if its p-value is below 0.05. These reactions are represented by substrate and product lipids, with red and blue text 
                        denoting the fold change of lipids. A comprehensive summary of all significant reactions is provided in the accompanying table. For a more detailed understanding of 
                        how we calculate reaction scores, please refer to the information below."),
                    tags$img(src = "react_img.png", align = "center", width = "100%"),

                    a(name = "reactscore"),  
                    h4("3.2.1 Reaction scoring method"),
                    p("To identify potential driver reactions or enzymes in the biosynthetic network, we have developed an algorithm to compute the perturbation score for each reaction (edge). 
                        The idea is to compare the ratio of product over reactant and assess its statistical significance across different conditions."),
                    p("Firstly, we calculate this ratio for each sample in a specific reaction. Then, the fold change is obtained by averaging the ratios from the control group and dividing it 
                        by the averaged ratios from the experimental group. Subsequently, we employ two-tailed Student's t-tests to compute the p-value of the ratios between the different 
                        conditions. Using these values, we calculate the final perturbation score for each reaction as: $-log_{10}(\\text{p-value}) \\times log_2(\\text{fold change})$. A score is considered significant 
                        when $\\text{p-value} < 0.05$."),
                    p("Genes/enzymes involved in the reactions are labeled mainly based on the LIPID MAPS and KEGG database. It's important to note that these databases and studies primarily 
                        provide information at the lipid class level, and the corresponding genes/enzymes in lipid molecular reactions highlighted by iLipidome do not take into account lipid 
                        species specificities."),

                    a(name = "network"),
                    h3("4 Lipid network"),
                    p("In the 'Lipid network' section, we constructed the Fatty Acid/Lipid Species/Lipid Class Network and highlighted the top 5 significantly increased/decreased representative 
                        pathways and reactions. In the network visualization, red and blue colors indicate increase and decrease, respectively. The line width and color depth reflect the importance 
                        of pathways, while the text size represents the significance of reactions. Additionally, the nodes in the figure are filled based on the $log_2(\\text{fold change})$ values, and their 
                        sizes represent $-log_{10}(\\text{adjusted p-value})$. If a node exhibits significant changes in abundance, its border will be highlighted in purple. It's important to note that for the 
                        Lipid Species Network, we only include the significant pathways that belong to the top 5 increased and decreased representative pathways to simplify the connections and 
                        enhance the clarity of the network visualization."),
                    tags$img(src = "net_img.png", align = "center", width = "100%"),
                ),
                column(2)
            )
        ),
        tabPanel("FAQ",
            fluidRow(
                column(2),
                column(8, 
                    tags$div(
                        tags$img(src = "ilipidome_title_logo.png", width = "400px", height = "100px"),
                        style = "text-align: center;"
                    ),
                    h3("Frequently Asked Questions"),
                    p("In this section, you will find answers to commonly asked questions. If there is an issue that you can't find the answer to here, 
                        please contact us using the info in the 'Contact Us' section of the website"),
                    br(),
                    p("Q1: How do I use the iLipidome website?"),
                    p("Answer: To get started with iLipidome, click on the 'Lipid Substructure Analysis' tab, select your dataset and parameters (if necessary), 
                        and click 'Run Analysis'. A status bar will appear in the bottom-right hand side of your screen to indcate the progress of the computation, 
                        and once it is done, you will be able to view the results in the five tabs directly below the 'Run Analysis' button."),
                    br(),
                    p("Q2: What types of files are accepted?"),
                    p("Answer: Currently, only .csv and .tsv files are accepted. "),
                    br(),
                    p("Q3: How should I format my data?"),
                    p("Answer: Guidelines on how to format your data can be found in section 2.2.1 of the 'Tutorial' tab of the website."),
                    br(),
                    p("Q4: I clicked the 'Run Analysis' button, but nothing is happening..."),
                    p(" "),
                    tags$div("Answer: If you clicked the 'Run Analysis' button and nothing is happening, please first make sure you have either selected the example dataset or have 
                        chosen your own dataset, and that there are no error messages pertaining to your own dataset. Otherwise, if the website grayed out, then an unexpected 
                        error occurred somewhere along the computational process. To continue using the website, please refresh the webpage. If this persists, then please open an issue", 
                        tags$a(href = "https://github.com/LewisLabUCSD/iLipidome-website/issues", "here.")
                    ),
                    br(),
                    p("Q5: The website suddenly turned gray, what happened?"),
                    tags$div("Answer: If the website suddenly turns gray, then an unexpected error occurred. To continue using the website, please refresh the webpage. 
                        If the issue persists, please open an issue ", 
                        tags$a(href = "https://github.com/LewisLabUCSD/iLipidome-website/issues", "here.")
                    ),
                    br(),
                    p("Q6: How can I cite iLipidome?"),
                    p("Answer: Please cite the following paper: Lin and Chiang et al. Nature Commun. 2023 XX(XX):XXXX. We would greatly appreciate it if you could also include our URL 
                        [http://www.bioinfomics.org/iLipidome/] in your text. If you find iLipidome useful, please help us spread the word."),
                    br(),
                    p("Q7: What types of lipid identifiers are supported?"),
                    p("Answer: iLipidome requires users to upload a lipidomics dataset in a specific format. For detailed instructions on how to prepare your dataset, please refer to the 
                        'How to prepare your dataset' section in the tutorial. It provides comprehensive information on the required format and steps to ensure compatibility with 
                        iLipidome's analysis pipeline."),
                    br(),
                    p("Q8: What species does iLipidome support?"),
                    tags$div("Answer: Currently, iLipidome website supports lipidomics data analysis for three species: human, mouse, and rat. We are continuously working to expand our support 
                        for additional species in the future. However, regarding the handling of species or reaction not currently covered in our repository, we provide advanced 
                        users with the flexibility to customize the reference fatty acid or lipid class pathways in the iLipidome package. For more detailed information, click ", 
                        tags$a(href = "https://github.com/LewisLabUCSD/iLipidome-website/issues", "here.")
                    ),
                    br(),
                    p("Q9: Is there an application programming interface (API) for iLipidome?"),
                    p("Answer: Please note that there is no API available for automated web-based submission in iLipidome. The server's backend is unable to handle the anticipated volume of API requests, 
                        as most iLipidome analyses require significant computational resources. We would like to reserve the iLipidome server for biologists.
                        However, we do provide a powerful R-based package called iLipidome for bioinformaticians [https://github.com/LewisLabUCSD/iLipidome-package]. iLipidome package allows 
                        bioinformaticians to perform unlimited offline analyses using their own hardware."),
                    br(),
                    p("Q10: Does iLipidome save my datasets and/or information?"),
                    p("Answer: Lipidome tracks standard analytics data for support purposes. Uploaded user files are automatically deleted after browser is closed. It is important to note that the 
                        maintainers of the website do not access or examine your data files, and we have no intention of adding a feature that saves your information and/or datasets on the website."),
                    br(),
                    p("Q11: Can I provide my own lipid dataset(s)?"),
                    p("Answer: Yes, any lipid dataset that follows the formatting guidelines detailed in the data formatting section of the tutorial can be used."),
                    br(),
                    p("Q12: How can I obtain publication quality images of results created by iLipidome?"),
                    tags$div("Answer: To obtain high quality images to use in your own publications, please make use of the R package ", 
                        tags$a(href = "https://github.com/LewisLabUCSD/iLipidome-website/issues", "here"),
                        " to create your own images using your own hardware."
                    ),
                    br(),
                ),
                column(2),
            ),
        ),
        tabPanel("Contact Us",
            fluidRow(
                column(2),
                column(8,
                    tags$div(
                        tags$img(src = "ilipidome_title_logo.png", width = "400px", height = "100px"),
                        style = "text-align: center;"
                    ),
                    h3("Contact Us"),
                    p("If you would like to reach out to us, please send an email to <waxlos987@gmail.com> and/or <austin.chiang@gmail.com> with 'iLipidome'
                    followed by the topic of the email in the subject line. "
                    ),
                    br(),
                    p("  "),
                    tags$div("If there is an issue with the website you would like to report, please create an issue ", 
                        tags$a(href = "https://github.com/LewisLabUCSD/iLipidome-website/issues", "here.")
                    ),
                    br(), br(),
                    h3("References:"),
                    tags$ul(
                        tags$li("Lin W.J.*, Chiang A.W.T.*‡, Zhou E.H., Liang C., Liu C.H., Ma W.L., Cheng W.C., Lewis N.E.‡ iLipidome: hidden biosynthetic interdependencies in the lipidome enhance statistical power and interpretability. submitted (2023). Under review at Nature Communications."),
                        tags$li("Bao, B. et al. Correcting for sparsity and interdependence in glycomics by accounting for glycan biosynthesis. Nature communications 12, 4988, doi:10.1038/s41467-021-25183-5 (2021)."),
                        tags$li("Lin, W.-J. et al. LipidSig: a web-based tool for lipidomic data analysis. Nucleic Acids Research 49, W336-W345, doi:10.1093/nar/gkab419 (2021).")
                    ),
                    br(),
                    p("Please cite Lin and Chiang et al. Nature Commun. 2023 XX(XX):XXXX within any publication that makes use of analyses inspired by iLipidome."),
                    br(),
                    h3("Acknowledgements"),
                    p("This work was supported by NIGMS (R35 GM119850) and generous support from the Novo Nordisk Foundation provided to the National Biologics Facility at the Technical University of Denmark (NNF20SA0066621). It was also supported by National Health Research Institute: grant # NHRI-EX111-11110BI, National Science and Technology Council: grant # MOST 111-2320-B-039-011, and China Medical University Hospital: grant # DMR-110-231, and # DMR-111-118. WJL was supported by the Elite Cultivation Scholarship Program Sponsored by China Medical University.")
                ),
                column(2)
            )
            
        ),
    ),
)

server <- function(input, output, session) {

    # output$logo <- renderImage({
    #     list(
    #         src = file.path("www/logo.png"),
    #         contentType = "image/png",
    #         width = 80,
    #         height = 30
    #     )
    # }, deleteFile = FALSE)

    # initialize shinyhelper package for help modals
    observe_helpers()

    FA_substructure_result <- NULL
    LS_substructure_result <- NULL
    LC_substructure_result <- NULL

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

        FA_substructure_result <<- NULL
        output$FA_error <- renderText("")
        output$FA_nosig_path <- renderText("")
        output$FA_nosig_reaction <- renderText("")

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
            if (is.null(input$FAfile$datapath)) {
                showModal(modalDialog(
                    title = "Please Upload a file",
                    "Please upload a dataset in the 'Choose file' section"
                ))
                return()
            }
            
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
                FA_substructure_result <<- FA_substructure_analysis(FA_exp_raw, method = 't.test',
                                                   ctrl = 1:7, exp = 8:13,
                                                   unmapped_FA = c('w9-18:2;0','w3-20:4;0'),
                                                   exo_lipid = 'w3-22:6;0', species = 'rat', 
                                                   progress = FA_progress, update_progress = update_progress)
            }
            if (input$FAData == "FACustom") {
                FA_substructure_result <<- FA_substructure_analysis(FA_exp_raw, method = input$FAMethod,
                                        ctrl = unlist(lapply(str_trim(unlist(strsplit(input$FActrl, ","))), function(x) eval(parse(text = x)))),
                                        exp = unlist(lapply(str_trim(unlist(strsplit(input$FAexp, ","))), function(x) eval(parse(text = x)))),
                                        unmapped_FA = input$FAunmappedFA,
                                        exo_lipid = input$FAexolipid, species = input$FAspecies,
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

            output$FAsubres <- DT::renderDataTable({
                FA_substructure_result[[6]] %>%
                    DT::datatable() %>%
                    DT::formatRound(which(sapply(FA_substructure_result[[6]], is.numeric)), digits = 3)
            })

            output$FAvolcano <- renderPlot(plot(FA_substructure_result[[9]]))

            output$FAresDownload <- renderUI(expr = if (!is.null(FA_substructure_result)) {
                downloadButton("FAresBTN", "Download Tables and Plots")
            } else {
                NULL
            })

            output$FANodeEdge <- renderUI(expr = if(!is.null(FA_substructure_result[[7]]) && !is.null(FA_substructure_result[[8]])) {
                downloadButton("FANodeEdgeBTN", "Download Network Node and Edge Tables")
            } else {
                NULL
            })
        }
        else { # display errors
            #figure out how to show FA_format error to user
            output$FA_error <- renderText({
                paste(unlist(strsplit(FA_format, split = "!")), sep = "\n")
            })
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

        LS_substructure_result <<- NULL
        output$LS_error <- renderText("")
        output$LS_nosig_path <- renderText("")
        output$LS_nosig_reaction <- renderText("")

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
            if (is.null(input$LSfile$datapath)) {
                showModal(modalDialog(
                    title = "Please Upload a file",
                    "Please upload a dataset in the 'Choose file' section"
                ))
                return()
            }

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
                LS_substructure_result <<- lipid_species_substructure_analysis(LS_exp_raw, method = 't.test',
                                                                         ctrl = 1:7, exp = 8:13,
                                                                         non_missing_pct = 0.3,
                                                                         exo_lipid = NULL, species = 'rat',
                                                                         progress = LS_progress, update_progress = update_progress)
            }
            if (input$LSData == "LSCustom") {
                if (input$LSexolipid == "") { LSexo <- NULL }
                else { LSexo <- str_trim(strsplit(input$LSexolipid, ",")[[1]]) }

                LS_substructure_result <<- lipid_species_substructure_analysis(LS_exp_raw, method = input$LSMethod,
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

            output$LSsubres <- DT::renderDataTable({
                LS_substructure_result[[6]] %>%
                    DT::datatable() %>%
                    DT::formatRound(which(sapply(LS_substructure_result[[6]], is.numeric)), digits = 3)
            })

            output$LSvolcano <- renderPlot(plot(LS_substructure_result[[9]]))

            output$LSresDownload <- renderUI(expr = if (!is.null(LS_substructure_result)) {
                downloadButton("LSresBTN", "Download Tables and Plots")
            } else {
                NULL
            })

            output$LSNodeEdge <- renderUI(expr = if(!is.null(LS_substructure_result[[7]]) && !is.null(LS_substructure_result[[8]])) {
                downloadButton("LSNodeEdgeBTN", "Download Network Node and Edge Tables")
            } else {
                NULL
            })
        }
        else { # display errors
            output$LS_error <- renderText({
                paste(unlist(strsplit(LS_format, split = "!")), sep = "\n")
            })
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

        LC_substructure_result <<- NULL
        output$LC_error <- renderText("")
        output$LC_nosig_path <- renderText("")
        output$LC_nosig_reaction <- renderText("")

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
            if (is.null(input$LCfile$datapath)) {
                showModal(modalDialog(
                    title = "Please Upload a file",
                    "Please upload a dataset in the 'Choose file' section"
                ))
                return()
            }
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
                LC_substructure_result <<- lipid_class_substructure_analysis(LC_exp_raw, method = 't.test',
                                                   ctrl = 1:7, exp = 8:13,
                                                   exo_lipid = NULL, species = 'rat',
                                                   progress = LC_progress, update_progress = update_progress)
            }
            if (input$LCData == "LCCustom") {
                if (input$LCexolipid == "") { LCexo <- NULL }
                else { LCexo <- str_trim(strsplit(input$LCexolipid, ",")[[1]]) }

                LC_substructure_result <<- lipid_class_substructure_analysis(LC_exp_raw, method = input$LCMethod,
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

            output$LCsubres <- DT::renderDataTable({
                LC_substructure_result[[6]] %>%
                    DT::datatable() %>%
                    DT::formatRound(which(sapply(LC_substructure_result[[6]], is.numeric)), digits = 3)
            })

            output$LCvolcano <- renderPlot(plot(LC_substructure_result[[9]]))

            output$LCresDownload <- renderUI(expr = if (!is.null(LC_substructure_result)) {
                downloadButton("LCresBTN", "Download Tables and Plots")
            } else {
                NULL
            })

            output$FANodeEdge <- renderUI(expr = if(!is.null(LC_substructure_result[[7]]) && !is.null(LC_substructure_result[[8]])) {
                downloadButton("LCNodeEdgeBTN", "Download Network Node and Edge Tables")
            } else {
                NULL
            })
        }
        else { # display errors
            output$LC_error <- renderText({
                paste(unlist(strsplit(LC_format, split = "!")), sep = "\n")
            })
        }
    })

    # HANDLE RESULT PLOT TABLE DOWNLOAD HERE
    output$FAresBTN <- downloadHandler(
        filename = function () { paste("FA_results", Sys.Date(), ".zip", sep = "") },
        content = function(file) {
            temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
            dir.create(temp_directory)

            pdf(file.path(temp_directory, "FA_plots.pdf"))
            print(FA_substructure_result[[2]])
            print(FA_substructure_result[[4]])
            print(FA_substructure_result[[9]])
            # FA_substructure_result[[5]] %>% View
            dev.off()

            write.csv(FA_substructure_result[[1]], file.path(temp_directory, "FA_sig_pathways.csv"))
            write.csv(FA_substructure_result[[3]], file.path(temp_directory, "FA_sig_reactions.csv"))
            write.csv(FA_substructure_result[[6]], file.path(temp_directory, "FA_diff_exp.csv"))

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        },
        contentType = "application/zip"
    )

    output$LSresBTN <- downloadHandler(
        filename = function () { paste("LS_results", Sys.Date(), ".zip", sep = "") },
        content = function(file) {
            temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
            dir.create(temp_directory)

            pdf(file.path(temp_directory, "LS_plots.pdf"))
            print(LS_substructure_result[[2]])
            print(LS_substructure_result[[4]])
            print(LS_substructure_result[[9]])
            # LS_substructure_result[[5]] %>% View
            dev.off()

            write.csv(LS_substructure_result[[1]], file.path(temp_directory, "LS_sig_pathways.csv"))
            write.csv(LS_substructure_result[[3]], file.path(temp_directory, "LS_sig_reactions.csv"))
            write.csv(LS_substructure_result[[6]], file.path(temp_directory, "LS_diff_exp.csv"))

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        },
        contentType = "application/zip"
    )

    output$LCresBTN <- downloadHandler(
        filename = function () { paste("LC_results", Sys.Date(), ".zip", sep = "") },
        content = function(file) {
            temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
            dir.create(temp_directory)

            pdf(file.path(temp_directory, "LC_plots.pdf"))
            print(LC_substructure_result[[2]])
            print(LC_substructure_result[[4]])
            print(LC_substructure_result[[9]])
            # LC_substructure_result[[5]]
            dev.off()

            write.csv(LC_substructure_result[[1]], file.path(temp_directory, "LC_sig_pathways.csv"))
            write.csv(LC_substructure_result[[3]], file.path(temp_directory, "LC_sig_reactions.csv"))
            write.csv(LC_substructure_result[[6]], file.path(temp_directory, "LC_diff_exp.csv"))

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        },
        contentType = "application/zip"
    )

    output$FANodeEdgeBTN <- downloadHandler(
        filename = function() { paste("FA_Node_Edge_Tables", Sys.Date(), ".zip", sep = "") },
        content = function(file) {
            temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
            dir.create(temp_directory)

            write.csv(FA_substructure_result[[7]], file.path(temp_directory, "FA_network_node.csv"))
            write.csv(FA_substructure_result[[8]], file.path(temp_directory, "FA_network_edge.csv"))

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        },
        contentType = "application/zip"
    )

    output$LSNodeEdgeBTN <- downloadHandler(
        filename = function() { paste("LS_Node_Edge_Tables", Sys.Date(), ".zip", sep = "") },
        content = function(file) {
            temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
            dir.create(temp_directory)

            write.csv(LS_substructure_result[[7]], file.path(temp_directory, "LS_network_node.csv"))
            write.csv(LS_substructure_result[[8]], file.path(temp_directory, "LS_network_edge.csv"))

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        },
        contentType = "application/zip"
    )

    output$LCNodeEdgeBTN <- downloadHandler(
        filename = function() { paste("LC_Node_Edge_Tables", Sys.Date(), ".zip", sep = "") },
        content = function(file) {
            temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
            dir.create(temp_directory)

            write.csv(LC_substructure_result[[7]], file.path(temp_directory, "LC_network_node.csv"))
            write.csv(LC_substructure_result[[8]], file.path(temp_directory, "LC_network_edge.csv"))

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        },
        contentType = "application/zip"
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

    output$Supported_lipid_class_download <- downloadHandler(
        filename = "supported_lipid_class.csv",
        content = function(fileDownload) {
            file.copy("supported_lipid_class.csv", fileDownload)
        }
    )
}

shinyApp(ui, server)
