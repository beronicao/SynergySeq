#App7



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### 1. LOADING ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
library(ggplot2)
library(shiny)
library(plotly)
library(DT)
library(visNetwork)
library(reshape2)
library(shinythemes)
SM_MOA <- read.csv(file="data/L1000_SM_MOA.csv", header=TRUE)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### 2. UI ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
ui <- navbarPage(  theme = shinytheme("flatly"),responsive = TRUE,
                   title = div(
                              img(src = "Logos_3.png",height="30",align="left",position="fixed",hspace="10"),
                              "SynergySeq"
                              ),
                  
                   ### UI-1 Introduction ####
                   source("tabs/introduction_UI.R",  local = TRUE)$value,
                  
                   ### UI-2 Tools1 ####
                   tabPanel("Tools", tabsetPanel( 
                            source("tabs/Tools1_synergy_plot_with_reference_UI.R",  local = TRUE)$value
                            )),
                  
                   ### UI-3 About ####
                   source("tabs/about.R",  local = TRUE)$value
              )
 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### 3. SERVER ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
server <- function(input, output)
                {
                  ### Tools1  ####
                  source("tabs/Tools1_synergy_plot_with_reference_SERVER.R",  local = TRUE)$value
  
  
  
                }






shinyApp(ui, server)