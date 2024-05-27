# * Author:    Yang.Xu 
# * Created:   11:40 PM Monday, 22 April 2024
# * Copyright: AS IS
#' @importFrom stats cor dist optim predict na.omit var 
#' @importFrom utils write.csv  
#' @importFrom pls plsr RMSEP MSEP mvrValstats
#' @importFrom BGLR BGLR
#' @importFrom glmnet cv.glmnet glmnet 
#' @importFrom xgboost xgboost 
#' @importFrom lightgbm lightgbm
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach %dopar% foreach
#' @importFrom graphics barplot text
#' @importFrom predhy cv convertgen predhy.predict predhy.predict_NCII infergen kin mixed crodesign
#' @importFrom data.table fread
#' @importFrom shiny fluidPage navbarPage tabPanel req navlistPanel fluidRow fileInput helpText column radioButtons conditionalPanel selectInput numericInput textOutput plotOutput downloadHandler downloadLink sliderInput checkboxInput reactive renderText renderPlot shinyApp 
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom htmltools h2 h3 h4 code 
NULL
