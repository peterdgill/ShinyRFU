
ShinyRFU<-function(){
library(shiny)
library(euroformix)
library(shinybusy)
# Define UI ----
ui <- fluidPage(
#Busy indicator---
add_busy_spinner(spin = "cube-grid",color = "#1D6FB3",position=c("full-page"),height = "120px", width = "120px"),

##App title---
  titlePanel("Average RFU input"),
  
  fluidRow(
    
    column(3, wellPanel(
           h3("Input Data (defaults shown)"),
  numericInput("Fstvalue",label= ("Fst:"),value=0.01),

  numericInput("ATvalue",label= ("Analytical Threshold (AT):"),value=50),

  numericInput("PCvalue",label= ("Drop-in probability (pC):"),value=0.05),

  numericInput("lambdavalue",label= ("Lambda (Drop-in parameter):"),value=0.01),

selectInput("kit", h3("Select kit"), 
                       choices = list("ESX16", "ESX17",
                                      "ESX17Fast", "ESI17Fast",
                                      "Fusion", "Fusion 6C",
						  "SGMPlus","Identifiler",
                                      "NGM","NGMSElect",
                                      "NGMDetect", "GlobalFiler",
                                      "PowerPlex16", "PowerPlex21",
                                      "24Plex", "ESSPlex",
                                      "ESSplexPlus", "ESSplexSEPlus",
                                      "ESSplexSEQS", "ForenSeq", "VFP"),selected = "Fusion 6C")                              
),

),

#####file input
fluidRow(
column(4,
h3("Input data files"),
# Input: Select a file ----
      fileInput("evidence", "Choose Evidence File",
                multiple = FALSE,
                accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
# Input: Select a file ----
      fileInput("reference", "Choose Reference File",
                multiple = FALSE,
                accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),

# Input: Select a file ----
      fileInput("hypotheses", "Choose Hypotheses File",
                multiple = FALSE,
                accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),


#Choose database file

# Input: Select a file ----
      fileInput("file1", "Choose population frequency database.CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),

# set download button
actionButton("goButton", "Press to Calculate", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
),

##Choose a kit
fluidRow(
    
column(4, offset=4, 
               downloadButton("downloadData", "Download Data",)),
 )))


# Define server logic ----
server <- function(input, output) {
observeEvent(input$goButton,{
#code here
kit<-input$kit
AT<-input$ATvalue
fst<-input$Fstvalue
pC<-input$PCvalue
lambda<-input$lambdavalue

popFreq_file=input$file1
evidData_file=input$evidence
refData_file=input$reference
hypData_file=input$hypotheses

#READ DATA:
  evidData = euroformix::tableReader(input$evidence$datapath) 
  refData = euroformix::tableReader(input$reference$datapath)
  hypData = euroformix::tableReader(input$hypotheses$datapath)
  popFreq <- freqImport(input$file1$datapath)[[1]] #import population freqs.

###Error catching##################
errorlogList <- list()

read = function(Pfile) {
  #stop("Something wrong")
  
  data = NA
  tryCatch({
    #data = read.table(Pfile)
  #evidData = euroformix::tableReader(input$evidence$datapath) 
data = euroformix::tableReader(Pfile) 

  },error = function(errorMsg) {
    errorlogList[[length(errorlogList)+1]] <<- errorMsg
    
  })
}
#Pfile="evids.csv"
#read(Pfile)
#errorlogList
####################################
##Validate function

validate(
need(evidData,'Error in file evidData')
)




###end error catch#################


##Run the function
outfile = "resultfile.csv"
resTable=NA
resTable<<-AveRFU3(popFreq,evidData,refData,hypData,kit,AT,fst,pC,lambda)

})###end of action button

output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$downloadData, input$filetype, sep = ";")
    },
    content = function(file) {
      write.table(resTable, file, sep=";", row.names=FALSE)
    }
) 
}

#Run the app
shinyApp(ui = ui, server = server)
####################################
}
