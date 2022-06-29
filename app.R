####################################
#   Artificial Intelligence for    #
#SARS-COV-2 Variant Prediction(AIS)#
#   Lue Ping Zhao & Michael Zhao   #
#         June 15, 2022            #
####################################
#source("ais.R")
#results=AIS(preddata, variant_prop, variant_corehap)

####################################
# User interface                   #
####################################
ui <- fluidPage(theme = shinytheme("united"),shinyjs::useShinyjs(),
								# Page header
								headerPanel('A.I. for SARS-COV-2 Variant Predictions'),
								tags$style(
									type = 'text/css',
									'table.dataTable td {white-space: nowrap;}'
								),
								
								
								# Input values
								sidebarPanel(
									HTML("<h3>Input Data</h3>"),
									
									wellPanel(
										tags$label("Step 1 - Enter your data",style="float: none; width: 100%;"),
										actionLink("addlink", "Insert example data"),
										tags$textarea(id="predata", rows=10, cols=100, style="float: none; width:100%;", "EPI_ISL_11657253,EPI_ISL_11657258,EPI_ISL_4579449,..."),
										uiOutput('file1_ui'),
										tags$label("Step 2 - Submit your job",style="float: none; width: 100%;"),
										actionButton("submitbutton", "Predict", class = "btn btn-primary"),
										actionButton("clear", "Clear", class = "btn btn-primary"),
									),
								),
								# Output panel
								mainPanel(
									verbatimTextOutput('contents'),
									navbarPage(title=("RESULTS"),
														 tabPanel("Summary",value="",
														 				 wellPanel(
														 				 	div(id='table-container', style='width: 100% !important; overflow-x: auto; overflow-y: auto; height: 70vh;',
														 				 			DT::dataTableOutput('summary'),  # Prediction results table
														 				 	),
														 				 	downloadButton("downloadData", "Download")
														 				 )
														 ),
														 tabPanel("Probability", value="",
														 				 wellPanel(
														 				 	div(id='table-container', style='width: 100% !important; overflow-x: auto; overflow-y: auto; height: 70vh;',
														 				 			DT::dataTableOutput('probability'),  # Prediction results table
														 				 	),
														 				 	downloadButton("downloadprob", "Download")
														 				 )
														 				 ),
														 tabPanel("Kappa", value="",
														 				 wellPanel(
														 				 	div(id='table-container', style='width: 100% !important; overflow-x: auto; overflow-y: auto; height: 70vh;',
														 				 			DT::dataTableOutput('kappa'),  # Kappa analysis with AIS prediction
														 				 	),
														 				 	downloadButton("downloadkappa", "Download")
														 				 )
														 ),
														 tabPanel("Kappapp", value="",
														 				 wellPanel(
														 				 	div(id='table-container', style='width: 100% !important; overflow-x: auto; overflow-y: auto; height: 70vh;',
														 				 			DT::dataTableOutput('kappapp'),  # Kappa analysis with post-prediction modification
														 				 	),
														 				 	downloadButton("downloadkappapp", "Download")
														 				 )
														 ),
														 
														 tabPanel("Your Data", value="",
														 				 wellPanel(
														 				 	div(id='table-container', style='width: 100% !important; overflow-x: auto; overflow-y: auto; height: 70vh;',
														 				 			DT::dataTableOutput('inputdata'),  # Prediction results table
														 				 	),
														 				 	downloadButton("downloadinput", "Download")
														 				 )
														 ),
									)	
								)
)

####################################
# Server                           #
####################################
server <- function(input, output, session) {
	# Default example data
	observe({
		FASTADATA <- ''
		fastaexample <- 'EPI_ISL_11657253,EPI_ISL_11657258,EPI_ISL_4579449,EPI_ISL_11657167,EPI_ISL_4579450,EPI_ISL_11524493,EPI_ISL_11657362,EPI_ISL_11657367,EPI_ISL_11657241,EPI_ISL_11657365'
		xdata=""
		m_parameters=c(prob_threshold=0.99, min_prob=0.1, min_AA=2)
		
		if(input$addlink>0) {
			isolate({
				FASTADATA <- fastaexample
				updateTextInput(session, inputId = "predata", value = FASTADATA)
			})
		}
	})
	
	# Input Data
	datasetInput <- reactive({  
		inFile <- input$file1 
		inTextbox <- input$predata
		shinyjs::logjs("hello12");
		
		if (inTextbox=="" & is.null(inFile)) {
			xdata=""
			return(data.frame(WARNING="Please insert/upload prediction data in GISAID format"))
		} else {
			if (is.null(inFile)) {
				xdata <- inTextbox
				xdata=gsub(" ","", xdata)
				xdata=strsplit2(xdata,",")
				xdata=xdata[substr(xdata,1,8)=="EPI_ISL_"]
				xdata=unique(xdata)
				xdata=intersect(xdata, rownames(metadata))
				xN=length(xdata)
				# Status/Output Text Box
				output$contents <- renderPrint({
					if (input$submitbutton>0 & xN>xN_limit) return(paste("Exceed the limit of", xN_limit, "viruses", sep=" "))
				})
				if (is.null(xdata)) return("invalid GISAID ID or no selected data") else {
					if (xN<=xN_limit) {
						xdata=intersect(xdata, rownames(metadata))
						results<-AIS(metadata[xdata,], variant_prop, variant_corehap,m_parameters)
						return(results)						
					}
				}
			} else {
				xdata <- read.csv(inFile$datapath, row.names = 1)
				xN=NROW(xdata)
				# Status/Output Text Box
				output$contents <- renderPrint({
					if (input$submitbutton>0 & xN>xN_limit) return(paste("Exceed the limit of", xN_limit, "viruses", sep=" "))
					if (input$submitbutton>0 & sum(colnames(xdata)=="AA.Substitutions")==0) return("Where is AA.Substitution column?")
				})
				if (xN<=xN_limit) {
					results<-AIS(xdata, variant_prop, variant_corehap,m_parameters)
					return(results)						
				}
			}
		}
	})
	
	observeEvent(input$clear, {
		updateTextInput(session, inputId='predata', value = "")
	})
	
	output$file1_ui <- renderUI({
		input$clear ## Create a dependency with the reset button
		fileInput('file1', 'or upload file',accept=c('text','.csv','.txt'))
	})
	
	# Prediction results
	output$summary <- DT::renderDataTable({
		if (input$submitbutton>0) { 
			isolate(datasetInput()[["summary"]])
		} 
	})
	
	# Prediction probabilities
	output$probability <- DT::renderDataTable({
		if (input$submitbutton>0) { 
			isolate(datasetInput()[["probability"]]) 
		} 
	})
	
	# Kappa
	output$kappa <- DT::renderDataTable({
		if (input$submitbutton>0) { 
			isolate(datasetInput()[["kappa"]]) 
		} 
	})
	
	# Kappa
	output$kappapp <- DT::renderDataTable({
		if (input$submitbutton>0) { 
			isolate(datasetInput()[["kappa_pp"]]) 
		} 
	})
	
	# Input data
	output$inputdata <- DT::renderDataTable({
		if (input$submitbutton>0) { 
			isolate(datasetInput()[["inputdata"]]) 
		} 
	})
	
	# Downloadable csv of key result
	output$downloadData <- downloadHandler(
		filename = function() {
			paste(input$dataset, "****.csv", sep = "")
		},
		content = function(file) {
			write.csv(datasetInput()[["summary"]], file, row.names = T)
		}
	)
	
	# Downloadable csv of posterior probabilities
	output$downloadprob <- downloadHandler(
		filename = function() {
			paste(input$dataset, "****.csv", sep = "")
		},
		content = function(file) {
			write.csv(datasetInput()[["probability"]], file, row.names = T)
		}
	)

		# Downloadable csv of kappa frequencies
	output$downloadkappa <- downloadHandler(
		filename = function() {
			paste(input$dataset, "****.csv", sep = "")
		},
		content = function(file) {
			write.csv(datasetInput()[["kappa"]], file, row.names = T)
		}
	)
	
	# Downloadable csv of kappa frequencies for post-prediction modifications
	output$downloadkappapp <- downloadHandler(
		filename = function() {
			paste(input$dataset, "****.csv", sep = "")
		},
		content = function(file) {
			write.csv(datasetInput()[["kappa_pp"]], file, row.names = T)
		}
	)
	# Downloadable csv of inputdata
	output$downloadinput <- downloadHandler(
		filename = function() {
			paste(input$dataset, "****.csv", sep = "")
		},
		content = function(file) {
			write.csv(datasetInput()[["inputdata"]], file, row.names = T)
		}
	)
	
}

####################################
# Create the shiny app             #
####################################
shinyApp(ui = ui, server = server)
