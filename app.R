library(shiny)
library(bslib)
library(shinyFiles)

source("helpers.R")

# Roots exposed to shinyFiles — shared between UI and server
ROOTS <- c(Home = path.expand("~"), Root = "/")

ui <- page_sidebar(
  title = "Flow Cytometry Analysis",
  sidebar = sidebar(
    width = 280,
    shinyDirButton(
      id    = "folder",
      label = "Select Folder",
      title = "Choose a folder containing .fcs files"
    ),
    uiOutput("file_status"),
    hr(),
    uiOutput("channel_ui"),
    hr(),
    p("Dot plot controls will appear here", class = "text-muted fst-italic")
  ),
  card(
    full_screen = TRUE,
    card_header("Dot Plot"),
    card_body(min_height = 300)
  ),
  card(
    full_screen = TRUE,
    card_header("Histogram"),
    card_body(min_height = 250)
  ),
  card(
    full_screen = TRUE,
    card_header("Results Table"),
    card_body(min_height = 200)
  )
)

server <- function(input, output, session) {

  shinyDirChoose(
    input   = input,
    id      = "folder",
    roots   = ROOTS,
    session = session
  )

  # Resolve the chosen directory to an absolute path string.
  # parseDirPath returns character(0) until the user picks a folder.
  folder_path <- reactive({
    req(input$folder)
    if (is.integer(input$folder)) return(NULL)   # not yet chosen
    path <- parseDirPath(ROOTS, input$folder)
    if (length(path) == 0) return(NULL)
    path
  })

  # Load all .fcs files in the chosen folder into a flowSet.
  flow_set <- reactive({
    req(folder_path())
    load_fcs_folder(folder_path())
  })

  # Small status line showing how many files were loaded.
  output$file_status <- renderUI({
    req(folder_path())          # show nothing until a folder is chosen
    fs <- flow_set()
    if (is.null(fs)) {
      p("No .fcs files found in this folder.",
        class = "text-warning small mt-1")
    } else {
      p(paste0("\u2713 ", length(fs), " .fcs file(s) loaded."),
        class = "text-success small mt-1")
    }
  })

  # Channel checkboxes — replaces the static placeholder.
  output$channel_ui <- renderUI({
    fs <- flow_set()
    if (is.null(fs)) {
      return(p("Channels will appear here", class = "text-muted fst-italic"))
    }
    checkboxGroupInput(
      inputId  = "channels",
      label    = "Channels",
      choices  = colnames(fs),
      selected = colnames(fs)   # all checked by default
    )
  })
}

shinyApp(ui, server)
