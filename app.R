library(shiny)
library(bslib)
library(shinyFiles)

source("helpers.R")

ui <- page_sidebar(
  title = "Flow Cytometry Analysis",
  sidebar = sidebar(
    width = 280,
    shinyDirButton(
      id    = "folder",
      label = "Select Folder",
      title = "Choose a folder containing .fcs files"
    ),
    hr(),
    p("Channels will appear here", class = "text-muted fst-italic"),
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
    roots   = c(Home = path.expand("~"), Root = "/"),
    session = session
  )
}

shinyApp(ui, server)
