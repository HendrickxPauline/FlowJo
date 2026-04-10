library(shiny)
library(bslib)
library(shinyFiles)
library(rhandsontable)

source("helpers.R")

# Auto-detect all available drives/volumes (home, external drives, network, etc.)
ROOTS <- getVolumes()()

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
    uiOutput("plate_layout_ui"),
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
    req(folder_path())
    fs <- flow_set()
    if (is.null(fs)) {
      p("No .fcs files found in this folder.",
        class = "text-warning small mt-1")
    } else {
      p(paste0("\u2713 ", length(fs), " .fcs file(s) loaded."),
        class = "text-success small mt-1")
    }
  })

  # ── Plate layout table ────────────────────────────────────────────────────

  # Wrap the table in a uiOutput so it is hidden until files are loaded.
  output$plate_layout_ui <- renderUI({
    req(flow_set())
    tagList(
      hr(),
      tags$label("Sample Labels", class = "form-label fw-semibold"),
      rHandsontableOutput("plate_layout", width = "100%")
    )
  })

  # Build the base data frame from filenames whenever a new folder is loaded.
  # Re-renders the table (and resets any edits) only when flow_set() changes.
  output$plate_layout <- renderRHandsontable({
    req(flow_set())
    df <- data.frame(
      Filename      = sampleNames(flow_set()),
      `Sample Name` = rep("", length(flow_set())),
      stringsAsFactors = FALSE,
      check.names   = FALSE
    )
    rhandsontable(df, stretchH = "all", rowHeaders = NULL) |>
      hot_col("Filename",    readOnly = TRUE)  |>
      hot_col("Sample Name", readOnly = FALSE)
  })

  # Reactive data frame of current table contents (used by plots later).
  # Falls back to the filename-only frame before the user makes any edits.
  sample_map <- reactive({
    if (!is.null(input$plate_layout)) {
      hot_to_r(input$plate_layout)
    } else {
      req(flow_set())
      data.frame(
        Filename      = sampleNames(flow_set()),
        `Sample Name` = rep("", length(flow_set())),
        stringsAsFactors = FALSE,
        check.names   = FALSE
      )
    }
  })

  # Print to console whenever the user edits the table.
  observeEvent(input$plate_layout, {
    df <- hot_to_r(input$plate_layout)
    cat("\n--- Sample map updated ---\n")
    print(df)
  })

  # ── Channel checkboxes ────────────────────────────────────────────────────

  output$channel_ui <- renderUI({
    fs <- flow_set()
    if (is.null(fs)) {
      return(p("Channels will appear here", class = "text-muted fst-italic"))
    }
    checkboxGroupInput(
      inputId  = "channels",
      label    = "Channels",
      choices  = colnames(fs),
      selected = colnames(fs)
    )
  })
}

shinyApp(ui, server)
