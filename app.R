library(shiny)
library(bslib)
library(shinyFiles)
library(rhandsontable)
library(flowCore)
library(ggplot2)
library(plotly)
library(scales)

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
    uiOutput("hist_channel_ui"),
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
    card_body(
      plotlyOutput("histogram"),
      uiOutput("threshold_slider_ui")
    )
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
  folder_path <- reactive({
    req(input$folder)
    if (is.integer(input$folder)) return(NULL)
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

  output$plate_layout_ui <- renderUI({
    req(flow_set())
    tagList(
      hr(),
      tags$label("Sample Labels", class = "form-label fw-semibold"),
      rHandsontableOutput("plate_layout", width = "100%")
    )
  })

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

  # Reactive data frame of current table contents — used by plots later.
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

  observeEvent(input$plate_layout, {
    cat("\n--- Sample map updated ---\n")
    print(hot_to_r(input$plate_layout))
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

  # ── Histogram ─────────────────────────────────────────────────────────────

  # Channel selector — appears once files are loaded.
  output$hist_channel_ui <- renderUI({
    req(flow_set())
    selectInput(
      inputId  = "hist_channel",
      label    = "Histogram Channel",
      choices  = colnames(flow_set()),
      selected = colnames(flow_set())[1]
    )
  })

  # Combined data for the selected channel across all samples.
  channel_data <- reactive({
    req(flow_set(), input$hist_channel)
    extract_channel_data(flow_set(), input$hist_channel)
  })

  # Threshold stored as a reactiveVal so other panels can read it later.
  threshold_val <- reactiveVal(NULL)
  observeEvent(input$threshold, threshold_val(input$threshold))

  # Slider range resets automatically whenever the selected channel changes.
  output$threshold_slider_ui <- renderUI({
    req(channel_data())
    vals <- channel_data()$value
    lo   <- floor(min(vals,   na.rm = TRUE))
    hi   <- ceiling(max(vals, na.rm = TRUE))
    sliderInput(
      inputId = "threshold",
      label   = "Threshold",
      min     = lo,
      max     = hi,
      value   = round((lo + hi) / 2),
      width   = "100%",
      step    = 1
    )
  })

  output$histogram <- renderPlotly({
    req(channel_data())
    df <- channel_data()

    p <- ggplot(df, aes(x = value)) +
      geom_histogram(bins = 100, fill = "#4C72B0", colour = NA, alpha = 0.85) +
      scale_x_continuous(
        trans  = pseudo_log_trans(sigma = 1),
        labels = label_number(scale_cut = cut_short_scale())
      ) +
      labs(x = input$hist_channel, y = "Count") +
      theme_bw(base_size = 13)

    # Add threshold line once the slider has rendered.
    if (!is.null(input$threshold)) {
      p <- p + geom_vline(xintercept = input$threshold,
                          colour = "red", linetype = "dashed", linewidth = 0.8)
    }

    ggplotly(p) |> layout(showlegend = FALSE)
  })
}

shinyApp(ui, server)
