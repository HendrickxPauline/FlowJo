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
  title    = "Flow Cytometry Analysis",
  fillable = TRUE,
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
    # View switcher — always visible so the user can pick before loading files
    radioButtons(
      inputId  = "active_view",
      label    = "Show",
      choices  = c("Histogram", "Dot Plot", "Results Table"),
      selected = "Histogram"
    ),
    # Context-sensitive controls for the active view
    uiOutput("view_controls_ui")
  ),
  # Single card that fills the main panel; full_screen = TRUE adds an expand button
  card(
    fill        = TRUE,
    full_screen = TRUE,
    card_header(textOutput("card_title", inline = TRUE)),
    card_body(
      fill    = TRUE,
      padding = "1rem",
      # ── Histogram ──────────────────────────────────────────────────────
      conditionalPanel(
        "input.active_view === 'Histogram'",
        plotlyOutput("histogram", width = "100%", height = "520px"),
        uiOutput("threshold_slider_ui")
      ),
      # ── Dot Plot (placeholder) ──────────────────────────────────────────
      conditionalPanel(
        "input.active_view === 'Dot Plot'",
        p("Dot plot will appear here.", class = "text-muted fst-italic p-2")
      ),
      # ── Results Table (placeholder) ────────────────────────────────────
      conditionalPanel(
        "input.active_view === 'Results Table'",
        p("Results table will appear here.", class = "text-muted fst-italic p-2")
      )
    )
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
      hot_col("Filename",    readOnly = TRUE) |>
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

  # ── View selector ─────────────────────────────────────────────────────────

  # Card header mirrors the active view name.
  output$card_title <- renderText({
    req(input$active_view)
    input$active_view
  })

  # Sidebar controls that change depending on the active view.
  # Shown only after files are loaded.
  output$view_controls_ui <- renderUI({
    req(flow_set())
    channels <- colnames(flow_set())

    if (input$active_view == "Histogram") {
      tagList(
        selectInput(
          inputId  = "hist_channel",
          label    = "Channel",
          choices  = channels,
          selected = channels[1]
        )
      )
    } else if (input$active_view == "Dot Plot") {
      p("Dot plot controls will appear here.",
        class = "text-muted fst-italic small")
    } else {
      NULL
    }
  })

  # ── Histogram ─────────────────────────────────────────────────────────────

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

    if (!is.null(input$threshold)) {
      p <- p + geom_vline(xintercept = input$threshold,
                          colour = "red", linetype = "dashed", linewidth = 0.8)
    }

    ggplotly(p) |> layout(showlegend = FALSE)
  })
}

shinyApp(ui, server)
