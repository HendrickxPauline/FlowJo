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

# Pseudo-log transform helpers (sigma = 1, matching the histogram x-axis).
# plotly stores data in transformed space, so shapes need transformed coordinates.
plog_fwd <- function(x) asinh(x)          # original  → plotly axis
plog_inv <- function(y) sinh(y)           # plotly axis → original

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
    radioButtons(
      inputId  = "active_view",
      label    = "Show",
      choices  = c("Histogram", "Dot Plot", "Results Table"),
      selected = "Histogram"
    ),
    uiOutput("view_controls_ui")
  ),
  card(
    fill        = TRUE,
    full_screen = TRUE,
    card_header(textOutput("card_title", inline = TRUE)),
    card_body(
      fill    = TRUE,
      padding = "1rem",
      conditionalPanel(
        "input.active_view === 'Histogram'",
        plotlyOutput("histogram", width = "100%", height = "520px"),
        uiOutput("threshold_slider_ui")
      ),
      conditionalPanel(
        "input.active_view === 'Dot Plot'",
        p("Dot plot will appear here.", class = "text-muted fst-italic p-2")
      ),
      conditionalPanel(
        "input.active_view === 'Results Table'",
        p("Results table will appear here.", class = "text-muted fst-italic p-2")
      )
    )
  )
)

server <- function(input, output, session) {

  shinyDirChoose(input = input, id = "folder", roots = ROOTS, session = session)

  folder_path <- reactive({
    req(input$folder)
    if (is.integer(input$folder)) return(NULL)
    path <- parseDirPath(ROOTS, input$folder)
    if (length(path) == 0) return(NULL)
    path
  })

  flow_set <- reactive({
    req(folder_path())
    load_fcs_folder(folder_path())
  })

  output$file_status <- renderUI({
    req(folder_path())
    fs <- flow_set()
    if (is.null(fs)) {
      p("No .fcs files found in this folder.", class = "text-warning small mt-1")
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
      tags$small("Click a row to filter the histogram to that sample.",
                 class = "text-muted d-block mb-1"),
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
    # selectCallback = TRUE fires input$plate_layout_select on every cell click
    rhandsontable(df, stretchH = "all", rowHeaders = NULL,
                  selectCallback = TRUE) |>
      hot_col("Filename",    readOnly = TRUE) |>
      hot_col("Sample Name", readOnly = FALSE)
  })

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

  # Which sample(s) to show in the histogram.
  # A single click on a row shows just that sample; nothing selected → all samples.
  selected_samples <- reactive({
    fs <- flow_set()
    req(fs)
    all_names <- sampleNames(fs)
    sel <- input$plate_layout_select
    if (is.null(sel) || is.null(sel$r) || sel$r < 0 || sel$r >= length(all_names)) {
      return(all_names)
    }
    # sel$r is 0-based; collect every row in the selection range
    rows <- seq(sel$r, sel$r2) + 1L
    rows <- rows[rows >= 1L & rows <= length(all_names)]
    if (length(rows) == 0L) all_names else all_names[rows]
  })

  # ── View selector ─────────────────────────────────────────────────────────

  output$card_title <- renderText({
    req(input$active_view)
    input$active_view
  })

  output$view_controls_ui <- renderUI({
    req(flow_set())
    channels <- colnames(flow_set())
    if (input$active_view == "Histogram") {
      tagList(
        selectInput("hist_channel", "Channel",
                    choices = channels, selected = channels[1])
      )
    } else if (input$active_view == "Dot Plot") {
      p("Dot plot controls will appear here.", class = "text-muted fst-italic small")
    } else {
      NULL
    }
  })

  # ── Histogram ─────────────────────────────────────────────────────────────

  # Data for the selected channel, filtered to the clicked sample (or all).
  channel_data <- reactive({
    req(flow_set(), input$hist_channel)
    df <- extract_channel_data(flow_set(), input$hist_channel)
    df[df$sample %in% selected_samples(), ]
  })

  # Min/max of the current channel data — shared by slider and drag handler.
  slider_range <- reactive({
    req(channel_data())
    vals <- channel_data()$value
    list(
      lo = floor(min(vals,   na.rm = TRUE)),
      hi = ceiling(max(vals, na.rm = TRUE))
    )
  })

  # Threshold stored for use by other panels later.
  threshold_val <- reactiveVal(NULL)

  # When the slider moves: store the value and update the shape via proxy
  # (no full histogram re-render needed).
  observeEvent(input$threshold, {
    threshold_val(input$threshold)
    plotlyProxy("histogram", session) |>
      plotlyProxyInvoke("relayout", list(
        shapes = list(list(
          type = "line",
          x0 = plog_fwd(input$threshold), x1 = plog_fwd(input$threshold),
          y0 = 0, y1 = 1, yref = "paper",
          line = list(color = "red", dash = "dot", width = 2)
        ))
      ))
  })

  # When the shape is dragged on the plot: sync the slider (which triggers
  # the observer above to also update threshold_val and the proxy shape).
  observeEvent(event_data("plotly_relayout", source = "hist"), {
    d <- event_data("plotly_relayout", source = "hist")
    if (is.null(d) || is.null(d[["shapes[0].x0"]])) return()

    # Shape x is in transformed space; convert back to original data units
    original_x <- plog_inv(as.numeric(d[["shapes[0].x0"]]))
    rng <- slider_range()
    original_x <- max(rng$lo, min(rng$hi, original_x))  # clamp to slider range

    updateSliderInput(session, "threshold", value = round(original_x))
  })

  # Slider UI — range resets when channel or sample selection changes.
  output$threshold_slider_ui <- renderUI({
    rng <- slider_range()
    sliderInput(
      inputId = "threshold",
      label   = "Threshold (drag the slider or the red line on the plot)",
      min     = rng$lo,
      max     = rng$hi,
      value   = round((rng$lo + rng$hi) / 2),
      width   = "100%",
      step    = 1
    )
  })

  # Histogram — only re-renders when channel data changes, NOT when the
  # threshold slider moves (the shape is updated via plotlyProxy instead).
  output$histogram <- renderPlotly({
    req(channel_data())
    df       <- channel_data()
    samples  <- selected_samples()
    subtitle <- if (length(samples) == 1L) basename(samples) else
                  paste(length(samples), "samples")

    p <- ggplot(df, aes(x = value)) +
      geom_histogram(bins = 100, fill = "#4C72B0", colour = NA, alpha = 0.85) +
      scale_x_continuous(
        trans  = pseudo_log_trans(sigma = 1),
        labels = label_number(scale_cut = cut_short_scale())
      ) +
      labs(subtitle = subtitle, x = input$hist_channel, y = "Count") +
      theme_bw(base_size = 13)

    # Use isolate() so the threshold slider does NOT trigger a full re-render;
    # the plotlyProxy observer keeps the shape in sync instead.
    thresh <- isolate(input$threshold)

    plt <- ggplotly(p, source = "hist") |> layout(showlegend = FALSE)

    if (!is.null(thresh)) {
      plt <- plt |> layout(shapes = list(list(
        type = "line",
        x0 = plog_fwd(thresh), x1 = plog_fwd(thresh),
        y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "red", dash = "dot", width = 2)
      )))
    }

    # shapePosition = TRUE makes the vertical line draggable
    plt |> config(edits = list(shapePosition = TRUE))
  })
}

shinyApp(ui, server)
