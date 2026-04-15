library(shiny)
library(bslib)
library(shinyFiles)
library(rhandsontable)
library(flowCore)
library(ggplot2)
library(plotly)
library(scales)
library(DT)

source("helpers.R")

# Auto-detect all available drives/volumes (home, external drives, network, etc.)
ROOTS <- getVolumes()()

# Pseudo-log transform helpers (sigma = 1, matching the histogram x-axis).
# plotly stores data in transformed space, so shapes need transformed coordinates.
plog_fwd <- function(x) asinh(x)   # original  → plotly axis
plog_inv <- function(y) sinh(y)    # plotly axis → original

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
        plotlyOutput("dotplot", width = "100%", height = "520px")
      ),
      conditionalPanel(
        "input.active_view === 'Results Table'",
        DT::dataTableOutput("results_table", width = "100%")
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
      tags$small("Edit the Sample Name column to label your samples.",
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

  # ── Sample selection ──────────────────────────────────────────────────────
  # Uses input$selected_sample (a selectInput rendered inside view_controls_ui).
  # Returns a character vector of sampleNames to include in channel_data().

  selected_samples <- reactive({
    req(flow_set())
    all_names <- sampleNames(flow_set())
    sel <- input$selected_sample
    if (is.null(sel) || sel == "__all__") return(all_names)
    if (sel %in% all_names) sel else all_names
  })

  # ── View selector ─────────────────────────────────────────────────────────

  output$card_title <- renderText({
    req(input$active_view)
    input$active_view
  })

  output$view_controls_ui <- renderUI({
    req(flow_set())
    channels  <- colnames(flow_set())
    all_names <- sampleNames(flow_set())

    if (input$active_view == "Histogram") {
      tagList(
        selectInput("hist_channel", "Channel",
                    choices = channels, selected = channels[1]),
        selectInput("selected_sample", "Sample",
                    choices  = c("All samples" = "__all__",
                                 setNames(all_names, basename(all_names))),
                    selected = "__all__")
      )
    } else if (input$active_view == "Dot Plot") {
      tagList(
        selectInput("dot_x", "X Axis",
                    choices = channels, selected = channels[1]),
        selectInput("dot_y", "Y Axis",
                    choices = channels, selected = channels[2])
      )
    } else {
      # Results Table view: show which channel and threshold the table is based on
      ch  <- hist_channel_val()
      thr <- threshold_val()
      if (is.null(ch) || is.null(thr)) {
        p("Open the Histogram view first to choose a channel and set a threshold.",
          class = "text-muted fst-italic small")
      } else {
        tagList(
          tags$p(tags$b("Channel: "), ch, class = "small mb-1"),
          tags$p(tags$b("Threshold: "),
                 format(round(thr), big.mark = ",", scientific = FALSE),
                 class = "small mb-0")
        )
      }
    }
  })

  # ── Histogram ─────────────────────────────────────────────────────────────

  # Data for the selected channel, filtered to the chosen sample (or all).
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

  # Persist the histogram channel across view switches so the Results Table
  # can read it even when the Histogram selectInput is not in the DOM.
  hist_channel_val <- reactiveVal(NULL)
  observeEvent(input$hist_channel, hist_channel_val(input$hist_channel))

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

    original_x <- plog_inv(as.numeric(d[["shapes[0].x0"]]))
    rng        <- slider_range()
    original_x <- max(rng$lo, min(rng$hi, original_x))

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

  # ── Results table ─────────────────────────────────────────────────────────

  # Recomputes whenever channel, threshold, or sample labels change.
  results_data <- reactive({
    req(flow_set(), hist_channel_val(), is.numeric(threshold_val()))
    compute_results_table(
      fs         = flow_set(),
      channel    = hist_channel_val(),
      threshold  = threshold_val(),
      sample_map = sample_map()
    )
  })

  output$results_table <- DT::renderDataTable({
    req(results_data())
    df <- results_data()
    DT::datatable(
      df,
      rownames  = FALSE,
      selection = "none",
      options   = list(
        dom        = "t",           # table only — no search bar or pagination chrome
        pageLength = nrow(df),
        ordering   = TRUE,
        scrollX    = TRUE
      )
    ) |>
      DT::formatRound("MFI", digits = 0) |>
      DT::formatStyle(
        "% Above Threshold",
        background         = DT::styleColorBar(c(0, 100), "#4C72B0"),
        backgroundSize     = "98% 88%",
        backgroundRepeat   = "no-repeat",
        backgroundPosition = "center"
      )
  })

  # ── Dot plot ──────────────────────────────────────────────────────────────

  # Subsample + density computed once per channel change, not per render.
  dot_data <- reactive({
    req(flow_set(), input$dot_x, input$dot_y)
    df <- extract_two_channels(flow_set(), input$dot_x, input$dot_y)
    df <- df[is.finite(df$x) & is.finite(df$y), ]

    # Cap at 20 000 points so plotly stays responsive
    if (nrow(df) > 20000L) {
      set.seed(42L)
      df <- df[sample(nrow(df), 20000L), ]
    }

    df$density <- point_density(df$x, df$y)
    df[order(df$density), ]   # low density first → dense points rendered on top
  })

  output$dotplot <- renderPlotly({
    req(dot_data())
    df <- dot_data()

    p <- ggplot(df, aes(x = x, y = y, colour = density)) +
      geom_point(size = 0.4, alpha = 0.6, stroke = 0) +
      scale_colour_gradientn(
        colours = c("#0000FF", "#00BFFF", "#00FF00", "#FFFF00", "#FF0000"),
        name    = "Density"
      ) +
      scale_x_continuous(
        trans  = pseudo_log_trans(sigma = 1),
        labels = label_number(scale_cut = cut_short_scale())
      ) +
      scale_y_continuous(
        trans  = pseudo_log_trans(sigma = 1),
        labels = label_number(scale_cut = cut_short_scale())
      ) +
      labs(x = input$dot_x, y = input$dot_y) +
      theme_bw(base_size = 13) +
      theme(legend.position = "right")

    ggplotly(p) |> layout(showlegend = TRUE)
  })

  # ── Histogram — only re-renders when channel data changes, NOT on slider moves.
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

    thresh <- isolate(input$threshold)
    plt    <- ggplotly(p, source = "hist") |> layout(showlegend = FALSE)

    if (!is.null(thresh)) {
      plt <- plt |> layout(shapes = list(list(
        type = "line",
        x0 = plog_fwd(thresh), x1 = plog_fwd(thresh),
        y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "red", dash = "dot", width = 2)
      )))
    }

    plt |> config(edits = list(shapePosition = TRUE))
  })
}

shinyApp(ui, server)
