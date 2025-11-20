library(shiny)
# your other libraries

## ---- UI ----
ui <- fluidPage(
  # ...
)

## ---- server ----
server <- function(input, output, session) {
  # ...
}

## ---- app object ----
app <- shinyApp(ui = ui, server = server)

## ---- return the app ----
app
