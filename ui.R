##################################################
## Project: Loratario - Shiny Visualisation Application
## Script purpose: Main UI Function
## Date: 06.11.2018
## Author: Ashleigh Myall
##################################################

# --- Main User Interface Function
ui <- dashboardPage(skin = "black",
  dashboardHeader(title = "Loratario",disable = F), 
  dashboardSidebar(
    useShinyalert(), # For generating custom alert pop ups 
    sidebarMenu(
      hr(),
      img(src="logo.png",height=220,width=220),
      hr(),
      menuItem("Data", tabName = "data", icon = icon("table")),
      menuItem("PSM Analysis", tabName = "analysis", icon = icon("bolt"),
               menuItem("Scatter", tabName = "scatter"),
               menuItem("Histograms", tabName = "histograms"),
               menuItem("Missed Cleavages", tabName = "cleavages"),
               menuItem("PTM", tabName = "ptm"),
               menuItem("Statistical", tabName = "statistical",icon=icon("line-chart"),
                        menuItem("Decoy", tabName = "decoy"),
                        menuItem("Score", tabName = "score"))),
      menuItem("Spectrum View", tabName = "specView", icon = icon("table")),
      hr(), # Line Break
      menuItem("Information", tabName = "info",icon=icon("info")),
      hr()
    )
  ),
  body # Callig the body for tabs 
)

