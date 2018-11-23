##################################################
## Project: Loratario - Shiny Visualisation Application
## Script purpose: Script managing each display tab
## Date: 06.11.2018
## Author: Ashleigh Myall
##################################################

# --- List of all Tabs
data_tab <- tabItem(tabName = "data",dataPageInput("datP"))
scatter_tab <- tabItem(tabName = "scatter",singleScatPageInput("scat"))
histograms_tab <- tabItem(tabName = "histograms",histInput("hist"))
cleavages_tab <- tabItem(tabName = "cleavages",cleavPageInput("cleav"))
ptm_tab <- tabItem(tabName = "ptm",ptmPageInput("ptmPage"))
decoy_tab <- tabItem(tabName = "decoy",decoyDisInput("dec"))
score_tab <- tabItem(tabName = "score",scorePageDisplayInput("scoreDis"))
specView_tab <- tabItem(tabName = "specView",specViewInput("specView"))
info_tab <- tabItem(tabName = "info",
                    box(width = 8, title = "Background",
                        h4("Loratario has been developed to support visualising the results of peptide search engines and offering a connection to their original MGF spectrum file."),
                        h4("Focused to work along side the Crowdsource search engine, Loratario also offers support for MS Amanda , MS-GF+ and others. Given a Csv of PSMs, Loratario can accept a wide range of datasets, requiring only a peptide and some form of scoring column."),
                        h4("The developed spectrum Viewer links the PSM ID to a scan ID in the MGF and allows for annotation; either by fragment generation through phpMs or by a direct string annotation which can be manually entered."),
                        a("https://github.com/ashm97/Standalone-Omics-Visualisation-App")
                        ),
                    box(width=4,title = "Profiles")
                    )