##################################################
## Project: Loratario - Shiny Visualisation Application
## Script purpose: Supporting functions for the app
## Date: 23.11.2018
## Author: Ashleigh Myall
##################################################



# -------------------------------------------------------------------

### Code that controls Global variable of the string for locations when passing by URL - see functions get URLs

# For example will be the form of http://pgb.liv.ac.uk/~andrew/crowdsource-server/src/public_html/results/  ***ID*** /psm.mzid
firstPartURL <- "http://pgb.liv.ac.uk/~andrew/crowdsource-server/src/public_html/results"
psmFileType <- "psm.mzid"
locationCsvName <- "locations.csv"


# -------------------------------------------------------------------

### Function to return a subsetted dataframe with new columns calculating

# for statistical anal. Takes only one column of the orignal data frame ('FP')
# which is used to determine between Decoys & hits and takes the string
# used to class the FP rows (default = 'DECOY')

get_df_stat <- function(df_FP_only, decoy_string){
  df_FP_only <- data.frame(df_FP_only)
  colnames(df_FP_only) <- c("Accession")
  df_FP_only <- transform(df_FP_only, FP = ifelse(grepl(decoy_string, Accession),1,0))  # Add an FP column which represents is a DECOY or hit, 1 for decoys
  df_FP_only <- transform(df_FP_only, Decoy = ifelse(grepl(decoy_string, Accession),"Decoy","Target")) # Add another column for Decoy or Target
  df_FP_only[,"Cum.FP"] <- cumsum(df_FP_only$FP) # Add cumulative column to the data frame
  df_FP_only$Hits.Above.Tresh <- seq(1,nrow(df_FP_only),1) # Add a column for the num of hits above a thresh
  df_FP_only <- transform(df_FP_only, TP = Hits.Above.Tresh - 2*Cum.FP)  # Add a column for TP which is = Hits.Above.Tresh - 2*Cum.FP
  df_FP_only <- transform(df_FP_only, FDR = Cum.FP / (TP  + Cum.FP)) # Add Column for FDR
  df_FP_only$Q.val <- return_Q_Val(df_FP_only$FDR,nrow(df_FP_only))  # Add a column for the Q value
  return(df_FP_only)
}

# -------------------------------------------------------------------

### Stats function where there exists an isDecoy Col

# Same as previous but uses the isDecoy column instead of regular expressions to identify decoys

getStatsDfExistingDecoy <- function(df_FP_only){
  
  df_FP_only <- data.frame(df_FP_only)
  colnames(df_FP_only) <- c("isDecoy")
  df_FP_only <- transform(df_FP_only, FP = ifelse(isDecoy,1,0))
  df_FP_only <- transform(df_FP_only, Decoy = ifelse(isDecoy,"Decoy","Target"))
  df_FP_only$isDecoy <- NULL
  df_FP_only[,"Cum.FP"] <- cumsum(df_FP_only$FP)
  df_FP_only$Hits.Above.Tresh <- seq(1,nrow(df_FP_only),1)
  df_FP_only <- transform(df_FP_only, TP = Hits.Above.Tresh - 2*Cum.FP)
  df_FP_only <- transform(df_FP_only, FDR = Cum.FP / (TP  + Cum.FP))
  df_FP_only$Q.val <- return_Q_Val(df_FP_only$FDR,nrow(df_FP_only))
  return(df_FP_only)
}

# -------------------------------------------------------------------

### Function to calculate the column for Q value

# Cycles through taking the current position to the end selecting the
# min value from FDR.                                                   

return_Q_Val <- function(dataCol,n){
  return_col <- rep(0,n)
  for(i in 1:n){
    return_col[i] <- min(dataCol[i:n])
  }
  return(return_col)
}


# -------------------------------------------------------------------

### Function to count the cleavages in each entry row

# Returns a column of the num of cleavages in the adjacent column

get_pep_cleav <- function(peptides_to_process){
  clev_df <- data.frame(peptides_to_process) #create the functions internal dataframe to handle 
  colnames(clev_df) <- c("Peptide")
  clev_df <- transform(clev_df, mod_p = gsub("^R|^K","A", Peptide)) #create a new column cutting the front off of any matching peptides
  clev_df <- transform(clev_df, mod_p_2 = gsub("R$|K$","A", mod_p)) #Create a new column cutting the end off of any matching peptides
  clev_df <- transform(clev_df, mod_p_3 = gsub("KP","A", mod_p_2)) #Create a col replacing where K followed by a P so not identified as a cleav
  clev_df <- transform(clev_df, pep_clev_count = str_count(mod_p_3,"R|K")) #Create a new col which has a count of the number of strings with a space seperating them
  return(clev_df$pep_clev_count)
}


# -------------------------------------------------------------------

### Function to get the current data Set

# Takes the file target and returns the first proccesed daataframe for use. It returns a list, first component is 
# the peptide list and second is the modication dataframe. This function autodetects file type and handles acordingly

get_current_dataSet <- function(in_File,passedUrlData,queryLen){
  # Selecting from the default or the uploaded
  if (is.null(in_File)&(length(queryLen)<1)){ # if both no file has been uploaded and no file passed by the url
    returnList <- list("pep" = NULL, "mod" = NULL, "ptmrsString" = NULL) #Return data as a list without mods
    return(returnList) # give the default dataSet
    
  }else{
    #if no file has been uploaded and there is a passed URL file, then set in_File to the passedUrl and validate
    if(is.null(in_File)&(length(queryLen)>=1)){
      
      
      ## If a CSV handle as a CSV and if mzIdentML handle acccordingly
      if(grepl("mzid",passedUrlData)|grepl(".gz",passedUrlData)){
        
        
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Getting URL Data", value = 0.99)
        
        #save the URL to www file
        download.file(url=passedUrlData,"./dat/urlMzIdentML.mzid" )
        
        #set the address ofthe file to be in the www folder
        
        #with the mzid psm list we filter the mod file by the same logiv using the spectrum ID of the filtered psm list
        filteredPsm <- handleFileMzid('./dat/urlMzIdentML.mzid')
        
        returnList <- list("pep" = filteredPsm, "mod" = getMzidMod('./dat/urlMzIdentML.mzid',filteredPsm), "ptmrsString" = getPtmrsString(filteredPsm))
        file.remove('./dat/urlMzIdentML.mzid') # delete the temp file
        return(returnList)
        
      }else{
        #else assume we are being passed a csv
        peptideDf <- get_url_dat(passedUrlData)
        returnList <- list("pep" = peptideDf, "mod" = getModCsv(peptideDf), "ptmrsString" = getPtmrsString(peptideDf))
        return(returnList)
      }
      
    }
    #check in File
    if(validate_file(in_File)){
      if(grepl("csv",in_File$datapath)){
        peptideDf <- handleFileCsv(in_File)
        returnList <- list("pep" = peptideDf, "mod" = getModCsv(peptideDf), "ptmrsString" = getPtmrsString(peptideDf))
        return(returnList)
        
      }else if(grepl("mzid",in_File$datapath)){ #IF .mzid
        
        filteredPsm <- handleFileMzid(in_File$datapath)
        
        returnList <- list("pep" = filteredPsm, "mod" = getMzidMod(in_File$datapath,filteredPsm), "ptmrsString" = getPtmrsString(filteredPsm))
        return(returnList)
      }else{
        #na (gzip currently which is handled the same as mzid)
        
        filteredPsm <- handleFileMzid(in_File$datapath)
        
        returnList <- list("pep" = filteredPsm, "mod" = getMzidMod(in_File$datapath,filteredPsm), "ptmrsString" = getPtmrsString(filteredPsm))
        return(returnList)
      }
      
    }else{
      shinyalert(title = "Bad upload file!",text = "Please upload a file of: .csv .mzid .gz", type = "warning")
      
      returnList <- list("pep" = NULL, "mod" = NULL, "ptmrsString" = NULL)
      return(returnList) # give the default dataSet
    }
  }
}


# -------------------------------------------------------------------

### Function to Extract PtmRS string from columns if present

# Takes the PSM list, checks for a column name - if exists collapse it and return a string for annotation

getPtmrsString <- function(psmDf){
  if("ptmRS.Result" %in% colnames(psmDf)){
    col <- psmDf$ptmRS.Result
    return(paste(col[col != ""], collapse = ", "))
  }else{
    return(NULL)
  }
}


# -------------------------------------------------------------------

### Function to create a new dataframe for behinds the scenes which has the score set, has been ordered by score

# This creates the servers version of the data frame with a few exra cols, and also ordered by score as needed by the
# stats function for calcualted FDR correctly. It takes the column to set to score as an arguement. An original version
# of the primary dataSet is passed each time to this function so that mutliple columns aren't renamed score.

get_dataSet_withScore <- function(df_to_use,selected_col){
  try(df_to_use <- transform(df_to_use, z = round(Mass / m.z))) #Create a new column of Charge
  try(df_to_use <- transform(df_to_use, Mass = round(z * m.z))) #Create a new column of Charge
  try(df_to_use <- transform(df_to_use, ppm = round(ppm,digits = 4))) #Round the ppm col to 4dp
  colnames(df_to_use)[colnames(df_to_use)==selected_col] <- "Score" #set the score column
  df_to_use <- df_to_use[rev(order(df_to_use$Score)),] #order by score
  df_to_use$counted_cleavages <- as.factor(get_pep_cleav(df_to_use$Peptide)) # Count the cleavages and add to the data set as a new column
  
  #replace spectrum ID to just numeric
  try(df_to_use <- transform(df_to_use, ID = as.numeric(gsub("index=","",spectrumID))))
  return(df_to_use)
}


# -------------------------------------------------------------------

### Function to return min, max and a default of the score range

get_score_range <- function(df_to_use){
  rang <- c(min(df_to_use$pep$Score),max(df_to_use$pep$Score))
  return(rang)
  
}


# -------------------------------------------------------------------

### Function for mzid and gz file input handler:

# Takes an mzid target, imports, merges across the score and peptide list, then filters by max rank=1 and subsets
# for only unique peptide sequences - selecting the top scoring entry - where score is by default always the
# second column of the score df. Then renames some columns and calculates ppm and gets the cleav count

handleFileMzid <- function(targetFile){
  print(targetFile)
  mzid <- openIDfile(targetFile)
  
  mzid_df <- cbind(psms(mzid), score(mzid)[-1]) # col bind the two dataframes (psms with scores leaving out the specID from the score df)
  mzid_df <- filter(mzid_df , rank == 1) #subset only rank 1
  mzid_df <- mzid_df[-which(duplicated(mzid_df$spectrumID)),] #filter ut duplicated rows using spectrumID
  
  #subet for unqiue peptides
  #colnames(mzid_df )[colnames(mzid_df )==colnames(score(mzid)[2])] <- "ScoreFilterCol"
  #mzid_df  <- mzid_df  %>% group_by(sequence) %>% slice(which.max(ScoreFilterCol))
  
  #Format the cols for the app
  colnames(mzid_df)[colnames(mzid_df)=="ScoreFilterCol"] <- colnames(score(mzid)[2])
  mzid_df <- getFormtedColDf(mzid_df)
  
  # Add extra columns
  mzid_df <- transform(mzid_df, ppm = ((m.z - calculatedMassToCharge)*1000000)/m.z) # Add a col for ppm
  
  return(mzid_df)
}


# -------------------------------------------------------------------

### Function to get the modification df from mzIdenML fies

getMzidMod <- function(targetFile,filteredPSM){
  mods <- modifications(openIDfile(targetFile))
  #filter mod df by finding join of spectrum ID
  mods <- subset(mods,spectrumID %in% filteredPSM$spectrumID)
  return(mods)
}


# -------------------------------------------------------------------

### Function to return the CSV dataset

# Firstly the entire peptide column is set to uppercase for comparison (some inpput files were seen to have mixed).
# Tries to detect a recognised colname for filtering the unqiue peptides by max score. If not then tries to use the
# first column containing the string 'Score'. If no match found we do not filter - could select a peptide with the 
# low score

handleFileCsv <- function(in_File){
  
  uploaded_df <- read.csv(in_File$datapath, header=TRUE)
  
  #If the df does not contains more than 5 rows than we assume the wrong format and proceed to check for tab del with skip=0 & skip =1
  if(!NCOL(uploaded_df) >= 5){
    tryCatch({
      uploaded_df <- suppressWarnings(read_delim(in_File$datapath, "\t", #This is especially for MSAmanda
                                                 escape_double = FALSE, trim_ws = TRUE, 
                                                 skip = 0))
    },error = function(e) {stop(safeError(e))})# return a safeError if a parsing error occurs
    
    
    if(!NCOL(uploaded_df) >= 5){
      tryCatch({
        uploaded_df <- suppressWarnings(read_delim(in_File$datapath, "\t", 
                                                   escape_double = FALSE, trim_ws = TRUE, 
                                                   skip = 1))
      },error = function(e) {stop(safeError(e))})# return a safeError if a parsing error occurs
      
      
      if(!NCOL(uploaded_df) >= 5){
        shinyalert(title = "Warning!",text = "Could not read in CSV. Please ensure comma or tab delimited", type = "warning")
        return(NULL) #could read in CSV
      }
    }
  }
  
  uploaded_df <- getFormtedColDf(uploaded_df) #Rename cols
  uploaded_df$Peptide = toupper(uploaded_df$Peptide) #upper case the enitre col
  
  
  #print(head(uploaded_df))
  ## Performing filtering if the columns exist otherwise dont filter
  
  if("Rank" %in% colnames(uploaded_df)){
    uploaded_df <- filter(uploaded_df , Rank == 1) #subset only rank 1
  }
  
  if("rank" %in% colnames(uploaded_df)){
    uploaded_df <- filter(uploaded_df , rank == 1)
  }
  
  if("spectrumID" %in% colnames(uploaded_df)){
    uploaded_df <- uploaded_df[-which(duplicated(uploaded_df$spectrumID)),]
  }
  
  return(uploaded_df)
}


# -------------------------------------------------------------------

### Fucntion to filter for unique peptides

# Filter for best scoring peptide sequence by score, if cannot get Score col then no filtering is done.

getFilteredUniqueDfPep <- function(df){
  if("X.10lgP" %in% colnames(df)){ #This is the general case for crowdsource data (the score col)
    df <- df %>% group_by(Peptide) %>% slice(which.max(X.10lgP))
  }else if(grep("Score", colnames(df))){
    ## R makes this problamatic as it's difficutl to select a colname by a variable
    # Soltuion: set variable column to fixed col name then return it to original name later
    colnames <- colnames(df)
    col_mat <- grepl("Score",colnames)
    selected_default <- colnames[which(col_mat)]
    colnames(df)[colnames(df)== selected_default[1] ] <- "ScoreFilterCol"
    df <- df %>% group_by(Peptide) %>% slice(which.max(ScoreFilterCol))
    colnames(df)[colnames(df)=="ScoreFilterCol"] <- selected_default[1]
    
  }else{
    #no match found, so rather than filter at random we choose to leave
  }
  return(df)
}


# -------------------------------------------------------------------

### Function to read in data from a web URL

get_url_dat <- function(webAddress){
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Reading in Data", value = 0.99)
  webDat <- read.csv(url(webAddress))
  webDat$counted_cleavages <- as.factor(get_pep_cleav(webDat$Peptide))
  return(webDat)
}


# -------------------------------------------------------------------

### Function to check the file uploaded is valid

validate_file <- function(inputFile){
  #First check to make sure of the right file format
  if(grepl(".csv",inputFile$datapath)|grepl(".mzid",inputFile$datapath)|grepl(".gz",inputFile$datapath)){
    #Contains one of the right formats
    return(TRUE)
  }else{
    return(FALSE)
  }
}


# -------------------------------------------------------------------

### function to check if score the score column is a numeric

checkScoreColNum <- function(df_to_use,selected_col){
  colnames(df_to_use)[colnames(df_to_use)==selected_col] <- "Check" #set the score column
  return(!is.numeric(df_to_use$Check))
}


# -------------------------------------------------------------------

### Function to return the url correct form

returnDataUrl <- function(string){
  string <- sub("[?]","",string)
  string <- sub("id=","",string)
  return(paste(firstPartURL,string,psmFileType,sep = "/"))
}


# -------------------------------------------------------------------

### Function to return the url for Server Data csv

returnServerDataCsvUrl <- function(string){
  string <- sub("[?]","",string)
  string <- sub("id=","",string)
  return(paste(firstPartURL,string,locationCsvName,sep = "/"))
}


# -------------------------------------------------------------------

### Functions to check the columns required exist and if not return an error message instead of plot

checkDataUploaded <- function(df){
  validate(
    need(is.null(df), erMessage)
  )
}

checkColExist <- function(dfCol,erMessage){
  validate(
    need(dfCol != "", erMessage)
  )
}


# -------------------------------------------------------------------

### Function to return the current server side version of the dataset

# cbinds the server dataframe with a calculated stats df for plots like FDR and Q curve.
# function of choice depends whether the column isDecoy exists.

returnCurrentServerDF <- function(current_dataSet,col,decoyString){
  returnPepDf <- get_dataSet_withScore(current_dataSet()$pep,col) #Set the score col
  if("isDecoy" %in% colnames(returnPepDf)){ #get the stats cols
    statsDf <-getStatsDfExistingDecoy(returnPepDf$isDecoy)
  }else{
    statsDf <- get_df_stat(returnPepDf$Accession,decoyString)
  }
  returnPepDf <- cbind.data.frame(returnPepDf,statsDf) #Combine pep and stats df
  returnList <- list("pep" = returnPepDf, "mod" = current_dataSet()$mod, "ptmrsString" = current_dataSet()$ptmrsString) #return all
  return(returnList)
  
}


# -------------------------------------------------------------------

### Function to return numeric colnames as a list without the stats DF elements

getColNames <- function(df){
  colNames <- names(Filter(is.numeric,df)) # Get the data set with the appropriate name
  colNames <- colNames[1:(length(colNames)-6)] # remove the last 6 cols which are ffrom the stats calc
  return(colNames)
}


# -------------------------------------------------------------------

### Function for v&h line on the fdr curve

vline <- function(x = 0, color = "grey") {
  list(type = "line", y0 = 0, y1 = 1, yref = "paper",x0 = x, x1 = x, line = list(color = color))
}

hline <- function(y = 0, color = "grey") {
  list(type = "line", x0 = 0, x1 = 1, xref = "paper",y0 = y, y1 = y, line = list(color = color))
}


# -------------------------------------------------------------------

### Function to get intercept

# Approximates the point where FDR percent meets the FDR curve, takes two arguemnts, the df containing the FDR curve
# and the input of FDR percentage

getIntercept <- function(df,fdr){
  intercepDF <- subset(df, Q.val > 0) #no zero interpolation
  interCep <- approx(x = as.numeric(intercepDF$Q.val), y = as.numeric(intercepDF$TP), xout = fdr/100)
  return(round(interCep$y))
}


# -------------------------------------------------------------------

### Function to return a dataframe with counts of mods, only pass it a dataframe of ID and name

# This counts how many mods are present per type of modification and returns a df for plotting a bar chart

getModCount <- function(modificationDf,numOfPep){
  colnames(modificationDf) <- c("spectrumID", "name")
  uniqueMods <- unique(modificationDf[ , 1:2 ] ) #unique 
  pepIdWithAnyMod <- data.frame(unique(modificationDf[ , 1] )) #peptide ids with any modification
  pepWithoutMod <- numOfPep-nrow(pepIdWithAnyMod) # num of peptides with no mods
  returnDf <- rbind.data.frame(data.frame("Var1" = "no mods", "Freq" = pepWithoutMod),as.data.frame(table(unlist(uniqueMods$name))))
  colnames(returnDf) <- c("Modification", "Frequency")
  return(returnDf)
}


# -------------------------------------------------------------------

### Function to return a modification dataframe for the non .mzid files (.csv) enables to treat it the same for ptm plots

# Takes the modification column and spectrum ID then creates a new df where each individual mod in the PTM col
# has its own entry with the relevant spectrum ID to be used as a key. Note: this is for processing csv into the same
# format as mzIdentML mod dataframes

getModCsv <- function(peptideDf){
  
  #Check the colnames exist otherwise return NULL
  if("PTM" %in% colnames(peptideDf)){
    ptmList <- as.data.frame(peptideDf$PTM) #First step is to isolate the ptm column from the peptide list
  }else if("post" %in% colnames(peptideDf)){
    ptmList <- as.data.frame(peptideDf$post) #First step is to isolate the ptm column from the peptide list
  }else{
    return(NULL)
  }
  
  
  ptmList$spectrumID <- seq.int(nrow(ptmList)) #add an ID col
  ptmList <- na.omit(ptmList) #remove na rows
  colnames(ptmList) <- c("name", "spectrumID")
  mods <- cSplit(ptmList, "name", sep = ";", direction = "long") #split so each mod is its own row
  return(mods)
}


# -------------------------------------------------------------------

### Function to return a dataframe with the count of modifications per peptide

getModCountPerPep <- function(modificationDf,numOfPep){
  colnames(modificationDf) <- c("spectrumID", "name")
  
  pepIdWithAnyMod <- data.frame(unique(modificationDf$spectrumID))
  pepWithoutMod <- numOfPep-nrow(pepIdWithAnyMod)
  returnDf <- data.frame("Var1" = 0, "Freq" = pepWithoutMod)
  
  countPerPep <- as.data.frame(table(unlist(modificationDf$spectrumID)))
  countPerNum <- as.data.frame(table(unlist(countPerPep$Freq)))
  
  #to delete repition remove count of 0 from the unlisted mod list
  countPerNum <- countPerNum[countPerNum$Var1 != 0, ]
  
  returnDf <- rbind.data.frame(returnDf,countPerNum)
  colnames(returnDf) <- c("Modifications", "Frequency")
  return(returnDf)
}


# -------------------------------------------------------------------

### Function to take a dataframe (input file) and rename the columns in a correct format

# It became problatic to take files with differeing column names so here is a function to sweep through ranges
# of potential column names and format that in an acceptable way

getFormtedColDf <- function(df){
  
  #General .mzid names
  try(colnames(df)[colnames(df)=="DatabaseAccess"] <- "Accession")
  try(colnames(df)[colnames(df)=="sequence"] <- "Peptide")
  try(colnames(df)[colnames(df)=="experimentalMassToCharge"] <- "m.z")
  try(colnames(df)[colnames(df)=="chargeState"] <- "z")
  try(colnames(df)[colnames(df)=="spectrumID"] <- "ID")
  try(df$ID <- as.numeric(numextract(df$ID))) # only numeric values
  
  
  #MS Amanda renames
  try(colnames(df)[colnames(df)== "Sequence" ] <- "Peptide")
  try(colnames(df)[colnames(df)== "Protein Accessions" ] <- "Accession")
  try(colnames(df)[colnames(df)== "Modifications" ] <- "PTM")
  try(colnames(df)[colnames(df)== "Charge" ] <- "z")
  
  return(df)
}

# -------------------------------------------------------------------

## Function to extract numbers from a string

numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 

# -------------------------------------------------------------------

### Function to return column

# Checks the columns for which should be chosen as the intial scoring column
# works from known commonly used score columns present in a range of file uploads

getInitScoreCol <- function(df){
  colnames <- names(df) # Get the data set with the appropriate name
  
  if(is.element("PSM.level.p.value", colnames)){
    selected_default_vec <- colnames[which(grepl("PSM.level.p.value",colnames))]
    return(selected_default_vec)
  }else if(is.element("peptide.sequence.level.p.value", colnames)){
    return(colnames[which(grepl("peptide.sequence.level.p.value",colnames))])
  }else if(is.element("mzid.Scoring", colnames)){
    return(colnames[which(grepl("mzid.Scoring",colnames))])
  }else if(is.element("X.10lgP", colnames)){
    return(colnames[which(grepl("X.10lgP",colnames))])
  }else if("scr.PEAKS.peptideScore" %in% colnames){
    return(colnames[which(grepl("scr.PEAKS.peptideScore",colnames))])
  }else if("Score" %in% colnames){
    selected_default_vec <- colnames[which(grepl("Score",colnames))]
    return(selected_default_vec[1])
  }else if(is.element("value", colnames)){
    selected_default_vec <- colnames[which(grepl("value",colnames))]
    return(selected_default_vec[1])
  }else{
    return(NULL)
  }
  return(NULL)
}


# -------------------------------------------------------------------

### Function to get spectrum file 

# read in the spectrum file if doesnt exist return NULL

getSpecFile <- function(inFile){
  
  if(is.null(inFile$datapath)){return(NULL)}
  
  # <- checking that an MGF has been referenced too
  if(identical(inFile$datapath, character(0))){
    return(NULL)
  }else{
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Uploading MGF", value = 0.99)
    
    
    ##
    ## Area to insert the call phpMS Command with spectrum
    ##
    
    
    ## Php Command
    phpCom <- "php /home/sgamyall/phpMs-CLI/src/Mgf2Csv.php"
    mgfPath <- inFile$datapath
    #locationOut <- paste("./dat/spec",inFile$size,sep = "") #name of file dependant on size
    locationOut <- paste("./dat/spec",MHmakeRandomString(),sep = "") #with random string after
    
    #delete any current files in that directory
    #unlink("./dat/spec", recursive = TRUE)
    
    #print(paste(phpCom,mgfPath,locationOut, sep = " "))
    
    
    try(system(paste(phpCom,mgfPath,locationOut, sep = " ")))
    
    
    
    ##
    ##  Delete this when uploading!!!
    ##
    
    #locationOut <- "./dat/spec"
    
    return(locationOut)
  }
  
}

# -------------------------------------------------------------------

#Make random string for temp file name 

MHmakeRandomString <- function(n=1, lenght=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    lenght, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}


# -------------------------------------------------------------------

### Function to get a return list

# Unpack the list then combine into one list to be returned by the datapage module

getDataReturnList <- function(current_dataSet_server_side, current_specDataSet){
  return(list("pep" = current_dataSet_server_side()$pep, "mod" = current_dataSet_server_side()$mod, "mgf" = current_specDataSet(), 
              "ptmrsString" = current_dataSet_server_side()$ptmrsString))
}


# -------------------------------------------------------------------

### Function to return a dataframe with select columns

# Use to provide a visual table for the user to select a PSM spectrum from takes input of the columns to include - then returns dataframe

getSelectCols <- function(df, selectCols){
  return(df[,selectCols])
}









# ---------------------------------------------------------------------------------

### Function to return Annotated Spectrum

# Takes Spectrum, fragment df and params to return a new df with an annoated col - currently only functional with b and y ions

returnAnnotatedSpectrum <- function(spectrumDf, fragListDfs, tolBin, binType){
  
  #spectrumDf <- spectrum
  
  #fragListDfs <- ionListString
  
  #tolBin <- 5
  
  #binType <- "Da"
  
  
  if(is.null(spectrumDf)){  #If no spectra return NULL
    return(NULL)
  }
  
  #create a return Df with 2 additional columns for the annotation and anno position
  returnSpectrumDf <- spectrumDf
  returnSpectrumDf$anno <- rep(NA,nrow(spectrumDf))
  returnSpectrumDf$annoPos <- rep(NA,nrow(spectrumDf))
  if(paste(names(fragListDfs[1]),"charge",sep = ".") %in% colnames(as.data.frame(fragListDfs[1]))){
    returnSpectrumDf$annoPosCharge <- rep(NA,nrow(spectrumDf))
  }
  
  if(is.null(fragListDfs)){return(returnSpectrumDf)}
  
  # For each item in the fragment list of dfs annotate the spectrum df
  for (j in 1:length(fragListDfs)) {
    #set the ion name of iteration
    ion <- names(fragListDfs)[j]
    ionDf <- as.data.frame(fragListDfs[j])
    if(paste(ion,"charge",sep = ".") %in% colnames(ionDf)){
      colnames(ionDf) <- c("mz","charge","position","ion")
    }else{
      colnames(ionDf) <- c("mz","position","ion")
    }
    
    
    
    #generate search bins for the ppm method
    if(binType == "ppm"){
      searchRange <- returnUpperLowerPpm(as.numeric(as.character(ionDf$mz)),tolBin,as.numeric(as.character(ionDf$position)))
    }else{# generrate an upper and a lower for the window for the simple daltons 
      searchRange <- returnUpperLower(as.numeric(as.character(ionDf$mz)),tolBin,as.numeric(as.character(ionDf$position)))
    }
    
    
    # for each row find if a value exists in the spectrum files m/z col
    highLightValueVec <- c(NA)
    for(i in 1:nrow(searchRange)){
      ## annotate - which values are true between the ranges returned by the search range function
      returnSpectrumDf$anno[which(returnSpectrumDf$mz == doesExist(searchRange$Upper[i],searchRange$Lower[i],spectrumDf$mz))] <- ion 
      returnSpectrumDf$annoPos[which(returnSpectrumDf$mz == doesExist(searchRange$Upper[i],searchRange$Lower[i],spectrumDf$mz))] <- searchRange$position[i]
      
      if("charge" %in% colnames(ionDf)){ # annotate for the charge if there 
        returnSpectrumDf$annoPosCharge[which(returnSpectrumDf$mz == doesExist(searchRange$Upper[i],searchRange$Lower[i],spectrumDf$mz))] <- "+"
      }
    }
    
  }
  
  return(returnSpectrumDf)
}




# ---------------------------------------------------------------------------------

### Function to create a dataframe of upper and lower values around a column

returnUpperLower <- function(column, tolBin, position){
  upperCol <- column + tolBin
  lowwerCol <- column - tolBin
  returnDB <- data.frame("Upper" = upperCol, "Lower" = lowwerCol, "position" = position)
  return(returnDB)
}


# ---------------------------------------------------------------------------------

### Funtion to create a dataframe of upper and lower around a column for ppm

returnUpperLowerPpm <- function(column, ppm, position){
  x <- data.frame("mz" = column)
  x <- transform(x, Upper = mz + mz* ppm/1000000)
  x <- transform(x, Lower = mz - mz* ppm/1000000)
  x$position <- position
  x$mz <- NULL
  return(x)
}


# ---------------------------------------------------------------------------------

### Function to return closet value within a range or NA is no value exists

doesExist <- function(up, low, col){
  
  #return rows from range
  col <- col[col > low & col < up]
  
  if(length(col) == 1){
    return(col)
  }else if(length(col) > 1){
    #code to select the closest value (if more than 1 row find the cloest to the midpoint of upper and lower)
    return(col[nearestValueSearch((up+low)/2,col)])
  }else{
    return(NA)
  }
}


# ---------------------------------------------------------------------------------

### Function to find nearest value

nearestValueSearch = function(x, w){
  ## A simple binary search algo
  ## Assume the w vector is sorted so we can use binary search
  left = 1
  right = length(w)
  while(right - left > 1){
    middle = floor((left + right) / 2)
    if(x < w[middle]){
      right = middle
    }
    else{
      left = middle
    }
  }
  if(abs(x - w[right]) < abs(x - w[left])){
    return(right)
  }
  else{
    return(left)
  }
}


# ---------------------------------------------------------------------------------

### Function to return a list of dataframes one for each ion present

returnIonDfList <- function(df){
  
  if(is.null(df)){
    return(NULL)
  }
  
  #rename cols <---- edit this later
  if("charge" %in% colnames(df)){ #accoutning for a charge column
    colnames(df) <- c("IonPos","mz","charge")
  }else{
    colnames(df) <- c("IonPos","mz")
  }
  
  
  df$position <- str_extract(df$IonPos, "[:digit:]+$")
  df$ion <- str_replace(df$IonPos, "[:digit:]+$", "")
  df$IonPos <- NULL
  #split the first column for position and ion
  #df <- separate(data = df, col = IonPos, into = c("ion", "position"), sep = 1)
  
  
  
  #lower case all ions
  df$ion <- tolower(df$ion)
  
  #return(df)
  
  # Find unique values in column ion and put into a vec
  ionListDf <- split(df, df$ion)
  
  return(ionListDf)
}


# ---------------------------------------------------------------------------------

### Function to pass a string and return a dataframe for annotations

getAnnoStringAsDf <- function(string){
  
  if(nchar(string) == 0){
    return(NULL) #error handling for empty string
  }
  
  # REmove later
  
  #string <- "b1+-Phospho(8): 734.29, b1+-Phospho(10): 902.38, b1+-Phospho(17): 1610.67, b2+-Phospho(3): 136.54, b2+-Phospho(5): 233.59, b2+-Phospho(8): 367.65, b2+-Phospho(26): 1301.54, b2+-Phospho(28): 1407.58, b2+-Phospho(34): 1758.75, b1+(2): 187.07"
  
  
  s <-  unlist(strsplit(string, ","))
  s <- as.data.frame(s)
  
  
  #seperate for m/z
  s <- separate(data = s, col = s, into = c("annotation", "mz"), sep = ":")
  #convert ti numeric
  s$mz <- as.numeric(s$mz)
  
  #column for charge
  s <- transform(s, charge = ifelse(grepl("[+]",annotation),"+","-"))
  
  #seperate for position (n) where n is the length of ion identified 
  s <- separate(data = s, col = annotation, into = c("anno", "position"), sep = "[(]")
  #remove remaining close brakcet and convert to numeric
  s <- transform(s, position = as.numeric(str_replace(position, "[)]", "")))
  
  #clean the -Phospho from the text
  s <- transform(s, anno = str_replace(anno, "-Phospho", ""))
  
  #remove whitespace
  s$anno <- str_replace_all(s$anno, " ", "")
  
  #remove the charge completely
  s <- separate(data = s, col = anno, into = c("ion", "chargee"), sep = "[:digit:][+-]")
  s$chargee <- NULL
  
  #final split after the first character
  #s <- separate(data = s, col = anno, into = c("ion", "charge"), sep = 1)
  
  s$ion <- paste(s$ion,s$position,sep = "")
  s <- s[c(1,3,4)]
  
  return(s)
  
}




# ---------------------------------------------------------------------------------

### Function to calculate the tolerance under ppm mode

# bin number of parts per million then divided by charge 

getPpmBin <- function(bin){
  return((bin/1000000))
}


# ---------------------------------------------------------------------------------

### Fucntions related to the spectrum plot


##Function to add the missing ion annotation

annoMissIon <- function(p,ionsAnnoRange,ionList,annotatedSpec,colVec,plotHeight){
  
  returnPlot <- p
  
  for (j in 1:length(ionList)){
    ionDf <- ionList[[j]]
    
    allIons <- data.frame("mz" = as.numeric(as.character(ionDf$mz)), "i" = rep(plotHeight,length(ionDf$mz)), "anno" = paste(ionDf$ion, ionDf$position, sep = ""))
    
    returnPlot <- returnPlot %>% 
      
      
      
      add_segments(x = allIons$mz, xend = allIons$mz, y = 0, yend = allIons$i,opacity = 0.25, line=list(color=colVec[j]), showlegend = FALSE) %>%
      
      
      
      
      #add_trace(x = allIons$mz, y = allIons$i,opacity = 0.25, mode = "bar", width = 0.1, marker = list(color = colVec[j]), showlegend = F)%>%
      
      add_text(x = allIons$mz, y = allIons$i,opacity = 0.25, text = allIons$anno, type = 'scatter', mode = 'text', color = colVec[j],
               marker = list(color = colVec[j],size = 0.01), showlegend = F, textposition = "top right",textfont = list(color = colVec[j])) 
    
  }
  
  return(returnPlot)
}





## GGplot rendition of the plotly function as above

annoMissIonGG <- function(p,ionsAnnoRange,ionList,annotatedSpec,colVec,plotHeight){
  
  returnPlot <- p
  
  for (j in 1:length(ionList)){
    ionDf <- ionList[[j]]
    
    allIons <- data.frame("mz" = as.numeric(as.character(ionDf$mz)), "i" = rep(plotHeight,length(ionDf$mz)), "anno" = paste(ionDf$ion, ionDf$position, sep = ""))
    
    
    returnPlot = returnPlot +
      geom_segment(data = allIons, aes(x = mz, y = 0, xend = mz, yend = i),alpha = 0.25,colour = colVec[j]) +
      geom_text(data = allIons, aes(x = mz, y = i, label = anno),alpha = 0.25,colour = colVec[j], nudge_y = 1000)
    
    
  }
  
  return(returnPlot)
  
}







## Function to check if the ions combination exists 

checkIonCombsExist <- function(ionsList){
  if("b" %in% ionsList & "y" %in% ionsList){
    return(TRUE)
  }else if("c" %in% ionsList & "z" %in% ionsList){
    return(TRUE)
  }else if("a" %in% ionsList & "x" %in% ionsList){
    return(TRUE)
  }else{
    return(FALSE)
  }
}




## Function to annoate for the ion ladder

annoSpectrumIonLad <- function(p,ionsList,pepCharVecOrg,plotHeight){
  
  peptideTraceFunc <- function(pObj,trace,colour,heightMult,pepCharVec,ionList){
    
    
    #find a vector of positions 
    ionPos <- as.numeric(as.character(ionList[[trace]]$mz))
    ionLowPos <- append(ionPos, 0, after = 0)[-(1+length(ionPos))]
    annotationPositionVec <- (ionPos + ionLowPos)/2
    
    #Check the num of frags == to length of pepchar vec
    print(length(pepCharVec))
    print(length(annotationPositionVec))
    
    if(length(pepCharVec) != length(annotationPositionVec)){
      return(pObj)
    }
      
    annotationStringDf <- data.frame("mz" = annotationPositionVec, "i" = rep(heightMult*plotHeight,length(pepCharVec)),"anno" = pepCharVec)
    
    #add trace for the peptide characters
    pObj <- pObj %>%
      
      add_annotations(
        x= annotationStringDf$mz,
        y= annotationStringDf$i,
        text = annotationStringDf$anno,
        showarrow = F,
        font = list(color = colour,
                    size = 14)
      )
    
    
    #add trace on the same line for the breaks
    breakDf <- data.frame("mz" = ionPos, "i" = rep(heightMult*plotHeight,length(pepCharVec)))
    
    pObj <- pObj %>% 
      
      add_annotations(
        x= breakDf$mz,
        y= breakDf$i,
        text = "X",
        showarrow = F,
        font = list(color = colour,
                    size = 8),
        opacity = 0.5
      )
  }
  
  
  
  p <- peptideTraceFunc(p,1,'red',1.05,pepCharVecOrg,ionsList)
  p <- peptideTraceFunc(p,2,'green',1,rev(pepCharVecOrg),ionsList)
  
  
  return(p)
  
}



## Function to annoate for the ion ladder --- GG Rendition

annoSpectrumIonLadGg <- function(p,ionsList,pepCharVecOrg,plotHeight){
  
  peptideTraceFunc <- function(pObj,trace,colour,heightMult,pepCharVec,ionList){
    
    #find a vector of positions 
    ionPos <- as.numeric(as.character(ionList[[trace]]$mz))
    ionLowPos <- append(ionPos, 0, after = 0)[-(1+length(ionPos))]
    annotationPositionVec <- (ionPos + ionLowPos)/2
    
    annotationStringDf <- data.frame("mz" = annotationPositionVec, "i" = rep(heightMult*plotHeight,length(pepCharVec)),"anno" = pepCharVec)
    
    #add trace for the peptide characters
    pObj <- pObj + geom_text(data = annotationStringDf, aes(x = mz, y = i, label = anno),colour = colour)
    
    
    #add trace on the same line for the breaks
    breakDf <- data.frame("mz" = ionPos, "i" = rep(heightMult*plotHeight,length(pepCharVec)))
    
    pObj <- pObj + geom_text(data = breakDf, aes(x = mz, y = i, label = "X"),colour = colour,size=1)
  }
  
  
  
  p <- peptideTraceFunc(p,1,'red',1.05,pepCharVecOrg,ionsList)
  p <- peptideTraceFunc(p,2,'green',1,rev(pepCharVecOrg),ionsList)
  
  
  return(p)
  
}


# ---------------------------------------------------------------------------------

### Fucntions to read in and return csv's from phpMS


## Function to read in the spectrum for a given peptide - first check the file exists. If doesn't then return NULL

getSpectrum <- function(rowSelected = 0,mgfLocation){
  if(is.null(mgfLocation)){
    return(NULL) #if no spectra been uploaded then return null
  }
  
  #
  # Read the csv in the location provided
  #
  
  # If no selection by default return the firt row
  if(is.null(rowSelected)){
    spectrumID = 0
  }else{
    spectrumID = rowSelected$ID
  }
  
  
  #Check that the file in qeustion exists
  if(!file.exists(paste(mgfLocation,"/spectrum_",spectrumID,".csv",sep = ""))){
    return(NULL)
  }else{
    
    
    
    csvSpec <- read_csv(paste(mgfLocation,"/spectrum_",spectrumID,".csv",sep = ""))
    colnames(csvSpec) <- c("mz","i")
    spectrum <- data.frame("mz" =  csvSpec$mz, "i" = csvSpec$i)
    
    return(spectrum)
  }
}




## Function to read in the fragment csv

# Pass it the args to call phpMS then read in the return from a destination file 

getFragmentDf <- function(pepSeq,fragMeth,charge){
  
  #return(NULL) # <------CHANGE THIS WHEN RETURNING TO SERVER
  
  ##
  ##  Need to add code so that the different charge states are accounted for
  ##
  
  systemOut <- NA
  
  
  # PhpMs Command
  com <- paste("php /home/sgamyall/phpMs-CLI/src/Fragment.php", toupper(pepSeq), fragMeth, charge, sep = " ")
  
  try(systemOut <- system(com, intern = T))
  
  if(!is.na(systemOut)){
    
    fragDf <- as.data.frame(systemOut)
    
    fragDf <- separate(data = fragDf, col = systemOut, into = c("ion", "mz","amino"), sep = ",")
    
    fragDf$amino <- NULL
    
    return(fragDf)
  }else{
    return(NULL)
  }
  
  
  
}




# ---------------------------------------------------------------------------------

## Code to get volumes for the shiny files button

function () 
{
  osSystem <- Sys.info()["sysname"]
  if (osSystem == "Darwin") {
    volumes <- list.files("/Volumes/", full.names = T)
    names(volumes) <- basename(volumes)
  }
  else if (osSystem == "Linux") {
    volumes <- c(Computer = "/")
    media <- list.files("/media/", full.names = T)
    names(media) <- basename(media)
    volumes <- c(volumes, media)
  }
  else if (osSystem == "Windows") {
    volumes <- system("wmic logicaldisk get Caption", intern = T)
    volumes <- sub(" *\\r$", "", volumes)
    keep <- !tolower(volumes) %in% c("caption", "")
    volumes <- volumes[keep]
    volNames <- system("wmic logicaldisk get VolumeName", 
                       intern = T)
    volNames <- sub(" *\\r$", "", volNames)
    volNames <- volNames[keep]
    volNames <- paste0(volNames, ifelse(volNames == "", "", 
                                        " "))
    volNames <- paste0(volNames, "(", volumes, ")")
    names(volumes) <- volNames
  }
  else {
    stop("unsupported OS")
  }
  if (!is.null(exclude)) {
    volumes <- volumes[!names(volumes) %in% exclude]
  }
  volumes
}
