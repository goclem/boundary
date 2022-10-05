# !/usr/bin/env Rscript
# Usage: Tracks administrative changes in EHESS records
# Author: Clement Gorin
# Contact: gorinclem@gmail.com
# Date: June 2021

pacman::p_load(data.table, dplyr, gtools, pbmcapply, readstata13, reshape2, rvest, sf, stringr, zoo, xml2)
setwd("~/Dropbox/research/arthisto/data_ehess")

# Functions ---------------------------------------------------------------

# Displays record
see_record <- function(idehess) {
  cmds <- paste0("open http://cassini.ehess.fr/cassini/fr/html/fiche.php?select_resultat=", idehess)
  sapply(cmds, system)
}

# Displays administrative table
see_table <- function(idehess) {
  url <- paste0("http://cassini.ehess.fr/cassini/fr/html/fiche.php?select_resultat=", idehess)
  tab <- url %>% read_html() %>% html_node("table.fiche:nth-child(27)") %>% htmlTable::htmlTable()
  return(tab)
}

# Reads an HTML list from RDS
read_rvest <- function(file, subset = NULL) {
  htmls <- readRDS(file)
  if(!is.null(subset))
    htmls <- htmls[subset]
  htmls <- lapply(htmls, function(html) tryCatch(read_html(html), error = function(err) NA))
  return(htmls)
}

# Removes spacial characters from string
clean_libcos <- function(x) {
  x <- x %>%
    str_replace_all("\u0092", "'") %>%
    str_replace_all("'", "@ ") %>%
    iconv(to = "ASCII//TRANSLIT") %>%
    str_remove_all("\\^|'|`") %>%
    str_remove_all('"') %>%
    str_replace_all("-", " ") %>%
    str_to_title(locale = "fr") %>%
    str_replace_all("@ ", "'") %>%
    str_remove_all("^\\s+|\\s+$")
  return(x)
}

# Encodes administrative changes (used in extract_changes)
encode_change <- function(actions_i, dates_i, mains_i, others_i, com_i) {
  change <- data.table(comA  = NULL, action= NULL, comB  = NULL, date = NULL)
  if(actions_i == "Absorbe") {
    change <- data.table(
      comA   = mains_i,
      action = "merges with",
      comB   = c(com_i, na.omit(others_i)),
      date   = dates_i)
  }
  if(actions_i == "Reunie") {
    change <- data.table(
      comA   = c(com_i, na.omit(others_i)),
      action = "merges with",
      comB   = mains_i,
      date   = dates_i)
  }
  if(actions_i == "Cede") {
    change <- data.table(
      comA   = mains_i,
      action = "created from",
      comB   = c(com_i, na.omit(others_i)),
      date   = dates_i)
  }
  if(actions_i == "Creee") {
    change <- data.table(
      comA   = c(com_i, na.omit(others_i)),
      action = "created from",
      comB   = mains_i,
      date   = dates_i)
  }
  if(actions_i == "Donne") {
    change <- data.table(
      comA   = tail(mains_i, -1),
      action = "merges with",
      comB   = head(mains_i, 1),
      date   = dates_i)
  }
  if(actions_i == "Recoit") {
    change <- data.table(
      comA   = tail(mains_i, -1),
      action = "merges with",
      comB   = com_i,
      date   = dates_i)
  }
  if(actions_i == "Transferee") {
    change <- data.table(
      comA   = c(com_i, na.omit(others_i)),
      action = "merges with",
      comB   = tail(mains_i, 1),
      date   = dates_i)
  }
  return(change)
}

# Extracts administrative changes from EHESS records
extract_changes <- function(id) {
  # Extracts data from HTML
  html <- htmls[[as.character(id)]]                       # Extracts html
  data <- html %>% html_node("table.fiche:nth-child(27)") # Extracts data
  com  <- html %>% html_node("h4:nth-child(1)") %>% html_text() # Extracts commune 
  # Skips records with no data
  if(is.na(data)) return(NULL)                   # Records with no table
  if(html_table(data)[1, 2] == "-") return(NULL) # Records with empty table
  # Extracts communes involved from HTML hyperlinks
  idcom <- data %>% html_nodes("td") %>% .[2] %>% html_nodes("a")
  idcom <- unique(data.table(
    com = c(com, idcom %>% html_text()), # Fix for records 7795, 60688
    id  = c(id,  idcom %>% html_attrs() %>% unlist() %>% str_remove_all("\\D"))))
  # Keeps breaks inside cells
  xml_find_all(data, ".//br") %>% xml_add_sibling("p", "\n")
  xml_find_all(data, ".//br") %>% xml_remove()
  # Separates action substrings
  data <- data %>% html_table() %>% .[1, 2] # Extracts string
  data <- data %>% str_split("\\\n+|\\s//\\s") %>% unlist() # Separates actions
  # Extracts actions and dates
  actions <- "(Absorbe|cède|Créée|donne|Reçoit|réunie|transférée)"
  actions <- data %>% str_extract(actions) %>% clean_libcos() # Extracts actions
  dates   <- data %>% str_extract("(\\d{4}-\\d{4}|\\d{4})") %>% zoo::na.locf() # Extract dates, corrects missing dates after "//"
  # Extracts communes
  coms    <- idcom$com %>% str_replace_all("\\(", "\\\\(") %>% str_replace_all("\\)", "\\\\)") %>% str_sort(decreasing = T) # Fix for com with parenthesis in the name and com with the same substring
  coms    <- paste0("(?<=\\s|^)(", paste(coms, collapse = "|"), ")(?=\\s|\\)|$)") # Regex for all communes
  mains   <- data %>% str_remove("\\(avec\\s.*\\)")  %>% str_extract_all(coms) # Extracts communes inside parenthesis
  others  <- data %>% str_extract("\\(avec\\s.*\\)") %>% str_extract_all(coms) # Extracts communes outside parenthesis
  # Encodes relationships
  changes <- Map(encode_change, actions, dates, mains, others, com)
  # Attributes identifiers
  changes <- lapply(changes, function(change) {
    change$idA <- idcom$id[match(change$comA, idcom$com)]
    change$idB <- idcom$id[match(change$comB, idcom$com)]
    return(change)
  })
  changes <- rbindlist(changes) # Binds changes
  return(changes)
}

# Formats extractd changes in EHESS records (used in compute_changes)
format_change <- function(change, start) {
  if(start) { # Computes the starting id
    if(all(change$action == "created from")) idchange <- str_c(unique(change$idB), collapse = " | ")
    if(all(change$action == "merges with"))  idchange <- unique(change$idA)
  }
  if(!start) { # Computes ids for all other changes
    if(all(change$action == "created from")) idchange <- unique(change$idA)
    if(all(change$action == "merges with"))  idchange <- str_c(unique(change$idB), collapse = " | ")
  }
  change <- data.table(idehess = unique(change$idA), idchange = idchange)
  return(change)
}

# Computes the ids for all changes
compute_changes <- function(record) {
  # Splits changes by date
  changes <- split(record, record$date)
  # Computes the starting id
  start   <- format_change(getElement(changes, 1), start = T)
  # Computes all other changes (the first element of changes must be included)
  others  <- lapply(changes, function(change) format_change(change, start = F))
  changes <- c(list(`1791` = start), others)
  # Merges changes recursively using idehess
  changes <- lapply(names(changes), function(year) setnames(changes[[year]], c("idehess", str_c("idehess", year))))
  changes <- Reduce(function(change0, change1) merge(change0, change1, by = "idehess", all.x = T, allow.cartesian = T), changes)
  # Computes the ids for all years
  changes[, (setdiff(paste0("idehess", 1791:2006), names(changes))) := NA]
  changes <- data.table(cbind(idehess = changes$idehess, t(na.locf(t(changes[, str_sort(names(changes)[-1]), with = F])))))
  return(changes)
}

# Extracts list of changes ------------------------------------------------
# Description: Extracts and encodes changes from EHESS records

# Data
ids     <- readRDS("ehess_scrapping/urls.RDS") %>% str_remove_all("\\D") %>% as.numeric()
htmls   <- setNames(read_rvest("ehess_scrapping/htmls.RDS"), ids)

# Changes
changes <- pbmclapply(ids, function(id) tryCatch(extract_changes(id), error = function(.) NA), mc.cores = 5L) # Computes changes
changes <- setNames(changes, ids)               # Assigns record identifier for checking
changes <- changes[!(sapply(changes, is.null))] # Removes no changes
changes <- unique(rbindlist(changes))           # Removes duplicated changes
changes <- changes[!is.na(comA)]                # Removes changes with missing communes (id: 12537, 30082)
changes[, `:=`(comA = clean_libcos(comA), comB = clean_libcos(comB))] # Cleans commune names
changes[, `:=`(idA  = as.numeric(idA), idB = as.numeric(idB))]        # Sets ids to numeric
changes[, date := sapply(str_split(changes$date, "-"), function(date) round(mean(as.numeric(date))))] # Replaces periods by mean year (!)
changes <- changes[, c("idA", "comA",  "date", "action", "idB", "comB"), with = F][order(idA, date)]  # Orders dataset

saveRDS(changes, "ehess_scrapping/changes.RDS")

# Computes table of changes -----------------------------------------------

pacman::p_load(data.table, dplyr, gtools, pbmcapply, readstata13, reshape2, rvest, sf, stringr, zoo, xml2)
setwd("~/Dropbox/research/arthisto/data_ehess")

# Data
ids     <- readRDS("ehess_scrapping/urls.RDS") %>% str_remove_all("\\D") %>% as.character()
changes <- readRDS("ehess_scrapping/changes.RDS")

# (!) Fixes changes -------------------------------------------------------

# Adds change & corrects date PPC
changes <- rbind(changes, data.table(idA = 20433, comA = "Lux", action = c("merges with", "created from"), date = c(1833, 1869), idB = 36265, comB = "Sevrey"))[order(idA)] # Missing changes in EHESS
changes[comA == "Mouaville" & date == 1847, date := 1857] # Mouaville splits from Bechamps in 1857 (i.e. PPC from population data)
# changes <- changes[!(comB %in% c("Paris"))] # Removes Paris changes

# Merges and splits the same year (no change recorded) Clément
changes <- changes[!(idA == "9148"  & date == 1833)]
changes <- changes[!(idA == "12436" & date == 1837)]
changes <- changes[!(idA == "21033" & date == 1982 & idB == "17099")]

# Changes table -----------------------------------------------------------

# Computes all changes
changes <- split(changes, changes$idA)
changes <- pbmclapply(changes, compute_changes, mc.cores = 5L)
changes <- rbindlist(changes)
changes[, anychange := 1] # Dummy

# Diagnostic
# visdat::vis_miss(changes, warn_large_data = F)

# Appends records with no changes
nochanges <- setdiff(ids, changes$idehess)
nochanges <- data.table(idehess = nochanges)
nochanges[, (names(changes)[-1]) := idehess] 
nochanges[, anychange := 0]

changes <- rbind(changes, nochanges)
changes <- changes[, c("idehess", "anychange", setdiff(names(changes), c("idehess", "anychange"))), with = F][mixedorder(idehess)]

save.dta13(changes, "ehess_data/changes.dta")

# Checks  -----------------------------------------------------------------

# check <- 34629
# check <- 33718
# check <- 16608
# check <- 15611
# t(changes[idehess == check])
# see_record(check)

# Notes -------------------------------------------------------------------

# Absorbe en 1835, (avec Grenois / Taconnay) Huban
# [com] absorbe en [date], (avec [others]) [mains]
# [mains] merge with [com] and [others]
# réunie en 1844, (avec Ayssènes-Broquiès / La Romiguière) à Le True
# [com] reunie en [date], (avec others) à [mains]
# [com] and [others] merge with [mains]
# cède en 1827, (avec Aoste / La Bâtie-Montgascon) Chimilin
# [com] cede en [date], (avec [others]) [mains]
# [mains] split from [com] and [others]
# Créée en 1827, à partir de Aoste / La Bâtie-Montgascon / Corbelin
# [com] creee en [date], (avec [others]) à partir de [mains]
# [com] and [others] created from [mains]
# donne en 1854, à Châteauneuf-de-Randon les h. Laubert / Montbel
# [com] donne en [date], à [mains_1] les [mains_-1]
# [mains_-1] merge with [mains_1] # Destination-based
# Reçoit en 1854, de Allenc les h. Laubert / Montbel
# [com] recoit en [date], de [mains_1] les [mains_-1]
# [mains_-1] merges with [com]
# transférée en 1854, (avec Montbel) de Allenc à Châteauneuf-de-Randon
# [com] tranferee en [date], (avec [others]) de [mains_1] à [mains_2]
# [com] and [others] merge with [mains_2]