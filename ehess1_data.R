#!/usr/bin/env Rscript
# Usage: Scrapes EHESS data
# Author: Clement Gorin
# Contact: gorinclem@gmail.com
# Date: June 2021

pacman::p_load(data.table, gtools, pbmcapply, readstata13, rvest, sf, tidyverse, zoo)
setwd("~/Dropbox/research/arthisto/data_ehess")

# Functions ---------------------------------------------------------------

# Displays record
see_record <- function(idehess) {
  cmds <- paste0("open http://cassini.ehess.fr/cassini/fr/html/fiche.php?select_resultat=", idehess)
  sapply(cmds, system)
}

# Extracts the urls of records from the html pages
get_urls <- function(page) {
  urls <- read_html(page) %>% html_nodes("div#liste_g,#liste_d") %>% html_children() 
  regx <- "fiche\\.php\\?select_resultat=\\d+"
  urls <- urls %>% as.character() %>% str_subset(regx) %>% str_extract(regx)
  urls <- urls %>% unique() %>% mixedsort()
  urls <- str_c("http://cassini.ehess.fr/cassini/fr/html/", urls)
  return(urls)
}

# Extracts commune name form record htmls
get_name <- function(html) {
  name <- html %>% html_node("h4:nth-child(1)") %>% html_text() 
  name <- setNames(name, "commune")
  return(name)
}

# Extracts general information form record htmls
get_info <- function(html) {
  info <- html %>% html_node("table.fiche:nth-child(3)") %>% html_table()
  info <- do.call(cbind, c(info))[, -3] %>% t
  info <- info[2, ] %>% str_remove_all("(?<=\\d)\\s(?=\\d)| ha| m") %>% str_replace_all(",", "\\.") %>% setNames(info[1, ])
  xy   <- info["coordonnées"] %>% str_split_fixed("(x|y)\\s+:\\s+", 3) %>% c() %>% setNames(c("crs", "x", "y"))
  info <- c(info[names(info) != "coordonnées"], xy)
  names(info) <- names(info) %>% str_replace_all("\\s", "_") %>% str_remove_all("\\(.*\\)")
  return(info)
}

# Extracts population data form record htmls
get_pop <- function(html) {
  pop <- html %>% html_node("table.recensement")
  if(!is.na(pop)) {
    pop <- pop %>% html_table()
    pop <- do.call(c, c(pop))[-1] %>% str_split_fixed(., "\\s+", 2) %>% t()
    pop <- setNames(pop[2, ], paste0("pop", pop[1, ]))
  }
  return(pop)
}

# Extracts administrative data record htmls
get_admin <- function(html) {
  admin <- html %>% html_node("table.fiche:nth-child(15)")    # Extracts administrative data
  if(!is.na(admin)) {
    admin <- admin %>% html_table() %>% setNames(c("key", "value"))
    admin <- split(admin, 1:nrow(admin))                      # Places all changes in a list
    admin <- lapply(admin, function(row) {
      years <- unlist(str_extract_all(row$value, "\\d{4}"))   # Extracts years
      labs  <- str_split(row$value, str_c(years, collapse = "|")) %>% unlist() %>% tail(-1)
      regx  <- str_c(c("\\[*Tcl\\.", "av\\.", "prov\\.", "\\(Ch\\.-l\\..*\\)", "[^[:alnum:]\\s-\\[']"), collapse = "|")
      labs  <- str_remove_all(labs, regx)                     # Removes annotations (Tcl = perd chef lieu, Ch-l = recoit chef lieu)
      labs  <- str_replace_all(labs, "-", " ")                # Removes dashes
      labs  <- str_trim(labs)                                 # Trims white spaces
      labs  <- str_replace_all(labs, "\\[(.+)\\s*$", "(\\1)") # Puts secondary name into parenthesis
      labs  <- str_remove_all(labs, "\\s*\\[")                # Remove is square bracket is not a secondary name
      labs  <- str_replace_all(labs, "\\s+", " ")             # Replaces multiple spaces with one
      labs  <- str_to_title(labs, locale = "fr")              # Formats names as proper noun
      labs[labs == ""] <- NA                                  # Labels empty values as NA
      if (length(labs) > 0)
        labs <- setNames(labs, paste(row$key, years, sep = "_"))
      return(labs)
    })
    admin <- Reduce(c, admin)
    names(admin) <- names(admin) %>% iconv(to = "ASCII//TRANSLIT") %>% str_remove_all("'")
  }
  return(admin)
}

# Wrapper for extraction utilities
extract_data <- function(html) {
  info <- c(get_name(html), get_info(html), get_admin(html), get_pop(html))
  info <- data.table(data.frame(as.list(info)))
  return(info)
}

# Writes an HTML list as RDS
write_rvest <- function(data, file, ...) {
  data <- lapply(data, as.character)
  saveRDS(data, file, ...)
}

# Reads an HTML list from RDS
read_rvest <- function(file) {
  data <- readRDS(file)
  data <- lapply(data, function(.) tryCatch(read_html(.), error = function(err) NA))
  return(data)
}

# Removes spacial characters from string
remove_special <- function(string) {
  string <- string %>%
    str_replace_all("'", "@") %>%
    iconv(to = "ASCII//TRANSLIT") %>%
    str_replace_all("\\^|'|`", "") %>%
    str_replace_all('"', "") %>%
    str_replace_all("@", "'") %>%
    str_replace_all("@", "'")
  return(string)
}

# Scrapping ---------------------------------------------------------------

if(!file.exists("ehess_scrapping/data_raw.RDS")) {
  # Extracts record urls
  pages <- paste0("http://cassini.ehess.fr/cassini/fr/html/6d_index.php?alpha=", LETTERS)
  urls  <- lapply(pages, get_urls) 
  urls  <- urls %>% unlist() %>% unique() %>% mixedsort()
  saveRDS(urls, "ehess_scrapping/urls.RDS")
  # Extracts record htmls
  if(!exists("htmls")) { # Safety
    htmls <- as.list(rep(NA, length(urls)))
  }
  nrecords <- length(urls)
  saveidx  <- c(seq(1, nrecords, 1000), nrecords)
  # Loops over records
  idehess  <- str_extract(urls, "\\d+")      # Printed
  progress <- seq(nrecords) / nrecords * 100 # Printed
  for(i in which(sapply(htmls, is.na))) {
    htmls[[i]] <- tryCatch(read_html(urls[i], encoding = "UTF-8"), error = function(.) NA)
    print(sprintf("i %i (%.4f%%) - idehess %s %s", i, progress[i], idehess[i], ifelse(is.na(htmls[[i]]), "failure", "success")))
    Sys.sleep(0.05)
    if(i %in% saveidx) {
      write_rvest(htmls, "ehess_scrapping/htmls.RDS")
    }
  }
  # Fix for records that cannot be read (idehess 4251)
  htmls[[4240]] <- read_html("~/Desktop/Le Bignon - Notice Communale.html")
  write_rvest(htmls, "ehess_scrapping/htmls.RDS")
}

data <- pbmclapply(htmls, function(html) tryCatch(extract_data(html), error = function(err) NA), mc.cores = 5L)
saveRDS(data, "ehess_scrapping/data.RDS")

# Loading data ------------------------------------------------------------

urls <- readRDS("ehess_scrapping/urls.RDS") 
ids  <- urls %>% str_remove_all("\\D") %>% as.numeric()
data <- readRDS("ehess_scrapping/data.RDS")
data <- rbindlist(data, fill = T)

# Cleaning information ----------------------------------------------------

# Data
info <- c("commune", "code_insee", "statut")
info <- data[, info, with = F]
info <- info %>% rename(insee_ehess = code_insee)

# Codes INSEE
info$insee_court <- str_remove_all(info$insee_ehess, "\\s*")
info$insee_court <- str_replace_all(info$insee_court, '\\s\\"\\?\\"', "XXX")
info$insee_court <- ifelse(
  str_length(info$insee_court) == 8,
  str_c(str_sub(info$insee_court, 1, 2), str_sub(info$insee_court, 6, 8)),
  info$insee_court
)

# Misc
info[, commune := remove_special(commune)] # Commune
info[, statut  := str_replace_all(remove_special(statut), "anII\\)", "anII) ")] # Statut
info[, idehess := ids]  # Unique id
info[, url     := urls] # URLS

# Order
info <- info[, c("idehess", "commune", "insee_ehess", "insee_court", "url", "statut"), with = F][order(idehess)]
save.dta13(info, "ehess_data/info.dta")

# Cleaning spatial --------------------------------------------------------

# Data
spatial <- c("superficie", "altitude", "crs", "x", "y")
spatial <- data[, spatial, with = F]

# Missing
spatial <- spatial[, lapply(.SD, function(.) str_replace_all(., "\\.{3}\\s/\\s\\.{3}|\\.{3}", ""))]
spatial <- spatial[, lapply(.SD, function(.) str_replace_all(., "(?<=\\d)\\s(?=\\d)", ""))]

# Altitude
spatial$altitude_min <- str_replace_all(spatial$altitude, "(\\d+)\\s/\\s(\\d+)", "\\1")
spatial$altitude_max <- str_replace_all(spatial$altitude, "(\\d+)\\s/\\s(\\d+)", "\\2")
spatial <- spatial[, -"altitude"]

# Numeric
vars <- c("superficie", "x", "y", "altitude_min", "altitude_max")
spatial[, (vars) := lapply(.SD, as.numeric), .SDcols = vars]

# Misc
spatial[, crs := remove_special(crs)]
spatial[nchar(crs) == 0, crs := NA]
spatial[, idehess := ids]
spatial <- spatial[, c("idehess", "crs", "x", "y", "superficie", "altitude_min", "altitude_max"), with = F][order(idehess)]

# Save
save.dta13(spatial, "ehess_data/spatial.dta")

# Creating geometries
spatial <- data.table(read.dta13("ehess_data/spatial.dta")) 
spatial <- merge(spatial, info) %>% select(idehess, commune, insee_court, x, y)
spatial <- spatial[!is.na(x) & !is.na(y)]
spatial <- st_as_sf(spatial, coords = c("x", "y"), crs = 27572) %>% st_transform(3035)

st_write(spatial, "ehess_data/spatial.gpkg", delete_dsn = T)

# Cleaning population -----------------------------------------------------

# Data
pop <- str_subset(names(data), "pop")
pop <- data[, pop, with = F]

# Diagnostics
table(str_subset(unlist(pop), "[0-9]+|\\s", negate = T))

# Recoding information
# - ...  commune n'existant pas à cette date
# - adm. commune recensée avec une autre
# - lac. commune oubliée sur la publication du recensement
# - ill. information illisible
# - empty information non disponible actuellement
# - abs. document (ou pages) disparu(es)

recode_values <- function(col) {
  col <- col %>%
    str_remove_all(",\\d+") %>%                 # Removes digits (2006 census)
    str_remove_all("\\s+") %>%                  # Removes space separators (thousands)
    str_to_lower() %>%                          # Harmonises labels (eg. Lac., lac.)
    str_remove("(?<=[a-z])\\.$") %>%            # Harmonises labels (eg. ill., ill)
    str_replace("avecbarcelonnette", "adm") %>% # Fix wrongly labeled observations (2)
    str_replace("\\.{3,4}", "ine") %>%          # Labels inexistant communes
    str_replace("^$", "ind")                    # Labels non-available yet populations
  return(col)
}

pop <- pop[, lapply(.SD, recode_values)]

# Fills missing values
fill_populations <- function(row) {
  idx <- str_detect(row, "^(ine|adm)$") # Index of values that should not be replaced
  row <- row %>% 
    str_replace_all("^(abs|adm|ill|ind|ine|lac)$", NA_character_) %>%
    as.numeric() %>%
    na.approx(na.rm = F) %>% 
    round()
  row[idx] <- NA
  return(row)
}

filled <- data.table(t(apply(pop, 1, fill_populations)))
pop    <- setnames(data.table(filled), names(pop))
rm(filled)

# Diagnostic
# visdat::vis_miss(pop, warn_large_data = F)

# Misc
pop[, idehess := ids]
pop <- pop[, c("idehess", setdiff(names(pop), "idehess")), with = F][order(idehess)]

# Add population series that are missing (idehess 6958)
# Carnetin
fix <- setNames(as.list(c(6958, 224, 233, 242, 208, 232, 225, 231, 217, 220, 177, 187, 180, 182, 183, 168, 185, 185, 193, 198, 169, 167, 188, 179, 220, 192, 213, 231, 215, 207, 311, 376, 411, 436, NA)), names(pop))
pop[idehess == 6958] <- fix
# Balaruc les Bains
fix <- setNames(as.list(c(2478, rep(NA, 15), c(775, 620, 1008, 1418, 1497,	1625,	1717,	1718,	2028,	1601, 1627,	1724, 1822,	1830,	2957,	4369,	5013,	5688,	6180))), names(pop))
pop[idehess == 2478] <- fix
rm(fix)

# Save
save.dta13(pop, "ehess_data/population.dta")
save.dta13(pop, "../shared_data/ehess/sources/population.dta")

# Cleaning administration -------------------------------------------------

# Data
admin <- str_subset(names(data), "souverainete|departement|district|arrondissement|canton|municipalite")
admin <- data[, admin, with = F]

# Extract duplicated columns
dupcols <- str_subset(names(admin), "\\.\\d+") %>% 
  str_remove("\\.\\d+") %>% 
  unique()

newvals <- lapply(dupcols, function(dupcol) {
  dupnams <- str_subset(names(admin), str_c(dupcol, "(\\.\\d+)?"))
  dupidx  <- !is.na(admin[[dupcol]]) 
  dupvals <- admin[dupidx, mixedsort(dupnams), with = F]
  unqval  <- apply(dupvals, 1, function(row) {
    unqval <- unlist(row) %>% unique() %>% na.omit %>% str_c(collapse = " | ")
    return(unqval)
  })
  newval <- admin$dupcol
  newval[dupidx] <- unqval
  return(newval)
})

admin[, (dupcols) := newvals]
admin <- admin[, str_detect(names(admin), "\\.\\d+", negate = T), with = F]
rm(dupcols, newvals)

# Misc
admin <- admin[, lapply(.SD, remove_special)] # Special characters
admin <- admin[, lapply(.SD, function(.) str_replace_all(., "\\b(\\w+)(?:\\s+\\1\\b)+", "\\1"))] # Removes duplicated word in parenthesis after removing accents
admin <- admin[, lapply(.SD, function(.) str_replace_all(., "^$", NA_character_))] # Labels empty strings as missing
admin[, idehess := ids]
admin <- admin[, c("idehess", setdiff(names(admin), "idehess")), with = F][mixedorder(idehess)]

# Cleaning arrondissement -------------------------------------------------

# Data
arrond <- str_subset(names(admin), "arrondissement")
arrond <- admin[, arrond, with = F]
arrond <- setnames(arrond, str_replace(names(arrond), "arrondissement_", "ar"))
arrond[, (setdiff(str_c("ar", 1791:2006), names(arrond))) := NA] # Adds all years to the changes table
arrond <- arrond[, mixedorder(names(arrond)), with = F]
filled <- t(apply(arrond, 1, function(row) na.locf(row, na.rm = F)))
arrond <- setnames(data.table(filled), names(arrond))
rm(filled)

# Select only the years of census
# surveys <- str_c("ar", c(1793, 1800, 1806, 1821, 1831, 1836, 1841, 1846, 1851, 1856, 1861, 1866, 1872, 1876, 1881, 1886, 1891, 1896, 1901, 1906, 1911, 1921, 1926, 1931, 1936, 1946, 1954, 1962, 1968, 1975, 1982, 1990, 1996, 1999, 2001, 2006))
# arrond  <- arrond[, surveys, with = F]

# Misc
arrond[, idehess := ids]
arrond <- arrond[, c("idehess", setdiff(names(arrond), "idehess")), with = F][mixedorder(idehess)]

save.dta13(arrond, "ehess_data/arrondissements.dta")