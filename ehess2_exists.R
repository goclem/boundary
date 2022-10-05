#!/usr/bin/env Rscript
# Usage: Scrapes EHESS data
# Author: Clement Gorin
# Contact: gorinclem@gmail.com
# Date: June 2021

pacman::p_load(data.table, gtools, pbmcapply, readstata13, rvest, sf, zoo, tidyverse)
setwd("~/Dropbox/research/arthisto/data_ehess")

# Functions ---------------------------------------------------------------

# Displays record
see_record <- function(idehess) {
  cmds <- paste0("open http://cassini.ehess.fr/cassini/fr/html/fiche.php?select_resultat=", idehess)
  sapply(cmds, system)
}

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
    str_replace("^$", "ind") %>%                # Labels non-available yet populations
    str_replace("\\d+", "exi")                  # Labels existing observations
  return(col)
}

# Computes status from population -----------------------------------------

# Data
urls <- readRDS("ehess_scrapping/urls.RDS") 
ids  <- urls %>% str_remove_all("\\D") %>% as.numeric()
data <- readRDS("ehess_scrapping/data.RDS")
data <- rbindlist(data, fill = T)

# Diagnostics
# table(str_subset(unlist(pop), "[0-9]+|\\s", negate = T))

status <- str_subset(names(data), "pop")
status <- data[, status, with = F]
status <- setnames(status, str_replace(names(status), "pop", "status"))
status <- status[, lapply(.SD, recode_values)]
status <- status[, lapply(.SD, function(col) str_replace_na(col, "ind"))]

status[, idehess := ids]
status <- status[, c("idehess", setdiff(names(status), "idehess")), with = F][order(idehess)]

# Population fix Carnetin
cols <- setdiff(names(status), c("idehess", "status2006"))
status[idehess == 6958, (cols) := "exi", .SDcols = cols]

# Computes status from changes --------------------------------------------

regx   <- str_c(c("^idehess$", str_extract(names(status)[-1], "\\d+")), collapse = "|")
exists <- data.table(read.dta13("ehess_data/changes.dta"))
exists <- exists[, str_subset(names(exists), regx), with = F]
cols   <- setdiff(names(exists), "idehess")
exists <- cbind(idehess = exists$idehess, exists[, lapply(.SD, function(col) col == exists$idehess), .SDcols = cols])

# Updates status
status  <- cbind(idehess = changes$idehess, data.table(mapply(function(col, exists) {
  col <- ifelse(!exists & col == "exi", "spe", col) # Population and inexistant
  col <- ifelse(!exists & col == "ind", "ine", col) # No population and inexistant
  return(col)
  }, status[, -1], exists[, -1])))

save.dta13(status, "~/Dropbox/research/arthisto/shared_data/ehess/sources/statut.dta")

# Diagnostic
stats <- table(unlist(status[, -1]))
data.table(label = names(stats), number = c(unname(stats)), percent = c(unname(round(prop.table(stats) * 100, 4))))


status[!changes$idehess1856]
