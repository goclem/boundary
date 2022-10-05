# !/usr/bin/env Rscript
# Usage: Commutes communes from EHESS records for the period 2006 - 2020
# Author: Clement Gorin
# Contact: gorinclem@gmail.com
# Date: June 2021

pacman::p_load(data.table, readstata13, pbmcapply, sf, tidyverse, zoo)
setwd("~/Dropbox/research/arthisto/data_ehess")
options(warn = -1)

# Functions ---------------------------------------------------------------

# Cleans communes label
clean_libcos <- function(libcos) {
  libcos <- libcos %>%
    stringi::stri_trans_general('Latin-ASCII') %>%
    str_replace_all("-", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_remove_all("^\\s+|\\s+$") %>%
    str_to_title(locale = "fr")
  return(libcos)
}

# Reads Stata data file as a data table
read_dta <- function(file, ...) {
  data <- read.dta13(file, ...)
  data <- data.table(data)
  return(data)
}

# Cleans changes ----------------------------------------------------------

# Mergers
merges <- data.table(readxl::read_xlsx("ehess_communes/table_passage_2003_2021.xlsx", 2, skip = 5))
setnames(merges, c("date", "coyear", "idco", "libcoyear", "libidco"))
merges[, (c("libcoyear", "libidco"))   := lapply(.SD, clean_libcos), .SDcols = c("libcoyear", "libidco")]
merges[, (c("date", "coyear", "idco")) := lapply(.SD, as.numeric),   .SDcols = c("date", "coyear", "idco")]
merges <- merges[order(date), c("idco", "coyear", "date", "libidco", "libcoyear"), with = F]
merges[, `:=`(action = "merge", origco = idco)] # Origco is the starting id
merges <- merges[, c("idco", "coyear", "date", "action", "origco")]

# Splits
splits <- data.table(readxl::read_xlsx("ehess_communes/table_passage_2003_2021.xlsx", 3, skip = 5))
splits <- setnames(splits, c("date", "idco", "origco", "libidco", "liborigco"))
splits[, (c("date", "idco", "origco")) := lapply(.SD, as.numeric),   .SDcols = c("date", "idco", "origco")]
splits[, (c("libidco", "liborigco"))   := lapply(.SD, clean_libcos), .SDcols = c("libidco", "liborigco")]
splits <- splits[order(date), c("origco", "idco", "date", "liborigco", "libidco"), with = F]
splits[, `:=`(coyear = idco, action = "split")] # A commune is mapped to itself
splits <- splits[, c("idco", "coyear", "date", "action", "origco")]

# Computes start
changes <- rbind(merges, splits)
start   <- changes[, .SD[which.min(date)], by = idco, .SDcols = c("origco")] # Keeps only the first change by idco
setnames(start, c("idco", "coyear"))

# Computes changes
changes <- split(changes[, c("idco", "coyear")], changes$date)
changes <- c(list(`2003` = start), changes)
changes <- lapply(names(changes), function(year) setnames(changes[[year]], c("idco", str_c("co", year))))
changes <- Reduce(function(changeyear0, changeyear1) merge(changeyear0, changeyear1, by = "idco", all.x = T, allow.cartesian = T), changes)
changes[, (setdiff(paste0("co", 2003:2021), names(changes))) := NA]
changes <- changes[, c("idco", str_sort(names(changes)[-1])), with = F]
changes <- data.table(cbind(idco = changes$idco, t(na.locf(t(changes[, -c("idco")])))))

# Checks ------------------------------------------------------------------

# Note: Two communes split in 2012 and are two separate 2015 geometries. To compute the geometry 2011, the two 2015 geometries should merge. The commune 52504 should have the geometries of the communes 52124 and 52504

# check <- c(14697, 14654, 14472)
# check <- c(55138, 55298)
# check <- c(71578, 71138)
# check <- 14666 # Merge and split
# check <- c(52124, 52504) # Split
# check <- c(49018, 49213, 49245, 49303, 49372)
# geo2015[geo2015$co2015 %in% check, ]
# merges[idco %in% check | origco %in% check, ]
# splits[idco %in% check | origco %in% check, ]
# changes[idco %in% check]
# idcogeom[idco %in% check]
# idcoyr[idco %in% check]
# idcoyr[co2008 %in% check]
# rm(merges, splits, start, check)

# Fixes geometries --------------------------------------------------------

geo2015 <- st_read("../data_project/co2015.gpkg", quiet = T)

# Fixes incorrectly mapped changes ----------------------------------------

# 55138 split from 55298 in 2015, pushing back the change to 2016 to match geometries
changes[idco == 55138, co2015 := 55298]

# 14472 merged with 14654 in 2017
# 14697 merged with 14472 in 2015
vars <- str_subset(names(changes), str_c(2017:2021, collapse = "|"))
changes[idco == 14697, (vars) := 14654]

# 28361 merged with 28015 in 2016
# 28042 merged with 28361 in 2012
vars <- str_subset(names(changes), str_c(2016:2021, collapse = "|"))
changes[idco == 28042, (vars) := 28015]

# 49065 merges with 49080 in 2019
# 49051, 49096, 49105, 49189, 49254, 49335 merged with 49065
vars <- str_subset(names(changes), str_c(2019:2021, collapse = "|"))
changes[idco %in% c(49051, 49096, 49105, 49189, 49254, 49335), (vars) := 49080]

# 49101 merges with 49018 in 2016
# 49380 merges with 49101 in 2013
vars <- str_subset(names(changes), str_c(2016:2021, collapse = "|"))
changes[idco == 49380, (vars) := 49018]

# 49149 merges with 49261 in 2020
# 49094 49154 49279 49346 merge with 49149 in 2016
vars <- str_subset(names(changes), str_c(2020:2021, collapse = "|"))
changes[idco %in% c(49094, 49154, 49279, 49346), (vars) := 49261]

# Computes sets -----------------------------------------------------------

idcogeom <- changes[, c("idco", "co2015")]
years    <- c(2008, 2014, 2020)

# Matches idco to geometries
idcoyrs <- lapply(years, function(year) {
  idcoyr <- changes[, c("idco", str_c("co", year)), with = F] # What idco was each commune in year
  idcoyr <- merge(idcoyr, idcogeom, by = "idco")              # Attributes a geometry to each idco
  idcoyr <- unique(dplyr::select(idcoyr, -idco))
  vars   <- c("coyear", "co2015")
  setnames(idcoyr, vars)
  idcoyr <- idcoyr[, (vars) := lapply(.SD, as.integer), .SDcols = vars]
  return(idcoyr)
})
idcoyrs <- setNames(idcoyrs, years)

# Computes sets of communes
dir.create("~/Desktop/tmp")

cosetsyrs <- lapply(years, function(year) {
  # Assigns sets identifiers
  cosets <- idcoyrs[[as.character(year)]]
  cosets[duplicated(coyear) | duplicated(coyear, fromLast = T), idcoset := coyear]
  cosets[duplicated(co2015) | duplicated(co2015, fromLast = T), idcoset := co2015]
  cosets[is.na(idcoset), idcoset := co2015]
  # Computes labels
  labels <- cosets[, .(
    co2015set = str_c(unique(co2015), collapse = " | "),
    coyearset = str_c(unique(coyear), collapse = " | ")), 
    by = idcoset]
  # Matches geometries
  cosets <- unique(select(cosets, -coyear))
  cosets <- merge(geo2015, cosets, by = "co2015", all.x = T)
  cosets <- mutate(cosets, idcoset = ifelse(!is.na(idcoset), idcoset, co2015))
  cosets <- select(cosets, idcoset)
  # Writes shapefile
  st_write(cosets, sprintf("~/Desktop/tmp/co%s.shp", year), delete_dsn = T, quiet = T)
  return(labels)
})
cosetsyrs <- setNames(cosetsyrs, years)

# Merges sets of geometries
commands <- sprintf('ogr2ogr -f "GPKG" ~/Desktop/tmp/co%s.gpkg ~/Desktop/tmp/co%s.shp -dialect sqlite -sql "SELECT idcoset, ST_Union(geometry), idcoset FROM co%s GROUP BY idcoset"', years, years, years)
sapply(commands, system)

# Formats sets geometries
lapply(years, function(year) {
  coyear <- st_read(sprintf("~/Desktop/tmp/co%s.gpkg", year), quiet = T)
  coyear <- rename(coyear, geometry = `ST_Union(geometry)`)
  coyear <- merge(coyear, cosetsyrs[[as.character(year)]], by = "idcoset", all.x = T)
  coyear <- mutate(coyear, 
    co2015set = ifelse(!is.na(co2015set), co2015set, idcoset),
    coyearset = ifelse(!is.na(coyearset), coyearset, idcoset)
    )
  coyear <- setNames(coyear, str_replace(names(coyear), "year", as.character(year)))
  st_write(coyear, sprintf("ehess_communes/co%s.gpkg", year), delete_dsn = T, quiet = T)
  st_write(coyear, sprintf("../shared_data/ehess/communes/co%s.gpkg", year), delete_dsn = T, quiet = T)
})

# Computes pco2015coyear
lapply(years, function(year) {
  pcocoyear <- st_read(sprintf("ehess_communes/co%s.gpkg", year), quiet = T)
  pcocoyear <- setNames(pcocoyear, str_replace(names(pcocoyear), as.character(year), "year"))
  pcocoyear <- data.table(st_drop_geometry(pcocoyear))
  tmp1 <- pcocoyear[, .(co2015 = as.integer(unlist(str_split(co2015set, " \\| ")))), by = idcoset]
  tmp2 <- pcocoyear[, .(coyear = as.integer(unlist(str_split(coyearset, " \\| ")))), by = idcoset]
  pcocoyear <- merge(tmp1, tmp2, by = "idcoset", all = T)
  pcocoyear <- setNames(pcocoyear, str_replace(names(pcocoyear), "year", as.character(year)))
  save.dta13(pcocoyear, sprintf("../shared_data/ehess/communes/pco2015co%i.dta", year))
})

unlink("~/Desktop/tmp", recursive = T)

# DEPECIATED --------------------------------------------------------------

# # Fix INSEE changes with no geometries
# fixes    <- data.table(
#   oldco = c(14472, 55138, 71578),
#   newco = c(14654, 55298, 71138)
# )
# 
# # Fixes communes in change
# for(i in seq(nrow(fixes))) {
#   vars    <- names(changes)
#   changes <- changes[, (vars) := lapply(.SD, function(var) ifelse(var == fixes$oldco[i], fixes$newco[i], var))]  
#   changes <- unique(changes)
# }
# rm(fixes)
