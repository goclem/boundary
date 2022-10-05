# !/usr/bin/env Rscript
# Usage: Commutes communes from EHESS records for the period 1861 - 2006
# Author: Clement Gorin
# Contact: gorinclem@gmail.com
# Date: June 2021

pacman::p_load(data.table, gtools, pbmcapply, readstata13, reshape2, sf, stringr, dplyr)
setwd("~/Dropbox/research/arthisto/data_ehess")

# Functions ---------------------------------------------------------------


# Displays record
see_record <- function(idehess) {
  cmds <- paste0("open http://cassini.ehess.fr/cassini/fr/html/fiche.php?select_resultat=", idehess)
  sapply(cmds, system)
}

# Displays idehess that belongs to different communes due to forcing the sets
get_conflicts <- function(year, cos, droponly = T) {
  idcoyear  <- compute_coyear(year)
  conflicts <- idcoyear[((duplicated(idehess) | duplicated(idehess, fromLast = T)) & co %in% cos)]$idehess
  conflicts <- idcoyear[idehess %in% conflicts]
  if(droponly) {
    conflicts <- conflicts[!(co %in% cos)] 
  }
  return(conflicts)
}

# Matches communes
compute_coyear <- function(year) {
  # Extracts variables
  idyear <- changes[, c("idehess", str_c("idehess", year)), with = F]
  idyear <- setnames(idyear, c("idehess", "idehessyear")) # Renames for simplicity
  # For each idehess the idehess it belonged to at [year]
  idyear <- idyear[, .(idehessyear = as.numeric(unlist(str_split(idehessyear, " \\| ")))), by = idehess] # Unpack idehessyear
  # Matches with communes
  coyear <- merge(idyear, pidco, by.x = "idehessyear", by.y = "idehess") # Matches with communes
  coyear <- unique(coyear[, -c("idehessyear")])
  coyear <- unique(rbind(coyear, pidco)) # Adds self idehess to compute the appropriate sets
  coyear <- setnames(coyear[order(idehess)], c("idehess", "co"))
  return(coyear)
}

# Computes sets of geometries strictly containing records
# Note: The problem is geometries containing multiple records
compute_cosets <- function(coyear) {
  cosets <- split(pull(coyear, -idehess), pull(coyear, idehess))
  cosets <- cosets[sapply(cosets, length) > 1] # Only geometries with multiple records
  cos    <- sort(unique(unlist(cosets)))
  # Merges sets recursively
  merge_sets <- function(sets, co) {
    set <- sets[sapply(sets, function(set) is.element(co, set))]
    set <- sort(Reduce(union, set))
    set <- unique(set)
    return(set)
  }
  cnd <- any(duplicated(unlist(cosets))) # Stopping rule
  while(cnd) {
    cosets <- pbmclapply(cos, function(co) merge_sets(cosets, co), mc.cores = 5L) %>% unique()
    cnd    <- any(duplicated(unlist(cosets)))
    print(length(cosets))
  }
  return(cosets)
}

# Gives each set an id with 99[dep][set]
compute_idcosets <- function(sets) {
  # Computes [dep] 99 corresponds to sets across multiple departement
  sets <- unique(sets)
  ids  <- sapply(sets, function(set) {
    ids <- unique(str_sub(str_pad(set, 5, pad = "0"), 1, 2))
    ids <- ifelse(length(ids) > 1, 99, ids)
    ids <- str_c(99, ids)
    return(ids)
  })
  sets <- setNames(sets, ids)
  sets <- sets[order(names(sets))]
  # Computes [set] attributes a number to the sets
  ids  <- duplicated(names(sets))
  ids  <- ave(ids, cumsum(!ids), FUN = cumsum) + 1
  ids  <- str_pad(ids, 2, pad = 0)
  sets <- setNames(sets, str_c(names(sets), ids))
  return(sets)
}

# Forces a particular set
force_coset <- function(cos, cosets) {
  cosets <- lapply(cosets, function(coset) setdiff(coset, cos)) # Removes the coset from all others
  cosets <- c(cosets, list(cos))                                # Appends the new coset
  cosets <- cosets[sapply(cosets, length) > 1]                  # Removes cosets of 1
  return(cosets)
}

# Fixes geometries --------------------------------------------------------

# Data
geo1861 <- st_read("../data_project/co2015.gpkg", quiet = T)
geo1861 <- rename(geo1861, co = co2015, libco = libco2015)

# Paris
paris1861 <- geo1861 %>%
  filter(co %in% c(75101:75120)) %>%
  mutate(co = 75056, libco = "Paris") %>%
  group_by(co, libco) %>%
  summarise(.groups = "drop")
  
# Lyon
lyon1861 <- geo1861 %>%
  filter(co %in% c(69381:69389)) %>%
  mutate(co = 69123, libco = "Lyon") %>%
  group_by(co, libco) %>%
  summarise(.groups = "drop")

# Marseille
marseille1861 <- geo1861 %>%
  filter(co %in% c(13201:13216)) %>%
  mutate(co = 13055, libco = "Marseille") %>%
  group_by(co, libco) %>%
  summarise(.groups = "drop")

# Merges changes
geo1861 <- geo1861 %>%
  filter(!(co %in% c(75101:75120, 69381:69389, 13201:13216))) %>%
  bind_rows(paris1861, lyon1861, marseille1861)

rm(paris1861, lyon1861, marseille1861)

# Check
# st_write(geo1861, "~/Desktop/geo1861_base.gpkg", delete_dsn = T)

# Fixes records -----------------------------------------------------------

# Cleans info
info <- data.table(read.dta13("ehess_data/info.dta", select.cols = c("idehess", "commune", "insee_court")))
info <- info %>% 
  filter(str_detect(insee_court, "2A|2B", negate = T)) %>%                     # Removes Corse
  mutate(coehess = str_replace_all(insee_court, "[^0-9]+", NA_character_)) %>% # Question marks in postcode
  mutate(coehess = as.numeric(coehess)) %>%
  dplyr::select(-insee_court) %>%
  rename(libehess = commune)

# Fixes records with no geometry
fixes <- rbindlist(list(
  # Unmatched records
  data.table(idehess = 61897, coehess = NA,    libehess = "Esquieu"),
  data.table(idehess = 61615, coehess = NA,    libehess = "Halmoy"),
  data.table(idehess = 37225, coehess = NA,    libehess = "Tersou"),
  data.table(idehess = 62016, coehess = NA,    libehess = "Trois Fontaines"),
  data.table(idehess = 61075, coehess = NA,    libehess = "Loucamp"),
  # Other changes (PPC)
  data.table(idehess = 9771,  coehess = 71578, libehess = "Clux"),
  data.table(idehess = 40327, coehess = 71578, libehess = "La Villeneuve")
))
info$coehess[match(fixes$idehess, info$idehess)] <- fixes$coehess
rm(fixes)

# L'Oudon Area was renamed in the 2015 geometries
info[coehess == 14697, coehess := 14472]

# Fixes geometries with no record
fixes <- rbindlist(list(
  data.table(coehess = 14482, idehess = 25911, labehess = "Ouezy"),
  data.table(coehess = 31300, idehess = 19460, labehess = "Lieoux"),
  data.table(coehess = 35317, idehess = 34808, labehess = "Saint Symphorien"),
  data.table(coehess = 51201, idehess = 11323, labehess = "Cuisles"),
  data.table(coehess = 52033, idehess = 2216,  labehess = "Avrecourt"),
  data.table(coehess = 52124, idehess = 9333,  labehess = "Chezeaux"),
  data.table(coehess = 52266, idehess = 24877, labehess = "Laneuville A Remy"),
  data.table(coehess = 52278, idehess = 40330, labehess = "Lavilleneuve Au Roi"),
  data.table(coehess = 52465, idehess = 35517, labehess = "Saulxures"),
  data.table(coehess = 62847, idehess = 39532, labehess = "Verquigneul"),
  data.table(coehess = 67057, idehess = 5010,  labehess = "Bosselshausen"),
  data.table(coehess = 89326, idehess = 29703, labehess = "Rosoy")
  ))
info$coehess[match(fixes$idehess, info$idehess)] <- fixes$coehess
rm(fixes)

# Fixes wrong information in bvarco (Paris, Lyon, Marseille)
bvarco <- data.table(read.dta13("../shared_data/project/bvarco.dta", select.cols = c("co", "co15")))
bvarco <- filter(bvarco, !(co %in% c(75056, 69123, 13055, 71578, 14472, NA))) # Paris, Lyon, Marseille, Clux
fixes  <- data.table(
  co   = c(75056, 69123, 13055, 71578, 14472),
  co15 = c(75056, 69123, 13055, 71578, 14472)
)
bvarco <- rbind(bvarco, fixes)
rm(fixes)

# Matches records to communes (pidco) -------------------------------------

# Matches records to 2015 communes
pidco <- merge(info, bvarco, by.x = "coehess", by.y = "co", all.x = T) # Corrects coehess using bvarco
pidco <- pidco[, -c("coehess"), with = F]
pidco <- rename(pidco, co = co15)                                      # Variable co is co15 using the 1861 definition
pidco <- merge(pidco, st_drop_geometry(geo1861), by = "co", all = T)  # Merges with communes geometries
pidco <- pidco[!is.na(co), c("co", "idehess")]                         # Removes records with no geometries: Camazes, Tersou, Loucamp, Halmoy, Esquieu, Trois Fontaines

saveRDS(rename(pidco, co2015 = co), "ehess_communes/pidco2015.RDS")

# Diagnostics 
# pidco[is.na(pidco$co), ] %>% data.frame()      # Records with no geometries
# pidco[is.na(pidco$idehess), ] %>% data.frame() # Geometries with no record

# Computes changes --------------------------------------------------------

# Data
changes <- data.table(read.dta13("ehess_data/changes.dta"))
changes <- changes[idehess %in% pidco$idehess] # Removes Corse & records with no geometry
changes[, idehess := as.numeric(idehess)]

# Computes cosets
years    <- c(1861, 1866, 1872, 1876, 1881, 1886, 1891, 1896, 1901, 1906, 1911, 1921, 1926, 1931, 1936, 1946, 1954, 1962, 1968, 1975, 1982, 1990, 1996, 1999, 2001, 2006)
idcoyrs  <- setNames(lapply(years, compute_coyear), years) # Matches records to a commune for a given year
cosetyrs <- lapply(idcoyrs, compute_cosets)                # Computes sets of communes

# Utilities
# save <- cosetyrs
# save -> cosetyrs 

# Fixes Paris coset
cosetyrs <- lapply(cosetyrs, function(cosetyr) {
  cosetyr <- force_coset(75056, cosetyr) # Paris
  cosetyr <- force_coset(92012, cosetyr) # Boulogne Billancourt
  cosetyr <- force_coset(94018, cosetyr) # Charenton le Pont
  cosetyr <- force_coset(93070, cosetyr) # Saint Ouen
  cosetyr <- force_coset(93066, cosetyr) # Saint Denis
  cosetyr <- force_coset(93048, cosetyr) # Montreuil
  cosetyr <- force_coset(c(92051, 92044, 92024), cosetyr) # Neuilly, Levallois, Clichy
  return(cosetyr)
})

# Fixes Paris idehess
fixes <- get_conflicts(1861, c(75056), droponly = T)
idcoyrs  <- lapply(idcoyrs, function(idcoyr) {
  idcoyr <- idcoyr[!(idehess %in% fixes$idehess & co %in% fixes$co)] # Drops conflicts which don't have Paris as communes (checked manually)
  return(idcoyr)
})

# Fixes Arles region (i.e. arrondissements)
idx <- years < 1904
cosetyrs[idx] <- lapply(cosetyrs[idx], function(cosetyr) {
  cosetyr <- force_coset(c(13004, 13078, 13097), cosetyr) # Arles, Saint Martin de Crau, Port Saint Louis (i.e. arrondissement)  
  return(cosetyr)
})

# Fixes Arles region idehess
fixes <- get_conflicts(1861, c(13004, 13078, 13097))
idcoyrs[idx] <- lapply(idcoyrs[idx], function(idcoyr) {
  idcoyr <- idcoyr[!(idehess %in% fixes$idehess & co %in% fixes$co)]
  return(idcoyr)
})

idx <- years > 1904 & years < 1925
cosetyrs[idx] <- lapply(cosetyrs[idx], function(cosetyr) {
  cosetyr <- force_coset(c(13004, 13097), cosetyr) # Arles, Saint Martin de Crau
  cosetyr <- force_coset(13078, cosetyr)           # Port Saint Louis
  return(cosetyr)
})
fixes <- get_conflicts(1905, c(13004, 13078, 13097)) # Nothing to fix

idx <- years > 1925
cosetyrs[idx] <- lapply(cosetyrs[idx], function(cosetyr) {
  cosetyr <- force_coset(13004, cosetyr) # Arles
  cosetyr <- force_coset(13078, cosetyr) # Port Saint Louis
  cosetyr <- force_coset(13097, cosetyr) # Saint Martin de Crau
  return(cosetyr)
})
fixes <- get_conflicts(1926, c(13004, 13078, 13097))  # Nothing to fix

saveRDS(cosetyrs, "ehess_communes/cosetyrs_1861_2006.RDS")

# Computes cosets identifiers ---------------------------------------------

# This section computes sets identifiers that are compatible across 1793-1856 and 1861-2001. Only run to create new sets identifiers
# allsets <- lapply(c("ehess_communes/cosetyrs_1793_1856.RDS", "ehess_communes/cosetyrs_1861_2006.RDS"), readRDS)
# allsets <- unlist(allsets, recursive = F)
# 
# idcosets <- compute_idcosets(unlist(allsets, recursive = F)) # Attributes sets identifiers
# idcosets <- lapply(idcosets, function(coset) str_c(sort(coset), collapse = " | "))
# idcosets <- setnames(data.table(melt(idcosets)), c("coset", "idcoset"))[order(idcoset)]
# ifelse(any(duplicated(idcosets$idcoset)), "Identifiers must be unique, add padding!", "All ok") # Security
# saveRDS(idcosets, "ehess_communes/cosetyrs_idcosets.RDS")
# rm(allsets)

idcosets <- readRDS("ehess_communes/cosetyrs_idcosets.RDS")

# Saves database PPC
# tmp <- idcosets[, .(co = as.integer(unlist(str_split(coset, " \\| ")))), by = list(idcoset)]
# tmp[, idcoset:= as.integer(idcoset)]
# tmp <- setnames(tmp, c("cogrp", "co"))
# tmp <- tmp[order(cogrp), ]
# save.dta13(tmp, "../shared_data/ehess/communes/bcogrp.dta")
# rm(tmp)

# Attributes cosets identifiers
cosetyrs <- setnames(data.table(melt(cosetyrs)), c("co", "idxcoset", "year"))
cosetyrs <- cosetyrs[, .(coset = str_c(sort(unlist(co)), collapse = " | ")), by = list(year, idxcoset)]
cosetyrs <- merge(cosetyrs[, -c("idxcoset"), with = F], idcosets, by = "coset", all.x = T)
naDiag(cosetyrs) # All ok!
# Unpacks coset
cosetyrs <- cosetyrs[, .(co = as.numeric(unlist(str_split(coset, " \\| ")))), by = list(year, idcoset)] # Unpack
cosetyrs <- split(cosetyrs[, c("co", "idcoset"), with = F], cosetyrs$year)

# Merges information
cosetyrs <- Map(function(cosetyr, idcoyr) {
  cosetyr <- merge(idcoyr, cosetyr, by = "co", all.x = T)
  cosetyr[is.na(idcoset), idcoset := co] # Replace communes outside of sets with co
  labset <- cosetyr[, .(
    idehessset = str_c(unique(idehess), collapse = " | "),
    coset      = str_c(unique(co),      collapse = " | ")
  ), by = idcoset]
  cosetyr <- unique(cosetyr[, -c("idehess"), with = F])
  cosetyr <- merge(cosetyr, labset, by = "idcoset")
  return(cosetyr)
}, cosetyrs, idcoyrs)

# Merges geometries -------------------------------------------------------

dir.create("~/Desktop/tmp", showWarnings = F)

# Adds geometries
pbmclapply(years, function(year) {
  coyear <- cosetyrs[[as.character(year)]]
  coyear <- merge(geo1861[, c("co")], coyear[, c("co", "idcoset")], by = "co") # Reference geometry
  coyear <- st_collection_extract(coyear, "POLYGON") # Data type fix
  coyear <- mutate(coyear, co = as.integer(co), idcoset = as.integer(idcoset))
  st_write(coyear, sprintf("~/Desktop/tmp/co%s.shp", year), delete_dsn = T, quiet = T)
}, mc.cores = 4L)

# Merges sets of geometries
commands <- sprintf('ogr2ogr -f "GPKG" ~/Desktop/tmp/co%s.gpkg ~/Desktop/tmp/co%s.shp -dialect sqlite -sql "SELECT idcoset, ST_Union(geometry), idcoset FROM co%s GROUP BY idcoset"', years, years, years)
sapply(commands, system)

# Formats output
cosetyrs <- lapply(cosetyrs, function(cosetyr) unique(cosetyr[, -c("co")]))

pbmclapply(years, function(year) {
  coyear <- st_read(sprintf("~/Desktop/tmp/co%s.gpkg", year), quiet = T)
  coyear <- rename(coyear, geometry = `ST_Union(geometry)`)
  coyear <- merge(coyear, cosetyrs[[as.character(year)]], by = "idcoset")
  coyear <- mutate(coyear, idcoset = as.integer(idcoset))
  coyear <- coyear[order(coyear$idcoset), ]
  st_write(coyear, sprintf("ehess_communes/co%s.gpkg", year), delete_dsn = T, quiet = T)
  st_write(coyear, sprintf("../shared_data/ehess/communes/co%s.gpkg", year), delete_dsn = T, quiet = T)
}, mc.cores = 4L)

# Confirms that each idehess has one geometry only
for(year in years) {
  cat(year)
  coyear <- st_read(sprintf("ehess_communes/co%s.gpkg", year), quiet = T)
  test   <- any(duplicated(unlist(str_split(coyear$idehessset, " \\| "))))
  print(ifelse(test, "Duplicate idehess!", "All ok."))
}

unlink("~/Desktop/tmp", recursive = T)