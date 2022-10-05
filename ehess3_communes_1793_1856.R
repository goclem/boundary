# !/usr/bin/env Rscript
# Usage: Commutes communes from EHESS records for the period 1793 - 1856
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
get_conflicts <- function(year, cos) {
  idcoyear  <- compute_coyear(year)
  conflicts <- idcoyear[((duplicated(idehess) | duplicated(idehess, fromLast = T)) & co %in% cos)]$idehess
  conflicts <- idcoyear[idehess %in% conflicts]
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
geo1793 <- st_read("../data_project/co2015.gpkg", quiet = T)
geo1793 <- rename(geo1793, co = co2015, libco = libco2015)

# Paris 
fixes <- rbindlist(list(
  data.table(co_old = 94018, co_new = 75112, libco_new = "12e Arrondissement"), # Charenton
  data.table(co_old = 94067, co_new = 75112, libco_new = "12e Arrondissement"), # Saint Mandé
  data.table(co_old = 94037, co_new = 75113, libco_new = "13e Arrondissement"), # Gentilly
  data.table(co_old = 94043, co_new = 75113, libco_new = "13e Arrondissement"), # Kremlin Bicêtre
  data.table(co_old = 92049, co_new = 75114, libco_new = "14e Arrondissement"), # Montrouge
  data.table(co_old = 92040, co_new = 75115, libco_new = "15e Arrondissement"), # Issy
  data.table(co_old = 92012, co_new = 75116, libco_new = "16e Arrondissement"), # Boulogne
  data.table(co_old = 92024, co_new = 75117, libco_new = "17e Arrondissement"), # Clichy
  data.table(co_old = 92044, co_new = 75117, libco_new = "17e Arrondissement"), # Levallois-Perret
  data.table(co_old = 92051, co_new = 75117, libco_new = "17e Arrondissement"), # Neuilly
  data.table(co_old = 93061, co_new = 75119, libco_new = "19e Arrondissement"), # Pré Saint Gervais    
  data.table(co_old = 93006, co_new = 75120, libco_new = "20e Arrondissement")  # Bagnolet
))

geo1793$libco[match(fixes$co_old, geo1793$co)] <- fixes$libco_new
geo1793$co[match(fixes$co_old, geo1793$co)]    <- fixes$co_new

# Paris
paris1793 <- st_read("ehess_arrondissements/ar1850_paris_corrected.gpkg", quiet = T) %>%
  rename(co = co15, libco = libco15) %>%
  mutate(co = 75900)

# Arrondissements
arronds1793 <- geo1793 %>%
  filter(co %in% c(75112:75120)) %>%
  group_by(co, libco) %>%
  summarise(.groups = "drop") %>%
  mutate(geom = st_difference(geom, paris1793)) %>%
  mutate(co    = as.numeric(str_replace(co, "(?<=75)1", "9")),
         libco = str_c(libco, " (1793)"))

# Lyon
lyon1793 <- geo1793 %>%
  filter(co %in% c(69381:69389)) %>%
  mutate(co = 69123, libco = "Lyon") %>%
  group_by(co, libco) %>%
  summarise(.groups = "drop")

# Marseille
marseille1793 <- geo1793 %>%
  filter(co %in% c(13201:13216)) %>%
  mutate(co = 13055, libco = "Marseille") %>%
  group_by(co, libco) %>%
  summarise(.groups = "drop")

# Merges changes
geo1793 <- geo1793 %>%
  filter(!(co %in% c(fixes$co_old, 75101:75120, 69381:69389, 13201:13216))) %>%
  bind_rows(arronds1793, paris1793, lyon1793, marseille1793)
  
rm(paris1793, arronds1793, lyon1793, marseille1793, fixes)

# Check
# st_write(geo1793, "~/Desktop/geo1793_base.gpkg", delete_dsn = T))

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
  # Paris
  data.table(idehess = 26207, coehess = 75900, libehess = "Paris"),
  data.table(idehess = 3719,  coehess = 75912, libehess = "Bercy"),
  data.table(idehess = 8423,  coehess = 75912, libehess = "Charenton Le Pont"),
  data.table(idehess = 33100, coehess = 75912, libehess = "Saint Mande"),
  data.table(idehess = 15323, coehess = 75913, libehess = "Gentilly"),
  data.table(idehess = 18163, coehess = 75913, libehess = "Le Kremlin Bicetre"),
  data.table(idehess = 23851, coehess = 75914, libehess = "Montrouge"),
  data.table(idehess = 16141, coehess = 75915, libehess = "Grenelle"),
  data.table(idehess = 17659, coehess = 75915, libehess = "Issy Les Moulineaux"),
  data.table(idehess = 39042, coehess = 75915, libehess = "Vaugirard"),
  data.table(idehess = 1960,  coehess = 75916, libehess = "Auteuil"),
  data.table(idehess = 5199,  coehess = 75916, libehess = "Boulogne Billancourt"),
  data.table(idehess = 62049, coehess = 75916, libehess = "Chaillot"),
  data.table(idehess = 26284, coehess = 75916, libehess = "Passy"),
  data.table(idehess = 2924,  coehess = 75917, libehess = "Batignolles Monceau"),
  data.table(idehess = 9727,  coehess = 75917, libehess = "Clichy"),
  data.table(idehess = 24784, coehess = 75917, libehess = "Neuilly Sur Seine"),
  data.table(idehess = 19362, coehess = 75917, libehess = "Levallois Perret"),
  data.table(idehess = 8165,  coehess = 75918, libehess = "La Chapelle"),
  data.table(idehess = 62051, coehess = 75918, libehess = "Le Roule"),
  data.table(idehess = 23671, coehess = 75918, libehess = "Montmartre"),
  data.table(idehess = 3564,  coehess = 75919, libehess = "Belleville"),
  data.table(idehess = 40600, coehess = 75919, libehess = "La Villette"),
  data.table(idehess = 27988, coehess = 75919, libehess = "Le Pre Saint Gervais"),
  data.table(idehess = 2378,  coehess = 75920, libehess = "Bagnolet"),
  data.table(idehess = 8509,  coehess = 75920, libehess = "Charonne"),
  data.table(idehess = 62050, coehess = 75920, libehess = "Menilmontant"),
  # Other changes (PPC)
  data.table(idehess = 9771,  coehess = 71578, libehess = "Clux"),
  data.table(idehess = 40327, coehess = 71578, libehess = "La Villeneuve"),
  # Unmatched records
  data.table(idehess = 61075, coehess = NA,    libehess = "Loucamp"),
  data.table(idehess = 61897, coehess = NA,    libehess = "Esquieu"),
  data.table(idehess = 61615, coehess = NA,    libehess = "Halmoy"),
  data.table(idehess = 37225, coehess = NA,    libehess = "Tersou"),
  data.table(idehess = 62016, coehess = NA,    libehess = "Trois Fontaines")
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
  co   = c(75900, 75912:75920, 69123, 13055, 71578, 14472),
  co15 = c(75900, 75912:75920, 69123, 13055, 71578, 14472)
)
bvarco <- rbind(bvarco, fixes)
rm(fixes)

# Matches records to communes (pidco) -------------------------------------

# Matches records to 2015 communes
pidco <- merge(info, bvarco, by.x = "coehess", by.y = "co", all.x = T) # Corrects coehess using bvarco
pidco <- pidco[, -c("coehess"), with = F]
pidco <- rename(pidco, co = co15)                                      # Variable co is co15 using the 1861 definition
pidco <- merge(pidco, st_drop_geometry(geo1793), by = "co", all = T)  # Merges with communes geometries
pidco <- pidco[!is.na(co), c("co", "idehess")]                         # Removes records with no geometries: Camazes, Tersou, Loucamp, Halmoy, Esquieu, Trois Fontaines
# Diagnostics 
# pidco[is.na(pidco$co), ] %>% data.frame()      # Records with no geometries
# pidco[is.na(pidco$idehess), ] %>% data.frame() # Geometries with no record

# Computes changes --------------------------------------------------------

# Data
changes <- data.table(read.dta13("ehess_data/changes.dta"))
changes <- changes[idehess %in% pidco$idehess] # Removes Corse & records with no geometry
changes[, idehess := as.numeric(idehess)]

# Computes cosets
years    <- c(1793, 1800, 1806, 1821, 1831, 1836, 1841, 1846, 1851, 1856)
idcoyrs  <- setNames(lapply(years, compute_coyear), years) # Matches records to a commune for a given year
cosetyrs <- lapply(idcoyrs, compute_cosets)                # Computes sets of communes

# Utilities
# save <- cosetyrs
# save -> cosetyrs 

# Fixes Paris coset
cosetyrs <- lapply(cosetyrs, function(cosetyr) {
  cosetyr <- force_coset(75900, cosetyr) # Paris
  cosetyr <- force_coset(75916, cosetyr) # Paris 16e arrondissement
  cosetyr <- force_coset(75918, cosetyr) # Paris 18e arrondissement
  cosetyr <- force_coset(75920, cosetyr) # Paris 20e arrondissement
  return(cosetyr)
})

# Fixes Paris idehess
# get_conflicts(1800, 75900)
# get_conflicts(1800, 75916)
# get_conflicts(1800, 75918)
# get_conflicts(1800, 75920)

idcoyrs <- lapply(idcoyrs, function(idcoyr) {
  idcoyr <- idcoyr[!(idehess %in% c(62049, 62050, 62051) & co == 75900)]
  idcoyr <- idcoyr[!(idehess == 19576 & co == 75920)]
  return(idcoyr)
})

saveRDS(cosetyrs, "ehess_communes/cosetyrs_1793_1856.RDS")

# Computes cosets identifiers ---------------------------------------------

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

# Attributes cosets identifiers
cosetyrs <- setnames(data.table(melt(cosetyrs)), c("co", "idxcoset", "year"))
cosetyrs <- cosetyrs[, .(coset = str_c(sort(unlist(co)), collapse = " | ")), by = list(year, idxcoset)]
cosetyrs <- merge(cosetyrs[, -c("idxcoset"), with = F], idcosets, by = "coset", all.x = T)

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

dir.create("~/Desktop/tmp")

# Adds geometries
pbmclapply(years, function(year) {
  coyear <- cosetyrs[[as.character(year)]]
  coyear <- merge(geo1793[, c("co")], coyear[, c("co", "idcoset")], by = "co") # Reference geometry
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
  check  <- unlist(str_split(coyear$idehessset, " \\| "))
  test   <- any(duplicated(check))
  # coyear[str_detect(coyear$idehessset, "^21 "), ]
  print(ifelse(test, "Duplicate idehess!", "All ok."))
}

unlink("~/Desktop/tmp", recursive = T)