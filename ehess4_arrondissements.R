# !/usr/bin/env Rscript
# Usage: Computes 1861 arrondissement from EHESS data
# Author: Clement Gorin
# Contact: gorinclem@gmail.com
# Date: June 2021

pacman::p_load(caret, data.table, lwgeom, future.apply, raster, sf, rmapshaper, readstata13, pbmcapply, tidyverse)
options(dplyr.summarise.inform = F, warn = -1)
setwd("~/Dropbox/research/arthisto/data_ehess")
plan(multiprocess, workers = 10L)

# Functions ---------------------------------------------------------------

# Displays record
see_record <- function(idehess) {
  cmds <- paste0("open http://cassini.ehess.fr/cassini/fr/html/fiche.php?select_resultat=", idehess)
  sapply(cmds, system)
}

# Matches arrondissements
match_arrond <- function(year, pidco) {
  regx   <- str_c(c("^idehess$", year), collapse = "|")
  # Extracts variables
  chyear <- changes[, str_subset(names(changes), regx), with = F]
  aryear <- arronds[, str_subset(names(arronds), regx), with = F]
  # For each record the records it belonged to
  index <- chyear %>%
    pull(str_c("idehess", year)) %>% 
    str_split(" \\| ")
  index   <- setNames(unlist(index), Vectorize(rep)(chyear$idehess, sapply(index, length)))
  newar <- arronds[match(index, idehess), str_c("ar", year), with = F]
  newar[, idehess := as.numeric(names(index))]
  newar <- newar[order(newar), c("idehess", str_c("ar", year)), with = F]
  return(newar)
}

# Approximates the border of a commune with multiple arrondissements
approximate_border <- function(dupar) {
  # Geometry to predict
  dupgeom <- dupar %>%
    dplyr::select(geometry) %>% 
    unique()
  # Geometries of neighbours
  ngbgeom <- ar1861 %>% 
    filter(arset1861 %in% unlist(str_split(dupar$arset1861, " \\| "))) %>%
    filter(st_intersects(., dupgeom, sparse=F))
  # Training points
  ngbpts <- ngbgeom %>%
    dplyr::select(arset1861) %>%
    st_segmentize(50) %>% 
    st_cast("POINT") %>% 
    st_as_sf() %>%
    filter(st_intersects(., dupgeom, sparse=F))
  ngbpts <- data.table(arset1861 = as.factor(ngbpts$arset1861), st_coordinates(ngbpts))
  knn_fh <- train(arset1861 ~ X + Y, data = ngbpts, method = "knn", tuneGrid = data.frame(k = 1))
  # Creates observation grid
  dupgrd <- dupgeom %>%
    st_buffer(500) %>%
    st_make_grid(cellsize = 100)
  # Creates observation points
  duppts <- dupgrd %>%
    st_centroid() %>%
    st_coordinates() %>%
    data.table()
  # Predicts class for grid cells
  dupgeom <- st_sf(arset1861 = predict(knn_fh, duppts), geometry = dupgrd) %>%
    group_by(arset1861) %>% 
    summarise() %>%
    st_intersection(dupgeom)
  # Formats output
  dupar <- data.table(st_drop_geometry(dupar))
  dupar <- dupar[, .(arset1861 = unlist(str_split(arset1861, " \\| ")), multiple = 0), by = c("idcoset", "idehessset", "coset")]
  dupar <- merge(dupgeom, dupar, by = "arset1861")
  return(dupar)
}

# Visualises approximate_border
test_border <- function(dupar) {
  dupar <- approximate_border(dupar)
  ngb   <- filter(ar1861, arset %in% dupar$arset & !(co1861 %in% unique(dupar$co1861)))
  plot(dupar[, "arset"], reset = F, border = "red")
  plot(ngb[, "arset"], add = T)
}

# Fixes arrondissement names ----------------------------------------------

co1861 <- st_read("ehess_communes/co1861.gpkg", quiet = T)
ar1861 <- data.table(read.dta13("ehess_data/arrondissements.dta", select.cols = c("idehess", "ar1861"))) %>% mutate(ar1861 = ifelse(ar1861 == "", NA, ar1861))
ch1861 <- data.table(read.dta13("ehess_data/changes.dta",         select.cols = c("idehess", "idehess1861"))) %>% mutate(idehess = as.numeric(idehess))

# Removes records with no communes
filter <- as.numeric(unlist(str_split(co1861$idehessset, " \\| ")))
ar1861 <- ar1861[idehess %in% filter]
ch1861 <- ch1861[idehess %in% filter]

# Match 1861 arrondissement to the idehess of the record in 1861
ch1861 <- ch1861[, .(idehess1861 = as.numeric(unlist(str_split(idehess1861, " \\| ")))), by = idehess]
ar1861 <- merge(ch1861, ar1861, by.x = "idehess1861", by.y = "idehess", all.x = T)
ar1861 <- unique(ar1861[, idehess1861 := NULL])
ar1861 <- ar1861[order(idehess)]
ar1861[str_detect(ar1861, "\\(.*\\)"), ar1861 := str_extract(ar1861, "(?<=\\(\\s?).*(?=\\))")] # Keep the text in parenthesis

fixes <- lapply(list(
  c("^Cassel \\| Hazebrouck$", "Hazebrouck"),
  c("^Tarascon \\| Foix$", "Foix"),
  c("^Nice \\| Puget Theniers$", "Puget Theniers"),
  c("^Ales$", "Alais"),
  c("^Bagneres De Bigorre$", "Bagneres"),
  c("^Bagneres De Bigorre Bagneres De Bigorre$", "Bagneres"),
  c("^Bar Le Duc Bar Le Duc$", "Bar Le Duc"),
  c("^Boulogne$", "Boulogne Sur Mer"),
  c("^Boulay$", "Metz"),
  c("^Chalons$", "Chalons Sur Marne"),
  c("^Chalons Sur Marne Chalons Sur Marne$", "Chalons Sur Marne"),
  c("^Chalon Sur Saone Chalon Sur Saone", "Chalon Sur Saone"),
  c("^Clermont Ferrand Clermont Ferrand$", "Clermont Ferrand"),
  c("^Clermont$", "Clermont En Beauvaisis"),
  c("^Fontenay Le Comte Fontenay Le Comte$", "Fontenay Le Comte"),
  c("^Napoleon Vendee$", "La Roche Sur Yon"),
  c("^Nogaro$", "Mirande"),
  c("^L'argentiere$", "Largentiere"),
  c("^Le Havre Le Havre$", "Le Havre"),
  c("^Le Mans Le Mans$", "Le Mans"),
  c("^Moulins$", "Moulins Sur Allier"),
  c("^Mauleon Licharre$", "Mauleon"),
  c("^Montreuil$", "Montreuil Sur Mer"),  
  c("^Mortagne Sur Huisne Mortagne Sur Huisne", "Mortagne"),
  c("^Mortagne Sur Huisne", "Mortagne"),
  c("^Mirepoix Pamiers$", "Pamiers"),
  c("^Neufchatel$", "Neufchatel En Bray"),
  c("^Oloron$", "Oloron Sainte Marie"),
  c("^Napoleonville$", "Pontivy"),
  c("^Saint Brieuc Saint Brieuc$", "Saint Brieuc"),
  c("^Saint Hipolite$", "Montbeliard"),
  c("^Saint Lo Saint Lo$", "Saint Lo"),
  c("^Saint Amand Montrond Saint Amand Montrond$", "Saint Amand Montrond"),
  c("^Semur En Auxois Semur En Auxois$", "Semur En Auxois"),
  c("^Tartas$", "Mont De Marsan"),
  c("^Verdun Sur Meuse$", "Verdun"),
  c("^Vitry Le Francois Vitry Le Francois$", "Vitry Le Francois"),
  c("^Villefranche Sur Saone Villefranche Sur Saone$", "Villefranche Sur Saone"),
  c("^\\s+|\\s+$", "")
), function(fix) setNames(fix, c("pattern", "replacement")))

for(fix in fixes) {
  ar1861 <- mutate(ar1861, ar1861 = str_replace(ar1861, fix["pattern"], fix["replacement"]))
}

# Fixes Nice, Puget Thenier
# ar1861[idehess %in% c(1480, 4066, 19495, 26412, 26771, 29153),  ar1861 := "Puget Theniers"]  # From geometries
# ar1861[idehess %in% c(4899, 15491, 18868, 29026, 37794, 37906), ar1861 := "Nice"] # From geometries

# Fixes Villefranche using commune information
pidco  <- readRDS("ehess_communes/pidco2015.RDS")
ar1861 <- merge(ar1861, pidco, by = "idehess", all.x = T)
ar1861[ar1861 == "Villefranche" & str_detect(co2015, "^12"), ar1861 := "Villefranche De Rouergue"]
ar1861[ar1861 == "Villefranche" & str_detect(co2015, "^31"), ar1861 := "Villefranche De Lauragais"]
ar1861 <- ar1861[, -c("co2015")]

fixes <- lapply(list(
  c(2377, "Domfront"),
  c(4251, "Nantes"),
  c(7582, "Bagneres"),
  c(20484, "Laon"),
  c(36598, "Nantes"),
  c(60525, "Saint Affrique"),
  c(60552, "Saint Affrique"),
  c(60665, "Tulle"),
  c(60670, "Tulle"),
  c(20783, "Vannes"),
  c(12621, "Lavaur"),
  c(62274, "Lavaur"),
  c(26074, "Lavaur"),
  c(62275, "Lavaur"),
  c(27396, "Arras"),
  c(22886, "Arras")
), function(fix) setNames(fix, c("idehess", "ar1861")))
fixes <- data.table(do.call(rbind, fixes))

ar1861$ar1861[match(fixes$idehess, ar1861$idehess)] <- fixes$ar1861
ar1861 <- unique(ar1861) # Corrections corrected names that are the same
ar1861 <- ar1861[!is.na(ar1861)] # 1 duplicate observation with a NA

# There are duplicates when
# - The record contains two arronsisement
# - The record belongs to two other records whose arrondissements are different

save.dta13(ar1861, "../shared_data/ehess/arrondissements/pidar1861.dta")

# Matches arrondissements -------------------------------------------------

idset  <- data.table(st_drop_geometry(co1861))
idset  <- idset[, .(idehess = as.numeric(unlist(str_split(idehessset, " \\| ")))), by = idcoset]
ar1861 <- merge(idset, ar1861, by = "idehess", all.x = T)
ar1861 <- ar1861[, .(
  idehessset = str_c(unique(na.omit(idehess)), collapse = " | "),
  arset1861  = str_c(unique(na.omit(ar1861)),  collapse = " | ")
), by = idcoset]
ar1861 <- merge(ar1861, dplyr::select(st_drop_geometry(co1861), idcoset, coset), by = "idcoset", all.x = T)

# Fixes Paris sets
fixes <- lapply(list(
  c(75056,  "Paris"),       # Auteuil
  c(92012,  "Saint Denis"), # Batignolles-Monceau
  c(93066,  "Saint Denis"), # Paris
  c(93070,  "Saint Denis"),
  #c(99924,  "Saint Denis"),
  c(93048,  "Sceaux")
  #c(99931,  "Saint Denis"),
  #c(991312, "Arles")
), function(fix) setNames(fix, c("idcoset", "arset1861")))
fixes <- data.table(do.call(rbind, fixes))
fixes[, idcoset := as.numeric(idcoset)]

ar1861$arset1861[match(fixes$idcoset, ar1861$idcoset)] <- fixes$arset1861

ar1861[, multiple := as.numeric(str_detect(arset1861, "\\|"))]

# Fixes multiple arrondissements ------------------------------------------

ar1861 <- merge(dplyr::select(co1861, idcoset, geom), ar1861, by = "idcoset") # Appends geometries

dupars <- ar1861 %>%
  filter(multiple == 1) %>% 
  group_by(idcoset) %>% 
  group_split()

dupars <- future_lapply(dupars, approximate_border)
dupars <- bind_rows(dupars)

ar1861 <- ar1861 %>% 
  filter(multiple == 0) %>%
  bind_rows(., dupars) %>%
  dplyr::select(-multiple) %>%
  rename(ar1861 = arset1861)

# Creates pcoar1861 -------------------------------------------------------

pcoar1861 <- data.table(st_drop_geometry(ar1861))
pcoar1861 <- pcoar1861[, .(co1861 = as.numeric(unlist(str_split(idcoset, " \\| ")))), by = ar1861]
pcoar1861 <- pcoar1861[, c("co1861", "ar1861"), with = F][order(co1861)]
save.dta13(pcoar1861, "~/Dropbox/research/arthisto/shared_data/ehess/pcoar1861.dta")

# Creates geometries ------------------------------------------------------

st_write(st_set_precision(ar1861, 100), "~/Desktop/ar1861.shp", delete_dsn = T, quiet = T)
system('ogr2ogr -f "GPKG" ~/Desktop/ar1861.gpkg ~/Desktop/ar1861.shp -dialect sqlite -sql "SELECT ST_Union(geometry), ar1861 FROM ar1861 GROUP BY ar1861"')

ar1861 <- st_read("~/Desktop/ar1861.gpkg", quiet = T) %>% 
  rename(geometry = `ST_Union(geometry)`) %>%
  arrange(ar1861)

st_write(ar1861, "ehess_arrondissements/ar1861.gpkg", delete_dsn = T)
st_write(ar1861, "../shared_data/ehess/arrondissements/ar1861.gpkg", delete_dsn = T)

# Computes pcaar1861.dta --------------------------------------------------

pacman::p_load(data.table, raster, readstata13)
setwd("~/Dropbox/research/arthisto")

# Data
ca        <- raster("data_project/ca.tif")
idar1861  <- raster("data_ehess/ehess_arrondissements/ar1861.tif")
libar1861 <- st_read("data_ehess/ehess_arrondissements/ar1861.gpkg", quiet = T) %>% st_drop_geometry()
libar1861 <- data.table(idar1861 = as.numeric(row.names(libar1861)), libar1861 = libar1861$ar1861)

# Dataset
pcaar1861 <- data.table(ca = getValues(ca), idar1861 = getValues(ar1861))
pcaar1861 <- pcaar1861[!is.na(ca)]
pcaar1861 <- merge(pcaar1861, libar1861, by = "idar1861", all.x = T)
pcaar1861 <- pcaar1861[, -"idar1861"]
pcaar1861 <- pcaar1861[order(ca)]

naDiag(pcaar1861)
save.dta13(pcaar1861, "shared_data/ehess/arrondissements/pcaar1861.dta")