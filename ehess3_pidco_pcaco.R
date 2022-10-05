# !/usr/bin/env Rscript
# Usage: Commutes communes from EHESS records for the period 1861 - 2006
# Author: Clement Gorin
# Contact: gorinclem@gmail.com
# Date: June 2021

pacman::p_load(data.table, gtools, pbmcapply, readstata13, reshape2, sf, stringr, dplyr, zoo)
setwd("~/Dropbox/research/arthisto/data_ehess")

# Functions ---------------------------------------------------------------

duplicates <- function(vec) duplicated(vec) | duplicated(vec, fromLast = T)

# Displays record
see_record <- function(idehess) {
  cmds <- paste0("open http://cassini.ehess.fr/cassini/fr/html/fiche.php?select_resultat=", idehess)
  sapply(cmds, system)
}

# Computes pidco ----------------------------------------------------------

years   <- c(1793, 1800, 1806, 1821, 1831, 1836, 1841, 1846, 1851, 1856, 1861, 1866, 1872, 1876, 1881, 1886, 1891, 1896, 1901, 1906, 1911, 1921, 1926, 1931, 1936, 1946, 1954, 1962, 1968, 1975, 1982, 1990, 1996, 1999, 2001, 2006)
pidco2015 <- readRDS("ehess_communes/pidco2015.RDS")

pidcoyrs <- pbmclapply(years, function(year) {
  idcoyear <- st_read(sprintf("ehess_communes/co%s.gpkg", year), quiet = T)
  idcoyear <- unique(data.table(st_drop_geometry(idcoyear))[, c("idehessset", "idcoset")])
  idcoyear <- idcoyear[, .(idehess = as.numeric(unlist(str_split(idehessset, " \\| ")))), by = idcoset] # Unpacks
  idcoyear <- idcoyear[order(idehess), c("idehess", "idcoset")]
  idcoyear <- setnames(idcoyear, c("idehess", str_c("co", year)))
  return(idcoyear)
}, mc.cores = 4L)

pidcoyrs <- c(pidcoyrs, list(pidco2015))
pidcoyrs <- Reduce(function(x, y) merge(x, y, by = "idehess", all = T), pidcoyrs)
naDiag(pidcoyrs)

save.dta13(pidcoyrs, "../shared_data/ehess/communes/pidco.dta")

# Compute rasters and pcaco -----------------------------------------------

pacman::p_load(data.table, dplyr, future.apply, raster, readstata13, sf, stringr)
setwd("~/Dropbox/research/arthisto/data_ehess")
plan(multiprocess, workers = 10L)

fr     <- st_read("../data_project/fr15.gpkg", quiet = T)
ca     <- raster("../data_project/ca.tif")

border <- st_cast(fr, "MULTILINESTRING")
xy     <- rasterToPoints(ca) %>% data.table()
xy     <- split(xy, sort(1:nrow(xy) %% 10))
xy     <- lapply(xy, function(chunk) st_set_crs(st_as_sf(chunk, coords = c(1, 2)), 3035))
years  <- c(1800, 1861, 1793, 1866) # Priority
years  <- c(years, setdiff(c(1793, 1800, 1806, 1821, 1831, 1836, 1841, 1846, 1851, 1856, 1861, 1866, 1872, 1876, 1881, 1886, 1891, 1896, 1901, 1906, 1911, 1921, 1926, 1931, 1936, 1946, 1954, 1962, 1968, 1975, 1982, 1990, 1996, 1999, 2001, 2006), years))
years  <- c(2008, 2014, 2020)
  
for(year in years) {
  print(year)
  coyear <- st_read(sprintf("ehess_communes/co%s.gpkg", year), quiet = T)
  coyear <- dplyr::select(coyear, idcoset) %>% mutate(idcoset = as.integer(idcoset)) # 1793 - 2001 version
  coyear <- st_transform(coyear, 3035)
  # Computes buffer
  index  <- st_intersects(coyear, border) # (!) time
  buffer <- coyear[sapply(index, length) > 0, ]
  buffer <- st_buffer(buffer, 250)
  buffer <- st_difference(buffer, fr)
  coyear <- bind_rows(coyear, buffer)
  rm(index, buffer)
  # Rasterize
  cayear <- future_lapply(xy, function(chunk) st_join(chunk, coyear, join = st_intersects), future.seed = T) # 2 mins  
  cayear <- data.table(do.call(rbind, lapply(cayear, st_drop_geometry)))
  cayear <- cayear[!(duplicated(cayear, by = "ca", fromLast = T))]
  coyear <- setValues(ca, cayear$idcoset[match(getValues(ca), cayear$ca)]) # 1793 - 2001 version
  writeRaster(coyear, sprintf("ehess_communes/co%s.tif", year), overwrite = T)
  # Pcaco
  pcaco <- data.table(ca = getValues(ca), coyear = getValues(coyear))
  pcaco <- na.omit(pcaco)
  pcaco <- setnames(pcaco, c("ca", str_c("co", year)))
  save.dta13(pcaco, sprintf("../shared_data/ehess/communes/pcaco%s.dta", year))
}

# lapply(c(1831), function(year) {
#   coyear <- raster(sprintf("ehess_communes/co%s.tif", year))
#   pcaco  <- data.table(ca = getValues(ca), coyear = getValues(coyear))
#   pcaco  <- na.omit(pcaco)
#   pcaco  <- setnames(pcaco, c("ca", str_c("co", year)))
#   save.dta13(pcaco, sprintf("../shared_data/ehess/communes/pcaco%s.dta", year))
# })