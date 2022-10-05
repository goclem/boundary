# !/usr/bin/env Rscript
# Usage: Commutes communes from EHESS records for the period 1861 - 2006
# Author: Clement Gorin
# Contact: gorinclem@gmail.com
# Date: June 2021

pacman::p_load(data.table, gtools, pbmcapply, readstata13, reshape2, sf, stringr, dplyr)
setwd("~/Dropbox/research/arthisto/data_ehess")

duplicates <- function(vec) duplicated(vec) | duplicated(vec, fromLast = T)

years  <- c(1793, 1800, 1806, 1821, 1831, 1836, 1841, 1846, 1851, 1856, 1861, 1866, 1872, 1876, 1881, 1886, 1891, 1896, 1901, 1906, 1911, 1921, 1926, 1931, 1936, 1946, 1954, 1962, 1968, 1975, 1982, 1990, 1996, 1999, 2001, 2006)

pidco <- lapply(years, function(year) {
  idcoyear <- st_read(sprintf("ehess_communes/co%s.gpkg", year), quiet = T)
  idcoyear <- unique(data.table(st_drop_geometry(idcoyear))[, c("idehessset", "idcoset")])
  idcoyear <- idcoyear[, .(idehess = as.numeric(unlist(str_split(idehessset, " \\| ")))), by = idcoset]
  # Fixes forced sets (we forced the communes but not the idehess, which creates duplicates)
  idcoyear <- idcoyear[!(duplicates(idehess) & idcoset %in% c(75900, 75056))] 
  idcoyear <- idcoyear[!(duplicates(idehess) & idcoset < 99000)]
  idcoyear <- idcoyear[!(duplicates(idehess) & idehess == 27686 & idcoset == 991313)]
  print(ifelse(any(duplicates(idcoyear$idehess)), "duplicates!", "All ok")) # Check
  idcoyear <- idcoyear[order(idehess), c("idehess", "idcoset")]
  idcoyear <- setnames(idcoyear, c("idehess", str_c("co", year)))
  return(idcoyear)
})

pidco <- Reduce(function(x, y) merge(x, y, by = "idehess", all = T), pidco)
pidco[idehess == 8509] <- as.data.table(as.list(na.locf(unlist(pidco[idehess == 8509]))))

save.dta13(pidco, "/Users/clementgorin/Dropbox/research/arthisto/shared_data/ehess/communes/pidco.dta")

pidco[idehess == 25298]


idcosets <- readRDS("/Users/clementgorin/Dropbox/research/arthisto/data_ehess/ehess_communes/cosetyrs_idcosets.RDS")
idcosets[idcoset == 999305]

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
years  <- c(1793, 1800, 1806, 1821, 1831, 1836, 1841, 1846, 1851, 1856, 1861, 1866, 1872, 1876, 1881, 1886, 1891, 1896, 1901, 1906, 1911, 1921, 1926, 1931, 1936, 1946, 1954, 1962, 1968, 1975, 1982, 1990, 1996, 1999, 2001, 2006)

for(year in years) {
  print(year)
  coyear <- st_read(sprintf("ehess_communes/co%s.gpkg", year), quiet = T)
  coyear <- dplyr::select(coyear, idcoset) %>% mutate(idcoset = as.numeric(idcoset))
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
  coyear <- setValues(ca, cayear$idcoset[match(getValues(ca), cayear$ca)])
  writeRaster(coyear, sprintf("ehess_communes/co%s.tif", year), overwrite = T)
  # Pcaco
  pcaco <- data.table(ca = getValues(ca), coyear = getValues(coyear))
  pcaco <- na.omit(pcaco)
  pcaco <- setnames(pcaco, c("ca", str_c("co", year)))
  save.dta13(pcaco, sprintf("../shared_data/ehess/communes/pcaco%s.dta", year))
}

# Computes pcaco ----------------------------------------------------------

pacman::p_load(raster, data.table, readstata13)
setwd("~/Dropbox/research/arthisto/data_ehess")

geom1793 <- raster("ehess_communes/geom1793.tif")
ca     <- raster("../data_project/ca.tif")

pcageom1793 <- data.table(ca = getValues(ca), geom1793 = getValues(geom1793))
pcageom1793 <- na.omit(pcageom1793)
save.dta13(pcageom1793, "../shared_data/ehess/communes/pcageom1793.dta")

co2015  <- raster("../data_project/co15.tif")
ca      <- raster("../data_project/ca.tif")
pop2015 <- raster("../data_project/pop2015s.tif")

pcaco2015 <- data.table(ca = getValues(ca), co2015 = getValues(co2015), pop2015s = getValues(pop2015))
pcaco2015 <- na.omit(pcaco2015)
save.dta13(pcaco2015, "../shared_data/ehess/communes/pcaco2015.dta")
