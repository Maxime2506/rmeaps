library(data.table)
library(matrixStats)
library(r3035)
library(tidyverse)
library(ofce)
library(rmeaps)
library(Matrix)
library(tictoc)

setwd("~/larochelle2")

load("baselayer.rda")

load("output/data_LR.rda")
load(file= "output/dist.com.srda")

nb_tirages  <-  64
mat_names <- dimnames(mat_rang_na)
mat_names <- list(r3035::contract_idINS(mat_names[[1]]), r3035::contract_idINS(mat_names[[2]]))
dimnames(mat_rang_na) <- dimnames(mat_distance_na) <- mat_names

time_triplet <- mat_distance_na |> 
  as_tibble() |> 
  mutate(fromidINS = rownames(mat_distance_na)) |>
  pivot_longer(-fromidINS, names_to = "toidINS", values_to = "temps") |> 
  drop_na(temps) 

froms <- distinct(time_triplet, fromidINS) |> pull(fromidINS) |> as.character()
tos <- distinct(time_triplet, toidINS) |> pull(toidINS) |> as.character()

c200 <- qs::qread(c200ze_file)
cc <- c200 |> mutate(idINS = r3035::contract_idINS(idINS)) |> pull(com21, name = idINS) 
coms <- cc[froms]
dclts <- cc[tos]

emplois  <-  les_emplois |> mutate(idINS = r3035::contract_idINS(idINS)) |>  pull(emplois_reserves_c, name = idINS)
emplois <- emplois[tos]
actifs <-  les_actifs |> mutate(idINS = r3035::contract_idINS(idINS)) |>  pull(act_18_64c, name = idINS)
actifs <- actifs[froms]
fuite <-  les_actifs |> mutate(idINS = r3035::contract_idINS(idINS)) |>  pull(fuite, name = idINS)
fuite <- fuite[froms]
N <- length(froms)
K <- length(tos)

noshuf <- matrix(1:N, ncol = N, nrow = 1, dimnames = list(NULL, froms))
shufs <- emiette(les_actifs = actifs, nshuf = 64, seuil = 200) # ordre correct !

pp <- list(
  rkdist = matrixStats::rowRanks(mat_distance_na[froms, tos], ties = "random"), 
  emplois = emplois, 
  actifs = actifs, 
  f = fuite, 
  shuf = noshuf, 
  nthreads = 16)

odds <- tibble(fromidINS = time_triplet$fromidINS, toidINS = time_triplet$toidINS, odd = 0)

# multishuf
modds <- matrix(1, ncol = K, nrow = N)
dimnames(modds) <- dimnames(pp$rkdist)

tic();ref_ms <- rmeaps::meaps_multishuf(
  rkdist = pp$rkdist,
  emplois = emplois,
  actifs = actifs,
  f = fuite,
  modds = modds, 
  shuf = noshuf) |> 
  communaliser(as.integer(coms), as.integer(dclts)) |> 
  as_tibble(rownames= "COMMUNE") |>
  pivot_longer(-COMMUNE, names_to = "DCLT", values_to = "flux") |> 
  filter(flux > 0) |>
  arrange(desc(flux)); toc()

tic();ref_ms_shufs <- rmeaps::meaps_multishuf(
  rkdist = pp$rkdist,
  emplois = emplois,
  actifs = actifs,
  f = fuite,
  modds = modds, 
  shuf = shufs) |> 
  communaliser(as.integer(coms), as.integer(dclts)) |> 
  as_tibble(rownames= "COMMUNE") |>
  pivot_longer(-COMMUNE, names_to = "DCLT", values_to = "flux") |>
  filter(flux > 0) |>
  arrange(desc(flux)); toc()

tic();MOD2 <- rmeaps::meaps_continu(
  dist = time_triplet,
  emplois = emplois,
  actifs = actifs,
  f = fuite,
  shuf = noshuf,
  nthreads = 16L) |>
  mutate(COMMUNE = coms[fromidINS], DCLT = dclts[toidINS]) |> 
  group_by(COMMUNE, DCLT) |>
  summarize(flux = sum(flux), .groups = "drop") |>
  filter(flux > 0) |>
  arrange(desc(flux));toc()

tic();MOD3 <- rmeaps::meaps_continu(
  dist = time_triplet,
  emplois = emplois,
  actifs = actifs,
  f = fuite,
  shuf = shufs,
  nthreads = 16L) |>
  mutate(COMMUNE = coms[fromidINS], DCLT = dclts[toidINS]) |> 
  group_by(COMMUNE, DCLT) |>
  summarize(flux = sum(flux), .groups = "drop") |>
  filter(flux > 0) |>
  arrange(desc(flux));toc()


dpns <- rmeaps::prep_meaps_dist(time_triplet, emplois, actifs, fuite, noshuf, coms, dclts)
dps <- rmeaps::prep_meaps_dist(time_triplet, emplois, actifs, fuite, shufs, coms, dclts)

tic();modons <- rmeaps::meaps_optim(prep = dpns);toc()
tic();modos <- rmeaps::meaps_optim(prep = dps);toc()
tic();modons <- rmeaps::meaps_optim(prep = dpns);toc()
tic();modos <- rmeaps::meaps_optim(prep = dps);toc()
tic();modosm <- rmeaps::meaps_optim(prep = dps, attraction="marche", param = c(10,1));toc()
tic();modonsm <- rmeaps::meaps_optim(prep = dpns, attraction="marche", param = c(10,1));toc()

fluxs <- ref_ms |> rename(flux.ms.ns = flux) |> 
  full_join(MOD2 |> rename(flux.c.ns = flux), by = c("COMMUNE", "DCLT")) |> 
  full_join(modons |> rename(flux.o.ns = flux), by = c("COMMUNE", "DCLT")) |>
  full_join(modonsm |> rename(flux.o.ns.m = flux), by = c("COMMUNE", "DCLT")) |>
  full_join(ref_ms_shufs |> rename(flux.ms.s= flux), by = c("COMMUNE", "DCLT")) |>
  full_join(MOD3 |> rename(flux.c.s = flux), by = c("COMMUNE", "DCLT")) |> 
  full_join(modos |> rename(flux.o.s = flux), by = c("COMMUNE", "DCLT")) |>
  full_join(modosm |> rename(flux.o.s.m = flux), by = c("COMMUNE", "DCLT")) |>
  full_join(mobpro, by = c("COMMUNE", "DCLT"))

fluxs |> summarize(across(starts_with("flux"), ~sum(.x, na.rm=TRUE)))

bootstrap <- map_dfr(0:48/4, ~{
  shuf <- emiette(les_actifs = actifs, nshuf = round(2^.x), seuil = 200) 
  tic(); 
  MOD <- rmeaps::meaps_continu(
    dist = time_triplet,
    emplois = emplois,
    actifs = actifs,
    f = fuite,
    shuf = shuf,
    nthreads = 16L, progress = FALSE, quiet = TRUE) |>
    mutate(COMMUNE = coms[fromidINS], DCLT = dclts[toidINS]) |> 
    group_by(COMMUNE, DCLT) |>
    summarize(flux = sum(flux), .groups = "drop") |>
    filter(flux > 0) |>
    arrange(desc(flux))
  time <- toc(quiet = TRUE)
  time <- time$toc - time$tic
  MOD |> mutate(nshuf = round(2^.x), time = as.numeric(time))
}, .progress=TRUE)
