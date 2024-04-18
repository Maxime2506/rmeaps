library(rmeaps)
library(Matrix)
library(matrixStats)
library(tidyverse)

actifs <- c(11, 11, 11, 11, 11)
fuite <- c(0.1, 0.1, 0.1, 0.1, 0.1)
names(actifs) <- str_c("a", 1:5)
names(fuite) <- names(actifs)
emplois <- c(5, 5, 5, 5, 30)
names(emplois) <- str_c("e", 1:5)
distances <- matrix(
  c(0, 1, 2, 3, NA,
    1, 0, 1, NA, NA,
    2, NA, 0, 1, 2,
    3, 2, NA, 0, 1,
    4, 3, 2, NA, 0),
  nrow = 5, dimnames = list(names(actifs), names(emplois)))
triplet <- distances |>
  as_tibble(rownames = "fromidINS") |>
  pivot_longer(cols = -fromidINS, names_to = "toidINS", values_to = "distance") |>
  filter(!is.na(distance))

N <- length(actifs)
K <- length(emplois)
shuf <- emiette(actifs, 16)
noshuf <- matrix(1:5, nrow =1, ncol=N, dimnames = list(NULL, names(actifs)))
modds <- matrix(1, ncol = K, nrow = N)
dimnames(modds) <- dimnames(distances)

c0 <- meaps_multishuf(
  rkdist = rowRanks(distances, ties="random"),
  emplois = emplois,
  actifs = actifs,
  modds = modds,
  f = fuite, 
  shuf = noshuf, 
  progress = FALSE)

c1 <- meaps_continu(
  dist = triplet,
  emplois = emplois,
  actifs = actifs,
  f = fuite, 
  shuf = noshuf, 
  progress = FALSE) |> 
  as_tibble() |> 
  arrange(desc(flux))

prep1 <- prep_meaps_dist(dist = triplet, emplois = emplois, actifs = actifs, fuite = fuite, shuf = noshuf, 
                         groups_from = set_names(1:5, names(actifs)),
                         groups_to = set_names(1:5, names(emplois)))
prep2 <- prep_meaps_dist(dist = triplet, emplois = emplois, actifs = actifs, fuite = fuite, shuf = noshuf, 
                         groups_from = set_names(c(1,1,2,2,2), names(actifs)),
                         groups_to = set_names(1:5, names(emplois)))
prep3 <- prep_meaps_dist(dist = triplet, emplois = emplois, actifs = actifs, fuite = fuite, shuf = noshuf, 
                         groups_from = set_names(1:5, names(actifs)),
                         groups_to = set_names(c(1,1,1,2,2), names(emplois)))

o1 <- meaps_optim(prep1)
o2 <- meaps_optim(prep2)
o3 <- meaps_optim(prep3)
sum(o1$flux)
sum(o2$flux)
sum(o3$flux)
sum(c1$flux)  
sum(c0)

# marche ----------------
modds[distances<=2] <- 10
modds[distances>2] <- 1
m0 <- meaps_multishuf(
  rkdist = rowRanks(distances, ties="random"),
  emplois = emplois,
  actifs = actifs,
  modds = modds,
  f = fuite, 
  shuf = noshuf, 
  progress = FALSE) |> 
  as_tibble(rownames = "COMMUNE") |> 
  pivot_longer(cols = -COMMUNE, names_to = "DCLT", values_to = "flux") |> 
  mutate(DCLT = str_sub(DCLT,2, -1)) |> 
  arrange(desc(flux)) |> 
  filter(flux>0)

m1 <- meaps_continu(
  dist = triplet,
  emplois = emplois,
  actifs = actifs,
  attraction = "marche",
  param = c(2, 10),
  f = fuite, 
  shuf = noshuf, 
  progress = FALSE) |> 
  as_tibble() |> 
  arrange(desc(flux)) |> 
  rename(COMMUNE = fromidINS, DCLT = toidINS) |> 
  filter(flux>0) |> 
  mutate(DCLT = str_sub(DCLT, 2, -1),
         COMMUNE = str_sub(COMMUNE, 2, -1))

m2 <- meaps_optim(prep1, attraction = "marche", param=c(2, 10))
odds <- triplet |> 
  mutate(odd = ifelse(distance <=2, log(10), 0)) |> 
  select(fromidINS, toidINS, odd)
prep_o3 <- prep_meaps_odds_on_dist(odds, prep1)
m3 <- meaps_optim(prep1, odds_prep = prep_o3)
m4 <- meaps_continu(
  dist = triplet,
  emplois = emplois,
  actifs = actifs,
  attraction = "odds",
  modds = odds,
  f = fuite, 
  shuf = noshuf, 
  progress = FALSE) |> 
  as_tibble() |> 
  arrange(desc(flux)) |> 
  rename(COMMUNE = fromidINS, DCLT = toidINS) |> 
  filter(flux>0) |> 
  mutate(DCLT = str_sub(DCLT, 2, -1),
         COMMUNE = str_sub(COMMUNE, 2, -1))

m0 |>
  left_join(m1 |> rename(flux.1 = flux), by=c("COMMUNE", "DCLT")) |> 
  left_join(m2 |> rename(flux.2 = flux), by=c("COMMUNE", "DCLT")) |> 
  left_join(m3 |> rename(flux.3 = flux), by=c("COMMUNE", "DCLT")) |>
  left_join(m4 |> rename(flux.4 = flux), by=c("COMMUNE", "DCLT")) |> print(n=100)
  
# matrice de odds

todds <- triplet |> 
  mutate(odd = ifelse(distance <=2, log(2), 0))
podds <- prep_meaps_odds(todds |> select(fromidINS, toidINS, odd), prep1$cle_from, prep1$cle_to)

o4 <- meaps_optim(prep1, odds_prep = podds)
