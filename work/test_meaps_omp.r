library(tidyverse)

# Ligne 5 cases
cref = 0.6
cabs = rep(cref, 5)
pabs = cabs / (1 + cabs)
passe = cumprod(1-pabs)
reste = lag(passe, default = 1) * pabs

la_fuite = 0

meaps_oneshuf(
  rkdist = matrix(1:5, nrow = 1),
  emplois = rep(.3,5),
  actifs = 1,
  modds = matrix(1, nrow = 1, ncol = 5),
  f = la_fuite,
  shuf = 1,
  normalisation = FALSE)

nb_shuf = 4
meaps_multishuf(rkdist = matrix(1:5, nrow = 1),
                emplois = rep(1,5),
                actifs = 1,
                modds = matrix(1, nrow = 1, ncol = 5),
                f = la_fuite,
                shuf = matrix(rep(1L, nb_shuf), ncol = 1))


# carré 2x2
meaps_oneshuf(
  rkdist = matrix(c(1:2, 2:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,2),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 2),
  f = rep(0, 2),
  shuf = 1:2)

# nb : p_abs théorique pour fuite_min = 1e-3 : 1 - sqrt(1e-3) = 0.9683772


meaps_multishuf(
  rkdist = matrix(c(1:2, 2:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,2),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 2),
  f = rep(0.1, 2),
  shuf = matrix(1:2, ncol = 2, byrow = TRUE))

meaps_alt(
  rkdist = t(matrix(c(1:2, 2:1), nrow = 2, byrow = TRUE)),
  emplois = rep(1,2),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 2),
  f = rep(0.1, 2),
  shuf = t(matrix(1:2, ncol = 2, byrow = TRUE)))



meaps_multishuf(
  rkdist = matrix(c(1:2, 2:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,2),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 2),
  f = rep(0, 2),
  shuf = matrix(1:2, ncol = 2, nrow = 5, byrow = TRUE))

# carré 2x3
meaps_oneshuf(
  rkdist = matrix(c(1:3, 3:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,3),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 3),
  f = rep(0, 2),
  shuf = 1:2)

meaps_multishuf(
  rkdist = matrix(c(1:3, 3:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,3),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 3),
  f = rep(0, 2),
  shuf = matrix(1:2, ncol = 2, byrow = TRUE)) |> rowSums()




# deux lignes
meaps_oneshuf(
  rkdist = matrix(c(1:5, 5:1), nrow = 2, byrow = TRUE),
  emplois = rep(2,5),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 5),
  f = rep(la_fuite, 2),
  shuf = 1:2)

meaps_multishuf(
  rkdist = matrix(c(1:5, 5:1), nrow = 2, byrow = TRUE),
  emplois = rep(2, 5),
  actifs = c(1, 1),
  modds = matrix(1, nrow = 2, ncol = 5),
  f = rep(la_fuite, 2),
  shuf = matrix(c(1:2, 2:1), ncol = 2))


# Grille 4x4
library(sf)
library(matrixStats)
library(microbenchmark)

residences <- expand_grid(x=1:4, y=1:4) |> 
  as.matrix() |> 
  st_multipoint(dim = "XY") |> 
  st_sfc() |> 
  st_cast(to = "POINT")

emplois <- expand_grid(x=2:3, y=2:3) |> 
  as.matrix() |> 
  st_multipoint(dim = "XY") |> 
  st_sfc() |> 
  st_cast(to = "POINT")

distance <- st_distance(residences, emplois)
rkdist <- rowRanks(distance, ties.method = "random")

NB_emplois = 500
marge_emplois <- tibble(position = emplois) |> 
  mutate(dense = 1 / (st_distance(position, st_point(c(2.5, 2.5)), by_element = FALSE))^2,
         emplois = NB_emplois * dense / sum(dense)) |> 
  pull(emplois)

NB_actifs = 500
la_fuite = 0

marge_actifs <- tibble(position = residences) |> 
  mutate(dense = 1 / (st_distance(position, st_point(c(2.5, 2.5)), by_element = FALSE)),
         actifs = NB_actifs * dense / sum(dense)) |> 
  pull(actifs)

mat_odds <- matrix(1, nrow = 16, ncol = 4)


chances_absorption(
  rkdist = rkdist,
  emplois = marge_emplois,
  modds = mat_odds,
  f = rep(la_fuite, 16)) -> A

A / (A+1)


meaps_oneshuf(
  rkdist = rkdist,
  emplois = marge_emplois,
  actifs = marge_actifs,
  modds = mat_odds,
  f = rep(la_fuite, 16),
  shuf = 1:16) -> B

communaliser(B, rep(1:4, 4), c(1,1,2,2))

shuf_mat = matrix(1:16, nrow = 1)
sm2 = rbind(shuf_mat, shuf_mat)
sm3 = rbind(sm2, shuf_mat)
sm12 = rbind(sm3, sm3, sm3, sm3)

meaps_multishuf(
  rkdist = rkdist,
  emplois = marge_emplois,
  actifs = marge_actifs,
  modds = mat_odds,
  f = rep(la_fuite, 16),
  shuf = shuf_mat)

meaps_multishuf(
  rkdist = rkdist,
  emplois = marge_emplois,
  actifs = marge_actifs,
  modds = mat_odds,
  f = rep(la_fuite, 16),
  shuf = sm2)


meaps_tension(
  rkdist = rkdist,
  emplois = marge_emplois,
  actifs = marge_actifs,
  modds = mat_odds,
  f = rep(la_fuite, 16),
  shuf = sm2)

altrk = t(rkdist)
altsm2 = t(sm2)

### zz <- microbenchmark("old" = meaps_multishuf())

# gros test
# 
library(tidyverse)
library(sf)
library(matrixStats)
library(rmeaps)
maxx <- 5
maxy <- 5
n <- 3
k <- 3
NB_actifs  <- 2
NB_emplois = 2
la_fuite <-  0.1

residences <- expand_grid(x=seq(0, maxx, length.out=n), y=seq(0, maxy, length.out=n)) |> 
  as.matrix() |> 
  st_multipoint(dim = "XY") |> 
  st_sfc() |> 
  st_cast(to = "POINT")

emplois <- expand_grid(x=seq(0, maxx, length.out=k), y=seq(0, maxy, length.out=k)) |> 
  as.matrix() |> 
  st_multipoint(dim = "XY") |> 
  st_sfc() |> 
  st_cast(to = "POINT")

distance <- st_distance(residences, emplois)
rkdist <- rowRanks(distance, ties.method = "random")

marge_emplois <- tibble(position = emplois) |> 
  mutate(dense = 1 / (st_distance(position, st_point(c(2.5, 2.5)), by_element = FALSE))^2,
         emplois = NB_emplois * dense / sum(dense)) |> 
  pull(emplois)

marge_actifs <- tibble(position = residences) |> 
  mutate(dense = 1 / (1+st_distance(position, st_point(c(2.5, 2.5)), by_element = FALSE)),
         actifs = NB_actifs * dense / sum(dense)) |> 
  pull(actifs)

mat_odds <- matrix(1, nrow = nrow(marge_actifs), ncol = nrow(marge_emplois))

shuf <- map(1:64, ~sample.int(nrow(marge_actifs), nrow(marge_actifs)))
shuf <- do.call(rbind, shuf)
tens <- meaps_oneshuf(
  rkdist = rkdist,
  emplois = marge_emplois,
  actifs = marge_actifs,
  modds = mat_odds,
  f = rep(la_fuite, nrow(marge_actifs)),
  shuf = shuf[1, , drop=FALSE])
sum(tens)
