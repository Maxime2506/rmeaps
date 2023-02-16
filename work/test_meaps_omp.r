library(tidyverse)

# Ligne 5 cases
cref = 0.6
cabs = rep(cref, 5)
pabs = cabs / (1 + cabs)
passe = cumprod(1-pabs)
reste = lag(passe, default = 1) * pabs

la_fuite = 0

meaps_single(
  rkdist = matrix(1:5, nrow = 1),
  emplois = rep(1,5),
  actifs = 1,
  modds = matrix(1, nrow = 1, ncol = 5),
  f = la_fuite,
  shuf = 1)

nb_shuf = 4
meaps_bootstrap2(rkdist = matrix(1:5, nrow = 1),
                emplois = rep(1,5),
                actifs = 1,
                modds = matrix(1, nrow = 1, ncol = 5),
                f = la_fuite,
                shuf = matrix(rep(1L, nb_shuf), ncol = 1))


# carré 2x2
meaps_single(
  rkdist = matrix(c(1:2, 2:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,2),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 2),
  f = rep(0, 2),
  shuf = 1:2)

# nb : p_abs théorique pour fuite_min = 1e-3 : 1 - sqrt(1e-3) = 0.9683772

meaps_bootstrap2(
  rkdist = matrix(c(1:2, 2:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,2),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 2),
  f = rep(0.1, 2),
  shuf = matrix(1:2, ncol = 2, byrow = TRUE))

meaps_bootstrap2(
  rkdist = matrix(c(1:2, 2:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,2),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 2),
  f = rep(0, 2),
  shuf = matrix(1:2, ncol = 2, nrow = 5, byrow = TRUE))

# carré 2x3
meaps_single(
  rkdist = matrix(c(1:3, 3:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,3),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 3),
  f = rep(0, 2),
  shuf = 1:2)

meaps_bootstrap2(
  rkdist = matrix(c(1:3, 3:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,3),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 3),
  f = rep(0, 2),
  shuf = matrix(1:2, ncol = 2, byrow = TRUE)) |> rowSums()




# deux lignes
meaps_single(
  rkdist = matrix(c(1:5, 5:1), nrow = 2, byrow = TRUE),
  emplois = rep(2,5),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 5),
  f = rep(la_fuite, 2),
  shuf = 1:2)

meaps_bootstrap2(
  rkdist = matrix(c(1:5, 5:1), nrow = 2, byrow = TRUE),
  emplois = rep(2, 5),
  actifs = c(1, 1),
  modds = matrix(1, nrow = 2, ncol = 5),
  f = rep(la_fuite, 2),
  shuf = matrix(c(1:2, 2:1), ncol = 2))


# Grille 4x4
library(sf)
library(matrixStats)
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
  mutate(dense = 1 / (st_distance(position, st_point(c(2.5, 2.5)), by_element = TRUE))^2,
         emplois = NB_emplois * dense / sum(dense)) |> 
  pull(emplois)

NB_actifs = 500
la_fuite = 0

marge_actifs <- tibble(position = residences) |> 
  mutate(dense = 1 / (st_distance(position, st_point(c(2.5, 2.5)), by_element = TRUE)),
         actifs = NB_actifs * dense / sum(dense)) |> 
  pull(actifs)

mat_odds <- matrix(1, nrow = 16, ncol = 4)


meaps_single(
  rkdist = rkdist,
  emplois = marge_emplois,
  actifs = marge_actifs,
  modds = mat_odds,
  f = rep(la_fuite, 16),
  shuf = 1:16)

shuf_mat = matrix(1:16, nrow = 1)
sm2 = rbind(shuf_mat, shuf_mat)
sm3 = rbind(sm2, shuf_mat)
sm12 = rbind(sm3, sm3, sm3, sm3)

meaps_bootstrap2(
  rkdist = rkdist,
  emplois = marge_emplois,
  actifs = marge_actifs,
  modds = mat_odds,
  f = rep(la_fuite, 16),
  shuf = shuf_mat)

meaps_bootstrap2(
  rkdist = rkdist,
  emplois = marge_emplois,
  actifs = marge_actifs,
  modds = mat_odds,
  f = rep(la_fuite, 16),
  shuf = sm2)

