library(tidyverse)

# carré 2x2
# nb : p_abs théorique pour fuite_min = 1e-3 : 1 - sqrt(1e-3) = 0.9683772
meaps_oneshuf(
  rkdist = matrix(c(1:2, 2:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,2),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 2),
  f = rep(1e-3, 2),
  shuf = 1:2)

meaps_multishuf(
  rkdist = matrix(c(1:2, 2:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,2),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 2),
  f = rep(0.1, 2),
  shuf = matrix(1:2, ncol = 2, byrow = TRUE),
  normalisation=TRUE)

meaps_alt(
  rkdist = matrix(c(1:2, 2:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,2),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 2),
  f = rep(0.1, 2),
  shuf = matrix(1:2, ncol = 2, byrow = TRUE),
  normalisation = TRUE)

# carré 2x3
meaps_oneshuf(
  rkdist = matrix(c(1:3, 3:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,3),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 3),
  f = rep(1e-3, 2),
  shuf = 1:2)

meaps_multishuf(
  rkdist = matrix(c(1:3, 3:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,3),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 3),
  f = rep(0, 2),
  shuf = matrix(1:2, ncol = 2, byrow = TRUE),
  normalisation=TRUE) 

meaps_alt(
  rkdist = matrix(c(1:3, 3:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,3),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 3),
  f = rep(0, 2),
  shuf = matrix(1:2, ncol = 2, byrow = TRUE),
  normalisation=TRUE, progress=FALSE) 





meaps_alt(
  rkdist = matrix(c(1:3, 3:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,3),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 3),
  mode ="discret",
  f = rep(0, 2),
  shuf = matrix(1:2, ncol = 2, byrow = TRUE)) 

meaps_alt(
  rkdist = matrix(c(1:3, 3:1), nrow = 2, byrow = TRUE),
  emplois = rep(1,3),
  actifs = c(1,1),
  modds = matrix(1, nrow = 2, ncol = 3),
  mode ="subjectif_c",
  oddssubjectifs = c(1,1,1),
  f = rep(0, 2),
  shuf = matrix(1:2, ncol = 2, byrow = TRUE)) 



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

# gros test ----------------------
# 
library(tidyverse)
library(sf)
library(matrixStats)
maxx <- 5
maxy <- 5
n <- 30
k <- 30
NB_actifs  <- 20000
NB_emplois <- 18000
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
  mutate(dense = 1 / (1+st_distance(position, st_point(c(2.5, 2.5)), by_element = FALSE))^2,
         emplois = NB_emplois * dense / sum(dense)) |> 
  pull(emplois)

marge_actifs <- tibble(position = residences) |> 
  mutate(dense = 1 / (1+st_distance(position, st_point(c(2.5, 2.5)), by_element = FALSE)),
         actifs = NB_actifs * dense / sum(dense)) |> 
  pull(actifs)

mat_odds <- matrix(1, nrow = nrow(marge_actifs), ncol = nrow(marge_emplois))


shuf <- map(1:128, ~sample.int(nrow(marge_actifs), nrow(marge_actifs)))
shuf <- do.call(rbind, shuf)
tens <- meaps_tension(
  rkdist = rkdist,
  emplois = marge_emplois,
  actifs = marge_actifs,
  modds = mat_odds,
  f = rep(la_fuite, nrow(marge_actifs)),
  shuf = shuf, seuil_dispo = 0.5)
sum(tens$flux)
tens$tension
ggplot(tibble(r = tens$tension))+geom_histogram(aes(x=r), bins=100)


# la version alt prend en input la transposée pour les matrices.
# t_rkdist <- t(rkdist)
# t_mat_odds <- t(mat_odds)
# t_shuf <- t(shuf)


bench::mark(
  multishuf = meaps_multishuf(rkdist, marge_emplois, marge_actifs, mat_odds, rep(la_fuite, nrow(marge_actifs)), shuf),
  alt = meaps_alt(rkdist, marge_emplois, marge_actifs, mat_odds, rep(la_fuite, nrow(marge_actifs)), shuf),
  )
# comparaison 
quantile(abs(multishuf - alt)/multishuf, probs = c(.8, .9, .95, .99, .999, .9999, 1))





microbenchmark::microbenchmark(
  multishuf = meaps_multishuf(rkdist, marge_emplois, marge_actifs, mat_odds, rep(la_fuite, nrow(marge_actifs)), shuf),
  alt = meaps_alt(rkdist, marge_emplois, marge_actifs, mat_odds, rep(la_fuite, nrow(marge_actifs)), shuf),
  times = 2
)



meaps_alt(t_rkdist, marge_emplois, marge_actifs, t_mat_odds, rep(la_fuite, nrow(marge_actifs)), t_shuf,
          mode = "continu")

meaps_alt(t_rkdist, marge_emplois, marge_actifs, t_mat_odds, rep(la_fuite, nrow(marge_actifs)), t_shuf,
          mode = "discret")


odsubj <- (900:1)/mean(900:1)

meaps_alt(t_rkdist, marge_emplois, marge_actifs, t_mat_odds, rep(la_fuite, nrow(marge_actifs)), t_shuf,
          mode = "subjectif", oddssubjectifs = odsubj)


