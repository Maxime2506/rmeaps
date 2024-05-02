library(tidyverse)
library(sf)
library(Matrix)
library(rmeaps)

#----- Définition des points -----
actifs_points <- tribble(
  ~fromidINS, ~x, ~y,
  1, -10, -1,
  2, -10, 0,
  3, -10, 1,
  4, 0, 1,
  5, 0, 0,
  6, 5, 1,
  7, 10, 0,
  8, 10, -1,
  9, 10, 1
) |> st_as_sf(coords = c("x", "y"))

emplois_points <- tribble(
  ~toidINS, ~x, ~y,
  1, -10, 0,
  2, -5, 0,
  3, 0, 0,
  4, 0, 1,
  5, 0, 2,
  6, 0, -4,
  7, 10, 0,
  8, 10, 2,
  9, -15, 5,
  10, 15, -5,
  11, 0, 8,
  12, 0, -8
  ) |> st_as_sf(coords = c("x", "y")) 

dist <- st_distance(actifs_points, emplois_points) 

# Introduction de NA sur des distances longues.
dist[dist > 22] <- NA
dist[dist == 0] <- .5

temp <- as(as(as(dist, "dMatrix"), "generalMatrix"), "TsparseMatrix")
dist_triplet <- tibble(fromidINS = as.character(temp@i), toidINS = as.character(temp@j), dist = temp@x) |> 
  drop_na() |> 
  arrange(fromidINS, dist, toidINS)

froms <- dist_triplet$fromidINS |> unique() |> sort()
tos <- dist_triplet$toidINS |> unique() |> sort()

actifs <- c(5, 5 , 5 , 10, 25, 10, 5, 4, 2)
emplois <- c(2, 2, 2, 8, 10, 15, 6, 4, 2, 4, 2, 3)
names(actifs) <- froms
names(emplois) <- tos

f <- 1 - sum(emplois) / sum(actifs)
fuites <- rep(f, 9)
names(fuites) <- froms

dist_triplet <- dist_triplet |>
  left_join(tibble(i=seq_along(froms)-1L, fromidINS = froms), by = "fromidINS") |> 
  left_join(tibble(j=seq_along(tos)-1L, toidINS = tos), by = "toidINS") |> 
  select(i, j, value = dist)
md <- new("MeapsData", dist_triplet,fromidINS = froms, toidINS = tos, actifs = actifs, emplois = emplois, fuite = fuites)

# p_dist <- aggregate(md@distances$i, by = list(md@distances$i), FUN = length)$x |> cumsum()
# p_dist <- c(0L, p_dist)
# 
# meaps_all_in(md@distances$j, p_dist, md@distances$value, md@emplois, md@actifs, md@fuite, parametres = 1.0, xr_odds = 1.0,
#              attraction = "constant", nthreads = 0, verbose = TRUE, normalisation = FALSE, fuite_min = 1e-3)
# 

all_in(md)
all_in(md, attraction = "logistique", parametres = c(1,1,.1))

les_communes <- c(1,1,1, 2,2,2, 3,3,3) |> as.character()
les_dclts <- c(rep(1L, 4), rep(2L,4), rep(3L, 4)) |> as.character()

mobpro <- expand.grid(dclt = 1:3, commune = 1:3) |> 
  mutate(value = 60/9)

mdg <- new("MeapsDataGroup", dist_triplet,fromidINS = froms, toidINS = tos, group_from = les_communes, group_to = les_dclts, 
           actifs = actifs, emplois = emplois, fuite = fuites, cible = mobpro)


all_in_grouped(mdg)


object = mdg
arg <- list(attraction = "logistique", parametres = c(1,1,.1))


meaps_optim(mdg, attraction = "logistique", parametres = c(1,1,.1), method = "L-BFGS-B")



shuf <- matrix(NA, ncol = 9, nrow = 64)
for (i in 1:64) shuf[i, ] <- sample(1:9)
colnames(shuf) <- actifs_points$fromidINS - 1L
names(marges_actifs) <- actifs_points$fromidINS -1L
#----

res_another <- another_meaps(dist_triplet, marges_emplois, marges_actifs, fuites, attraction = "constant", nthreads = 1L)
res_continu <- meaps_continu(dist_triplet, marges_emplois, marges_actifs, fuites, shuf, attraction = "constant")


res_another |> summarise(flux = sum(flux))
sum(marges_emplois)
sum(marges_actifs)
sum(marges_actifs * (1-fuites))

res_another |> group_by(fromidINS) |> summarise(flux = sum(flux)) |> mutate(actifs_occ = marges_actifs * (1-fuites) ) 
res_another |> group_by(toidINS) |> summarise(flux = sum(flux)) |> mutate(emplois = marges_emplois)


m_another <- sparseMatrix(i = res_another$fromidINS + 1L, j = res_another$toidINS + 1L, x = res_another$flux)
m_another |> as.matrix() |> rowSums()
m_another |> as.matrix() |> colSums()

another_meaps(dist_triplet, marges_emplois, marges_actifs, fuites, attraction = "logistique", param = c(1, 1, .1)) |> summarise(sum(flux))

