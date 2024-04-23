library(tidyverse)
library(sf)
library(Matrix)

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
#dist[dist > 22] <- NA
dist[dist == 0] <- .5

temp <- as(as(as(dist, "dMatrix"), "generalMatrix"), "TsparseMatrix")
dist_triplet <- tibble(fromidINS = temp@i, toidINS = temp@j, dist = temp@x) |> 
  drop_na()


marges_actifs <- c(5, 5 , 5 , 10, 25, 10, 5, 4, 2)
marges_emplois <- c(2, 2, 2, 8, 10, 15, 6, 4, 2, 4, 2, 3)
f <- 1 - sum(marges_emplois) / sum(marges_actifs)
fuites <- rep(f, 9)

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

