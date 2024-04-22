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
dist[dist > 22] <- NA
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

#----

res_another <- another_meaps(dist_triplet, marges_emplois, marges_actifs, fuites, attraction = "constant")
res_continu <- meaps_continu(dist_triplet, marges_emplois, marges_actifs, fuites, shuf, attraction = "constant")




x <- c(1.2, 2.9, 2.0) ; xr <- c(1.2, 2.0, 2.9) 
j <- c(1L, 0L, 2L) ; jr <- c(1L, 2L, 0L)
p <- c(0L, 1L, 3L)
dim <- c(2L,3L)

mat <- new(RankedRSMatrix, xr, jr, p, dim) # constructeur sur des données rangées.
dg <- mat$unrank()


tripl <- data.frame(fromidINS = c("a", "b", "b"), toidINS = c("1", "2", "3"), dist = xr)

odds <- matrix(0, nrow = 2, ncol = 3)
odds[1, 2] <- 1
odds[2,1] <- 3
odds[2, 3] <- 5
triodds <- data.frame(fromidINS = c("a", "b", "b"), toidINS = c("1", "2", "3"), dist = c(0,2,2))


modds$rankby(mat)

meaps_continu(tripl, emplois = c(4,9,5), actifs = c(10, 10), f = c(.1, .1), shuf = matrix(1:2, nrow = 1), attraction = "odds", 
            modds = triodds)

gfrom <- 1:2
gto <- c(1L,1L,2L)

meaps_optim(mat, emplois = c(4,9,5), actifs = c(10, 10), fuite = c(.1, .1), shuf = matrix(1:2, nrow = 1), groups_from = gfrom, groups_to = gto, attraction = "odds", 
            odds_prep = list(attraction = "odds", lodds = mat))



matrix(0:5, nrow = 2) ->z
z
as(as(as(z, "dMatrix"), "generalMatrix"), "RsparseMatrix") -> w

new(RankedRSMatrix, w) -> ww 
ww$at(0,0)
ww$unrank()
ww$rankby(w)

