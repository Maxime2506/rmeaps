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
dist_triplet <- tibble(fromidINS = as.character(temp@i), toidINS = as.character(temp@j), metric = temp@x) |> 
  drop_na() |> 
  arrange(fromidINS, metric, toidINS)

froms <- dist_triplet$fromidINS |> unique() |> sort()
tos <- dist_triplet$toidINS |> unique() |> sort()

actifs <- c(5, 5 , 5 , 10, 25, 10, 5, 4, 2)
emplois <- c(2, 2, 2, 8, 10, 15, 6, 4, 2, 4, 2, 3)
names(actifs) <- froms
names(emplois) <- tos

f <- 1 - sum(emplois) / sum(actifs)
fuites <- rep(f, 9)
names(fuites) <- froms

MeapsData <- new("MeapsData", dist_triplet, actifs = actifs, emplois = emplois, fuite = fuites)

froms <- MeapsData@froms
tos <- MeapsData@tos

jlab <- seq_along(tos) - 1L
names(jlab) <- tos

les_j <- jlab[MeapsData@triplet$toidINS]

p_dist <- MeapsData@triplet |> group_by(fromidINS) |> summarize(n()) |> pull() |> cumsum()
p_dist <- c(0L, p_dist)

arg <- list(jr_dist = les_j,
             p_dist = p_dist,
             xr_dist = MeapsData@triplet$metric,
             emplois = MeapsData@emplois,
             actifs = MeapsData@actifs,
             fuites = MeapsData@fuites,
             parametres = 1.0,
             xr_odds = 1.0,
            nthreads = 1L,
            verbose = TRUE)





