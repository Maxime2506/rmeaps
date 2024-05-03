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

md <- new("MeapsData", dist_triplet, actifs = actifs, emplois = emplois, fuite = fuites)

# p_dist <- aggregate(md@distances$i, by = list(md@distances$i), FUN = length)$x |> cumsum()
# p_dist <- c(0L, p_dist)
# 
# meaps_all_in(md@distances$j, p_dist, md@distances$value, md@emplois, md@actifs, md@fuite, parametres = 1.0, xr_odds = 1.0,
#              attraction = "constant", nthreads = 0, verbose = TRUE, normalisation = FALSE, fuite_min = 1e-3)
# 

all_in(md)
all_in(md, attraction = "logistique", parametres = c(1,1,.1))

les_communes <- c(1,1,1, 2,2,2, 3,3,3) |> as.character()
names(les_communes) <- froms
les_dclts <- c(rep(1L, 4), rep(2L,4), rep(3L, 4)) |> as.character()
names(les_dclts) <- tos

mobpro <- expand.grid(group_to = 1:3, group_from = 1:3) |> 
  mutate(value = 60/9)

mdg <- new("MeapsDataGroup", dist_triplet, group_from = les_communes, group_to = les_dclts, 
           actifs = actifs, emplois = emplois, fuite = fuites, cible = mobpro)


mdg2 <- meapsdatagroup(md, group_from = les_communes, group_to = les_dclts, cible = mobpro)


all_in_grouped(mdg)

meaps_optim(mdg, attraction = "logistique", parametres = c(1,1,.1), method = "L-BFGS-B")





# test faux objet
faux_triplet <- dist_triplet |> arrange(metric)
new("MeapsData")
new("MeapsData", faux_triplet, actifs = actifs, emplois = emplois, fuite = fuites) 
new("MeapsData", faux_triplet, actifs = actifs, emplois = emplois, fuite = fuites) |> check_meapsdata()

meapsdatagroup(new("MeapsData", faux_triplet, actifs = actifs, emplois = emplois, fuite = fuites),  group_from = les_communes, group_to = les_dclts, cible = mobpro) |> check_meapsdatagroup(
  
)


