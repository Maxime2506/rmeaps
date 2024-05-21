devtools::load_all()
library(tidyverse)

load("~/Documents/LaRochelle.RData")
load("~/Documents/LaRochelleGroup.RData")


# Sélection des DLCT de cible
dclt <- cible |> distinct(group_to) |> pull()
group_to <- group_to[group_to %in% dclt]
emplois <- emplois[names(group_to)]
dist <- dist |> filter(toidINS %in% names(group_to))

dcltbis <- intersect(cible$group_to, group_to)
cible <- cible |> filter(group_to %in% dcltbis)
dist <- dist |> filter(toidINS %in% names(group_to))


LR <- new("MeapsData", dist, actifs = actifs, emplois = emplois, fuite = fuites)
LRG <- meapsdatagroup(LR, group_from, group_to, cible)