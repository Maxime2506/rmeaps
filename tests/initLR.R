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
cible <- cible |> filter(group_to %in% dcltbis) |> as.data.frame()
dist <- dist |> filter(toidINS %in% names(group_to)) |> as.data.frame()


LR <- new("MeapsData", dist, actifs = actifs, emplois = emplois, fuite = fuites)
LRG <- meapsdatagroup(LR, group_from, group_to, cible)


#cat("all_in 1 thread")
#all_in(LR, nthreads = 1)

#cat("all_in multithread")
#all_in(LR)

#cat("all_in grouped multi thread")
#all_in_grouped(LR)

multishuf(LR, nshuf = 8)
