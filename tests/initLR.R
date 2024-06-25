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

gfrom <- tibble(fromidINS = names(group_from), gf = group_from) |>
  mutate(gf = factor(gf) |> as.integer())

gto <- tibble(toidINS = names(group_to), gt = group_to) |>
  mutate(gt = factor(gt) |> as.integer())

strug <- dist |> left_join(gfrom, by = "fromidINS") |> left_join(gto, by = "toidINS")
strug <- strug |> group_by(gf, gt) |> mutate(cc = cur_group_id()-1) |> ungroup()
cible <- cible |> mutate(gf = factor(group_from) |> as.integer(), gt = factor(group_to) |> as.integer()) 
v1 <- strug |> group_by(gf, gt) |> summarise(.groups = "drop") |> left_join(cible, by = c("gf", "gt"))   |> mutate(value = replace_na(value, 0)) |>pull(value)


strug <- strug |> ungroup() |> transmute(fromidINS, toidINS, cc, qq = cut_number(metric, n = 10, labels = FALSE)-1)


v2 <- sum(cible$value)/10 |> rep(x=_, 10)


LRMC <- meapsdatamultichamps(LR, strug, cibles = list(v1, v2))

all_in_multichamps(MeapsDataMultiChamps = LRMC)

ponderer_champs(LRMC)

meaps_multiopt(LRMC, attraction = "marche", parametres = c(5, 40),
               control = list(maxit = 100))


