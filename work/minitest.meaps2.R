library(rmeaps)
library(Matrix)
library(matrixStats)
library(tidyverse)

actifs <- c(11, 11, 11, 11, 11)
fuite <- c(0.1, 0.1, 0.1, 0.1, 0.1)
names(actifs) <- str_c("a", 1:5)
names(fuite) <- names(actifs)
emplois <- c(5, 5, 5, 5, 30)
names(emplois) <- str_c("e", 1:5)
distances <- matrix(
  c(0, 1, 2, 3, 5,
    1, 0, 1, 2, 6,
    2, 4, 0, 1, 2,
    3, 2, 5, 0, 1,
    4, 3, 2, 3, 0),
  nrow = 5, dimnames = list(names(actifs), names(emplois)))
triplet <- distances |>
  as_tibble(rownames = "fromidINS") |>
  pivot_longer(cols = -fromidINS, names_to = "toidINS", values_to = "metric") |>
  filter(!is.na(metric)) |> 
  arrange(fromidINS, metric, toidINS)

data <- meapsdata(
  triplet = triplet, actifs = actifs, emplois = emplois, fuites=fuite, 
  nshuf=6, seed = 42)
cible <- all_in(data) |> arrange(desc(flux)) |> 
  rename(group_from = fromidINS, group_to = toidINS, value = flux)

multishuf_oc(data, verbose = TRUE)
all_in(data) |> pull(flux) |> sum()
data_g <- meapsdatagroup(data, group_from = set_names(names(actifs)), group_to = set_names(names(emplois)), cible = cible)

all_in_grouped(data_g)

multishuf_oc_grouped(data_g)
