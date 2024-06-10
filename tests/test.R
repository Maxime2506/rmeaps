# Tests de base.
library(tidyverse)
devtools::load_all()

# ----- 1. Cas simple -----
triplet <- expand.grid(toidINS = paste("to", 0:3), fromidINS = paste("from", 0:1)) |> 
  mutate(fromidINS = as.character(fromidINS), toidINS = as.character(toidINS), metric = c(1:4, 4:1)) |>
  arrange(fromidINS, metric)
emplois <- rep(20,4)
names(emplois) <- paste("to", 0:3)
actifs <- c(50, 50)
fuites <- c(.2, .2)
names(actifs)  <- names(fuites) <- paste("from", 0:1)

# ----- 3. Construction des objets MeapsData -----
md <- meapsdata(triplet, actifs, emplois, fuites)
mds <- meapsdata(triplet, actifs, emplois, fuites, nshuf = 16, seuil = 100)

# ----- 4. Tests des différentes methodes -----
all_in(md, attraction = "marche", parametres = c(2.5, .1))
multishuf_task(mds, attraction = "marche", parametres = c(2.5, .1))
multishuf_oc(mds, attraction = "marche", parametres = c(2.5, .1))$flux |> arrange(fromidINS)




