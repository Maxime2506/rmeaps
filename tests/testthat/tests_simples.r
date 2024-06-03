# Tests de base.
library(tidyverse)
# ---- 1. Premier cas ----
triplet <- expand.grid(toidINS = paste("to", 0:3), fromidINS = paste("from", 0:1)) |> 
  mutate(fromidINS = as.character(fromidINS), toidINS = as.character(toidINS), metric = c(1:4, 4:1)) |>
  arrange(fromidINS, metric)
emplois <- rep(20,4)
names(emplois) <- paste("to", 0:3)
actifs <- c(50, 50)
fuites <- c(.2, .2)
names(actifs)  <- names(fuites) <- paste("from", 0:1)

test_that("Création d'un MeapsData", {
  expect_no_warning(md <- meapsdata(triplet, actifs, emplois, fuites))
})

# Calcul direct des résultats
.f_distrib <- function(emp, f, act, attract = 1L) {
  accessibiliy <- cumsum(emp * attract)
  absorption <- -log(f)/max(accessibiliy)
  passants <- act * exp(-absorption * accessibiliy)
  lag(passants, default = act) - passants
}

ligne1 <- .f_distrib(emplois, fuites[1], actifs[1])
ligne2 <- .f_distrib(emplois, fuites[2], actifs[2])
exces <- emplois - ligne1 - ligne2
emploislibres <- ifelse(exces > 0, exces, 0)
actifslibres <- actifs/sum(actifs)*sum(emploislibres-exces)/(1-fuites)

ligne1.1 <- .f_distrib(emploislibres, fuites[1], actifslibres[1])
ligne2.1 <- .f_distrib(emploislibres, fuites[2], actifslibres[2])
exces.1 <- emploislibres - ligne1.1 - ligne2.1
emploislibres.1 <- ifelse(exces.1 > 0, exces.1, 0)
actifslibres.1 <- actifslibres/sum(actifslibres)*sum(emploislibres.1-exces.1)/(1-fuites)

ligne1.2 <- .f_distrib(emploislibres.1, fuites[1], actifslibres.1[1])
ligne2.2 <- .f_distrib(emploislibres.1, fuites[2], actifslibres.1[2])

emplois_pris <- emplois - emploislibres
ligne1 <- ligne1

test_that("Meaps simple", {
  all_in(md)

})
