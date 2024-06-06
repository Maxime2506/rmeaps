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

#' Fonctions de calcul de la solution R.

.f_distrib_subj <- function(emp, f, act, attract = 1L) {
  accessibiliy <- cumsum(emp * attract)
  absorption <- -log(f)/max(accessibiliy)
  passants <- exp(-absorption * accessibiliy) * as.numeric(act)
  res <- lag(passants, default = act) - passants
  names(res) <- names(emp)
  res
}

.f_distrib <- function(i, emp = emplois, act = actifs, attract = 1L) {
  z <- triplet$fromidINS |> unique() |> sort()
  tos <- triplet$toidINS |>unique()  |> sort()
  ord <- triplet |> filter(fromidINS == z[i]) |> arrange(metric) |> pull(toidINS)
  .f_distrib_subj(emp[ord], fuites[z[i]], act[z[i]], attract = attract)[tos]
}

.retourner <- function(list_distrib, places) {
  tot <- reduce(list_distrib, `+`)
  dis <- map2(list_distrib, names(actifs), \(x,y) as_tibble_col(x, column_name = y)) |> 
    list_cbind() |> 
    mutate(tot = tot)
  
  post <- dis |>
    mutate(across(.col = all_of(names(actifs)), 
                  .fns = \(x) if_else(tot > places,  x * places / tot, x)))

  retours <- dis |>
    mutate(across(.col = all_of(names(actifs)), 
                  .fns = \(x) if_else(tot > places,  x - x * places / tot, 0))) |>
    summarise(across(.col = all_of(names(actifs)), .fns = \(x) sum(x)))

  libres <- post |> select(-tot) |> rowSums()
  libres <- places - libres
  libres <- ifelse(libres < 0, 0, libres)
  retours <- (retours / (1 - fuites[names(actifs)])) |> as.numeric()
  names(retours) <- names(actifs)
  names(libres) <- names(emplois)
  list(distrib = post |>select(-tot) |>as.matrix() |> t(), 
       retours = retours,
       libres = libres)
}


# ----- 2. Solutions de référence ------
# Cas all_in
list_distrib <- map(1:2, .f_distrib)
retours <- .retourner(list_distrib, emplois)

list2 <- map(1:2, \(x) .f_distrib(x, emp = retours$libres, act = retours$retours))
retours2 <- .retourner(list2, retours$libres)

ref_allin <- retours$distrib + retours2$distrib
ref_allin <- Matrix::mat2triplet(ref_allin) |> as_tibble()
ref_allin <- ref_allin |>
  left_join(tibble(fromidINS = names(actifs), i = seq_along(actifs)), by = "i") |>
  left_join(tibble(toidINS = names(emplois), j = seq_along(emplois)), by = "j") |>
  select(-c(i,j)) 

ref_allin <- triplet |> 
  relocate(fromidINS, toidINS) |> 
  inner_join(ref_allin, by = c("fromidINS", "toidINS")) |>
  select(-metric) |> rename(flux = x)

attr(ref_allin, 'out.attrs') <- NULL

# Cas multishuf
ligne1 <- .f_distrib(1)
emplois1 <- emplois - ligne1
ligne2 <- .f_distrib(2, emp = emplois1)
prise2 <- ifelse(ligne2 > emplois1, emplois1, ligne2)
ret2 <- sum(ligne2-prise2)
emplois1_1 <- emplois1 - prise2

# En fait les actifs de la ligne 2 prennent tout ce qui reste.
ref_multishuf <- c(ligne1, emplois1) |> 
  matrix(byrow = TRUE, nrow = 2) |>
  Matrix::mat2triplet() |>
  as_tibble() |>
  transmute(fromidINS = names(actifs)[i], toidINS = names(emplois)[j], flux = x) |>
  arrange(fromidINS, toidINS, flux)

# ----- 3. Construction des objets MeapsData -----
md <- meapsdata(triplet, actifs, emplois, fuites)
mds <- meapsdata(triplet, actifs, emplois, fuites, nshuf = 1)

# ----- 4. Tests des différentes methodes -----
estim_allin <- all_in(md)

test_that("Meaps all_in cas simple", { 
  expect_equal(estim_allin, ref_allin, tolerance = 1e-3)
  }
)

estim_multi_origin <- multishuf_origin(mds) |> arrange(fromidINS, toidINS, flux)

test_that("Meaps multishuf origin cas simple", {
  expect_equal(estim_multi_origin, ref_multishuf, tolerance = 1e-3)
  }
)

estim_multi_task <- multishuf_task(mds) |> arrange(fromidINS, toidINS, flux)

test_that("Meaps multishuf task cas simple",
  expect_equal(estim_multi_task, ref_multishuf, tolerance = 1e-3)
)

estim_multi_oc <- multishuf_oc(mds)$flux |> arrange(fromidINS, toidINS, flux)

test_that("Meaps multishuf oc cas simple",
  expect_equal(estim_multi_oc, ref_multishuf, tolerance = 1e-3)
)

test_that("Meaps multishuf oc = task ?", {
  expect_equal(estim_multi_oc, estim_multi_task, tolerance = 1e-3)
})
