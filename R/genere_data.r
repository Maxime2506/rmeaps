#' Genere des données pour les tests
#'
#' @param n nombre de points de résidents (sur une grille 0:1x0:1) défaut : 2
#' @param k nombre de points d'emplois (sur une grille 0:1x0:1) défaut : 2
#' @param nshuf nombre de tirages défaut : 64
#' @param nb_actifs nombre d'actifs défaut : 10xn^2
#' @param nb_emplois nombre d'emplois défaut : 9xk^2
#' @param fuite fuite défaut : 10%
#' @param densite forme de la répartition des actifs et des emplois. Si "uniforme" alors c'est uniforme, "radial" c'ets en 1/r^2 ou r est la distance au centre
#'
#' @return une liste pour de données meaps (rkdist les rangs, dist les distances, actifs les actifs, emplois les emplois, shuf un shuf)
#' @export
#'
#' @examples
#' genere_data(densite="radiale")
genere_data <- function(n = 3, k = 3, nshuf = 64, 
                        nb_actifs = n*n*10, nb_emplois = nb_actifs*(1-fuite), 
                        fuite = 0.1,
                        densite = "uniforme",
                        ties = "random") {
  
  maxx <- 1
  maxy <- 1
  
  residences <- tidyr::expand_grid(
    x=seq(0, maxx, length.out=n),
    y=seq(0, maxy, length.out=n)) |> 
    as.matrix() |> 
    sf::st_multipoint(dim = "XY") |> 
    sf::st_sfc() |> 
    sf::st_cast(to = "POINT")
  
  emplois <- tidyr::expand_grid(x=seq(0, maxx, length.out=k), y=seq(0, maxy, length.out=k)) |> 
    as.matrix() |> 
    sf::st_multipoint(dim = "XY") |> 
    sf::st_sfc() |> 
    sf::st_cast(to = "POINT")
  
  distance <- sf::st_distance(residences, emplois)
  rkdist <- matrixStats::rowRanks(distance, ties.method = ties)
  fn_dens <- switch(
    densite,
    "uniforme" = \(x) 1,
    "radiale" =\(x) 1 / (0.5+sf::st_distance(x, sf::st_point(c(.5, .5)), by_element = FALSE))^2)
  
  marge_emplois <- tibble::tibble(position = emplois) |> 
    dplyr::mutate(dense = fn_dens(position),
           emplois = nb_emplois * dense / sum(dense)) |> 
    dplyr::pull(emplois)
  
  marge_actifs <- tibble::tibble(position = residences) |> 
    dplyr::mutate(dense = fn_dens(position),
           actifs = nb_actifs * dense / sum(dense)) |> 
    dplyr::pull(actifs)
  
  mat_odds <- matrix(1, nrow = length(marge_actifs), ncol = length(marge_emplois))
  
  shuf <- purrr::map(1:nshuf, ~sample.int(length(marge_actifs), length(marge_actifs)))
  shuf <- do.call(rbind, shuf)

  return(list(rkdist = rkdist,
              dist = distance,
              actifs = marge_actifs,
              emplois = marge_emplois,
              fuite = rep(fuite, length(marge_actifs)),
              modds = mat_odds,
              shuf = shuf))
}