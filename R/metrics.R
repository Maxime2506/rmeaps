#' Calcul de la béta-divergence Kullback-Leibler, soit l'entropie relative.
#' @param nb_mod Vecteur dont on veut évaluer la divergence.
#' @param nb_ref Vecteur de référence.
#' @param seuil_collapse Seuil en dessous duquel on aggrège les effectifs. Défaut : 0,1.
#' 
#' @return Valeur de Kullback-Leibler.
kullback_leibler <- function(nb_mod, nb_ref, seuil_collapse = .1) {
  
  tib <- tibble(x = nb_mod, y = nb_ref) |> 
    mutate(
      g = cur_group_rows(),
      g = ifelse( pmin(nb_ref, nb_mod) >= seuil_collapse, g, 0)
    ) |>
    group_by(g) |> 
    summarise(x = sum(x, na.rm=TRUE), y = sum(y, na.rm=TRUE)) |> 
    mutate(px = x / sum(x, na.rm=TRUE), py = y / sum(y, na.rm=TRUE))
  
  sum(tib$py * log(tib$py / tib$px), na.rm=TRUE)
}

#' Calcul de l'entropie de Shannon.
#' @param nb_ref Vecteur dont on veut mesurer l'entropie.
#' @param seuil_collapse Seuil en dessous duquel on aggrège les effectifs. Défaut : 1.
#' 
#' @return Entropie.
shannon <- function(nb_ref, seuil_collapse = 1) {
  tib <- tibble(y = nb_ref) |> 
    mutate(
      g = cur_group_rows(),
      g = ifelse(nb_ref >= seuil_collapse, g, 0)
    ) |>
    group_by(g) |> 
    summarise( y = sum(y, na.rm=TRUE)) |> 
    mutate(py = y / sum(y, na.rm=TRUE))
  
  sum(tib$py * log(tib$py), na.rm=TRUE)
}

entropie_relative <- function(x, y, floor) {
  py <- y / sum(y)
  px <- x / sum(x)
  py <- pmax(py, floor)
  px <- pmax(px, floor)
  sum( py * log(py/px) )
}

#' R2 sur la base de l'entropie relative
#' On prend comme référence une distribution uniforme
#' On calcule le KL de cette distribution et de la distributiopn observée (référence)
#' Puis le KL de la distribution prposée et de la distribution de référence
#' on défini alors R2 = 1-kl(est, obs)/kl(unif, obs)
#' Attention ! On exclue les valeurs nulles des deux distributions
#' @param estime distribution empirique estimée
#' @param observe distribution empirique observée
#' @param seuil Seuil de la distribution observée cumulée au dessus duquel on aggrège les effectifs. Défaut : 99%.
#' @param bruit remplace les valeurs inférieures au bruit par ce bruit. Défaut 1e-6
#' @return une valeur de R2 filtré ($r2kl) et une non filtrée
r2kl2 <- function(estime, observe, seuil = .99, bruit=1e-6) {
  tib <- tibble(x=estime, y=observe) |> 
    filter(x!=0, y!=0)
  r2kln0 <- 1-kl(tib$x, tib$y)/kl(rep(mean(tib$y), nrow(tib)), tib$y)
  tib <- tib |> 
    arrange(desc(y)) |> 
    mutate(cum = cumsum(y)/sum(y)) |> 
    filter(cum<=seuil)
  r2klnb <- 1-kl(tib$x, tib$y)/kl(rep(mean(tib$y), nrow(tib)), tib$y)
  return(list(r2kl = r2kln0, r2kl_l = r2klnb))
}

#' Un bête khi2
#' @param estime distribution empirique estimée
#' @param observe distribution empirique observée
#' @return khi2
khi2 <- function(estime, observe) {
  sum((estime - observe)^2/observe, na.rm=TRUE)
}