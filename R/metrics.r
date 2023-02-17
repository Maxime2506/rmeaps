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