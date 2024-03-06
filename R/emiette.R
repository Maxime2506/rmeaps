#' Emiette les actifs et produit un jeu d'arrangements (shuffles) des actifs
#'
#' @param les_actifs le tibble ou le vecteur des actifs
#' @param nshuf le nombre de tirages
#' @param seuil le seuil d'emmiettage
#' @param var si les actifs sont un tibble, le nom de la variable 'actifs' ('actifs' par défaut)
#' @param weighted si TRUE, les actifs sont pondérés par leur fréquence (TRUE par défaut)
#'
#' @return une matrice d'arrangement (autant de lignes que nshuf, autant de colonnes que d'actifs)
#' @export
#'

emiette <- function(les_actifs, nshuf = 256, seuil=40, var = "actifs", weighted=TRUE) {
  if(tibble::is_tibble(les_actifs))
    act <- les_actifs |> 
      pull(var)
  if(rlang::is_vector(les_actifs))
    act <- les_actifs
  freq <- act %/% seuil + 1
  set <- NULL
  w <- NULL
  for(i in 1:length(freq)) {
    set <- c(set, rep(i, freq[i]))
    w <- c(w, rep(act[i]/freq[i], freq[i]))
  }
  shuf <- matrix(NA, ncol = sum(freq), nrow = nshuf)
  if(weighted==FALSE)
    w <- NULL
  for(i in 1:nshuf) 
    shuf[i, ] <- sample(set, sum(freq), prob = w, replace=FALSE)
  return(shuf)
}