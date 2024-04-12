#' Emiette les actifs et produit un jeu d'arrangements (shuffles) des actifs
#'
#' @param les_actifs le tibble ou le vecteur des actifs (nommé par idINS)
#' @param nshuf le nombre de tirages
#' @param seuil le seuil d'emmiettage
#' @param var si les actifs sont un tibble, le nom de la variable 'actifs' ('actifs' par défaut)
#' @param weighted si TRUE, les actifs sont pondérés par leur fréquence (TRUE par défaut)
#'
#' @return une matrice d'arrangement (autant de lignes que nshuf, autant de colonnes que d'actifs, dont les colonnes sont l'idINS)
#' @export
#'

emiette <- function(les_actifs, nshuf = 256, seuil=40, var = "actifs", weighted=TRUE) {
  if(tibble::is_tibble(les_actifs))
    act <- les_actifs |> 
      pull(var, name = idINS)
  if(rlang::is_vector(les_actifs))
    act <- les_actifs
  freq <- act %/% seuil + 1
  ll <- sum(freq)
  k <- 0
  set <- numeric(ll)
  w <- numeric(ll)
  noms <- character(ll)
  for(i in 1:length(freq)) {
      j <- freq[i]
      set[(k+1):(k+j)] <-  rep(i, j)
      w[(k+1):(k+j)] <- rep(act[i]/j, j)
      noms[(k+1):(k+j)] <- rep(names(act)[i], j)
      k <- k+j
  }
  shuf <- matrix(NA, ncol = sum(freq), nrow = nshuf)
  if(weighted==FALSE)
    w <- NULL
  for(i in 1:nshuf) 
    shuf[i, ] <- sample(set, sum(freq), prob = w, replace=FALSE)
  colnames(shuf) <- noms
  return(shuf)
}