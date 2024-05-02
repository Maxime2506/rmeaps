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
      dplyr::pull(var, name = idINS)
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

#' Recalcule la matrice shuf sur un nouvel ordre d'actifs (l'initial est dans les colnames)
#'
#' @param shuf le tri initial, avec les noms des actifs ordonnés dans les colonnes
#' @param actifs Un vecteur d'actifs dans un autre ordre
#'
#' @return une matrice d'arrangement (autant de lignes que nshuf, autant de colonnes que d'actifs, dont les colonnes sont l'idINS)
#' @export
#'
reordonne_shuf <- function(shuf, actifs) {
  ori <- unique(colnames(shuf))
  new <- unique(names(actifs))
  if(length(setdiff(new, ori))!=0)
    stop("Des actifs ne sont pas dans la matrice shuf.")
  ori_inter <- intersect(ori, new)
  
  col_names <- colnames(shuf)
  col_names <- col_names[col_names%in%new]
  col_names <- factor(col_names, new) |>
    sort() |>
    as.character()
  n_act <- length(new) 
  new_index <- rlang::set_names(1:n_act, new)
  shuf_names <- ori[t(shuf)]
  new_shuf <- new_index[shuf_names[shuf_names%in%new]] |> 
    matrix(nrow = nrow(shuf), ncol = length(col_names), byrow=TRUE)
  colnames(new_shuf) <- col_names
  return(new_shuf)
}
