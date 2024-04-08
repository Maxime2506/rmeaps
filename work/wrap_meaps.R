#' La fonction meaps sur plusieurs shufs sur un format Sparse Row Matrix et un random des rangs sur les distances égales mélangées pour chacun des boot.
#' @param dist Matrice des distances où des résidents en ligne i rejoignent des opportunités en colonnes j. Les NA's sont évacués au format sparse. Les distances zeros ont le 1er rang. Au format matrix ou sparse (R). Dans ce dernier cas, la matrice est passée telle quelle.
#' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes).
#' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
#' @param modds La matrice des odds modifiant la chance d'absorption de chacun des sites j pour des résidents en i.
#' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude.
#' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs.
#' @param plafond_odd détermine la valeur maximale (et minimale) des odds, ramenés entre 1/plafond et plafond.
#' @param mode Choix du rôle de l'emploi disponible au cours du processus. Default : continu. Autre choix : discret, subjectif_c ou _d...
#' @param odds_subjectifs Attractivité des sites proches des actifs, pour le mode defini. default : null.
#' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto.
#' @param progress Ajoute une barre de progression. Default : true.
#' @param normalisation Calage des emplois disponibles sur le nombre d'actifs travaillant sur la zone. Default : false.
#' @param fuite_min Seuil minimal pour la fuite d'un actif. Doit être supérieur à 0. Défault = 1e-3.
#' 
#' @return renvoie une matrice avec les estimations du nombre de trajets de i vers j.
#' @export
meaps_sparse <- function(dist, emplois, actifs, modds, f, shuf, plafond_odd = 1000,
                         mode = "continu",
                         oddssubjectifs = NULL,
                         nthreads = 0,
                         progress = TRUE,
                         normalisation = FALSE,
                         fuite_min = 1e-3,
                         seuil_newton = 1e-6) {
  
  if (inherits(dist, "matrix")) {
    # On remplace les valeurs nulles par une petite valeur quelconque qui sera toujours au 1er rang.
    seuil_minimal <- min(as.numeric(dist[dist != 0]), na.rm = TRUE) / 2
    dist[dist == 0] <- seuil_minimal
    dist[is.na(dist)] <- 0

    dist_dgr <- as(dist, "RsparseMatrix")
  } else if (inherits(dist, "dsRMatrix")) {
    dist_dgr <- dist # On suppose qu'au format sparse la matrice est préparée comme voulue.
  } else {
    stop("Le format de dist n'est pas géré.")
  }
  
  if (inherits(modds, "matrix")) {
    if (anyNA(modds)) stop("Il y a des NA's dans la matrice des odds.")
    if (any(modds <= 0)) stop("Il y a des valeurs invalides dans la matrice des odds.")
    # On remplace les valeurs extrêmes par le plafond.
    plancher_odd <- 1/plafond_odd
    modds[modds < plancher_odd] <- plancher_odd
    modds[modds > plafond_odd] <- plafond_odd
    
    modds_dgr <- as(log(modds), "RsparseMatrix") # Passage en log pour profiter du sparse sur les odds=1. L'opération inverse est gérée en c++.
  } else if (inherits(modds, "dsRMatrix")) {
    modds_dgr <- modds # On suppose qu'au format sparse la matrice est préparée comme voulue au format LOG.
  } else {
    stop("Le format de modds n'est pas géré.")
  }
  
  if (any(dim(dist) != dim(modds))) stop("dist et modds ont des dimensions différentes.")
  
  dist_dgr@x <- meaps_dgr(j_dist = dist_dgr@j,
                          p_dist = dist_dgr@p,
                          x_dist = dist_dgr@x,
                          N = dist_dgr@Dim[1],
                          K = dist_dgr@Dim[2],
                          emplois = emplois,
                          actifs = actifs,
                          j_modds = modds_dgr@j,
                          p_modds = modds_dgr@p,
                          x_modds = modds_dgr@x,
                          f = f,
                          shuf = shuf,
                          mode = mode,
                          oddssubjectifs = oddssubjectifs,
                          nthreads = nthreads,
                          progress = progress,
                          normalisation = normalisation,
                          fuite_min = fuite_min,
                          seuil_newton = seuil_newton)
  
  dist_dgr
}