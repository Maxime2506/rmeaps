#' La fonction prep_meaps_dist prépare les données de distances pour traitement par meaps_optim
#' @param dist triplet des distances où des résidents en ligne i rejoignent des opportunités en colonnes j. 
#' 
#' @return renvoie une RankedRSMatrix des distances.
#' @import Matrix
#' @export
prep_meaps_dist <- function(dist) {
  
  if (!is_triplet(dist)) stop("Ce n'est pas un triplet valide.")
  dist <- triplet2listij(dist)
  
  list(RankedMat = new(RankedRSMatrix, dist$dgr), cle_from = dist$cle_from, cle_to = dist$cle_to) 
}

#' La fonction prep_meaps_odds prépare la matrice des odds pour traitement par meaps_optim.
#' @param modds un triplet des log odds.
#' @param cle_from la clé_from issue de distance.
#' @param cle_to la clé to 
#' 
#' @return renvoie une RankedRSMatrix des odds, rangés selon le rang des distances.
#' @import Matrix
#' @export
prep_meaps_odds <- function(modds, cle_from, cle_to) {
  
  lodds <- tripletlodds2dgr(modds, cle_from, cle_to)
  if (length(lodds) == 0) stop("La matrice est vide.")
  
  new(RankedRSMatrix, lodds)
}

#' La fonction meaps_optim sert pour la recherche des meilleurs paramètres ou odds lorsqu'un certain regroupement des flux est connu.
#' @param dist_prep Matrice des distances au format RankedRSMatrix (sortie de la fonction prep_meaps_dist).
#' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes).
#' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
#' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude.
#' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs.
#' @param group_from Le vecteur de regroupement des communes de départs (lignes).
#' @param group_to Le vecteur de regroupement des communes de destinations (colonnes).
#' @param attraction choix de la fonction d'attraction. Default = "constant". "marche", "logistique" ou "odds".
#' @param param paramètres pour mode = "marche" (le 1er est la distance de la marche, le 2nd est le plancher après la marche).
#' pour mode = "logistique" (le 1er est la distance du point de symétrie, le 2nd la raideur de la bascule, le 3ème le plancher).
#' @param odds_prep les log(odds) au format RankedRSMatrix, mais rangé selon dist (résultat de prep_meaps_odds).
#' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto.
#' @param progress Ajoute une barre de progression. Default : true.
#' @param normalisation Calage des emplois disponibles sur le nombre d'actifs travaillant sur la zone. Default : false.
#' @param fuite_min Seuil minimal pour la fuite d'un actif. Doit être supérieur à 0. Défault = 1e-3.
#' 
#' @return renvoie une matrice avec les estimations des flux regroupés.
#' @export
meaps_optim <- function(dist_prep, emplois, actifs, fuite, shuf, 
                        groups_from, groups_to,
                        attraction = "constant",
                        param = 0.0,
                        odds_prep = NULL,
                        nthreads = 0,
                        progress = TRUE,
                        normalisation = FALSE,
                        fuite_min = 1e-3) {
  
  if (sum(actifs * (1 - fuite)) != sum(emplois)) warning("Les actifs et les emplois ne correspondent pas, à la fuite près.")
  
  # contraintes sur les paramètres.
  if (attraction == "marche") {
    if (param[1] <= 0) stop("Le 1er paramètre doit indiquer la distance où se situe la marche.")
    if (param[2] <= 0) stop("Le 2nd paramètre doit représenter le facteur d'attractivité après la marche, (réf avant : 1")
  }
  
  if (attraction == "logistique") {
    if (param[1] <= 0) stop("Le 1er paramètre doit indiquer la distance où la logistique bascule (le point de symétrie).")
    if (param[2] <= 0) stop("Le 2nd paramètre indique la raideur de la bascule.")
    if (param[3] < 0) stop("Le 3ème paramètre indique un plancher (=limite pour x infini).")
  }
  
  if (!is.null(odds_prep)) {
    attraction <- odds_prep$attraction
    jr_odds <- odds_prep$lodds$jr
    p_odds <- odds_prep$lodds$p
    xr_odds <- odds_prep$lodds$xr
  } else {
    jr_odds <- 0L
    p_odds <- 0L
    xr_odds <- 0.0
  }
  
  row_group = groups_from - 1L # les indices c++ commencent à zéros.
  col_group = groups_to - 1L

  .meaps_optim(jr_dist = mat$jr,
                  p_dist = mat$p,
                  xr_dist = mat$xr,
                  emplois = emplois,
                  actifs = actifs,
                  f = fuite,
                  shuf = shuf,
                  row_group = row_group,
                  col_group = col_group,
                  param = param,
                  attraction = attraction,
                  jr_odds = jr_odds,
                  p_odds = p_odds,
                  xr_odds = xr_odds,
                  nthreads = nthreads,
                  progress = progress,
                  normalisation = normalisation,
                  fuite_min = fuite_min)
  
}


