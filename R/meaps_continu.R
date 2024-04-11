#' La fonction meaps sur plusieurs shufs selon le mode de calcul continu.
#' @param dist Soit un matrice des distances où des résidents en ligne i rejoignent des opportunités en colonnes j (en ce cas les NAs deviendront sparses).
#' Ou bien un objet de class dgRMatrix (Row Sparse Matrix du pkg Matrix).
#' Ou bien un triplet au format list ou data.frame. Dans ce cas, les variables sont désignés par i, j et x pour les valeurs.
#' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes).
#' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
#' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude.
#' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs.
#' @param attraction Choix de la fonction d'attraction des différents sites, appliquée à l'accessibilité. 
#' Par défaut, "constant" où aucun site n'a plus d'attrait qu'un autre. 
#' "marche" où l'attrait vaut 1 jusqu'à une certaine distance (param alpha) puis moins (param béta). f(x) = 1 si x < alpha, = béta si x > alpha.
#' "logistique" où l'attrait décroît selon une fonction logistique avec une distance de bascule (param alpha), une vitesse de bascule (param béta) 
#' et un seuil (param gamma). Si h(x) = exp( (x-alpha)/béta), f(x) = gamma + h(x) / (1 + h(x)).
#' "odds" où chaque flux (from, to) se voit attribuer un odds. Dans ce cas, on entre un Row Sparse Matrix des log(odds) selon ses éléments.
#' @param param Dans les cas de "marche" et "logistique", un vecteur avec dans l'ordre les valeurs des paramètres alpha, béta et, si nécessaire, gamma.
#' @param modds Dans le cas de "odds", une matrice des odds, un triplet au format list ou data.frame des odds ou une dgRMatrix des log(odds).
#' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto.
#' @param progress Ajoute une barre de progression. Default : true.
#' @param normalisation Calage des emplois disponibles sur le nombre d'actifs travaillant sur la zone. Default : false.
#' @param fuite_min Seuil minimal pour la fuite d'un actif. Doit être supérieur à 0. Défault = 1e-3.
#' 
#' @return renvoie un triplet au format data.frame des estimations du nombre de trajets de i vers j.
#' @import Matrix
#' @export
meaps_continu <- function(dist, emplois, actifs, f, shuf, 
                          attraction = "constant",
                          param = numeric(),
                          modds = NULL,
                          nthreads = 0,
                          progress = TRUE,
                          normalisation = FALSE,
                          fuite_min = 1e-3,
                          seuil_newton = 1e-6) {
  
  
  if (sum(actifs * (1 - f)) != sum(emplois)) warning("Les actifs restant dans la zone et les emplois ne correspondent pas.")
  
  dist <- if (inherits(dist, "matrix")) {
    .matrix2dgR(dist) 
  } else if (is_triplet(dist) ) {
    .triplet2dgR(dist)
  } else if (inherits(dist, "dgRMatrix")) {
    dist
  } else {
    stop("Format pour dist non reconnu.")
  }
  
  odds <- new("dgRMatrix")
  if (attraction == "odds") {
    odds <- if (inherits(modds, "matrix")) {
      .matrix2dgR(log(modds), na2zero = FALSE)
    } else if (is_triplet(modds)) {
      colnames(modds)[colnames(modds) == "fromidINS"] <- "i"
      colnames(modds)[colnames(modds) == "toidINS"] <- "j"
      colnames(modds)[!colnames(modds) %in% c("i", "j")] <- "x"
      sparseMatrix(i = modds$i, j = modds$j, x = modds$x, dims = dim(dist))
    } else if (inherits(modds, "dgRMatrix")){
      modds
    } else {
    stop("Format pour modds non reconnu.")
    }
    if (length(odds@x) == 0) attraction <- "constant"
  }
  
  dist@x <- meaps_continu_cpp(j_dist = dist@j,
                              p_dist = dist@p,
                              x_dist = dist@x,
                              emplois = emplois,
                              actifs = actifs,
                              f = f,
                              shuf = shuf,
                              attraction = attraction,
                              param = param,
                              j_odds = odds@j,
                              p_odds = odds@p,
                              x_odds = odds@x,
                              nthreads = nthreads,
                              progress = progress,
                              normalisation = normalisation,
                              fuite_min = fuite_min)
  
  dist <- as(as(as(dist, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  data.frame(i = dist@i, j = dist@j, flux = dist@x)
}


