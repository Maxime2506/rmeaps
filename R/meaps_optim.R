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
#' @import Matrix
#' @export
meaps_optim <- function(dist, emplois, actifs, f, shuf, 
                        groups_from, groups_to,
                        objectif = NULL,
                        attraction = "constant",
                        alpha = 10, beta = 1,
                        nthreads = 0,
                        progress = TRUE,
                        normalisation = FALSE,
                        fuite_min = 1e-3) {
  
  
  mat <- if (inherits(dist, "matrix")) {
          .transfom_matrix(dist) 
         } else if (inherits(dist, "dgRMatrix")) {
          .transform_triplet(dist)
         } else {
             stop("Format non reconnu")
           }
  # contraintes sur les paramètres.
  if (attraction == "marche") {
    if (alpha <= 0) stop("alpha doit indiquer la distance où se situe la marche.")
    if (beta <= 0) stop("beta doit représenter le facteur d'attractivité après la marche, (réf avant : 1")
  }
  
  if (attraction == "logistique") {
    if (alpha <= 0) stop("alpha doit indiquer la distance où la logistique bascule (le point de symétrie).")
    if (beta <= 0 | beta >= 1) stop("beta indique la force de l'effet (entre 0 et 1)")
  }
  
  
  
  mat <- new(RankedRSMatrix, mat)
  
  row_group = groups_from - 1L # les indices c++ commencent à zéros.
  col_group = groups_to - 1L
  
  meaps_optim_cpp(jr_dist = mat$jr,
                  p_dist = mat$p,
                  xr_dist = mat$xr,
                  emplois = emplois,
                  actifs = actifs,
                  f = f,
                  shuf = shuf,
                  row_group = row_group,
                  col_group = col_group,
                  attraction = attraction,
                  alpha = alpha, beta = beta,
                  nthreads = nthreads,
                  progress = progress,
                  normalisation = normalisation,
                  fuite_min = fuite_min)
  

}

#' Fonction de choix d'une petite distance pour remplacer les zéros éventuels.
#' Normalement il ne devrait pas y avoir de zéros. 
#' Au sein d'un même carreau de dim a, la distance entre deux points est d'environ a/2.
#' La distance retenue importe peu tant que les rangs sont respectés.
.petite_distance <- function(d, facteur = 10) {
  min(as.numeric(d[d != 0]), na.rm = TRUE) / facteur
}


#' Fonction de retraitement de la distance si elle est sous forme d'une matrice classique.
#' Les NA doivent devenir sparses. Les zéros deviennent simplement des petites valeurs.
#' 
#' @import Matrix
.transfom_matrix <- function(dist) {
  dist[dist == 0] <- .petite_distance(dist)
  dist[is.na(dist)] <- 0
  as(dist, "RsparseMatrix")
}

#' Fonction de retraitement de la distance si elle est sous forme de triplet dans une liste.
#' C'est en général la résultante de la fonction Matrix::mat2triplet.
#' On suppose que max(j) donne le nombre de colonnes, ce qui est attendu dans les analyses meaps.
#' Ici encore les valeurs x nulles deviennent des petites distances.
#' Liste de trois vecteurs : i, j et x.
#' @import Matrix
.transform_triplet <- function(dist) {
  if ( { setdiff(names(dist), c("i", "j", "x")) != character(0) } ||
       { length(dist$i) != length(dist$j) } ||
       { lenght(dist$i) != length(dist$x) } ||
       { !is.integer(dist$i) || !is.integer(dist$j) || !is.numeric(dist$x)}
  ) stop("Dist n'est pas une liste de triplets valide.")
  
  dist$x[dist$x == 0] <- .petite_distance(dist$x)
  
  spMatrix(nrow = length(dist$i), ncol = max(dist$j), i = dist$i, j = dist$j, x = dist$x)
  
}

