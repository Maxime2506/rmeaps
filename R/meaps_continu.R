#' La fonction meaps sur plusieurs shufs selon le mode de calcul continu.
#' @param dist un triplet sous forme de data.frame avec des colonnes fromidINS et toidINS.
#' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes).
#' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
#' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude.
#' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs.
#' @param attraction Choix de la fonction d'attraction des différents sites, appliquée à l'accessibilité. 
#' Par défaut, "constant" où aucun site n'a plus d'attrait qu'un autre. 
#' "marche" où l'attrait vaut béta jusqu'à une certaine distance (param alpha) puis 1. f(x) = béta si x < alpha, = 1 si x > alpha.
#' "logistique" où l'attrait décroît selon une fonction logistique avec une distance de bascule (param alpha), une vitesse de bascule (param béta) 
#' et un seuil (param gamma). Si h(x) = exp( (x-alpha)/béta), f(x) = gamma + h(x) / (1 + h(x)).
#' "odds" où chaque flux (from, to) se voit attribuer un odds. Dans ce cas, on entre un Row Sparse Matrix des log(odds) selon ses éléments.
#' @param param Dans les cas de "marche" et "logistique", un vecteur avec dans l'ordre les valeurs des paramètres alpha, béta et, si nécessaire, gamma.
#' @param modds un triplet au format data.frame des log(odds).
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
                          seuil_newton = 1e-6, 
                          quiet = FALSE) {
  
  delta <- (sum(actifs * (1 - f)) - sum(emplois))/sum(emplois)
  if (delta>10^(-5)&!quiet)
    warning(glue("Les actifs restant dans la zone et les emplois ne correspondent pas à {round(delta*100,1)}% près."))
  
  if (!is_triplet(dist)) stop("Format pour dist non reconnu.")
  dist <- triplet2listij(dist)
  cle_from <- dist$cle_from
  froms <- as.character(cle_from)
  cle_to <- dist$cle_to
  dist <- dist$dgr
  actifs <- actifs[froms]
  shuf <- shuf[, froms]
  f <- f[froms]
  emplois <- emplois[as.character(cle_to)]
  jodds <- podds <- xodds <- 1L
  if (is.null(modds)) { 
    if (attraction == "odds") stop("Il n'y a pas de odds définis.")
  }
  
  if(!is.null(modds)) {
    if (!is_triplet(modds)) {
      stop("Format pour modds non reconnu.")
    } 
    if (attraction != "odds"&!quiet) {
      warning("Attention: modds est défini mais ne correspond pas à attraction.")
    } 
    if(attraction == "odds") {
      lodds <- tripletlodds2dgr(modds, cle_from, cle_to)
      jodds <- lodds@j
      podds <- lodds@p
      xodds <- lodds@x
      if (length(xodds) == 0) { attraction <- "constant" }
    }
  }
  
  dist@x <- .meaps_continu(j_dist = dist@j,
                           p_dist = dist@p,
                           x_dist = dist@x,
                           emplois = emplois,
                           actifs = actifs,
                           f = f,
                           shuf = shuf,
                           attraction = attraction,
                           param = param,
                           j_odds = jodds,
                           p_odds = podds,
                           x_odds = xodds,
                           nthreads = nthreads,
                           progress = progress,
                           normalisation = normalisation,
                           fuite_min = fuite_min)
  
  dist <- as(as(as(dist, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  data.frame(i = dist@i + 1L, j = dist@j + 1L, flux = dist@x) |> 
    dplyr::left_join(data.frame(fromidINS = cle_from, i = seq_along(cle_from)), by = "i") |> 
    dplyr::left_join(data.frame(toidINS = cle_to, j = seq_along(cle_to)), by = "j") |> 
    dplyr::select(fromidINS, toidINS, flux)
  
}


