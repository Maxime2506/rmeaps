#' La fonction prep_meaps_dist prépare les données de distances pour traitement par meaps_optim
#' @param dist triplet des distances où des résidents en ligne i rejoignent des opportunités en colonnes j. 
#' 
#' @return renvoie une RankedRSMatrix des distances.
#' @import Matrix
#' @export
prep_meaps_dist <- function(dist, emplois, actifs, fuite, shuf, groups_from, groups_to) {
  
  if (!is_triplet(dist)) stop("Ce n'est pas un triplet valide.")
  dist <- triplet2listij(dist)
  froms <- as.character(dist$cle_from)
  tos <- as.character(dist$cle_to)
  
  row_group <- groups_from[froms]
  gg <- unique(row_group)
  gg <- rlang::set_names(1:length(gg), gg)
  row_group <- gg[as.character(row_group)] - 1L
  col_group <- groups_to[tos]
  gg <- unique(col_group)
  gg <- rlang::set_names(1:length(gg), gg)
  col_group <- gg[as.character(col_group)] -1L
  actifs <- actifs[froms]
  fuite <- fuite[froms]
  emplois <- emplois[tos]
  list(RankedMat = new(RankedRSMatrix, dist$dgr), 
       cle_from = dist$cle_from, 
       cle_to = dist$cle_to,
       emplois = emplois,
       actifs = actifs,
       fuite = fuite,
       shuf = reordonne_shuf(shuf, actifs),
       row_group = row_group,
       col_group = col_group) 
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

#' La fonction prep_meaps_odds_on_dist prépare la matrice des odds pour traitement par meaps_optim.
#' c'est une matrice ranké, mais selon les distances
#' @param odds un triplet des odds.
#' @param dist un triplet des distances.
#' 
#' @return renvoie une RankedRSMatrix des odds, rangés selon le rang des distances.
#' @import Matrix, data.table
#' @export
prep_meaps_odds_on_dist <- function(odds, prep) {
  if (!is_triplet(odds)) stop("Ce n'est pas un triplet d'odds valide.")
  colnames(odds)[!colnames(odds) %in% c("fromidINS", "toidINS")] <- "o"
  setDT(odds)
  lodds <- data.table(toidINS = prep$cle_to, j = 1:length(prep$cle_to)) |> 
    merge(odds, by = c("toidINS"), all.x = TRUE, all.y = FALSE)  
  lodds <- data.table(fromidINS = prep$cle_from, i = 1:length(prep$cle_from)) |> 
    merge(lodds, by = c("fromidINS"), all.x = TRUE, all.y = FALSE)
  oo <- set_names(lodds$o, lodds$j)
  p <- prep$RankedMat$p
  roo <- vector("numeric", length(prep$RankedMat$xr))
  for (i_p in 1:(length(p)-1)) {
    d <- p[i_p]+1
    f <- p[i_p+1]
    js <- as.character(prep$RankedMat$jr[d:f]+1)
    roo[d:f] <- oo[d:f][js]
  }
  res <- new(RankedRSMatrix, roo, prep$RankedMat$jr, p, prep$RankedMat$dim)
  return(res)
}


#' La fonction prep_meaps_odds_on_dist prépare la matrice des odds pour traitement par meaps_optim.
#' c'est une matrice ranké, mais selon les distances
#' @param prep une liste issue de prep_meaps_dist.
#' 
#' @return renvoie une RankedRSMatrix des odds, rangés selon le rang des distances.
#' @import Matrix, data.table
#' @export
prep_0lodds_on_dist <- function(prep) {
  roo <- numeric(length(prep$RankedMat$xr))
  res <- new(RankedRSMatrix, roo, prep$RankedMat$jr, p, prep$RankedMat$dim)
  return(res)
}


#' La fonction meaps_optim sert pour la recherche des meilleurs paramètres ou odds lorsqu'un certain regroupement des flux est connu.
#' @param prep liste, issue de prep_meaps_dist, qui contient plusieurs éléments précalculés.
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
meaps_optim <- function(prep,
                        attraction = "constant",
                        param = 0.0,
                        odds_prep = NULL,
                        nthreads = 0,
                        progress = TRUE,
                        normalisation = FALSE,
                        fuite_min = 1e-3) {
  actifs <- prep$actifs
  fuite <- prep$fuite
  emplois <- prep$emplois
  shuf <- prep$shuf
  mat <- prep$RankedMat
  jr_odds <- p_odds <- xr_odds <- 1L
  
  # contraintes sur les paramètres.
  if (attraction %in% c("marche", "marche_liss")) {
    if (param[1] <= 0) param[[1]] <- 0
    if (param[2] <= 0) param[[2]] <- 0
  }
  
  if (attraction == "logistique") {
    if (param[1] <= 0) param[[1]] <- 0
    if (param[2] <= 0) param[[2]] <- 0
    if (param[3] <= 0) param[[3]] <- 0
  }
  
  if (!is.null(odds_prep)) {
    attraction <- "odds"
    jr_odds <- odds_prep$jr
    p_odds <- odds_prep$p
    xr_odds <- odds_prep$xr
  } else {
    jr_odds <- 0L
    p_odds <- 0L
    xr_odds <- 0.0
  }
  res <- .meaps_optim(jr_dist = mat$jr,
                      p_dist = mat$p,
                      xr_dist = mat$xr,
                      emplois = emplois,
                      actifs = actifs,
                      f = fuite,
                      shuf = shuf,
                      row_group = prep$row_group,
                      col_group = prep$col_group,
                      param = param,
                      attraction = attraction,
                      jr_odds = jr_odds,
                      p_odds = p_odds,
                      xr_odds = xr_odds,
                      nthreads = nthreads,
                      progress = FALSE,
                      normalisation = normalisation,
                      fuite_min = fuite_min)
  
  coms <- tibble::tibble(COMMUNE = names(prep$row_group), ic = prep$row_group) |> 
    dplyr::group_by(COMMUNE) |> 
    dplyr::summarise(ic = as.character(dplyr::first(ic)+1))
  dclts <- tibble(DCLT = names(prep$col_group), id = prep$col_group) |> 
    dplyr::group_by(DCLT) |> 
    dplyr::summarise(id = as.character(dplyr::first(id)+1))
  colnames(res) <- stringr::str_c("id", 1:ncol(res))
  res |> 
    tibble::as_tibble(.name_repair = "unique") |> 
    dplyr::mutate(ic = as.character(1:nrow(res))) |> 
    tidyr::pivot_longer(cols = -ic, names_to = "id", values_to = "flux") |> 
    dplyr::filter(flux>0) |> 
    dplyr::mutate(id = stringr::str_sub(id,3,-1)) |>
    dplyr::left_join(coms, by = "ic") |>
    dplyr::left_join(dclts, by = "id") |>
    dplyr::select(-id,-ic) |>
    dplyr::arrange(dplyr::desc(flux))
}


