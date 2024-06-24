#' Définition d'une classe d'objet rmeaps de préparation des données.
#' MeapsData permet d'unifier un jeu de données en effectuant un examen de validation, garantissant que les données sont dans les bons formats et le bon ordre avant d'être invoqué dans la procédure rmeaps.
#' @slot triplet Matrice des distances (ou temps ou encore coûts généralisés) entre des départs et des arrivées. La matrice est donnée sous forme de triplet (une data.frame avec les colonnes fromidINS, toidINS et metric). Un tri lexicographique sur fromidINS, puis metric, est attendu.
#' @slot actifs Vecteur des résidents actifs au départ des fromidINS. Le vecteur doit être labélisé par les fromidINS.
#' @slot emplois Vecteurs des emplois à chacune des destinations. Le vecteur doit être labélisé par les toidINS.
#' @slot fuites Vecteurs des proportions d'actifs qui travaillent hors de la zone d'étude. Le vecteur doit être labélisé par les fromidINS.
#' @slot froms Vecteur des fromidINS dans l'ordre attendu par les différentes méthodes (automatique).
#' @slot tos Vecteur des toidINS dans l'ordre attendu par les différentes méthodes (automatique).
#' @slot shuf Une matrice de shuf si besoin de définir des priorités.
#' @slot jr_dist Vecteur des colonnes rangées par distance dans une représentation sparse matrix en ligne.
#' @slot p_dist Vecteur du nombre d'éléments cumulés par ligne dans une représentation sparse matrix en ligne.
setClass("MeapsData",
  representation = list(
    triplet = "data.frame",
    actifs = "numeric",
    emplois = "numeric",
    fuites = "numeric",
    froms = "character",
    tos = "character",
    shuf = "matrix",
    jr_dist = "integer",
    p_dist = "integer"
  ),
  prototype = list(
    triplet = data.frame(fromidINS = character(), toidINS = character(), metric = numeric()),
    actifs = numeric(),
    emplois = numeric(),
    fuites = numeric(),
    froms = character(),
    tos = character(),
    shuf = matrix(),
    jr_dist = integer(),
    p_dist = integer()
  )
)

#' Méthode d'initialisation pour MeapsData.
setMethod(
  f = "initialize", signature = "MeapsData",
  definition = function(
      .Object,
      triplet, actifs, emplois, fuites,
      nshuf = NULL, seuil = 40, 
      seed = NULL, quiet = FALSE) {
    fromidINS <- dplyr::distinct(triplet, fromidINS) |>
      dplyr::pull()
    # le tri est vérifié par ailleurs.
    toidINS <- dplyr::distinct(triplet, toidINS) |>
      dplyr::arrange(toidINS) |>
      dplyr::pull()

    .Object@triplet <- triplet
    .Object@actifs <- actifs[fromidINS]
    .Object@emplois <- emplois[toidINS]
    .Object@fuites <- fuites[fromidINS]
    .Object@froms <- fromidINS
    .Object@tos <- toidINS

    jlab <- seq_along(.Object@tos) - 1L
    names(jlab) <- .Object@tos

    les_j <- jlab[triplet$toidINS]
    names(les_j) <- NULL
    .Object@jr_dist <- les_j

    p_dist <- triplet |>
      group_by(fromidINS) |>
      summarize(n()) |>
      pull() |>
      cumsum()
    .Object@p_dist <- c(0L, p_dist)

    if (!is.null(nshuf)) {
      .Object@shuf <- emiette(
        les_actifs = .Object@actifs,
        nshuf = nshuf, seuil = seuil, seed = seed
      )
    }
    check_meapsdata(.Object, abort = TRUE, quiet = quiet)

    return(.Object)
  }
)

#' Constructeur de MeapsData.
#' @param triplet Un triplet (fromidINS, toidINS, metric) trié lexicographiquement selon fromidINS et metric.
#' @param actifs Un vecteur du nombre d'actifs, labélisé par fromidINS.
#' @param emplois Un vecteur du nombre d'emplois, labélisé par toidINS.
#' @param fuites Un vecteur de la proportion de fuite, labélisé selon fromidINS.
meapsdata <- function(triplet, actifs, emplois, fuites,
                      nshuf = NULL, seuil = 40, seed = NULL, quiet = FALSE) {
  new("MeapsData",
    triplet, actifs, emplois, fuites,
    nshuf = nshuf, seuil = seuil, seed = seed, quiet = quiet
  )
}

#' Méthode implicite show pour MeapsData.
setMethod("show", "MeapsData", function(object) {
  N <- length(object@froms)
  K <- length(object@tos)
  tx <- round(nrow(object@triplet) / N / K * 100, digits = 1)
  cat("MeapsData :\n")
  cat("Matrice des triplet", N, "x", K, "remplie à", tx, "%\n")
  cat("Nombre d'actifs =", sum(object@actifs), "\n")
  cat("Nombre d'emplois =", sum(object@emplois), "\n")
  cat("Nombre de fuyards =", sum(object@fuites * object@actifs))
})

#' Définition d'une sous-classe MeapsDataGroup qui ajoute à MeapsData des données pour regrouper les zones de départs et d'arrivées (notamment par IRIS ou par commune).
#' @slot group_from Le vecteur des codes de regroupement des actifs (à l'iris ou à la commune) labélisé selon les fromidINS.
#' @slot group_to Le vecteur des codes de regroupement des emplois (à l'iris ou à la commune) labélisé selon les toidINS.
#' @slot cible Un triplet (avec les colonnes group_from, group_to et value) des flux attendus.
#' @slot index_gfrom Vecteur des indices correspondant à group_from.
#' @slot index_gto Vecteur des indices correpsondant à group_to.
setClass("MeapsDataGroup",
  representation = list(
    group_from = "character",
    group_to = "character",
    cible = "data.frame",
    index_gfrom = "integer",
    index_gto = "integer"
  ),
  prototype = list(
    group_from = character(),
    group_to = character(),
    cible = data.frame(),
    index_gfrom = integer(),
    index_gto = integer()
  ),
  contains = "MeapsData"
)

#' Méthode d'initialisation pour MeapsDataGroup.
setMethod(
  f = "initialize", signature = "MeapsDataGroup",
  definition = function(.Object, triplet, actifs, emplois, fuites, shuf,
                        froms, tos, group_from, group_to, cible, jr_dist, p_dist, quiet = FALSE) {
    if (!is.null(cible)) {
      cible <- cible |>
        complete(group_from, group_to, fill = list(value = 0)) |>
        arrange(group_from, group_to)
    }
    .Object@triplet <- triplet
    .Object@actifs <- actifs[froms]
    .Object@emplois <- emplois[tos]
    .Object@fuites <- fuites[froms]
    .Object@froms <- froms
    .Object@tos <- tos
    .Object@group_from <- group_from[froms]
    .Object@group_to <- group_to[tos]
    .Object@cible <- cible
    .Object@shuf <- shuf
    .Object@jr_dist <- jr_dist
    .Object@p_dist <- p_dist

    g_froms <- unique(group_from) |> sort()
    g_tos <- unique(group_to) |> sort()

    gf_lab <- seq_along(g_froms) - 1L
    names(gf_lab) <- g_froms
    gt_lab <- seq_along(g_tos) - 1L
    names(gt_lab) <- g_tos

    .Object@index_gfrom <- gf_lab[group_from]
    .Object@index_gto <- gt_lab[group_to]

    check_meapsdata(.Object, abort = TRUE, quiet = quiet )
    check_meapsdatagroup(.Object, abort = TRUE, quiet = quiet)

    return(.Object)
  }
)

#' Méthode implicite show pour MeapsDataGroup.
setMethod("show", "MeapsDataGroup", function(object) {
  N <- length(object@froms)
  K <- length(object@tos)
  Ng <- unique(object@group_from) |> length()
  Kg <- unique(object@group_to) |> length()
  tx <- round(nrow(object@triplet) / N / K * 100, digits = 1)
  cat("MeapsDataGroup :\n")
  cat("Matrice des triplet", N, "x", K, "remplie à", tx, "%\n")
  cat("Nombre d'actifs =", sum(object@actifs), "\n")
  cat("Nombre d'emplois =", sum(object@emplois), "\n")
  cat("Nombre de fuyards =", sum(object@fuites * object@actifs), "\n")
  cat("Nombre de groupes de départ =", Ng, "\n")
  cat("Nombre de groupes d'arrivée =", Kg, "\n")
  cat("Flux observés (cible) =", sum(object@cible$value))
})

#' Constructeur de MeapsDataGroup.
#' @param MeapsData Un objet MeapsData qui décrit la zone d'étude.
#' @param group_from Un vecteur des codes de regroupement, labélisé par fromidINS.
#' @param group_to Un vecteur des codes de regroupement, labélisé par toidINS.
#' @param cible Un triplet (group_from, group_to, value) décrivant les flux groupés de référence.
#'
#' @import dplyr
meapsdatagroup <- function(MeapsData, group_from, group_to, cible, quiet=FALSE) {
  new(
    "MeapsDataGroup", MeapsData@triplet, MeapsData@actifs,
    MeapsData@emplois, MeapsData@fuites, MeapsData@shuf,
    MeapsData@froms, MeapsData@tos,
    group_from[MeapsData@froms], group_to[MeapsData@tos],
    cible,
    MeapsData@jr_dist, MeapsData@p_dist, quiet
  )
}


