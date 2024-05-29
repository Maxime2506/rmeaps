#' Définition d'une classe d'objet rmeaps de préparation des données.
#' MeapsData permet d'unifier un jeu de données en effectuant un examen de validation, garantissant que les données sont dans les bons formats et le bon ordre avant d'être invoqué dans la procédure rmeaps.
#' @slot triplet Matrice des distances (ou temps ou encore coûts généralisés) entre des départs et des arrivées. La matrice est donnée sous forme de triplet (une data.frame avec les colonnes fromidINS, toidINS et metric). Un tri lexicographique sur fromidINS, puis metric, est attendu.
#' @slot actifs Vecteur des résidents actifs au départ des fromidINS. Le vecteur doit être labélisé par les fromidINS.
#' @slot emplois Vecteurs des emplois à chacune des destinations. Le vecteur doit être labélisé par les toidINS.
#' @slot fuites Vecteurs des proportions d'actifs qui travaillent hors de la zone d'étude. Le vecteur doit être labélisé par les fromidINS.
#' @slot froms Vecteur des fromidINS dans l'ordre attendu par les différentes méthodes (automatique).
#' @slot tos Vecteur des toidINS dans l'ordre attendu par les différentes méthodes (automatique).
setClass("MeapsData", 
         representation = list(
           triplet = "data.frame",
           actifs = "numeric",
           emplois = "numeric",
           fuites = "numeric", 
           froms = "character",
           tos = "character"
         ),
         prototype = list(
           triplet = data.frame(fromidINS = character(), toidINS = character(), metric = numeric()),
           actifs = numeric(),
           emplois = numeric(),
           fuites = numeric(),
           froms = character(),
           tos = character()
         )
)

#' Méthode d'initialisation pour MeapsData.
setMethod(f = "initialize", signature = "MeapsData",
          definition = function(.Object, triplet, actifs, emplois, fuites) {
            fromidINS <- dplyr::distinct(triplet, fromidINS) |> 
              dplyr::pull()
            toidINS <- dplyr::distinct(triplet, toidINS) |>
              dplyr::arrange(toidINS) |> 
              dplyr::pull() |> 
              sort()
            .Object@triplet <- triplet
            .Object@actifs <- actifs[fromidINS]
            .Object@emplois <- emplois[toidINS]
            .Object@fuites <- fuites[fromidINS]
            .Object@froms <- fromidINS
            .Object@tos <- toidINS
            
            check_meapsdata(.Object, abort = TRUE)
            
            return(.Object)
          })

#' Constructeur de MeapsData.
#' @param triplet Un triplet (fromidINS, toidINS, metric) trié lexicographiquement selon fromidINS et metric.
#' @param actifs Un vecteur du nombre d'actifs, labélisé par fromidINS.
#' @param emplois Un vecteur du nombre d'emplois, labélisé par toidINS.
#' @param fuites Un vecteur de la proportion de fuite, labélisé selon fromidINS.
meapsdata <- function(triplet, actifs, emplois, fuites) {
  new("MeapsData", triplet, actifs, emplois, fuites)
}

#' Méthode implicite show pour MeapsData.
setMethod("show", "MeapsData", function(object) {
  N <- length(object@froms)
  K <- length(object@tos)
  tx <- round(nrow(object@triplet) / N / K * 100, digits = 1)
  cat("MeapsData :\n")
  cat("Matrice des triplet", N, "x", K, "remplie à", tx, "%\n")
  cat("Nombre d'actifs =", sum(object@actifs),"\n")
  cat("Nombre d'emplois =", sum(object@emplois), "\n")
  cat("Nombre de fuyards =", sum(object@fuites * object@actifs))
})

#' Fonction all_in. Applique rmeaps dans sa version "tous les actifs partent en même temps".
#' @param MeapsData Un jeu de données MeapsData décrivant la zone d'étude.
#' @param attraction Choix du mode d'attraction. Par défaut, "constant". Les autres modes possibles sont : "marche", "marche_liss", "decay" et "logistique". Chacune applique une pénalité selon la distance (au-delà de la simple logique des rangs, qui correspond au mode "constant").
#' @param parametres Un vecteur avec les valeurs initiales des paramètres. Dépend du mode d'attraction retenu.
#' @param nthreads Le nombre de threads pour mener les calculs parallélisés. Par defaut, le nombre maximal.
#' @param verbose TRUE par defaut.
#' 
#' @import dplyr
all_in <- function(MeapsData, attraction = "constant", parametres = 0, nthreads = 0L, verbose = TRUE) {
  
  if (!inherits(MeapsData, "MeapsData")) cli::cli_abort("Ce n'est pas un objet MeapsData.") 
  # Validation des paramètres
  if (!attraction %in% c("constant", "marche", "marche_liss", "decay", "logistique")) cli::cli_abort("Fonction attraction inconnue")
  
  if (attraction == "marche" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour marche invalide.")
  if (attraction == "marche_liss" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour marche_liss invalide.")
  if (attraction == "decay" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour decay invalide.")
  if (attraction == "logistique" && (length(parametres) != 3 || !is.numeric(parametres))) cli::cli_abort("Parametres pour logistique invalide.")
  
  froms <- MeapsData@froms
  tos <- MeapsData@tos
  
  jlab <- seq_along(tos) - 1L
  names(jlab) <- tos
  
  les_j <- jlab[MeapsData@triplet$toidINS]
  
  p_dist <- MeapsData@triplet |> group_by(fromidINS) |> summarize(n()) |> pull() |> cumsum()
  p_dist <- c(0L, p_dist)
  
  meapsclass(jr_dist = les_j,
               p_dist = p_dist,
               xr_dist = MeapsData@triplet$metric,
               emplois = MeapsData@emplois,
               actifs = MeapsData@actifs,
               fuites = MeapsData@fuites,
               parametres = parametres,
               attraction = attraction,
               nthreads = nthreads, verbose = verbose) |> 
    as.data.frame() |> 
    dplyr::left_join(data.frame(i = seq_along(froms) - 1L, fromidINS = froms), by = "i") |>
    dplyr::left_join(data.frame(j = seq_along(tos) - 1L, toidINS = tos), by = "j") |>
    dplyr::select(fromidINS, toidINS, flux)
  
}

#' Définition d'une sous-classe MeapsDataGroup qui ajoute à MeapsData des données pour regrouper les zones de départs et d'arrivées (notamment par IRIS ou par commune). 
#' @slot group_from Le vecteur des codes de regroupement des actifs (à l'iris ou à la commune) labélisé selon les fromidINS.
#' @slot group_to Le vecteur des codes de regroupement des emplois (à l'iris ou à la commune) labélisé selon les toidINS.
#' @slot cible Un triplet (avec les colonnes group_from, group_to et value) des flux attendus.
setClass("MeapsDataGroup", 
         representation = list(
           group_from = "character",
           group_to = "character",
           cible = "data.frame"
         ),
         prototype = list(
           group_from = character(),
           group_to = character(),
           cible = data.frame()
         ),
         contains = "MeapsData")

#' Méthode d'initialisation pour MeapsDataGroup.
setMethod(f = "initialize", signature = "MeapsDataGroup",
          definition = function(.Object, triplet, actifs, emplois, fuites,
                                froms, tos, group_from, group_to, cible) {
            froms <- sort(froms)
            tos <- sort(tos)
            .Object@triplet <- triplet
            .Object@actifs <- actifs[froms]
            .Object@emplois <- emplois[tos]
            .Object@fuites <- fuites[froms]
            .Object@froms <- froms
            .Object@tos <- tos
            .Object@group_from <- group_from[froms]
            .Object@group_to <- group_to[tos]
            .Object@cible <- cible |> dplyr::arrange(group_from, group_to)
            
            check_meapsdata(.Object, abort = TRUE)
            check_meapsdatagroup(.Object, abort = TRUE)
            
            return(.Object)
          })

#' Methode implicite show pour MeapsDataGroup.
setMethod("show", "MeapsDataGroup", function(object) {
  N <- length(object@froms)
  K <- length(object@tos)
  Ng <- unique(object@group_from) |> length()
  Kg <- unique(object@group_to) |> length()
  tx <- round(nrow(object@triplet) / N / K * 100, digits = 1)
  cat("MeapsDataGroup :\n")
  cat("Matrice des triplet", N, "x", K, "remplie à", tx, "%\n")
  cat("Nombre d'actifs =", sum(object@actifs),"\n")
  cat("Nombre d'emplois =", sum(object@emplois), "\n")
  cat("Nombre de fuyards =", sum(object@fuites * object@actifs), "\n")
  cat("Nombre de groupes de départ =", Ng,"\n")
  cat("Nombre de groupes d'arrivée =", Kg, "\n")
  cat("Flux observés (cible) =", sum(object@cible$value))
})

#' Constructeur de MeapsDataGroup.
#' @param MeapsData Un objet MeapsData qui décrit la zone d'étude.
#' @param group_from Un vecteur des codes de regroupement, labélisé par fromidINS.
#' @param group_to Un vecteur des codes de regroupement, labélisé par toidINS.
#' @param cible Un triplet (group_from, group_to, value) décrivant les flux groupés de référence.
#' 
#'@import dplyr
meapsdatagroup <- function(MeapsData, group_from, group_to, cible) {
  
  new("MeapsDataGroup", MeapsData@triplet, MeapsData@actifs, MeapsData@emplois, MeapsData@fuites,
      MeapsData@froms, MeapsData@tos,
      group_from[MeapsData@froms], group_to[MeapsData@tos],
      cible |> dplyr::arrange(group_from, group_to)
  )
}

#' Méthode de préparation des données de MeapsDataGroup pour rmeaps.
.prep_grouped <- function(MeapsDataGroup,  attraction = "constant", parametres = 0, nthreads = 0L, verbose = TRUE) {
  
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.") 
  # Validation des paramètres
  if (!attraction %in% c("constant", "marche", "marche_liss", "decay", "logistique")) cli::cli_abort("Fonction attraction inconnue")
  
  if (attraction == "marche" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour marche invalide.")
  if (attraction == "double_marche_liss" && (length(parametres) != 4 || !is.numeric(parametres))) cli::cli_abort("Parametres pour double_marche_liss invalide.")
  if (attraction == "decay" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour decay invalide.")
  if (attraction == "logistique" && (length(parametres) != 3 || !is.numeric(parametres))) cli::cli_abort("Parametres pour logistique invalide.")
  
  tos <- MeapsDataGroup@tos
  
  jlab <- seq_along(tos) - 1L
  names(jlab) <- tos
  
  les_j <- jlab[MeapsDataGroup@triplet$toidINS]
  
  p_dist <- MeapsDataGroup@triplet |> group_by(fromidINS) |> summarize(n()) |> pull() |> cumsum()
  p_dist <- c(0L, p_dist)
  
  g_froms <- unique(MeapsDataGroup@group_from) |> sort()
  g_tos <- unique(MeapsDataGroup@group_to) |> sort()
  
  gf_lab <- seq_along(g_froms) - 1L
  names(gf_lab) <- g_froms
  gt_lab <- seq_along(g_tos) - 1L
  names(gt_lab) <- g_tos
  
  index_gfrom <- gf_lab[MeapsDataGroup@group_from]
  index_gto <- gt_lab[MeapsDataGroup@group_to]
  
  index_value <- expand_grid(data.frame(index_from = gf_lab, group_from = names(gf_lab)),
                             data.frame(index_to = gt_lab, group_to = names(gt_lab))) |> 
    left_join(cible, by = c("group_from", "group_to")) |> 
    transmute(value = replace_na(value, 0)) |> 
    pull(value)
  
  return(list(jr_dist = les_j,
              p_dist = p_dist,
              xr_dist = MeapsDataGroup@triplet$metric,
              emplois = MeapsDataGroup@emplois,
              actifs = MeapsDataGroup@actifs,
              fuites = MeapsDataGroup@fuites,
              parametres = parametres,
              attraction = attraction,
              group_from = index_gfrom,
              group_to = index_gto,
              cible = index_value,
              nthreads = nthreads, verbose = verbose))
}

#' Fonction all_in_grouped applique meaps à un objet MeapsDataGroup et renvoie les flux estimés regroupés.
#' @import dplyr
all_in_grouped <- function(MeapsDataGroup,  attraction = "constant", parametres = 0, nthreads = 0L, verbose = TRUE) {
  
  les_args <- .prep_grouped(MeapsDataGroup = MeapsDataGroup, attraction = attraction, parametres = parametres, 
                nthreads = nthreads, verbose = verbose)
  
  do.call(meapsclass, les_args)
  
}

meaps_optim <- function(MeapsDataGroup,  attraction, parametres, 
                        method = "L-BFGS-B", objective = "KL", lower = NULL, upper = NULL,
                        nthreads = 0L, progress = TRUE, control = NULL) { 
  
  les_args <- .prep_grouped(MeapsDataGroup = MeapsDataGroup, attraction = attraction, parametres = parametres, 
                        nthreads = nthreads, verbose = FALSE)
  
  les_args$parametres <- NULL
  env <- environment()
  fn <- switch(
    objective, 
    "KL" = function(par) {
      if (progress) cli::cli_progress_update(.envir = env, status = par)
      do.call(meapsclass, args = append(les_args, list(parametres = par)))$kl
    }
  )
  
  if (is.null(fn)) stop("Fonction objective non définie !")
  
  nb_par <- length(parametres)
  if (is.null(lower)) lower <- rep(0, nb_par)
  if (is.null(upper)) upper <- rep(Inf, nb_par)
  
  cli::cli_progress_bar(.envir = env, clear = FALSE)
  stats::optim(par = parametres, fn = fn, method = method, lower = lower, upper = upper, control = control)
  cli::cli_progress_done(.envir = env)
}


multishuf <- function(MeapsData, nshuf = 64, seuil = 40, attraction = "constant", parametres = 0, nthreads = 0L, verbose = TRUE) {
  
  shuf <- emiette(MeapsData@actifs, nshuf = nshuf, seuil = seuil)
  paste("seed du shuf =", attr(shuf, "seed")) |> print()
  
  if (!inherits(MeapsData, "MeapsData")) cli::cli_abort("Ce n'est pas un objet MeapsData.") 
  # Validation des paramètres
  if (!attraction %in% c("constant", "marche", "marche_liss", "decay", "logistique")) cli::cli_abort("Fonction attraction inconnue")
  
  if (attraction == "marche" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour marche invalide.")
  if (attraction == "marche_liss" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour marche_liss invalide.")
  if (attraction == "decay" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour decay invalide.")
  if (attraction == "logistique" && (length(parametres) != 3 || !is.numeric(parametres))) cli::cli_abort("Parametres pour logistique invalide.")
  
  froms <- MeapsData@froms
  tos <- MeapsData@tos
  
  jlab <- seq_along(tos) - 1L
  names(jlab) <- tos
  
  les_j <- jlab[MeapsData@triplet$toidINS]
  
  p_dist <- MeapsData@triplet |> group_by(fromidINS) |> summarize(n()) |> pull() |> cumsum()
  p_dist <- c(0L, p_dist)
  
  newmultishuf(jr_dist = les_j,
             p_dist = p_dist,
             xr_dist = MeapsData@triplet$metric,
             emplois = MeapsData@emplois,
             actifs = MeapsData@actifs,
             fuites = MeapsData@fuites,
             parametres = parametres,
             shuf = shuf,
             attraction = attraction,
             nthreads = nthreads, verbose = verbose) |> 
    as.data.frame() |> 
    dplyr::left_join(data.frame(i = seq_along(froms) - 1L, fromidINS = froms), by = "i") |>
    dplyr::left_join(data.frame(j = seq_along(tos) - 1L, toidINS = tos), by = "j") |>
    dplyr::select(fromidINS, toidINS, flux)
  
}


