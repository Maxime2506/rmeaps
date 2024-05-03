#' Définition d'une classe d'objet rmeaps de préparation des données.
#' 
#' 
#' 
#' 
is_triplet_meaps <- function(object, quiet = TRUE) {

  if (!inherits(object, "data.frame") || length(object) != 3) {
    status <- 2L
  } else if (!setequal(names(object), c("fromidINS", "toidINS", "metric"))) {
      status <- 3L
  } else if (any(is.na(object))) { 
      status <- 4L 
  } else {
      status <- 1L
  }
  
  if (!quiet) { 
    switch(status,
         cli::cli_alert_success("Le triplet meaps est un triplet valide."),
         cli::cli_alert_warning("Le triplet meaps n'est pas un triplet."),
         cli::cli_alert_warning("Les noms des colonnes devraient être fromidINS, toidINS et metric."),
         cli::cli_alert_warning("Le triplet meaps contient des valeurs manquantes.")
    )
  }
  
  invisible(status==1)
}

is_triplet_ordered <- function(object, quiet = TRUE) {
  if ( is.unsorted(object$fromidINS) ) { 
    status <- 2L
  } else {
    sort_v <- by(object, object$fromidINS, \(x) is.unsorted(x$metric)) |> any()
    if (sort_v) {
      status <- 3L
    } else {
      sort_j <- by(object, list(object$fromidINS, object$metric), \(x) { is.unsorted(x$toidINS, strictly = TRUE) }) |> any(na.rm = TRUE)
      if (sort_j) { 
        status <- 4L 
      } else {
          status <- 1L
        }
    }
  }
  
  if (!quiet) { 
    switch(status,
         cli::cli_alert_success("Le triplet meaps est bien ordonné."),
         cli::cli_alert_warning("Le triplet meaps n'est pas ordonné selon fromidINS."),
         cli::cli_alert_warning("Le triplet meaps n'est pas ordonné selon metric"),
         cli::cli_alert_warning("Le triplet meaps n'est pas ordonné selon toidINS.")
         )
  }
  
  invisible(status==1)
}

is_meapsdata_valid <- function(object, quiet = TRUE) {
  status <- TRUE
  if (is.null(names(object@actifs))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("actifs n'a pas de labels") }
    }
  if (is.null(names(object@emplois))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("emplois n'a pas de labels") }
  }
  if (is.null(names(object@fuites))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("fuites n'a pas de labels") }
  }
  
  N <- unique(object@triplet$fromidINS) |> length()
  K <- unique(object@triplet$toidINS) |> length()
  
  if (!setequal(object@triplet$fromidINS, names(object@actifs)) || length(object@actifs) != N) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("triplet et actifs n'ont pas la même liste de fromidINS") }
  }
  if (!setequal(object@triplet$toidINS, names(object@emplois)) || length(object@emplois) != K) { 
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("triplet et emplois n'ont pas la même liste de toidINS") }
  }
  if (!setequal(object@triplet$fromidINS, names(object@fuites)) || length(object@fuites) != N) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("triplet et fuites n'ont pas la même liste de fromidINS") }
  }
  
  if (any(object@fuites <= 0 | object@fuites >= 1)) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("Le vecteur fuite a des valeurs invalides.") }
  }
  
  if (status) { cli::cli_alert_success("Le MeapsData est valide.") }
  
  invisible(status)
}

is_meaps_normalized <- function(object, just_a_warning = TRUE, seuil = 1e-4) {
  delta <- (sum(object@actifs * (1 - object@fuites)) - sum(object@emplois))/sum(object@emplois)
  if (delta > seuil) {
    if (just_a_warning) { 
      cli::cli_alert_warning("Les actifs restant dans la zone et les emplois divergent de {round(100 * delta, digits = 2)}%.")
    } else {
      cli::cli_abort("Les actifs restant dans la zone et les emplois divergent de {round(100 * delta, digits = 2)}%.")
    }
  }
  invisible(delta <= seuil)
}

check_meapsdata <- function(object, abort = FALSE) {
  status <- is_triplet_meaps(object@triplet, quiet = FALSE) && 
    is_triplet_ordered(object@triplet, quiet = FALSE) && 
    is_meapsdata_valid(object, quiet = FALSE) &&
    is_meaps_normalized(object)
  
  if (abort && !status) {cli::cli_abort("Invalide.")} else { return(status)}
}

is_meapsdatagroup_valid <- function(object, quiet = FALSE)  {
  
  if (!inherits(object, "data.frame") || length(object) != 3) {
    status <- 2L
  } else if (!setequal(names(object), c("group_from", "group_to", "value"))) {
    status <- 3L
  } else if (any(is.na(object))) { 
    status <- 4L 
  } else {
    status <- 1L
  }
  
  if (!quiet) { 
    switch(status,
           cli::cli_alert_success("Le triplet cible est un triplet valide."),
           cli::cli_alert_warning("La cible n'est pas un triplet."),
           cli::cli_alert_warning("Les noms des colonnes devraient être group_from, group_to et value."),
           cli::cli_alert_warning("La cible contient des valeurs manquantes.")
    )
  }
  
  invisible(status==1)
}

check_meapsdatagroup <- function(object, abort = FALSE){
  status <- TRUE
  if (length(object@group_from) != length(object@actifs)) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("group_from et actifs ne font pas la même longueur.") }
  }
  if (length(object@group_to) != length(object@emplois)) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("group_to et emplois ne font pas la même longueur.") }
  }
  
  if (is.null(names(object@group_from))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("group_from n'est pas labelisé") }
  }
  if (is.null(names(object@group_to))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("group_to n'est pas labelisé.") }
  }
  
  if (!setequal(names(object@group_from), names(object@actifs))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("group_from et actifs n'ont pas la même liste de fromidINS.") }
  }
  if (!setequal(names(object@group_to), names(object@emplois))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("group_to et emplois n'ont pas la même liste de toidINS.") }
  }
  
  if (!setequal(object@cible$group_from, object@group_from)) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("group_from et cible ne correspondent pas.") }
  }
  
  if (!setequal(object@cible$group_to, object@group_to)) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("group_to et cible ne correspondent pas.") }
  }
  
  check_meapsdata(object)
}


# Class de base pour les données meaps.
setClass("MeapsData", 
         representation = list(
           triplet = "data.frame",
           actifs = "numeric",
           emplois = "numeric",
           fuites = "numeric"
         ),
         prototype = list(
           triplet = data.frame(fromidINS = character(), toidINS = character(), metric = numeric()),
           actifs = numeric(),
           emplois = numeric(),
           fuites = numeric()
         )
)

setMethod(f = "initialize", signature = "MeapsData",
          definition = function(.Object, triplet, actifs, emplois, fuites) {
            
            fromidINS <- unique(triplet$fromidINS) # le tri est vérifié par ailleurs.
            toidINS <- unique(triplet$toidINS) |> sort()
            .Object@triplet <- triplet
            .Object@actifs <- actifs[fromidINS]
            .Object@emplois <- emplois[toidINS]
            .Object@fuites <- fuites[fromidINS]
            
            check_meapsdata(.Object, abort = TRUE)
            
            return(.Object)
          })

setMethod("show", "MeapsData", function(object) {
  N <- unique(object@triplet$fromidINS) |> length()
  K <- unique(object@triplet$toidINS) |> length()
  tx <- round(nrow(object@triplet) / N / K * 100, digits = 1)
  cat("MeapsData :\n")
  cat("Matrice des triplet", N, "x", K, "remplie à", tx, "%\n")
  cat("Nombre d'actifs =", sum(object@actifs),"\n")
  cat("Nombre d'emplois =", sum(object@emplois), "\n")
  cat("Nombre de fuyards =", sum(object@fuites * object@actifs))
})

#' @import dplyr
all_in <- function(MeapsData, attraction = "constant", parametres = 0, odds = 1, 
                   nthreads = 0L, verbose = TRUE) {
  
  if (!inherits(MeapsData, "MeapsData")) cli::cli_abort("Ce n'est pas un objet MeapsData.") 
  # Validation des paramètres
  if (!attraction %in% c("constant", "marche", "marche_liss", "double_marche_liss","decay", "logistique")) cli::cli_abort("Fonction attraction inconnue")
 
  # RQ : pas de méthode pour vérifier le bon ordre des odds.
  if (attraction == "odds" && length(odds) != nrow(object@triplet) ) cli::cli_abort("vecteur odds invalide.") 
  if (attraction == "marche" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour marche invalide.")
  if (attraction == "marche_liss" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour marche_liss invalide.")
  if (attraction == "double_marche_liss" && (length(parametres) != 4 || !is.numeric(parametres))) cli::cli_abort("Parametres pour double_marche_liss invalide.")
  if (attraction == "decay" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour decay invalide.")
  if (attraction == "logistique" && (length(parametres) != 3 || !is.numeric(parametres))) cli::cli_abort("Parametres pour logistique invalide.")
  
  froms <- unique(MeapsData@triplet$fromidINS)
  tos <- unique(MeapsData@triplet$toidINS) |> sort()
  
  jlab <- seq_along(tos) - 1L
  names(jlab) <- tos
  
  les_j <- jlab[MeapsData@triplet$toidINS]
  
  p_dist <- MeapsData@triplet |> group_by(fromidINS) |> summarize(n()) |> pull() |> cumsum()
  p_dist <- c(0L, p_dist)
  
  meaps_all_in(jr_dist = les_j,
               p_dist = p_dist,
               xr_dist = MeapsData@triplet$metric,
               emplois = MeapsData@emplois,
               actifs = MeapsData@actifs,
               fuites = MeapsData@fuites,
               parametres = parametres,
               xr_odds = odds,
               attraction = attraction,
               nthreads = nthreads, verbose = verbose) |>
    dplyr::left_join(data.frame(i = seq_along(froms) - 1L, fromidINS = froms), by = "i") |>
    dplyr::left_join(data.frame(j = seq_along(tos) - 1L, toidINS = tos), by = "j") |>
    dplyr::select(fromidINS, toidINS, flux)
  
}

# Class pour les regroupements
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


setMethod(f = "initialize", signature = "MeapsDataGroup",
          definition = function(.Object, triplet, actifs, emplois, fuites, group_from, group_to, cible) {
            
            fromidINS <- unique(triplet$fromidINS) # le tri est vérifié par ailleurs.
            toidINS <- unique(triplet$toidINS) |> sort()
            
            .Object@triplet <- triplet
            .Object@actifs <- actifs[fromidINS]
            .Object@emplois <- emplois[toidINS]
            .Object@fuites <- fuites[fromidINS]
            .Object@group_from <- group_from[fromidINS]
            .Object@group_to <- group_to[toidINS]
            .Object@cible <- cible |> arrange(group_from, group_to)
            
            check_meapsdatagroup(.Object)
            
            return(.Object)
          })

setMethod("show", "MeapsDataGroup", function(object) {
  N <- unique(object@triplet$fromidINS) |> length()
  K <- unique(object@triplet$toidINS) |> length()
  Ng <- unique(object@group_from) |> length()
  Kg <- unique(object@group_to) |> length()
  tx <- round(nrow(object@triplet) / N / K * 100, digits = 1)
  cat("MeapsData :\n")
  cat("Matrice des triplet", N, "x", K, "remplie à", tx, "%\n")
  cat("Nombre d'actifs =", sum(object@actifs),"\n")
  cat("Nombre d'emplois =", sum(object@emplois), "\n")
  cat("Nombre de fuyards =", sum(object@fuites * object@actifs), "\n")
  cat("Nombre de groupes de départ =", Ng,"\n")
  cat("Nombre de groupes d'arrivée =", Kg, "\n")
  cat("Flux observés (cible) =", sum(object@cible$value))
})

# Constructeur
meapsdatagroup <- function(MeapsData, group_from, group_to, cible) {
  
  fromidINS <- unique(MeapsData@triplet$fromidINS) # le tri est vérifié par ailleurs.
  toidINS <- unique(MeapsData@triplet$toidINS) |> sort()
  
  new("MeapsDataGroup", MeapsData@triplet, MeapsData@actifs, MeapsData@emplois, MeapsData@fuites,
      group_from[fromidINS], group_to[toidINS], cible |> arrange(group_from, group_to)
      )
}

#' @import dplyr
all_in_grouped <- function(MeapsDataGroup,  attraction = "constant", parametres = 0, odds = 1, 
                           nthreads = 0L, verbose = TRUE) {
  
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.") 
  # Validation des paramètres
  if (!attraction %in% c("constant", "marche", "marche_liss", "double_marche_liss","decay", "logistique")) cli::cli_abort("Fonction attraction inconnue")
  
  # RQ : pas de méthode pour vérifier le bon ordre des odds.
  if (attraction == "odds" && length(odds) != nrow(MeapsDataGroup@triplet) ) cli::cli_abort("vecteur odds invalide.") 
  if (attraction == "marche" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour marche invalide.")
  if (attraction == "marche_liss" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour marche_liss invalide.")
  if (attraction == "double_marche_liss" && (length(parametres) != 4 || !is.numeric(parametres))) cli::cli_abort("Parametres pour double_marche_liss invalide.")
  if (attraction == "decay" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour decay invalide.")
  if (attraction == "logistique" && (length(parametres) != 3 || !is.numeric(parametres))) cli::cli_abort("Parametres pour logistique invalide.")
  
  froms <- unique(MeapsDataGroup@triplet$fromidINS)
  tos <- unique(MeapsDataGroup@triplet$toidINS) |> sort()
  
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
  
  all_in_optim(jr_dist = les_j,
               p_dist = p_dist,
               xr_dist = MeapsDataGroup@triplet$metric,
               group_from = index_gfrom,
               group_to = index_gto,
               emplois = MeapsDataGroup@emplois,
               actifs = MeapsDataGroup@actifs,
               fuites = MeapsDataGroup@fuites,
               parametres = parametres,
               xr_odds = odds,
               attraction = attraction,
               nthreads = nthreads, verbose = verbose)
}
            
meaps_optim <- function(MeapsDataGroup,  attraction, parametres, odds = 1, 
                        method = "L-BFGS-B", objective = "KL", lower = NULL, upper = NULL,
                        nthreads = 0L, progress = TRUE) { 
          
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.") 
  if (!attraction %in% c("marche", "marche_liss", "double_marche_liss","decay", "logistique")) cli::cli_abort("Pas de fonction choisi ou fonction non paramètrique.")

  arg <- list(MeapsDataGroup, attraction = attraction, odds = odds, nthreads = nthreads, verbose = FALSE)
  
  fn <- switch(objective, 
               "KL" = function(par) {
                 if (progress) cat("*")
                 estim <- do.call(all_in_grouped, args = append(arg, list(parametres = par)))
                 entropie_relative( estim, object@cible$value, floor = 1e-3 / length(estim) )
                 
               }
  )
  
  if (is.null(fn)) stop("Fonction objective non définie !")
  
  nb_par <- length(parametres)
  if (is.null(lower)) lower <- rep(0, nb_par)
  if (is.null(upper)) upper <- rep(Inf, nb_par)
  
  stats::optim(par = parametres, fn = fn, method = method, lower = lower, upper = upper)
}


