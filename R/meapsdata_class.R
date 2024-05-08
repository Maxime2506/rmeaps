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
    sort_v <- object |> 
      dplyr::group_by(fromidINS) |> 
      dplyr::summarize(test = is.unsorted(metric)) |> 
      dplyr::pull() |> 
      any()
    if (sort_v) {
      status <- 3L
    } else {
      sort_j <- FALSE
      #object |> group_by(fromidINS, metric) |> summarize(test = is.unsorted(j)) |> pull() |> any()
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
  froms <- object@froms
  tos <- object@tos
  N <- froms |> length()
  K <- tos |> length()
  
  if (!setequal(froms, names(object@actifs)) || length(object@actifs) != N) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("triplet et actifs n'ont pas la même liste de fromidINS") }
  }
  if (!setequal(tos, names(object@emplois)) || length(object@emplois) != K) { 
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("triplet et emplois n'ont pas la même liste de toidINS") }
  }
  if (!setequal(froms, names(object@fuites)) || length(object@fuites) != N) {
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

meaps_has_shuf <- function(object, quiet = TRUE) {
  res <- !is.null(object@shuf)
  if(!quiet) {
    if(res) cli::cli_alert_info("un shuf est présent ({nrow(object@shuf)} tirages)")
  }
  res
}

check_meapsdata <- function(object, abort = FALSE, quiet = FALSE) {
  status <- is_triplet_meaps(object@triplet, quiet = quiet) && 
    is_triplet_ordered(object@triplet, quiet = quiet) && 
    is_meapsdata_valid(object, quiet = quiet) &&
    is_meaps_normalized(object) 
  meaps_has_shuf(object, quiet = quiet)
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

check_meapsdatagroup <- function(object, abort = FALSE, quiet = FALSE){
  
  status <- is_meapsdatagroup_valid(object@cible)
  
  if (length(object@group_from) != length(object@actifs)) {
    status <- FALSE
    cli::cli_alert_warning("group_from et actifs ne font pas la même longueur.")
  }
  if (length(object@group_to) != length(object@emplois)) {
    status <- FALSE
    cli::cli_alert_warning("group_to et emplois ne font pas la même longueur.")
  }
  
  if (is.null(names(object@group_from))) {
    status <- FALSE
    cli::cli_alert_warning("group_from n'est pas labelisé")
  }
  if (is.null(names(object@group_to))) {
    status <- FALSE
    cli::cli_alert_warning("group_to n'est pas labelisé.")
  }
  
  if (!setequal(names(object@group_from), names(object@actifs))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("group_from et actifs n'ont pas la même liste de fromidINS.") }
  }
  if (!setequal(names(object@group_to), names(object@emplois))) {
    status <- FALSE
    cli::cli_alert_warning("group_to et emplois n'ont pas la même liste de toidINS.")
  }
  
  if (!setequal(object@cible$group_from, object@group_from)) {
    status <- FALSE
    cli::cli_alert_warning("group_from et cible ne correspondent pas.")
  }
  
  if (!setequal(object@cible$group_to, object@group_to)) {
    status <- FALSE
    cli::cli_alert_warning("group_to et cible ne correspondent pas.")
  }
  
  if (status) { 
    cli::cli_alert_success("Le MeapsDataGroup est valide.") 
  } else if (abort) {
    cli::cli_abort("Le MeapsDataGroup n'est pas valide.")
  }
  
  invisible(status)
}


# Class de base pour les données meaps.
setClass("MeapsData", 
         representation = list(
           triplet = "data.frame",
           actifs = "numeric",
           emplois = "numeric",
           fuites = "numeric", 
           froms = "character",
           tos = "character",
           shuf = "matrix"
         ),
         prototype = list(
           triplet = data.frame(fromidINS = character(), toidINS = character(), metric = numeric()),
           actifs = numeric(),
           emplois = numeric(),
           fuites = numeric(),
           froms = character(),
           tos = character(),
           shuf = matrix()
         )
)

setMethod(f = "initialize", signature = "MeapsData",
          definition = function(.Object, triplet, actifs, emplois, fuites, nshuf = NULL, seuil = 40) {
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
            
            if(!is.null(nshuf))
              .Object@shuf <- emiette(les_actifs = .Object@actifs, nshuf = nshuf, seuil=seuil)
            
            check_meapsdata(.Object, abort = TRUE)
            
            return(.Object)
          })

# Constructeur
meapsdata <- function(triplet, actifs, emplois, fuites, nshuf = NULL, seuil = 40) {
  new("MeapsData", triplet, actifs, emplois, fuites, nshuf = nshuf, seuil = seuil)
}

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
  
  froms <- MeapsData@froms
  tos <- MeapsData@tos
  
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

#' Multishuf, l'original
#' algorithme original meaps, avec des odds discrets et quand même bien présents,
#' et un shuf, avec openMP
#' 
#' @param MeapsData 
#'
#' @param attraction 
#' @param parametres 
#' @param odds 
#' @param nshuf 
#' @param nthreads 
#' @param verbose 
#'
#' @import dplyr
multishuf <- function(MeapsData, attraction = "constant", parametres = 0, odds = 1, nshuf = 16,
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
  
  froms <- MeapsData@froms
  tos <- MeapsData@tos
  
  if(is.null(MeapsData@shuf)) {
    if(nshuf<=0) cli::cli_abort("pas de shufs et pas d'envie de shuf") 
    shuf <- emiette(MeapsData@actifs, nshuf = nshuf)
  } else 
    shuf <- MeapsData@shuf
  
  jlab <- seq_along(tos) - 1L
  names(jlab) <- tos
  
  les_j <- jlab[MeapsData@triplet$toidINS]
  
  p_dist <- MeapsData@triplet |>
    dplyr::group_by(fromidINS) |> 
    dplyr::summarize(n()) |> 
    dplyr::pull() |> 
    cumsum()
  
  p_dist <- c(0L, p_dist)
  res <- meaps_multishuf(
    jr_dist = les_j,
    p_dist = p_dist,
    xr_dist = MeapsData@triplet$metric,
    emplois = MeapsData@emplois,
    actifs = MeapsData@actifs,
    fuites = MeapsData@fuites,
    shuf = shuf,
    attraction = attraction,
    parametres = parametres,
    xr_odds = odds,
    nthreads = nthreads, verbose = verbose) 
  colnames(res) <- MeapsData@tos
  
  res |>
    tibble::as_tibble() |> 
    dplyr::mutate(fromidINS = MeapsData@froms) |> 
    tidyr::pivot_longer(cols = -fromidINS, names_to = "toidINS", values_to = "flux") |> 
    dplyr::filter(flux>0) |> 
    arrange(desc(flux))
}

#' Multishuf, la version continue
#' algorithme original meaps, 2eme version, avec des odds continus,
#' et un shuf, avec openMP
#' 
#' @param MeapsData 
#'
#' @param attraction 
#' @param parametres 
#' @param odds 
#' @param nshuf 
#' @param nthreads 
#' @param verbose 
#'
#' @import dplyr
multishuf_oc <- function(MeapsData, attraction = "constant", parametres = 0, odds = 1, nshuf = 16,
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
  
  froms <- MeapsData@froms
  tos <- MeapsData@tos
  
  if(is.null(MeapsData@shuf)) {
    if(nshuf<=0) cli::cli_abort("pas de shufs et pas d'envie de shuf") 
    shuf <- emiette(MeapsData@actifs, nshuf = nshuf)
  } else 
    shuf <- MeapsData@shuf
  
  jlab <- seq_along(tos) - 1L
  names(jlab) <- tos
  
  les_j <- jlab[MeapsData@triplet$toidINS]
  
  p_dist <- MeapsData@triplet |>
    dplyr::group_by(fromidINS) |> 
    dplyr::summarize(n()) |> 
    dplyr::pull() |> 
    cumsum()
  
  p_dist <- c(0L, p_dist)
  res <- multishuf_optim_cpp(
    jr_dist = les_j,
    p_dist = p_dist,
    xr_dist = MeapsData@triplet$metric,
    emplois = MeapsData@emplois,
    actifs = MeapsData@actifs,
    fuites = MeapsData@fuites,
    group_from = 0:(length(MeapsData@froms)-1),
    group_to = 0:(length(MeapsData@tos)-1),
    shuf = shuf,
    attraction = attraction,
    parametres = parametres,
    xr_odds = odds,
    nthreads = nthreads, verbose = verbose) 
  colnames(res) <- MeapsData@tos
  
  res |>
    tibble::as_tibble() |> 
    dplyr::mutate(fromidINS = MeapsData@froms) |> 
    tidyr::pivot_longer(cols = -fromidINS, names_to = "toidINS", values_to = "flux") |> 
    dplyr::filter(flux>0) |> 
    arrange(desc(flux))
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
          definition = function(.Object, triplet, actifs, emplois, fuites, shuf,
                                froms, tos, group_from, group_to, cible) {
            
            .Object@triplet <- triplet
            .Object@actifs <- actifs[froms]
            .Object@emplois <- emplois[tos]
            .Object@fuites <- fuites[froms]
            .Object@froms <- froms
            .Object@tos <- tos
            .Object@group_from <- group_from[froms]
            .Object@group_to <- group_to[tos]
            .Object@cible <- cible |> arrange(group_from, group_to)
            .Object@shuf <- shuf
            check_meapsdata(.Object, abort = TRUE,)
            check_meapsdatagroup(.Object, abort = TRUE)
            
            return(.Object)
          })

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

# Constructeur
meapsdatagroup <- function(MeapsData, group_from, group_to, cible) {
  
  new("MeapsDataGroup", MeapsData@triplet, MeapsData@actifs, MeapsData@emplois, MeapsData@fuites, MeapsData@shuf,
      MeapsData@froms, MeapsData@tos,
      group_from[MeapsData@froms], group_to[MeapsData@tos],
      cible |> arrange(group_from, group_to)
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
  
  froms <- MeapsDataGroup@froms
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
  
  res <- all_in_optim(jr_dist = les_j,
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
  MeapsDataGroup@cible |> 
    mutate(cible = value, value = res[res>0]) |> 
    arrange(desc(cible))
}



#' Multishuf continu groupé
#' 
#' Calcule les flux meaps à partir de l'algorithme multishuf et
#' retourne un résultat agrégé selon les groupes envoyés
#'
#' @param MeapsDataGroup object avec l'information de groupe
#' @param attraction fonction pour calculer les odds
#' @param parametres paramètre de la fonction "attraction"
#' @param odds possibilité de passer un vecteur d'odds complet
#' @param nthreads nombre de threads
#' @param verbose affiche des messages
#'
#' @return
#' @export
#'
multishuf_grouped <- function(MeapsDataGroup,  attraction = "constant", parametres = 0, odds = 1, 
                              nthreads = 0L, verbose = TRUE) {
  
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) 
    cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.") 
  # Validation des paramètres
  if (!attraction %in% c("constant", "marche", "marche_liss", "double_marche_liss","decay", "logistique")) 
    cli::cli_abort("Fonction attraction inconnue")
  
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
  
  res <- multishuf_optim_cpp(
    jr_dist = les_j,
    p_dist = p_dist,
    xr_dist = MeapsDataGroup@triplet$metric,
    group_from = index_gfrom,
    group_to = index_gto,
    emplois = MeapsDataGroup@emplois,
    actifs = MeapsDataGroup@actifs,
    fuites = MeapsDataGroup@fuites,
    shuf = MeapsDataGroup@shuf,
    parametres = parametres,
    xr_odds = odds,
    attraction = attraction,
    nthreads = nthreads, verbose = verbose)
  colnames(res) <- g_tos
  res |>
    tibble::as_tibble() |> 
    dplyr::mutate(group_from = g_froms) |> 
    tidyr::pivot_longer(cols = -group_from, names_to = "group_to", values_to = "value") |>  
    filter(value>0) |> 
    left_join(MeapsDataGroup@cible |> rename(target = value), by = c("group_from", "group_to")) |> 
    arrange(desc(value))
}


meaps_optim <- function(MeapsDataGroup,  attraction, parametres, odds = 1, 
                        method = "L-BFGS-B", objective = "KL", lower = NULL, upper = NULL,
                        nthreads = 0L, progress = TRUE) { 
  
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.") 
  if (!attraction %in% c("marche", "marche_liss", "double_marche_liss","decay", "logistique")) cli::cli_abort("Pas de fonction choisi ou fonction non paramètrique.")
  
  arg <- list(MeapsDataGroup, attraction = attraction, odds = odds, nthreads = nthreads, verbose = FALSE)
  
  env <- environment()
  fn <- switch(
    objective, 
    "KL" = function(par) {
      if (progress) cli::cli_progress_update(.envir = env)
      estim <- do.call(all_in_grouped, args = append(arg, list(parametres = par)))
      entropie_relative( estim, MeapsDataGroup@cible$value, floor = 1e-3 / length(estim) )
      
    }
  )
  
  if (is.null(fn)) stop("Fonction objective non définie !")
  
  nb_par <- length(parametres)
  if (is.null(lower)) lower <- rep(0, nb_par)
  if (is.null(upper)) upper <- rep(Inf, nb_par)
  
  cli::cli_progress_bar(.envir = env, clear = FALSE)
  stats::optim(par = parametres, fn = fn, method = method, lower = lower, upper = upper)
  cli::cli_progress_done(.envir = env)
}


