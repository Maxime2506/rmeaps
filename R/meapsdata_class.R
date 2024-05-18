#' Définition d'une classe d'objet rmeaps de préparation des données.
#' 
#' 
#' 
#' 

# Class de base pour les données meaps.
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

setMethod(f = "initialize", signature = "MeapsData",
          definition = function(.Object, triplet, actifs, emplois, fuites) {
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
            
            check_meapsdata(.Object, abort = TRUE)
            
            return(.Object)
          })

# Constructeur
meapsdata <- function(triplet, actifs, emplois, fuites) {
  new("MeapsData", triplet, actifs, emplois, fuites)
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
          definition = function(.Object, triplet, actifs, emplois, fuites,
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
            
            check_meapsdata(.Object, abort = TRUE)
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
  
  new("MeapsDataGroup", MeapsData@triplet, MeapsData@actifs, MeapsData@emplois, MeapsData@fuites,
      MeapsData@froms, MeapsData@tos,
      group_from[froms], group_to[tos],
      cible |> arrange(group_from, group_to)
  )
}

#' @import dplyr
all_in_grouped <- function(MeapsDataGroup,  attraction = "constant", parametres = 0, nthreads = 0L, verbose = TRUE) {
  
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.") 
  # Validation des paramètres
  if (!attraction %in% c("constant", "marche", "marche_liss", "decay", "logistique")) cli::cli_abort("Fonction attraction inconnue")
  
  if (attraction == "marche" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour marche invalide.")
  if (attraction == "double_marche_liss" && (length(parametres) != 4 || !is.numeric(parametres))) cli::cli_abort("Parametres pour double_marche_liss invalide.")
  if (attraction == "decay" && (length(parametres) != 2 || !is.numeric(parametres))) cli::cli_abort("Parametres pour decay invalide.")
  if (attraction == "logistique" && (length(parametres) != 3 || !is.numeric(parametres))) cli::cli_abort("Parametres pour logistique invalide.")
  
  #froms <- unique(MeapsDataGroup@triplet$fromidINS)
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
  
  meapsclass(jr_dist = les_j,
             p_dist = p_dist,
             xr_dist = MeapsDataGroup@triplet$metric,
             emplois = MeapsDataGroup@emplois,
             actifs = MeapsDataGroup@actifs,
             fuites = MeapsDataGroup@fuites,
             parametres = parametres,
             attraction = attraction,
             group_from = index_gfrom,
             group_to = index_gto,
             cible = MeapsDataGroup@cible$value,
             nthreads = nthreads, verbose = verbose)
}

meaps_optim <- function(MeapsDataGroup,  attraction, parametres, 
                        method = "L-BFGS-B", objective = "KL", lower = NULL, upper = NULL,
                        nthreads = 0L, progress = TRUE) { 
  
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.") 
  if (!attraction %in% c("marche", "marche_liss", "decay", "logistique")) cli::cli_abort("Pas de fonction choisi ou fonction non paramètrique.")
  
  arg <- list(MeapsDataGroup, attraction = attraction, nthreads = nthreads, verbose = FALSE)
  
  env <- environment()
  fn <- switch(
    objective, 
    "KL" = function(par) {
      if (progress) cli::cli_progress_update(.envir = env)
      do.call(all_in_grouped, args = append(arg, list(parametres = par)))$kl
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


