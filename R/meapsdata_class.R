#' Définition d'une classe d'objet rmeaps de préparation des données.
#' 
#' 
#' 
#' 
is_triplet_meaps <- function(object, quiet = FALSE) {
  if (!inherits(object, "data.frame")) return(FALSE)
  if (!setequal(names(object), c("i", "j", "value"))) {
    if (quiet) { return(FALSE) } else stop("Les variables doivent être i, j et value.")
  }
  if (!is.integer(object$i) | !is.integer(object$j) | !is.numeric(object$value)) {
    if (quiet) { return(FALSE) } else stop("Le type des variables n'est pas valide.")
  }
  if (any(is.na(object))) {
    if (quiet) { return(FALSE) } else stop("Il y a des valeurs manquantes.")
  }
  TRUE
}

is_meapsdata_ordered <- function(object) {
  if ( is.unsorted(object$i) ) return(FALSE)
  sort_v <- by(object, object$i, \(x) is.unsorted(x$value)) |> any()
  if (sort_v) return(FALSE)
  sort_j <- by(object, list(object$i, object$value), \(x) {
    is.unsorted(x$j, strictly = TRUE)
    }) |> any(na.rm = TRUE)
  return(!sort_j)
}

is_idINS_valid <- function(id) {
  if (typeof(id) != "character") return(FALSE)
  !is.unsorted(id)
}

is_meapsdata_valid <- function(object) {
  if (min(object@distances$i) != 0) stop("i doit commencer à 0 (convention C++)")
  if (min(object@distances$j) != 0) stop("j doit commencer à 0 (convention C++)")
  
  N <- unique(object@distances$i) |> length()
  K <- unique(object@distances$j) |> length()
  
  if (length(object@fromidINS) != N) stop("fromidINS et distances ne correspondent pas.")
  if (length(object@toidINS) != K) stop("toidINS et distances ne correspondent pas.")
  if (nrow(object@actifs) != N) stop("habitants et distances ne correspondent pas.")
  if (nrow(object@emplois) != K) stop("opportunities et distances ne correspondent pas.")
  if (nrow(object@fuite) != N) stop("habitants et fuite n'ont pas la même longueurs.")
  
  delta <- (sum(object@actifs * (1 - object@fuite)) - sum(object@emplois))/sum(object@emplois)
  if (delta > 1e-5)  warning("Les actifs restant dans la zone et les emplois ne correspondent pas.")
  
  if (is.null(names(object@actifs))) stop("actifs n'est pas labelisée.")
  if (is.null(names(object@emplois))) stop("emplois n'est pas labelisée.")
  if (is.null(names(object@fuite))) stop("fuite n'est pas labelisée.")
  
  if (!setequal(names(object@actifs), object@fromidINS)) stop("Les labels d'actifs ne correspondent pas à fromidINS.")
  if (!setequal(names(object@emplois), object@toidINS)) stop("Les labels d'emplois ne correspondent pas à toidINS.")
  if (!setequal(names(object@fuite), object@fromidINS)) stop("Les labels de fuite ne correspondent pas à fromidINS.")
  TRUE
}

check_meapsdata <- function(object) {
  is_triplet_meaps(object) && 
    is_meapsdata_ordered(object) && 
    is_idINS_valid(object@fromidINS) &&
    is_idINS_valid(object@toidINS) &&
    is_meapsdata_valid(object)
}

formate_ciblemeaps <- function(x) {
  if (!inherits(x, "data.frame") || length(x) != 3) return(FALSE)
  filter_from <- grepl("commune", names(x))
  if (sum(filter_from) != 1) stop("Impossible d'identifier le groupe commune")
  filter_to <- grepl("dclt", names(x))
  if (sum(filter_to) != 1) stop("Impossible d'identifier le groupe dclt")
  names(x)[filter_from] <- "commune"
  names(x)[filter_to] <- "dclt"
  names(x)[!filter_from & !filter_to] <- "value"
  
  x <- x[order(x$commune, x$dclt), ]
  
  return(x)
}

check_meapsdatagroup <- function(object){
  
  if (length(object@group_from) != length(object@actifs)) stop("group_from et actifs ne font pas la même longueur.")
  if (length(object@group_to) != length(object@emplois)) stop("group_to et emplois ne font pas la même longueur.")
  
  if (is.null(names(object@group_from))) stop("group_from n'est pas labelisé.")
  if (is.null(names(object@group_to))) stop("group_to n'est pas labelisé.")
  
  if(!setequal(names(object@group_from), object@fromidINS)) stop("Les labels de group_from ne correspondent pas à fromidINS")
  if(!setequal(names(object@group_to), object@toidINS)) stop("Les labels de group_to ne correspondent pas à toidINS")
  
  # La cible (MOBPRO) est une dataframe avec les colonnes contenant commune et dclt, d'un côté, et autre chose de l'autre.
  cible <- formate_ciblemeaps(object@cible)
  if (!setequal(cible$commune, object@group_from)) stop("group_from et cible ne correspondent pas.")
  if (!setequal(cible$dclt, object@group_to)) stop("group_to et cible ne correspondent pas.")
  
  check_meapsdata(object)
}


# Class de base pour les données meaps.
setClass("MeapsData", 
         representation = list(
           distances = "data.frame",
           fromidINS = "character",
           toidINS = "character",
           actifs = "numeric",
           emplois = "numeric",
           fuite = "numeric"
         ),
         prototype = list(
           distances = NULL,
           fromidINS = character(),
           toidINS = character(),
           actifs = numeric(),
           emplois = numeric(),
           fuite = numeric()
         ),
         validity = check_meapsdata
)

setMethod(f = "initialize", signature = "MeapsData",
          definition = function(.Object, distances, fromidINS, toidINS, actifs, emplois, fuite) {
            
            .Object@distances <- distances
            .Object@fromidINS <- fromidINS
            .Object@toidINS <- toidINS
            .Object@actifs <- actifs[fromidINS]
            .Object@emplois <- emplois[toidINS]
            .Object@fuite <- fuite[fromidINS]
            
            return(.Object)
          })

setMethod("show", "MeapsData", function(object) {
  N <- unique(object@distances$i) |> length()
  K <- unique(object@distances$j) |> length()
  tx <- round(nrow(object@distances) / N / K * 100, digits = 1)
  cat("MeapsData :\n")
  cat("Matrice des distances", N, "x", K, "remplie à", tx, "%\n")
  cat("Nombre d'actifs =", sum(object@actifs),"\n")
  cat("Nombre d'emplois =", sum(object@emplois), "\n")
  cat("Nombre de fuyards =", sum(object@fuite * object@actifs))
})

setGeneric("all_in", function(object, ...) {standardGeneric("all_in")})

setMethod("all_in", 
          signature = "MeapsData",
          function(object, ...){
 
  arg <- list(...)
  # Définition de valeurs par défaut.
  if (is.null(arg$attraction)) { attraction <- "constant" } else { attraction <- arg$attraction }
  if (is.null(arg$parametres)) { parametres <- 0.0 } else { parametres <- arg$parametres }
  if (is.null(arg$odds)) { odds <- 1.0 } else { odds <- arg$odds }
  if (is.null(arg$nthreads)) { nthreads <- 0L} else { nthreads <- arg$nthreads }
  if (is.null(arg$verbose)) { verbose <- TRUE } else { verbose <- arg$verbose }
  if (is.null(arg$normalisation)) { normalisation <- FALSE } else { normalisation <- arg$normalisation }
  if (is.null(arg$fuite_min)) { fuite_min <- 1e-3 } else { fuite_min <- arg$fuite_min }
  
  # Validation des paramètres
  if (!attraction %in% c("constant", "marche", "marche_liss", "double_marche_liss","decay", "logistique")) stop("Fonction attraction inconnue")
  
  # RQ : pas de méthode pour vérifier le bon ordre des odds.
  if (attraction == "odds" && length(odds) != nrow(object@distances) ) stop("vecteur odds invalide.")
  
  if (attraction == "marche" && (length(parametres) != 2 || !is.numeric(parametres))) stop("Parametres pour marche invalide.")
  if (attraction == "marche_liss" && (length(parametres) != 2 || !is.numeric(parametres))) stop("Parametres pour marche_liss invalide.")
  if (attraction == "double_marche_liss" && (length(parametres) != 4 || !is.numeric(parametres))) stop("Parametres pour double_marche_liss invalide.")
  if (attraction == "decay" && (length(parametres) != 2 || !is.numeric(parametres))) stop("Parametres pour decay invalide.")
  if (attraction == "logistique" && (length(parametres) != 3 || !is.numeric(parametres))) stop("Parametres pour logistique invalide.")
  # p_dist <- aggregate(object@distances$i, by = list(object@distances$i), FUN = length)$x |> cumsum()
  p_dist <- object@distances |> group_by(i) |> summarize(n()) |> pull()
  p_dist <- c(0L, p_dist)
  

  meaps_all_in(jr_dist = object@distances$j,
               p_dist = p_dist,
               xr_dist = object@distances$value,
               emplois = object@emplois,
               actifs = object@actifs,
               fuite = object@fuite,
               parametres = parametres,
               xr_odds = odds,
               attraction = attraction,
               nthreads = nthreads, verbose = verbose, normalisation = normalisation, fuite_min = fuite_min) |>
    dplyr::left_join(data.frame(i = seq_along(object@fromidINS) - 1L, fromidINS = object@fromidINS), by = "i") |>
    dplyr::left_join(data.frame(j = seq_along(object@toidINS) - 1L, toidINS = object@toidINS), by = "j") |>
    dplyr::select(fromidINS, toidINS, flux)
})


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
         contains = "MeapsData",
         validity = check_meapsdatagroup)


setMethod(f = "initialize", signature = "MeapsDataGroup",
          definition = function(.Object, distances, fromidINS, toidINS, group_from, group_to, actifs, emplois, fuite, cible) {
            
            .Object@distances <- distances
            .Object@fromidINS <- fromidINS
            .Object@toidINS <- toidINS
            .Object@group_from <- sort(group_from)
            .Object@group_to <- sort(group_to)
            .Object@actifs <- actifs[fromidINS]
            .Object@emplois <- emplois[toidINS]
            .Object@fuite <- fuite[fromidINS]
            .Object@cible <- formate_ciblemeaps(cible)
            
            return(.Object)
          })

setMethod("show", "MeapsDataGroup", function(object) {
  N <- unique(object@distances$i) |> length()
  K <- unique(object@distances$j) |> length()
  Ng <- unique(object@group_from) |> length()
  Kg <- unique(object@group_to) |> length()
  tx <- round(nrow(object@distances) / N / K * 100, digits = 1)
  cat("MeapsData :\n")
  cat("Matrice des distances", N, "x", K, "remplie à", tx, "%\n")
  cat("Nombre d'actifs =", sum(object@actifs),"\n")
  cat("Nombre d'emplois =", sum(object@emplois), "\n")
  cat("Nombre de fuyards =", sum(object@fuite * object@actifs), "\n")
  cat("Nombre de groupes de départ =", Ng,"\n")
  cat("Nombre de groupes d'arrivée =", Kg, "\n")
})

# Constructeur
meapsdatagroup <- function(MeapsData, group_from, group_to) {

  new("MeapsDataGroup", MeapsData@distances, MeapsData@fromidINS, MeapsData@toidINS, 
      group_from[MeapsData@fromidINS], group_to[MeapsData@toidINS],
      MeapsData@actifs, MeapsData@emplois, MeapsData@fuite)
}

setGeneric("all_in_grouped", function(object, ...) {standardGeneric("all_in_grouped")})

setMethod("all_in_grouped", 
          signature = "MeapsDataGroup",
          function(object, ...){
            
            arg <- list(...)
            # Définition de valeurs par défaut.
            if (is.null(arg$attraction)) { attraction <- "constant" } else { attraction <- arg$attraction }
            if (is.null(arg$parametres)) { parametres <- 0.0 } else { parametres <- arg$parametres }
            if (is.null(arg$odds)) { odds <- 1.0 } else { odds <- arg$odds }
            if (is.null(arg$nthreads)) { nthreads <- 0L} else { nthreads <- arg$nthreads }
            if (is.null(arg$verbose)) { verbose <- TRUE } else { verbose <- arg$verbose }
            if (is.null(arg$normalisation)) { normalisation <- FALSE } else { normalisation <- arg$normalisation }
            if (is.null(arg$fuite_min)) { fuite_min <- 1e-3 } else { fuite_min <- arg$fuite_min }
            
            # Validation des paramètres
            if (!attraction %in% c("constant", "marche", "marche_liss", "double_marche_liss","decay", "logistique")) stop("Fonction attraction inconnue")
            
            # RQ : pas de méthode pour vérifier le bon ordre des odds.
            if (attraction == "odds" && length(odds) != nrow(object@distances) ) stop("vecteur odds invalide.")
            
            if (attraction == "marche" && (length(parametres) != 2 || !is.numeric(parametres))) stop("Parametres pour marche invalide.")
            if (attraction == "marche_liss" && (length(parametres) != 2 || !is.numeric(parametres))) stop("Parametres pour marche_liss invalide.")
            if (attraction == "double_marche_liss" && (length(parametres) != 4 || !is.numeric(parametres))) stop("Parametres pour double_marche_liss invalide.")
            if (attraction == "decay" && (length(parametres) != 2 || !is.numeric(parametres))) stop("Parametres pour decay invalide.")
            if (attraction == "logistique" && (length(parametres) != 3 || !is.numeric(parametres))) stop("Parametres pour logistique invalide.")
            
            p_dist <- aggregate(object@distances$i, by = list(object@distances$i), FUN = length)$x |> cumsum()
            p_dist <- c(0L, p_dist)
            
            lab_groupfrom <- unique(object@group_from)
            index_from <- sapply(object@group_from, FUN = \(x) which(x == lab_groupfrom) - 1L)
            lab_groupto <- unique(object@group_to)
            index_to <- sapply(object@group_to, FUN = \(x) which(x == lab_groupto) - 1L)
            
            all_in_optim(jr_dist = object@distances$j,
                         p_dist = p_dist,
                         xr_dist = object@distances$value,
                         group_from = index_from,
                         group_to = index_to,
                         emplois = object@emplois,
                         actifs = object@actifs,
                         fuite = object@fuite,
                         parametres = parametres,
                         xr_odds = odds,
                         attraction = attraction,
                         nthreads = nthreads, verbose = verbose, normalisation = normalisation, fuite_min = fuite_min) 
            
          })

setGeneric("meaps_optim", function(object, ...) {standardGeneric("meaps_optim")})

setMethod("meaps_optim", 
          signature = "MeapsDataGroup",
          function(object, ...) {
            arg <- list(...)
            
            if (is.null(arg$attraction) || arg$attraction == "constant") stop("Pas de fonction choisi ou fonction non paramètrique.")
            if (is.null(arg$parametres)) stop("Pas de paramètres définis")
            if (is.null(arg$method)) { arg$method <- "Nelder-Mead" }
            
            par0 <- arg$parametres
            arg$parametres <- NULL
            
            # Définition des valeurs par défaut
            if (is.null(arg$objective)) { objective <- "KL" } else { objective <- arg$objective}
            
            fn <- switch(objective, 
                         "KL" = function(par) {
                           les_arguments <- append(object, list(parametres = par)) |> append(arg)
                           estim <- do.call(all_in_grouped, args = les_arguments)
                           entropie_relative(estim, object@cible$value, floor = 1e-3 / length(estim) )
                           }
                         )
            
            if (is.null(fn)) stop("Fonction objective non définie !")
            
            stats::optim(par = par0, fn = fn, method = arg$method)

            
          })


