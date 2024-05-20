

# Class de base pour les données meaps.
setClass("MeapsData", 
         representation = list(
           triplet = "data.frame",
           actifs = "numeric",
           emplois = "numeric",
           fuites = "numeric", 
           froms = "character",
           tos = "character",
           shuf = "matrix",
           j_dist = "integer",
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
           j_dist = integer(),
           p_dist = integer()
         )
)

setMethod(f = "initialize", signature = "MeapsData",
          definition = function(
    .Object,
    triplet, actifs, emplois, fuites,
    nshuf = NULL, seuil = 40, seed = NULL) {
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
            .Object@j_dist <- les_j
            
            p_dist <- triplet |> 
              group_by(fromidINS) |> 
              summarize(n()) |>
              pull() |>
              cumsum()
            .Object@p_dist <- c(0L, p_dist)
            
            if(!is.null(nshuf))
              .Object@shuf <- emiette(
                les_actifs = .Object@actifs, 
                nshuf = nshuf, seuil=seuil, seed = seed)
            
            check_meapsdata(.Object, abort = TRUE)
            
            return(.Object)
          })

# Constructeur
meapsdata <- function(triplet, actifs, emplois, fuites,
                      nshuf = NULL, seuil = 40, seed = NULL) {
  new("MeapsData", 
      triplet, actifs, emplois, fuites, 
      nshuf = nshuf, seuil = seuil, seed = seed)
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

# Class pour les regroupements
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
            
            g_froms <- unique(group_from) |> sort()
            g_tos <- unique(group_to) |> sort()
            
            gf_lab <- seq_along(g_froms) - 1L
            names(gf_lab) <- g_froms
            gt_lab <- seq_along(g_tos) - 1L
            names(gt_lab) <- g_tos
            
            .Object@index_gfrom <- gf_lab[group_from]
            .Object@index_gto <- gt_lab[group_to]
            
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
  
  new("MeapsDataGroup", MeapsData@triplet, MeapsData@actifs,
      MeapsData@emplois, MeapsData@fuites, MeapsData@shuf,
      MeapsData@froms, MeapsData@tos,
      group_from[MeapsData@froms], group_to[MeapsData@tos],
      cible |> arrange(group_from, group_to)
  )
}


