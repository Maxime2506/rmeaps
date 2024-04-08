#' Définition d'une classe d'objet rmeaps de préparation des données.
#' 
#' 
#' 
#' 
is_triplet_meaps <- function(obj) {
  if (!inherits(obj, "data.frame")) return(FALSE)
  if (!setequal(names(obj), c("fromidINS", "toidINS", "distances"))) stop("Les variables doivent être fromidINS, toidINS et distances.")
  if (!is.integer(obj$fromidINS) | !is.integer(obj$toidINS) | !is.numeric(obj$distance)) return(FALSE)
  TRUE
}

# Class de base pour les données meaps.
setClass("MeapsData", 
         representation = list(
           distances = "data.frame",
           habitants = "data.frame",
           opportunites = "data.frame",
           fuite = "data.frame"
         ),
         prototype = list(
           distances = NULL,
           habitants = data.frame(),
           opportunites = data.frame(),
           fuite = data.frame()
         )
        )

setMethod(f = "initialize", signature = "MeapsData",
          definition = function(.Object, distances, habitants, opportunites, fuite) {
            
            if (!is_triplet_meaps(distances)) stop("distances n'est pas un triplet valide.")
            N <- unique(distances$fromidINS) |> length()
            K <- unique(distances$toidINS) |> length()
            
            if (nrow(habitants) != N) stop("habitants et distances ne correspondent pas.")
            if (nrow(opportunites) != K) stop("opportunities et distances ne correspondent pas.")
            if (nrow(fuite) != N) stop("habitants et fuite n'ont pas la même longueurs.")
                
            .Object@distances <- distances
            .Object@habitants <- habitants
            .Object@opportunites <- opportunites
            .Object@fuite <- fuite
            
            return(.Object)
          })

setGeneric(name = "distances", def = function(x) standardGeneric("distances"))
setGeneric(name = "distances<-", def = function(x, d) standardGeneric("distances<-"))

setMethod("distances<-", "MeapsData", function(x, d) {
  if (inherits("list")) { d <- as.data.frame(d) }
  if (is_triplet_meaps(d)) {
    x@distances <- d[order(d$fromidINS, d$toidINS),]
  } else stop("Objet invalide.")
})

setGeneric("habitants", function(x) standardGeneric("habitants"))
setGeneric("opportunites", function(x) standardGeneric("opportunites"))
setGeneric("habitants<-", function(x, df) standardGeneric("habitants<-"))
setGeneric("opportunites<-", function(x, df) standardGeneric("opportunites<-"))
setGeneric("fuite", function(x) standardGeneric("fuite"))
setGeneric("fuite<-", function(x, df) standardGeneric("fuite<-"))

setMethod("habitants", "MeapsData", function(x) x@habitants)
setMethod("habitants<-", "MeapsData", function(x, df) {
  if (is.null(x@distances)) stop("distances doit être défini préalablement.")
  if (!inherits(df, "data.frame")) stop("df n'est pas une data.frame.")
  if (!"fromidINS" %in% names(df)) stop("variable fromidINS manquante.")
  if (!"habitants" %in% names(df)) stop("variable habitants manquante")
  if (!is.integer(df$fromidINS)) stop("fromidINS n'est pas en integer.")
  if (!is.numeric(df$habitants)) stop("habitants n'est pas en numeric.")
  
  if (length(setdiff(x@distances$fromidINS, df$fromidINS)) != 0) stop("Certains carreaux de départs d'habitants manquent.")
  if (length(setdiff(df$fromidINS, x@distances$fromidINS)) != 0) stop("Certains carreaux de départs de distances manquent.")
  
  x@habitants <- df[order(df$fromidINS),]
  return(x)    
  }
)
  
setMethod("opportunites", "MeapsData", function(x) x@opportunites)
setMethod("opportunites<-", "MeapsData", function(x, df) {
  if (is.null(x@distances)) stop("distances doit être défini préalablement.")
  if (!inherits(df, "data.frame")) stop("df n'est pas une data.frame.")
  if (!"toidINS" %in% names(df)) stop("variable toidINS manquante.")
  if (!"opportunites" %in% names(df)) stop("variable opportunites manquante")
  if (!is.integer(df$toidINS)) stop("toidINS n'est pas en integer.")
  if (!is.numeric(df$opportunites)) stop("opportunites n'est pas en numeric.")
  
  if (length(setdiff(x@distances$toidINS, df$toidINS)) != 0) stop("Certains carreaux de départs d'habitants manquent.")
  if (length(setdiff(df$toidINS, x@distances$toidINS)) != 0) stop("Certains carreaux de départs de distances manquent.")
  
  x@opportunites <- df[order(df$toidINS),]
  return(x)    
  }
)

setMethod("fuite", "MeapsData", function(x) x@fuite)
setMethod("fuite<-", "MeapsData", function(x, df) {
  if (is.null(x@distances)) stop("distances doit être défini préalablement.")
  if (!inherits(df, "data.frame")) stop("df n'est pas une data.frame.")
  if (!"fromidINS" %in% names(df)) stop("variable fromidINS manquante.")
  if (!"fuite" %in% names(df)) stop("variable fuite manquante")
  if (!is.integer(df$fromidINS)) stop("fromidINS n'est pas en integer.")
  if (!is.numeric(df$fuite)) stop("fuite n'est pas en numeric.")
  
  if (length(setdiff(x@distances$fromidINS, df$fromidINS)) != 0) stop("Certains carreaux de départs de fuite manquent.")
  if (length(setdiff(df$fromidINS, x@distances$fromidINS)) != 0) stop("Certains carreaux de départs de distances manquent.")
  
  x@fuite <- df[order(df$fromidINS),]
  return(x)    
  }
)
  
  



