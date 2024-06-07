#' Définition d'une classe d'objet rmeaps de préparation des données.
#' 
#' 
#' 
#' 

# Class de base pour les données meaps.

is.meapsdata.ordered <- function(data) {
  if (is.unsorted(data$i)) return(FALSE)
  sortx <- by(data, INDICES = data$i, \(x) is.unsorted(x$x)) |> any(na.rm = TRUE)
  if (sortx) return(FALSE)
  sortj <- by(data, INDICES = list(data$i, data$x), \(x) is.unsorted(x$j)) |> any(na.rm = TRUE)
  if (sortj) return(FALSE)
  return(TRUE)
}

is.meapsdata.triplet <- function(data) {
  if ( !inherits(data, "data.frame") ) return(FALSE)
  if ( !setequal(c("i", "j", "x"), names(data)) ) return(FALSE)
  if ( !is.integer(data$i) | !is.integer(data$j) | !is.numeric(data$x)) return(FALSE)
  if ( min(data$i) != 0 | min(data$j) != 0 | min(data$x) <= 0) return(FALSE)
  return(TRUE)
}

setClass("MeapsData", 
         representation = list(
           distances = "data.frame",
           actifs = "numeric",
           emplois = "numeric",
           fuites = "numeric"
         ),
         prototype = list(
           distances = NULL,
           actifs = numeric(),
           emplois = numeric(),
           fuites = numeric()
         )
        )

setMethod(f = "initialize", signature = "MeapsData",
          definition = function(.Object, distances, actifs, emplois, fuites) {
            
            if (!is.meapsdata.triplet(distances)) stop("distances n'est pas un triplet valide.")
            if (!is.meapsdata.ordered(distances)) stop("distances n'est pas triés comme attendu.")
            N <- unique(distances$i) |> length()
            K <- unique(distances$j) |> length()
            
            if (length(actifs) != N) stop("actifs et distances ne correspondent pas.")
            if (length(emplois) != K) stop("emplois et distances ne correspondent pas.")
            if (length(fuites) != N) stop("fuite  et distances ne correspondent pas.")
            
            if ( abs(sum(actifs * (1 - fuites)) / sum(emplois) - 1) > 1e-5 ) warning("actifs non fuyants ne correspondent pas aux emplois.")
                
            .Object@distances <- distances
            .Object@actifs <- actifs
            .Object@emplois <- emplois
            .Object@fuites <- fuites
            
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
  
  



