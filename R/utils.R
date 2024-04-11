#' utils
#' 

#' Fonction de test de la validité d'un triplet sous de data.frame.
#' Nom de variables admises : fromidINS, toidINS.
#' Le nom de la troisième variable est libre.
#' 
is_triplet <- function(obj) {
  if ( !inherits(obj, "data.frame") & length(obj) != 3) return(FALSE)
  if ( length(setdiff(c("fromidINS", "toidINS"), names(obj))) != 0 ) return(FALSE)
  if ( !is.numeric(obj[[setdiff(names(obj), c("fromidINS", "toidINS"))]]) ) return(FALSE)
  return(TRUE)
}

#' Fonction de transformation d'un triplet fromidINS toidINS value à une liste avec i j x et les clés de correspondances.
triplet2listij <- function(data) {
  if (!is_triplet(data)) stop("Ce n'est pas un triplet valide.")
  colnames(data)[!colnames(data) %in% c("fromidINS", "toidINS")] <- "x"
  
  cle_from <- unique(data$fromidINS) |> sort()
  cle_to <- unique(data$toidINS) |> sort()
  
  data <- data[order(data$fromidINS, data$toidINS), ]
  
  data <- dplyr::inner_join(data, data.frame(fromidINS = cle_from, i = seq_along(cle_from)), by = "fromidINS") |> 
    dplyr::inner_join(data.frame(toidINS = cle_to, j = seq_along(cle_to)), by = "toidINS") |> 
    dplyr::select(i, j, x)
  
  data$x[data$x == 0] <- .petite_distance(data$x)
  
  list(dgr = sparseMatrix(i = data$i, j = data$j, x = data$x, dims = c(max(data$i), max(data$j)), repr = "R"), 
       cle_from = cle_from, 
       cle_to = cle_to)
}

#' Fonction de transformation d'un triplet de log(odds) avec fromidINS et toidINS en une matrice row sparse.
#' Les cle_from et cle_to sont celles tirées du triplet des distances.
#' @import Matrix
tripletlodds2dgr <- function(data, cle_from, cle_to) {
  if (!is_triplet(data)) stop("Ce n'est pas un triplet valide.")
  colnames(data)[!colnames(data) %in% c("fromidINS", "toidINS")] <- "x"
  
  data <- data[data$x != 0, ]
  data <- dplyr::inner_join(data, data.frame(fromidINS = cle_from, i = seq_along(cle_from)), by = "fromidINS") |> 
    dplyr::inner_join(data.frame(toidINS = cle_to, j = seq_along(cle_to)), by = "toidINS")
  
 sparseMatrix(i = data$i, j = data$j, x = data$x, dims = c(length(cle_from), length(cle_to)), repr = "R")
}


#' Fonction de choix d'une petite distance pour remplacer les zéros éventuels.
#' Normalement il ne devrait pas y avoir de zéros. 
#' Au sein d'un même carreau de dim a, la distance entre deux points est d'environ a/2.
#' La distance retenue importe peu tant que les rangs sont respectés.
.petite_distance <- function(d, facteur = 10) {
  min(as.numeric(d[d != 0]), na.rm = TRUE) / facteur
}

#' Fonction de retraitement de la distance si elle est sous forme d'une matrice classique.
#' Les NA deviennent sparses. Les zéros deviennent simplement des petites valeurs.
#' Le résultat est une dgRMatrix.
#' 
#' @import Matrix
.matrix2dgR <- function(x, na2zero = TRUE) {
  if (na2zero) {
    x[x == 0] <- .petite_distance(x)
    x[is.na(x)] <- 0
  }
  as(as(as(x, "dMatrix"), "generalMatrix"), "RsparseMatrix")
}


