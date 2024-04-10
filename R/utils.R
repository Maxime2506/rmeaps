#' utils
#' 
#' Fonction de test de la validité d'un triplet sous forme de liste ou de data.frame.
#' Nom de variables admises : i ou fromidINS, j ou toidINS.
#' Le nom de la troisième variable est libre.
#' 
is_triplet <- function(obj) {
  if ( !inherits(obj, "data.frame") | !inherits(obj, "list") | length(obj) != 3) return(FALSE)
  colnames(obj)[colnames(obj) == "fromidINS"] <- "i"
  colnames(obj)[colnames(obj) == "toidINS"] <- "j"
  colnames(obj)[!colnames(obj) %in% c("i", "j")] <- "x"
  
  if ( !setequal(names(obj), c("i", "j", "x")) ) return(FALSE)
  if ( !(is.integer(obj$i) & is.integer(obj$j) & is.numeric(obj$x)) ) return(FALSE)
  if ( length(obj$i) != length(obj$j) | length(obj$i) != length(obj$x)) return(FALSE)
  return(TRUE)
}

.make_triplet <- function(obj) {
  colnames(obj)[colnames(obj) == "fromidINS"] <- "i"
  colnames(obj)[colnames(obj) == "toidINS"] <- "j"
  colnames(obj)[!colnames(obj) %in% c("i", "j")] <- "x"
  
  if (inherits(obj, "list")) {
    as.data.frame(obj)
  } else {
    obj
  }
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

#' Fonction de retraitement de la distance si elle est sous forme de triplet dans une liste.
#' C'est en général la résultante de la fonction Matrix::mat2triplet.
#' On suppose que max(j) donne le nombre de colonnes, ce qui est attendu dans les analyses meaps.
#' On ne teste pas la validité du triplet.
#' Ici encore les valeurs x nulles deviennent des petites distances.
#' Liste de trois vecteurs : i, j et x.
#' @import Matrix
.triplet2dgR <- function(obj) {
  
  if ( setequal(names(obj), c("fromidINS", "toidINS", "distances")) ) {
    obj <- obj[, c("fromidINS", "toidINS", "distances")]
    names(obj) <- c("i", "j", "x")
  }
  
  obj$x[obj$x == 0] <- .petite_distance(obj$x)
  
  spMatrix(nrow = length(dist$i), ncol = max(dist$j), i = dist$i, j = dist$j, x = dist$x)
  
}

