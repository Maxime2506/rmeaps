#' Calcule la kl-divergence entre des flux "aux carreaux" et une matrice de flux "communaux" de référence.
#' @param flux Matrice de flux, typiquement issue de meaps.
#' @param flux_ref Matrice de flux regroupés par catégories (typiquement des communes), servant de référence.
#' @param group_orig Vecteur désignant la catégorie des "carreaux" de départ. Si null, les rownames de flux sont utilisées.
#' @param group_dest Vecteur désignant la catégorie des "carreaux" d'arrivée. Si null, les colnames de flux sont utilisées.
#' 
#' @return Retourne la valeur de Kullback-Leibler.
kl_fluxes <- function(flux, flux_ref, group_orig = NULL, group_dest = NULL, seuil_collapse = .1) {
  
  if (is.null(group_orig)) {
    group_orig <- rownames(flux)
  }
  if (is.null(group_dest)) {
    group_dest <- colnames(flux)
  }
  
  if (!is.integer(group_orig) | !is.integer(group_dest)) stop("Les catégories de regroupements ne sont pas en integer.")
  
  if (length(group_orig) != nrow(flux)) stop("Le nombre de lignes de flux ne correspond pas à la longueur de group_orig.")
  if (length(group_dest) != ncol(flux)) stop("Le nombre de colonnes de flux ne correspond pas à la longueur de group_dest.")
  
  orig_reordered <- unique(group_orig) |> sort()
  dest_reordered <- unique(group_dest) |> sort()
  
  if (length(orig_reordered) != nrow(flux_ref)) stop("Le nombre de lignes de flux_ref ne correspond pas au nombre de catégories de group_orig.")
  if (length(dest_reordered) != ncol(flux_ref)) stop("Le nombre de colonnes de flux_ref ne correspond pas au nombre de catégories de group_dest.")
  
  flux_ref <- flux_ref[orig_reordered, dest_reordered] # On réordonne, au cas où, la matrice de référence pour coller à la sortie de communaliser.
  
  res <- communaliser(flux, group_orig = group_orig, group_dest = group_dest)
  
  kullback_leibler(res, flux_ref, seuil_collapse = seuil_collapse)
  
  }