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
                      nthreads = 0L, verbose = TRUE, gbperthreads=4) {
  
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
  
  # check mem
  size <- length(MeapsData@triplet$metric)/1024^3
  large <- 4*size>gbperthreads/4
  ntr <- nthreads
  if(large) {
    gc()
    memused <- as.numeric(lobstr::mem_used())/1024^3
    if(nthreads==0) ntr <- max_threads()
    memleft <- gbperthreads*max_threads() - 4*size*ntr -20*size - memused
    if(memleft < ntr) {
      ntr <- min(max_threads(), max(1, round((gbperthreads*max_threads() - 20*size - memused)/(1.5+4*size))))
      cli::cli_warn("le nombre de threads est réduit à {ntr}")
    }
  }
  
  res <- meaps_multishuf_cpp(
    jr_dist = MeapsData@j_dist,
    p_dist = MeapsData@p_dist,
    xr_dist = MeapsData@triplet$metric,
    emplois = MeapsData@emplois,
    actifs = MeapsData@actifs,
    fuites = MeapsData@fuites,
    shuf = shuf,
    attraction = attraction,
    parametres = parametres,
    xr_odds = odds,
    nthreads = ntr, verbose = verbose) 
  
  res <- list(
    flux = tibble::tibble(
      fromidINS = MeapsData@triplet$fromidINS,
      toidINS = MeapsData@triplet$toidINS,
      flux = res$flux) |> 
      filter(flux>0) |> 
      arrange(desc(flux)))
  if(large) gc()
  return(res)
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
multishuf_grouped <- function(MeapsDataGroup, attraction = "constant", 
                              parametres = 0, odds = 1, nshuf = 16,
                              nthreads = 0L, verbose = TRUE) {
  
  res <- multishuf(MeapsDataGroup, attraction = attraction, 
                   parametres = parametres, odds = odds, nshuf = nshuf,
                   nthreads = nthreads, verbose = verbose) 

  res_grouped <- res$flux |>
    dplyr::mutate(group_from = MeapsDataGroup@group_from[fromidINS],
           group_to = MeapsDataGroup@group_to[toidINS]) |> 
    dplyr::group_by(group_from, group_to) |> 
    dplyr::summarize(flux = sum(flux)) |> 
    dplyr::left_join(MeapsDataGroup@cible, by = c("group_from", "group_to")) |> 
    dplyr::arrange(desc(flux))
  
  return(list(flux = res_grouped,
              kl = kl(res_grouped$flux, res_grouped$value)))
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
multishuf_oc <- function(MeapsData, attraction = "constant", 
                         parametres = 0, odds = 1, nshuf = 16,
                         nthreads = 0L, verbose = TRUE, gbperthreads=4) {
  
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
  
  # check mem
  size <- length(MeapsData@triplet$metric)/1024^3
  large <- 4*size>gbperthreads/4
  ntr <- nthreads
  if(large) {
    gc()
    memused <- as.numeric(lobstr::mem_used())/1024^3
    if(nthreads==0) ntr <- max_threads()
    memleft <- gbperthreads*max_threads() - 4*size*ntr -20*size - memused
    if(memleft < ntr) {
      ntr <- min(max_threads(), max(1, round((gbperthreads*max_threads() - 20*size - memused)/(1.5+4*size))))
      cli::cli_warn("le nombre de threads est réduit à {ntr}")
    }
  }
  
  res <- multishuf_oc_cpp(
    jr_dist = MeapsData@j_dist,
    p_dist = MeapsData@p_dist,
    xr_dist = MeapsData@triplet$metric,
    emplois = MeapsData@emplois,
    actifs = MeapsData@actifs,
    fuites = MeapsData@fuites,
    shuf = shuf,
    attraction = attraction,
    parametres = parametres,
    xr_odds = odds,
    nthreads = ntr, verbose = verbose)
  
  res <- list(
    flux = tibble::tibble(
      fromidINS = MeapsData@triplet$fromidINS,
      toidINS = MeapsData@triplet$toidINS,
      flux = res$flux) |> 
      dplyr::filter(flux>0) |> 
      dplyr::arrange(desc(flux)))
  if(large) gc()
  return(res)
}

#' @import dplyr
all_in_grouped <- function(MeapsDataGroup,  attraction = "constant", 
                           parametres = 0, odds = 1, 
                           nthreads = 0L, verbose = TRUE) {
  
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.") 
  # Validation des paramètres
  if (!attraction %in% c("constant", "marche", "marche_liss", "double_marche_liss","decay", "logistique")) cli::cli_abort("Fonction attraction inconnue")
  
  # RQ : pas de méthode pour vérifier le bon ordre des odds.
  if (attraction == "odds" && length(odds) != nrow(MeapsDataGroup@triplet) ) 
    cli::cli_abort("vecteur odds invalide.") 
  if (attraction == "marche" && (length(parametres) != 2 || !is.numeric(parametres))) 
    cli::cli_abort("Parametres pour marche invalide.")
  if (attraction == "marche_liss" && (length(parametres) != 2 || !is.numeric(parametres))) 
    cli::cli_abort("Parametres pour marche_liss invalide.")
  if (attraction == "double_marche_liss" && (length(parametres) != 4 || !is.numeric(parametres))) 
    cli::cli_abort("Parametres pour double_marche_liss invalide.")
  if (attraction == "decay" && (length(parametres) != 2 || !is.numeric(parametres))) 
    cli::cli_abort("Parametres pour decay invalide.")
  if (attraction == "logistique" && (length(parametres) != 3 || !is.numeric(parametres))) 
    cli::cli_abort("Parametres pour logistique invalide.")
  
  cible <- MeapsDataGroup@cible |> 
    tidyr::complete(group_from, group_to, fill = list(value = 0)) |> 
    dplyr::arrange(group_from, group_to) |> 
    dplyr::pull(value)
  
  res <- all_in_optim(jr_dist = MeapsDataGroup@j_dist,
                      p_dist = MeapsDataGroup@p_dist,
                      xr_dist = MeapsDataGroup@triplet$metric,
                      group_from = MeapsDataGroup@index_gfrom,
                      group_to = MeapsDataGroup@index_gto,
                      emplois = MeapsDataGroup@emplois,
                      actifs = MeapsDataGroup@actifs,
                      fuites = MeapsDataGroup@fuites,
                      cible = cible,
                      parametres = parametres,
                      xr_odds = odds,
                      attraction = attraction,
                      nthreads = nthreads, verbose = verbose)
  
  g_froms <- unique(MeapsDataGroup@group_from) |> sort()
  g_tos <- unique(MeapsDataGroup@group_to) |> sort()
  flux <- tibble::tibble(
    group_from = g_froms[res$i+1L],
    group_to = g_tos[res$j+1L],
    flux = res$flux) |> 
    dplyr::filter(flux>0)
  if( !is.null(MeapsDataGroup@cible) )
    flux <- flux |> left_join(MeapsDataGroup@cible |> rename(cible=value),
                              by=c("group_from", "group_to")) |>
    dplyr::mutate(flux = replace_na(flux, 0),
           cible = replace_na(cible, 0)) |> 
    dplyr::filter(flux>0|cible>0) 
  
  flux <- flux |> dplyr::arrange(dplyr::desc(flux))
  
  return( list(flux = flux, kl = res$kl) )
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
multishuf_oc_grouped <- function(
    MeapsDataGroup,  attraction = "constant", parametres = 0, odds = 1, 
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
  
  cible <- MeapsDataGroup@cible |> 
    tidyr::complete(group_from, group_to, fill = list(value = 0)) |> 
    dplyr::arrange(group_from, group_to) |> 
    dplyr::pull(value)
  
  res <- multishuf_oc_group_cpp(
    jr_dist = MeapsDataGroup@j_dist,
    p_dist = MeapsDataGroup@p_dist,
    xr_dist = MeapsDataGroup@triplet$metric,
    group_from = MeapsDataGroup@index_gfrom,
    group_to = MeapsDataGroup@index_gto,
    emplois = MeapsDataGroup@emplois,
    actifs = MeapsDataGroup@actifs,
    fuites = MeapsDataGroup@fuites,
    shuf = MeapsDataGroup@shuf,
    cible = cible,
    parametres = parametres,
    xr_odds = odds,
    attraction = attraction,
    nthreads = nthreads, verbose = verbose)
  
  g_froms <- unique(MeapsDataGroup@group_from) |> sort()
  g_tos <- unique(MeapsDataGroup@group_to) |> sort()
  flux <- tibble::tibble(
    group_from = g_froms[res$i+1L],
    group_to = g_tos[res$j+1L],
    flux = res$flux) |> 
    dplyr::filter(flux>0)
  if(!is.null(MeapsDataGroup@cible))
    flux <- flux |> dplyr::left_join(MeapsDataGroup@cible |> dplyr::rename(cible=value),
                              by=c("group_from", "group_to")) |>
    dplyr::mutate(flux = tidyr::replace_na(flux, 0),
           cible = tidyr::replace_na(cible, 0)) |> 
    dplyr::filter(flux>0|cible>0) 
  
  flux <- flux |> dplyr::arrange(dplyr::desc(flux))
  
  return( list(flux = flux, kl = res$kl) )
}


meaps_optim <- function(MeapsDataGroup,  attraction, parametres, odds = 1,
                        meaps_fun = "all_in",
                        method = "L-BFGS-B", objective = "KL",
                        lower = NULL, upper = NULL, control = NULL,
                        nthreads = 0L, progress = TRUE,
                        quiet = TRUE) { 
  
  if (!inherits(MeapsDataGroup, "MeapsDataGroup"))
    cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.") 
  if (!attraction %in% c("marche", "marche_liss", "double_marche_liss","decay", "logistique")) 
    cli::cli_abort("Pas de fonction choisi ou fonction non paramètrique.")
  if (!meaps_fun %in% c("all_in", "multishuf", "multishuf_discret")) 
    cli::cli_abort("meaps_fun doit être soit 'all_in', soit 'multishuf'")
  
  arg <- list(
    MeapsDataGroup,
    attraction = attraction,
    odds = odds, nthreads = nthreads, verbose = FALSE)
  meaps_fun_ <- switch(
    meaps_fun,
    "all_in" = all_in_grouped,
    "multishuf" = multishuf_oc_grouped,
    "multishuf_discret" = multishuf_grouped)
  env <- environment()
  
  fn <- switch(
    objective, 
    "KL" = function(par) {
      estim <- do.call(
        meaps_fun_, 
        args = append(arg, list(parametres = par)))
      kl <- estim$kl
      mes <- glue::glue("kl:{signif(kl, 4)} ; {str_c(signif(par,4), collapse=', ')}")
      if (progress) 
        cli::cli_progress_update(.envir = env)
      if(!quiet) cli::cli_progress_output(mes, .envir = env)
      return(kl)
    }
  )
  
  if (is.null(fn)) 
    cli::cli_abort("Erreur : la fonction objective est non définie")
  
  if(attraction == "marche") {
    nb_par <- length(parametres)-1
    
    if (is.null(lower)) lower <- rep(0, nb_par)
    if (is.null(upper)) upper <- rep(Inf, nb_par)
    d_min <- lower[[1]]
    d_max <- upper[[1]]
    if(d_min==0) d_min <- 5
    if(d_max>100) d_max <- 15
    cli::cli_progress_bar(.envir = env, clear = FALSE)
    bf <- purrr::map_dfr(d_min:d_max, \(d) { 
      opt <- stats::optim(
        par = parametres[[2]],
        fn = \(x) fn(c(d, x)),
        method = "Brent", lower = lower[[2]], upper = upper[[2]])
      tibble::tibble(d = d, x = opt$par, kl = opt$value, 
             convergence = opt$convergence, mes = opt$message)})
    cli::cli_progress_done(.envir = env)
    best <- bf |> dplyr::filter(kl == min(kl)) |> dplyr::slice(1)  
    res <- list(
      par = c(best$d, best$x),
      value = best$kl,
      counts = NA, 
      convergence = best$convergence,
      message = best$message,
      all_iter = bf
    )
    return(res)  
  }
  
  nb_par <- length(parametres)
  if (is.null(lower)) lower <- rep(0, nb_par)
  if (is.null(upper)) upper <- rep(Inf, nb_par)
  cli::cli_progress_bar(.envir = env, clear = FALSE)
  res <- stats::optim(
    par = parametres,
    fn = fn,
    method = method, lower = lower, upper = upper, control = control)
  cli::cli_progress_done(.envir = env)
  res
}

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
  
  p_dist <- MeapsData@triplet |> dplyr::group_by(fromidINS) |> dplyr::summarize(n()) |> dplyr::pull() |> cumsum()
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
