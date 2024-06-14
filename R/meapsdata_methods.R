## LES METHODES DE BASES

#' Multishuf, l'original dans sa version discrête
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
multishuf_origin <- function(MeapsData, attraction = "constant", parametres = 0, odds = 1, nshuf = 16,
                             nthreads = 0L, verbose = TRUE) {
  if (!inherits(MeapsData, "MeapsData")) cli::cli_abort("Ce n'est pas un objet MeapsData.")
  
  # RQ : pas de méthode pour vérifier le bon ordre des odds.
  if (attraction == "odds") {
    if (length(odds) != nrow(object@triplet)) cli::cli_abort("vecteur odds invalide.")
  } else {
    check_fct_attraction(attraction, parametres)
  }
  
  froms <- MeapsData@froms
  tos <- MeapsData@tos
  
  if (is.null(MeapsData@shuf)) {
    if (nshuf <= 0) cli::cli_abort("pas de shufs et pas d'envie de shuf")
    shuf <- emiette(MeapsData@actifs, nshuf = nshuf)
  } else {
    shuf <- MeapsData@shuf
  }
  
  res <- multishuf_origin_cpp(
    jr_dist = MeapsData@jr_dist,
    p_dist = MeapsData@p_dist,
    xr_dist = MeapsData@triplet$metric,
    emplois = MeapsData@emplois,
    actifs = MeapsData@actifs,
    fuites = MeapsData@fuites,
    shuf = shuf,
    attraction = attraction,
    parametres = parametres,
    xr_odds = odds,
    nthreads = nthreads, verbose = verbose) |>
    tibble::as_tibble() |>
    dplyr::bind_cols(MeapsData@triplet |> select(-metric)) |>
    #   dplyr::filter(flux > 0) |>
    dplyr::relocate(fromidINS, toidINS, flux) 
}

#' Multishuf, réécrit en version continue avec une clause omp task depend.
multishuf_task <- function(MeapsData, attraction = "constant", parametres = 0, nshuf = 16, seuil = 40, nthreads = 0L, verbose = TRUE) {
  if (!inherits(MeapsData, "MeapsData")) cli::cli_abort("Ce n'est pas un objet MeapsData.")
  
  check_fct_attraction(attraction, parametres)
  
  multishuf_task_cpp(
    jr_dist = MeapsData@jr_dist,
    p_dist = MeapsData@p_dist,
    xr_dist = MeapsData@triplet$metric,
    emplois = MeapsData@emplois,
    actifs = MeapsData@actifs,
    fuites = MeapsData@fuites,
    parametres = parametres,
    shuf = MeapsData@shuf,
    attraction = attraction,
    nthreads = nthreads, verbose = verbose) |>
    tibble::as_tibble() |>
    dplyr::bind_cols(MeapsData@triplet |> select(-metric)) |>
    #   dplyr::filter(flux > 0) |>
    dplyr::relocate(fromidINS, toidINS, flux) 
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
                         nthreads = 0L, verbose = TRUE, gbperthreads = 4) {
  if (!inherits(MeapsData, "MeapsData")) cli::cli_abort("Ce n'est pas un objet MeapsData.")
  
  # RQ : pas de méthode pour vérifier le bon ordre des odds.
  if (attraction == "odds") {
    if (length(odds) != nrow(object@triplet)) cli::cli_abort("vecteur odds invalide.")
  } else {
    check_fct_attraction(attraction, parametres)
  }
  
  froms <- MeapsData@froms
  tos <- MeapsData@tos
  
  if (is.null(MeapsData@shuf)) {
    if (nshuf <= 0) cli::cli_abort("pas de shufs et pas d'envie de shuf")
    shuf <- emiette(MeapsData@actifs, nshuf = nshuf)
  } else {
    shuf <- MeapsData@shuf
  }
  
  # check mem
  size <- length(MeapsData@triplet$metric) / 1024^3
  large <- 4 * size > gbperthreads / 4
  ntr <- nthreads
  if (large) {
    gc()
    memused <- as.numeric(lobstr::mem_used()) / 1024^3
    if (nthreads == 0) ntr <- max_threads()
    memleft <- gbperthreads * max_threads() - 4 * size * ntr - 20 * size - memused
    if (memleft < ntr) {
      ntr <- min(max_threads(), max(1, round((gbperthreads * max_threads() - 20 * size - memused) / (1.5 + 4 * size))))
      cli::cli_warn("le nombre de threads est réduit à {ntr}")
    }
  }
  
  res <- multishuf_oc_cpp(
    jr_dist = MeapsData@jr_dist,
    p_dist = MeapsData@p_dist,
    xr_dist = MeapsData@triplet$metric,
    emplois = MeapsData@emplois,
    actifs = MeapsData@actifs,
    fuites = MeapsData@fuites,
    shuf = shuf,
    parametres = parametres,
    xr_odds = odds,
    attraction = attraction,
    nthreads = ntr, verbose = verbose
  )
  
  res <- list(
    flux = tibble::tibble(
      fromidINS = MeapsData@triplet$fromidINS,
      toidINS = MeapsData@triplet$toidINS,
      flux = res$flux) |>
      #     dplyr::filter(flux > 0) |>
      dplyr::arrange(dplyr::desc(flux))
  )
  if (large) gc()
  return(res)
}

#' Fonction all_in. Applique rmeaps dans sa version "tous les actifs partent en même temps".
#' @param MeapsData Un jeu de données MeapsData décrivant la zone d'étude.
#' @param attraction Choix du mode d'attraction. Par défaut, "constant". Les autres modes possibles sont : "marche", "marche_liss", "decay" et "logistique". Chacune applique une pénalité selon la distance (au-delà de la simple logique des rangs, qui correspond au mode "constant").
#' @param parametres Un vecteur avec les valeurs initiales des paramètres. Dépend du mode d'attraction retenu.
#' @param nthreads Le nombre de threads pour mener les calculs parallélisés. Par defaut, le nombre maximal.
#' @param verbose TRUE par defaut.
#'
#' @import dplyr
all_in <- function(MeapsData, attraction = "constant", parametres = 0, nthreads = 0L, verbose = TRUE) {
  if (!inherits(MeapsData, "MeapsData")) cli::cli_abort("Ce n'est pas un objet MeapsData.")
  
  check_fct_attraction(attraction, parametres)
  
  meaps_all_in_cpp(
    jr_dist = MeapsData@jr_dist,
    p_dist = MeapsData@p_dist,
    xr_dist = MeapsData@triplet$metric,
    emplois = MeapsData@emplois,
    actifs = MeapsData@actifs,
    fuites = MeapsData@fuites,
    parametres = parametres,
    attraction = attraction,
    nthreads = nthreads,
    verbose = verbose) |>
    tibble::as_tibble() |>
    dplyr::bind_cols(MeapsData@triplet |> dplyr::select(-metric)) |>
    dplyr::relocate(fromidINS, toidINS, flux)
}

## UN WRAPPER : "one method fits all"
meaps <- function(
    MeapsData, version = "all_in", attraction = "constant", parametres = 0, nthreads = 0L, verbose = TRUE,
    nshuf = NULL, odds = NULL, seuil = NULL, gbperthreads = NULL) {
  
  if (!inherits(MeapsData, "MeapsData")) cli::cli_abort("Ce n'est pas un objet MeapsData.")
  
  # RQ : le check des odds est sommaire. Ordre non vérifié.
  if (attraction == "odds") {
    if (length(odds) != nrow(object@triplet)) cli::cli_abort("vecteur odds invalide.")
  } else {
    check_fct_attraction(attraction, parametres)
  }
  
  arg <- list(
    jr_dist = MeapsData@jr_dist,
    p_dist = MeapsData@p_dist,
    xr_dist = MeapsData@triplet$metric,
    emplois = MeapsData@emplois,
    actifs = MeapsData@actifs,
    fuites = MeapsData@fuites,
    parametres = parametres,
    attraction = attraction,
    nthreads = nthreads, 
    verbose = verbose
  )
  if (!is.null(odds)) arg <- append(arg, list(odds = odds))
  if (!is.null(nshuf)) arg <- append(arg, list(nshuf = nshuf))
  if (!is.null(gbperthreads)) arg <- append(arg, list(gbperthreads = gbperthreads))
  
  meaps_fun <- switch(version,
                      "all_in" = all_in_grouped,
                      "multishuf_oc" = multishuf_oc_grouped,
                      "multishuf_task" = multishuf_task,
                      "multishuf_origin" = multishuf_origin
  )
  
  do.call(meaps_fun, arg)
}

## METHODES POUR LES GROUPES
meaps_grouped <- function(MeapsDataGroup, version = "all_in", attraction = "constant",
                          parametres = 0, odds = NULL,
                          nthreads = 0L, verbose = TRUE) {
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.")
  if (attraction == "odds") {
    if (length(odds) != nrow(object@triplet)) cli::cli_abort("vecteur odds invalide.")
  } else {
    check_fct_attraction(attraction, parametres)
  }
  
  if (!is.null(MeapsDataGroup@cible)) {
    cible <- MeapsDataGroup@cible |> pull(value)
  } else {
    cible <- NULL
  }
  
  arg <- list(
    jr_dist = MeapsDataGroup@jr_dist,
    p_dist = MeapsDataGroup@p_dist,
    xr_dist = MeapsDataGroup@triplet$metric,
    group_from = MeapsDataGroup@index_gfrom,
    group_to = MeapsDataGroup@index_gto,
    emplois = MeapsDataGroup@emplois,
    actifs = MeapsDataGroup@actifs,
    fuites = MeapsDataGroup@fuites,
    cible = cible,
    parametres = parametres,
    attraction = attraction,
    nthreads = nthreads, verbose = verbose
  )
  if (!is.null(odds)) arg <- append(arg, list(odds = odds))
  if (!is.null(nshuf)) arg <- append(arg, list(nshuf = nshuf))
  if (!is.null(gbperthreads)) arg <- append(arg, list(gbperthreads = gbperthreads))
  
  meaps_fun <- switch(version,
                      "all_in" = all_in_grouped,
                      "multishuf_oc" = multishuf_oc_grouped,
                      "multishuf_task" = multishuf_task
  )
  
  res <- do.call(meaps_fun, arg)
  
  g_froms <- unique(MeapsDataGroup@group_from) |> sort()
  g_tos <- unique(MeapsDataGroup@group_to) |> sort()
  
  flux <- tibble::tibble(
    group_from = g_froms[res$i + 1L],
    group_to = g_tos[res$j + 1L],
    flux = res$flux,
    cible = cible
  ) |>
    dplyr::arrange(dplyr::desc(flux))
  
  if (is.null(cible)) {
    return(list(flux = flux))
  } else {
    return(list(flux = flux, kl = res$kl))
  }
}
#' Multishuf oc groupé
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
    MeapsDataGroup, attraction = "constant", parametres = 0, odds = 1,
    nthreads = 0L, verbose = TRUE) {
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) {
    cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.")
  }
  # RQ : pas de méthode pour vérifier le bon ordre des odds.
  
  if (attraction == "odds") {
    if (length(odds) != nrow(object@triplet)) cli::cli_abort("vecteur odds invalide.")
  } else {
    check_fct_attraction(attraction, parametres)
  }
  
  cible <- MeapsDataGroup@cible |> 
    tidyr::complete(group_from, group_to, fill = list(value = 0)) |> 
    dplyr::arrange(group_from, group_to) |> 
    dplyr::pull(value)
  
  res <- multishuf_oc_group_cpp(
    jr_dist = MeapsDataGroup@jr_dist,
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
    nthreads = nthreads, verbose = verbose
  )
  
  g_froms <- unique(MeapsDataGroup@group_from) |> sort()
  g_tos <- unique(MeapsDataGroup@group_to) |> sort()
  
  flux <- tibble::tibble(
    group_from = g_froms[res$i + 1L],
    group_to = g_tos[res$j + 1L],
    flux = res$flux) |>
    dplyr::filter(flux > 0)
  if (!is.null(MeapsDataGroup@cible)) {
    flux <- flux |>
      dplyr::left_join(MeapsDataGroup@cible |> dplyr::rename(cible = value),
                       by = c("group_from", "group_to")) |>
      dplyr::mutate(
        flux = tidyr::replace_na(flux, 0),
        cible = tidyr::replace_na(cible, 0)) |>
      dplyr::filter(flux > 0 | cible > 0)
  }
  
  flux <- flux |> dplyr::arrange(dplyr::desc(flux))
  
  return(list(flux = flux, kl = res$kl, lk=res$lk))
}

#' Fonction all_in groupé.
#' Si MeapsDataGroup a défini une cible, renvoie également l'entropie relative (kl).
#' @import dplyr
all_in_grouped <- function(MeapsDataGroup, attraction = "constant",
                           parametres = 0, odds = 1,
                           nthreads = 0L, verbose = TRUE) {
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.")
  check_fct_attraction(attraction, parametres)
  
  if (!is.null(MeapsDataGroup@cible)) {
    cible <- MeapsDataGroup@cible |> dplyr::pull(value)
  } else {
    cible <- NULL
  }
  
  res <- meaps_all_in_cpp(
    jr_dist = MeapsDataGroup@jr_dist,
    p_dist = MeapsDataGroup@p_dist,
    xr_dist = MeapsDataGroup@triplet$metric,
    group_from = MeapsDataGroup@index_gfrom,
    group_to = MeapsDataGroup@index_gto,
    emplois = MeapsDataGroup@emplois,
    actifs = MeapsDataGroup@actifs,
    fuites = MeapsDataGroup@fuites,
    cible = cible,
    parametres = parametres,
    attraction = attraction,
    nthreads = nthreads, verbose = verbose
  )
  
  g_froms <- unique(MeapsDataGroup@group_from) |> sort()
  g_tos <- unique(MeapsDataGroup@group_to) |> sort()
  
  flux <- tibble::tibble(
    group_from = g_froms[res$i + 1L],
    group_to = g_tos[res$j + 1L],
    flux = res$flux,
    cible = cible) |>
    dplyr::arrange(dplyr::desc(flux))
  
  if (is.null(cible)) {
    return(list(flux = flux))
  } else {
    return(list(flux = flux, kl = res$kl))
  }
}


#' Fonction all_in groupé.
#' Si MeapsDataGroup a défini une cible, renvoie également l'entropie relative (kl).
#' @import dplyr
multishuf_task_grouped <- function(MeapsDataGroup, attraction = "constant",
                                   parametres = 0, odds = 1,
                                   nthreads = 0L, verbose = TRUE) {
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.")
  check_fct_attraction(attraction, parametres)
  
  if (!is.null(MeapsDataGroup@cible)) {
    cible <- MeapsDataGroup@cible |> 
      tidyr::complete(group_from, group_to, fill = list(value = 0)) |> 
      dplyr::arrange(group_from, group_to) |> 
      dplyr::pull(value)
  } else {
    cible <- NULL
  }
  
  res <- multishuf_task_cpp(
    jr_dist = MeapsDataGroup@jr_dist,
    p_dist = MeapsDataGroup@p_dist,
    xr_dist = MeapsDataGroup@triplet$metric,
    group_from = MeapsDataGroup@index_gfrom,
    group_to = MeapsDataGroup@index_gto,
    emplois = MeapsDataGroup@emplois,
    actifs = MeapsDataGroup@actifs,
    fuites = MeapsDataGroup@fuites,
    cible = cible,
    parametres = parametres,
    shuf = MeapsDataGroup@shuf,
    attraction = attraction,
    nthreads = nthreads, verbose = verbose
  )
  
  g_froms <- unique(MeapsDataGroup@group_from) |> sort()
  g_tos <- unique(MeapsDataGroup@group_to) |> sort()
  
  flux <- tibble::tibble(
    group_from = g_froms[res$i + 1L],
    group_to = g_tos[res$j + 1L],
    flux = res$flux,
    cible = cible
  ) |>
    dplyr::arrange(dplyr::desc(flux))
  
  if (is.null(cible)) {
    return(list(flux = flux))
  } else {
    return(list(flux = flux, kl = res$kl))
  }
}

## METHODE OPTIM

#' Wrapper pour optimisation des différentes fonctions Meaps.
meaps_optim <- function(MeapsDataGroup, attraction, parametres, odds = NULL,
                        version = "all_in",
                        method = "L-BFGS-B", objective = "KL",
                        lower = NULL, upper = NULL, control = NULL,
                        discret = NULL,
                        nthreads = 0L, progress = TRUE,
                        quiet = TRUE) {
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) {
    cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.")
  }
  check_fct_attraction(attraction, parametres)
  
  if (!is.null(MeapsDataGroup@cible)) {
    cible <- MeapsDataGroup@cible |> dplyr::pull(value)
  } else {
    cible <- NULL
  }
  
  arg <- list(
    MeapsDataGroup,
    attraction = attraction,
    nthreads = nthreads, verbose = FALSE
  )
  if (!is.null(odds)) arg <- append(arg, list(odds = odds))
  
  if(!version %in%c("all_in", "multishuf_oc", "multishuf_task"))
    cli::cli_abort("moteur MEAPS inconnu")
  
  meaps_fun_ <- switch(version,
                       "all_in" = all_in_grouped,
                       "multishuf_oc" = multishuf_oc_grouped,
                       "multishuf_task" = multishuf_task_grouped
  )
  
  env <- environment()
  
  fn <- switch(
    objective,
    "KL" = function(par) {
      if(any(par<lower)|any(par>upper)) {
        kl <- Inf} else {
          estim <- do.call(
            meaps_fun_,
            args = append(arg, list(parametres = par))
          )
          kl <- estim$kl}
      mes <- glue("kl:{signif(kl, 4)} ; {str_c(signif(par,4), collapse=', ')}")
      if (progress) {
        cli::cli_progress_update(.envir = env, extra = list(mes=mes))
      }
      return(kl)
    }
  )
  
  if (is.null(fn)) {
    cli::cli_abort("Erreur : la fonction objective est non définie")
  }
  
  if (!is.null(discret)) {
    nb_par <- length(parametres)
    
    if (is.null(lower)) lower <- rep(0, nb_par)
    if (is.null(upper)) upper <- rep(Inf, nb_par)
    lower <- tail(lower, -1)
    upper <- tail(upper, -1)
    tp <- progress
    progress <- FALSE
    bf <- purrr::map_dfr(discret, \(d) {
      opt <- stats::optim(
        par = tail(parametres, -1),
        fn = \(x) fn(c(d, x)),
        method = "Brent", 
        lower = lower, upper = upper)
      tibble::tibble(
        d = d, x = opt$par, kl = opt$value,
        convergence = opt$convergence, mes = opt$message)
    }, .progress= tp)
    
    best <- bf |>
      dplyr::filter(kl == min(kl)) |>
      dplyr::slice(1)
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
  
  cli::cli_progress_bar(
    format = "{cli::pb_spin} Estimation de MEAPS {cli::pb_current} en {cli::pb_elapsed} ({round(1/cli::pb_rate_raw)} s/iter) ; {cli::pb_extra$mes}", 
    extra = list(mes=""),
    .envir = env, clear = FALSE)
  
  res <- stats::optim(
    par = parametres,
    fn = fn,
    method = method, lower = lower, upper = upper, control = control
  )
  cli::cli_progress_done(.envir = env)
  
  return(res)
}

#' version optimx.
#' Beaucoup plus de méthodes, dont certaines recommandées : nvm(), ncg(), ucminf::ucminf(), lbfgs::lbfgs(), nloptr::lbfgs()
#' @import optimx
#' @import glue
#' @import cli
#' @import dplyr
meaps_optimx <- function(MeapsDataGroup, attraction, parametres, 
                         version = "all_in",
                         method = "lbfgsb3c", objective = "KL",
                         amplitude_max = 1000, rayon_max = 10,
                         lower = NULL, upper = NULL, control = NULL,
                         discret = NULL,
                         nthreads = 0L, progress = TRUE,
                         quiet = TRUE) {
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) {
    cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.")
  }
  check_fct_attraction(attraction, parametres)
  
  if (!is.null(MeapsDataGroup@cible)) {
    cible <- MeapsDataGroup@cible |> dplyr::pull(value)
  } else {
    cible <- NULL
  }
  
  arg <- list(
    MeapsDataGroup,
    attraction = attraction,
    nthreads = nthreads, verbose = FALSE
  )
  
  if(!version %in%c("all_in", "multishuf_oc", "multishuf_task"))
    cli::cli_abort("moteur MEAPS inconnu")
  
  meaps_fun_ <- switch(
    version,
    "all_in" = all_in_grouped,
    "multishuf_oc" = multishuf_oc_grouped,
    "multishuf_task" = multishuf_task_grouped
  )
  
  env <- environment()
  
  fn <- switch(
    objective,
    "KL" = function(par) {
      if(any(par<lower)|any(par>upper)) {
        kl <- Inf} else {
          estim <- do.call(
            meaps_fun_,
            args = append(arg, list(parametres = par))
          )
          kl <- estim$kl}
      mes <- glue("kl:{signif(kl, 4)} ; {str_c(signif(par,4), collapse=', ')}")
      if (progress) {
        cli::cli_progress_update(.envir = env, extra = list(mes=mes))
      }
      return(kl)
    }
  )
  
  if (is.null(fn)) {
    cli::cli_abort("Erreur : la fonction objective est non définie")
  }
  
  if (!is.null(discret)) {
    nb_par <- length(parametres)
    
    if (is.null(lower)) lower <- rep(0, nb_par)
    if (is.null(upper)) upper <- rep(Inf, nb_par)
    tp <- progress
    progress <- FALSE
    bf <- purrr::map_dfr(discret, \(d) {
      opt <- stats::optim(
        par = tail(parametres, -1),
        fn = \(x) fn(c(d, x)),
        method = "Brent", 
        lower = tail(lower, -1), upper = tail(upper, -1))
      tibble::tibble(
        d = d, x = opt$par, kl = opt$value,
        convergence = opt$convergence, mes = opt$message)
    }, .progress= tp)
    
    best <- bf |>
      dplyr::filter(kl == min(kl)) |>
      dplyr::slice(1)
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
  if (is.null(upper)) {
    upper <- switch(
      attraction,
      "marche" = c(rayon_max, amplitude_max),
      "rampe" = c(rayon_max, amplitude_max),
      "grav_exp" = c(10 + log(1 + amplitude_max), amplitude_max),
      "grav_puiss" = c(5, .2, 2*amplitude_max)) # le paramètre p1 ne devrait pas dépasser 200m. On est en km.
  }
  
  cli::cli_progress_bar(
    format = "{cli::pb_spin} Estimation de MEAPS {cli::pb_current} en {cli::pb_elapsed} ({round(1/cli::pb_rate_raw)} s/iter) ; {cli::pb_extra$mes}", 
    extra = list(mes=""),
    .envir = env, clear = FALSE)
  
  res <- optimx::optimr(
    par = parametres,
    fn = fn,
    method = method, control = control
  )
  cli::cli_progress_done(.envir = env)
  
  return(res)
}


#' Refonte de l'optim en wrapper
#' STRATEGIE :
#' 1 = cmaes::cma_es = Covariance Matrix Adapting Evolutionary Strategy.
#' 2 = dfoptim::nmkb = Nelder-Mead derivative free algorithm.
#' 3 = dfoptim::mads = mesh adaptative direct search.
#' 
#' @import dfoptim
#' @import cmaes
#' @import glue
#' @import cli
#' @import dplyr
meaps_opt <- function(MeapsDataGroup, attraction, parametres, 
                         fct_meaps = "all_in", 
                         strategie = 1L, 
                         amplitude_max = 1000, rayon_max = 10,
                         lower = NULL, upper = NULL, control = NULL,
                         nthreads = 0L, progress = TRUE,
                         quiet = TRUE) {
  if (!inherits(MeapsDataGroup, "MeapsDataGroup")) cli::cli_abort("Ce n'est pas un objet MeapsDataGroup.")
  if (is.null(MeapsDataGroup@cible)) cli::cli_abort("Il n'y a pas de cible.")
  check_fct_attraction(attraction, parametres)
  
  arg <- list(
    MeapsDataGroup,
    attraction = attraction,
    nthreads = nthreads, verbose = FALSE
  )

  if(!fct_meaps %in%c("all_in", "multishuf_oc", "multishuf_task")) cli::cli_abort("moteur MEAPS inconnu")
  meaps_fun_ <- switch(fct_meaps,
                       "all_in" = all_in_grouped,
                       "multishuf_oc" = multishuf_oc_grouped,
                       "multishuf_task" = multishuf_task_grouped
  )
  
  env <- environment()
  kl <- 1

  fn <- function(par) {
    estim <- do.call(
      meaps_fun_,
      args = append(arg, list(parametres = par))
    )
    kl <- estim$kl
    min_kl <- min(get("kl", envir = env), kl)
    assign("kl", min_kl, envir = env)
    mes <- glue::glue("kl:{signif(kl, 4)} min:{signif(min_kl, 4)}% -> {str_c(signif(par,4), collapse=', ')}")
    if (progress) {
      cli::cli_progress_update(.envir = env, extra = list(mes=mes))
      }
    return(kl)
  }

  nb_par <- length(parametres)
  if (is.null(lower)) lower <- rep(0, nb_par)
  if (is.null(upper)) {
    upper <- switch(attraction,
      "marche" = c(rayon_max, amplitude_max),
      "rampe" = c(rayon_max, amplitude_max),
      "grav_exp" = c(10 + log(1 + amplitude_max), amplitude_max),
      "grav_puiss" = c(5, .2, 2 * amplitude_max)) # le paramètre p1 ne devrait pas dépasser 200m. On est en km.
    }
  
  optim_method <- switch(strategie,
    function(par) {
      cmaes::cma_es(par = par, fn = fn, lower = lower, upper = upper, control = control) },
    function(par)  {
      dfoptim::nmkb(par = par, fn = fn, lower = lower, upper = upper, control = control) },
    function(par)  {
      dfoptim::mads(par = par, fn = fn, lower = lower, upper = upper, control = control) }
  )

  cli::cli_progress_bar(
    format = "{cli::pb_spin} Estimation avec {attraction} n°{cli::pb_current} en {cli::pb_elapsed} : {cli::pb_extra$mes}", 
    extra = list(mes=""),
    .envir = env, clear = FALSE)
  
  res <- optim_method(par = parametres)

  cli::cli_progress_done(.envir = env)
  
  return(res)
}