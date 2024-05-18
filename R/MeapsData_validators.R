is_triplet_meaps <- function(object, quiet = TRUE) {
  
  if (!inherits(object, "data.frame") || length(object) != 3) {
    status <- 2L
  } else if (!setequal(names(object), c("fromidINS", "toidINS", "metric"))) {
    status <- 3L
  } else if (any(is.na(object))) { 
    status <- 4L 
  } else {
    status <- 1L
  }
  
  if (!quiet) { 
    switch(status,
           cli::cli_alert_success("Le triplet meaps est un triplet valide."),
           cli::cli_alert_warning("Le triplet meaps n'est pas un triplet."),
           cli::cli_alert_warning("Les noms des colonnes devraient être fromidINS, toidINS et metric."),
           cli::cli_alert_warning("Le triplet meaps contient des valeurs manquantes.")
    )
  }
  
  invisible(status==1)
}

is_triplet_ordered <- function(object, quiet = TRUE) {
  if ( is.unsorted(object$fromidINS) ) { 
    status <- 2L
  } else {
    sort_v <- object |> 
      dplyr::group_by(fromidINS) |> 
      dplyr::summarize(test = is.unsorted(metric)) |> 
      dplyr::pull() |> 
      any()
    if (sort_v) {
      status <- 3L
    } else {
      sort_j <- FALSE
      #object |> group_by(fromidINS, metric) |> summarize(test = is.unsorted(j)) |> pull() |> any()
      if (sort_j) { 
        status <- 4L 
      } else {
        status <- 1L
      }
    }
  }
  
  if (!quiet) { 
    switch(status,
           cli::cli_alert_success("Le triplet meaps est bien ordonné."),
           cli::cli_alert_warning("Le triplet meaps n'est pas ordonné selon fromidINS."),
           cli::cli_alert_warning("Le triplet meaps n'est pas ordonné selon metric"),
           cli::cli_alert_warning("Le triplet meaps n'est pas ordonné selon toidINS.")
    )
  }
  
  invisible(status==1)
}

is_meapsdata_valid <- function(object, quiet = TRUE) {
  status <- TRUE
  if (is.null(names(object@actifs))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("actifs n'a pas de labels") }
  }
  if (is.null(names(object@emplois))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("emplois n'a pas de labels") }
  }
  if (is.null(names(object@fuites))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("fuites n'a pas de labels") }
  }
  froms <- object@froms
  tos <- object@tos
  N <- froms |> length()
  K <- tos |> length()
  
  if (!setequal(froms, names(object@actifs)) || length(object@actifs) != N) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("triplet et actifs n'ont pas la même liste de fromidINS") }
  }
  if (!setequal(tos, names(object@emplois)) || length(object@emplois) != K) { 
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("triplet et emplois n'ont pas la même liste de toidINS") }
  }
  if (!setequal(froms, names(object@fuites)) || length(object@fuites) != N) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("triplet et fuites n'ont pas la même liste de fromidINS") }
  }
  
  if (any(object@fuites <= 0 | object@fuites >= 1)) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("Le vecteur fuite a des valeurs invalides.") }
  }
  
  if (status) { cli::cli_alert_success("Le MeapsData est valide.") }
  
  invisible(status)
}

is_meaps_normalized <- function(object, just_a_warning = TRUE, seuil = 1e-4) {
  delta <- (sum(object@actifs * (1 - object@fuites)) - sum(object@emplois))/sum(object@emplois)
  if (delta > seuil) {
    if (just_a_warning) { 
      cli::cli_alert_warning("Les actifs restant dans la zone et les emplois divergent de {round(100 * delta, digits = 2)}%.")
    } else {
      cli::cli_abort("Les actifs restant dans la zone et les emplois divergent de {round(100 * delta, digits = 2)}%.")
    }
  }
  invisible(delta <= seuil)
}

check_meapsdata <- function(object, abort = FALSE) {
  status <- is_triplet_meaps(object@triplet, quiet = FALSE) && 
    is_triplet_ordered(object@triplet, quiet = FALSE) && 
    is_meapsdata_valid(object, quiet = FALSE) &&
    is_meaps_normalized(object)
  
  if (abort && !status) {cli::cli_abort("Invalide.")} else { return(status)}
}

is_meapsdatagroup_valid <- function(object, quiet = FALSE)  {
  
  if (!inherits(object, "data.frame") || length(object) != 3) {
    status <- 2L
  } else if (!setequal(names(object), c("group_from", "group_to", "value"))) {
    status <- 3L
  } else if (any(is.na(object))) { 
    status <- 4L 
  } else {
    status <- 1L
  }
  
  if (!quiet) { 
    switch(status,
           cli::cli_alert_success("Le triplet cible est un triplet valide."),
           cli::cli_alert_warning("La cible n'est pas un triplet."),
           cli::cli_alert_warning("Les noms des colonnes devraient être group_from, group_to et value."),
           cli::cli_alert_warning("La cible contient des valeurs manquantes.")
    )
  }
  
  invisible(status==1)
}

check_meapsdatagroup <- function(object, abort = FALSE){
  
  status <- is_meapsdatagroup_valid(object@cible)
  
  if (length(object@group_from) != length(object@actifs)) {
    status <- FALSE
    cli::cli_alert_warning("group_from et actifs ne font pas la même longueur.")
  }
  if (length(object@group_to) != length(object@emplois)) {
    status <- FALSE
    cli::cli_alert_warning("group_to et emplois ne font pas la même longueur.")
  }
  
  if (is.null(names(object@group_from))) {
    status <- FALSE
    cli::cli_alert_warning("group_from n'est pas labelisé")
  }
  if (is.null(names(object@group_to))) {
    status <- FALSE
    cli::cli_alert_warning("group_to n'est pas labelisé.")
  }
  
  if (!setequal(names(object@group_from), names(object@actifs))) {
    status <- FALSE
    if (!quiet) { cli::cli_alert_warning("group_from et actifs n'ont pas la même liste de fromidINS.") }
  }
  if (!setequal(names(object@group_to), names(object@emplois))) {
    status <- FALSE
    cli::cli_alert_warning("group_to et emplois n'ont pas la même liste de toidINS.")
  }
  
  if (!setequal(object@cible$group_from, object@group_from)) {
    status <- FALSE
    cli::cli_alert_warning("group_from et cible ne correspondent pas.")
  }
  
  if (!setequal(object@cible$group_to, object@group_to)) {
    status <- FALSE
    cli::cli_alert_warning("group_to et cible ne correspondent pas.")
  }
  
  if (status) { 
    cli::cli_alert_success("Le MeapsDataGroup est valide.") 
  } else if (abort) {
    cli::cli_abort("Le MeapsDataGroup n'est pas valide.")
  }
  
  invisible(status)
}

