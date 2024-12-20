% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meapsdata_methods.R
\name{meaps_optim}
\alias{meaps_optim}
\title{Wrapper pour optimisation des différentes fonctions Meaps.}
\usage{
meaps_optim(
  MeapsDataGroup,
  attraction,
  parametres,
  odds = NULL,
  version = "all_in",
  method = "L-BFGS-B",
  objective = "KL",
  lower = NULL,
  upper = NULL,
  control = NULL,
  discret = NULL,
  nthreads = 0L,
  progress = TRUE,
  quiet = TRUE,
  klway = "kl",
  weights = NULL
)
}
\description{
Wrapper pour optimisation des différentes fonctions Meaps.
}
