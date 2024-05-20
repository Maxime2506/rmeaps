% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{kullback_leibler}
\alias{kullback_leibler}
\title{Calcul de la béta-divergence Kullback-Leibler, soit l'entropie relative.}
\usage{
kullback_leibler(nb_mod, nb_ref, seuil = 1e-06)
}
\arguments{
\item{nb_mod}{Vecteur dont on veut évaluer la divergence.}

\item{nb_ref}{Vecteur de référence.}

\item{seuil}{Seuil de bruit. Défaut : 1e-6.}
}
\value{
Valeur de Kullback-Leibler.
}
\description{
Calcul de la béta-divergence Kullback-Leibler, soit l'entropie relative.
}