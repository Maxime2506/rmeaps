% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emiette.R
\name{reordonne_shuf}
\alias{reordonne_shuf}
\title{Recalcule la matrice shuf sur un nouvel ordre d'actifs (l'initial est dans les colnames)}
\usage{
reordonne_shuf(shuf, actifs)
}
\arguments{
\item{shuf}{le tri initial, avec les noms des actifs ordonnés dans les colonnes}

\item{actifs}{Un vecteur d'actifs dans un autre ordre}
}
\value{
une matrice d'arrangement (autant de lignes que nshuf, autant de colonnes que d'actifs, dont les colonnes sont l'idINS)
}
\description{
Recalcule la matrice shuf sur un nouvel ordre d'actifs (l'initial est dans les colnames)
}
