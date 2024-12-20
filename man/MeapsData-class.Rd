% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meapsdata_class.R
\docType{class}
\name{MeapsData-class}
\alias{MeapsData-class}
\title{Définition d'une classe d'objet rmeaps de préparation des données.
MeapsData permet d'unifier un jeu de données en effectuant un examen de validation, garantissant que les données sont dans les bons formats et le bon ordre avant d'être invoqué dans la procédure rmeaps.}
\description{
Définition d'une classe d'objet rmeaps de préparation des données.
MeapsData permet d'unifier un jeu de données en effectuant un examen de validation, garantissant que les données sont dans les bons formats et le bon ordre avant d'être invoqué dans la procédure rmeaps.
}
\section{Slots}{

\describe{
\item{\code{triplet}}{Matrice des distances (ou temps ou encore coûts généralisés) entre des départs et des arrivées. La matrice est donnée sous forme de triplet (une data.frame avec les colonnes fromidINS, toidINS et metric). Un tri lexicographique sur fromidINS, puis metric, est attendu.}

\item{\code{actifs}}{Vecteur des résidents actifs au départ des fromidINS. Le vecteur doit être labélisé par les fromidINS.}

\item{\code{emplois}}{Vecteurs des emplois à chacune des destinations. Le vecteur doit être labélisé par les toidINS.}

\item{\code{fuites}}{Vecteurs des proportions d'actifs qui travaillent hors de la zone d'étude. Le vecteur doit être labélisé par les fromidINS.}

\item{\code{froms}}{Vecteur des fromidINS dans l'ordre attendu par les différentes méthodes (automatique).}

\item{\code{tos}}{Vecteur des toidINS dans l'ordre attendu par les différentes méthodes (automatique).}

\item{\code{shuf}}{Une matrice de shuf si besoin de définir des priorités.}

\item{\code{jr_dist}}{Vecteur des colonnes rangées par distance dans une représentation sparse matrix en ligne.}

\item{\code{p_dist}}{Vecteur du nombre d'éléments cumulés par ligne dans une représentation sparse matrix en ligne.}
}}

