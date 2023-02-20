#include <Rcpp.h>
#include "deborder.h"
using namespace Rcpp;
using namespace std;

//// ' La fonction déborder remplit une série de conteneurs dans l'ordre de présentation. Il passe au conteneur suivant si le précédent est rempli.
//// ' @param conteneurs Un vecteur donnant la capacité de chacun des conteneurs dans l'ordre voulu de remplissage.
//// ' @param quantité La quantité totale à verser dans les conteneurs.
//// '
//// ' @return Une liste avec : part = la répartition finale de la quantité dans les conteneurs et reste = la quantité qui n'a pas pu être versée si tout est rempli.
vector<double> deborder(
  const vector<double>& conteneurs,
  const double& quantite) {
  
  int N = conteneurs.size();
  vector<double> contenus(N);
  
  // Calcul du volume disponible sur les conteneurs positifs uniquement. 
  double volume_disp = 0.0;
  for (int i = 0; i < N; ++i) {
    contenus[i] = 0.0;
    if (conteneurs[i] > 0.0) {
      volume_disp += conteneurs[i];
      }
  }
  // Cas où aucun conteneur n'a de place disponible.
  // Cas où la quantité est nulle ou négative : on ne déborde pas.
  if (volume_disp == 0 || quantite <= 0) {
    return contenus; // Les valeurs sont nulles.
  }
  
  // On vérifie qu'il y a suffisamment de place, ou pas.
  double excedent = quantite - volume_disp;
  
  // Cas simple où ça déborde de partout.
  if (excedent >= 0.0) {
    return conteneurs;
  };
  
  // S'il n'y a pas de quoi tout remplir, on distribue la quantité dans l'ordre.
  excedent = quantite;
  int pos = 0;
  
  while (excedent > 0.0) {
    if (excedent <= conteneurs[pos]) {
      contenus[pos] = excedent;
    } else {
      if (conteneurs[pos] > 0.0) {
        contenus[pos] = conteneurs[pos];
        }
    } 
    excedent -= contenus[pos];
    ++pos;
  }
  
  return contenus;
}

//' La fonction déborder remplit une série de conteneurs dans l'ordre de présentation. Il passe au conteneur suivant si le précédent est rempli.
//' Celle-ci est exportée et donne le reste non distribué.
//' @param conteneurs Un vecteur donnant la capacité de chacun des conteneurs dans l'ordre voulu de remplissage.
//' @param quantité La quantité totale à verser dans les conteneurs.
//' 
//' @return Une liste avec : part = la répartition finale de la quantité dans les conteneurs et reste = la quantité qui n'a pas pu être versée si tout est rempli.
// [[Rcpp::export]]
List deborder(
    NumericVector conteneurs,
    const double& quantite) {
  
  int N = conteneurs.size();
  NumericVector contenus(N);
  
  // Calcul du volume disponible sur les conteneurs positifs uniquement. 
  double volume_disp = 0.0;
  for (int i = 0; i < N; ++i) {
    if (conteneurs[i] > 0.0) {
      volume_disp += conteneurs[i];
    }
  }
  // Cas où aucun conteneur n'a de place disponible.
  // Cas où la quantité est nulle ou négative : on ne déborde pas.
  if (volume_disp == 0 || quantite <= 0) {
    return List::create(Named("part", contenus), Named("reste", quantite));
  }
  
  // On vérifie qu'il y a suffisamment de place, ou pas.
  double excedent = quantite - volume_disp;
  
  // Cas simple où ça déborde de partout.
  if (excedent >= 0.0) {
    return List::create(Named("part", conteneurs), Named("reste", excedent));
  }
  
  // S'il n'y a pas de quoi tout remplir, on distribue la quantité dans l'ordre.
  excedent = quantite;
  int pos = 0;
  
  while (excedent > 0.0) {
    if (excedent <= conteneurs[pos]) {
      contenus[pos] = excedent;
    } else {
      if (conteneurs[pos] > 0.0) {
        contenus[pos] = conteneurs[pos];
      }
    } 
    excedent -= contenus[pos];
    ++pos;
  }
  
  return List::create(Named("part", contenus), Named("reste", 0.0));
}