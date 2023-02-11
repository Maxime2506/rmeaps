#include <Rcpp.h>
#include <algorithm>
#include "distribuer.h"
using namespace Rcpp;
using namespace std;

//' La fonction distribuer remplit un ensemble de conteneurs.
//' Cette version interne ne se préoccupe pas du reste non versé.
//' @param conteneurs Un vecteur donnant la capacité de chacun des conteneurs dans l'ordre voulu de remplissage.
//' @param proportion Un vecteur définissant les proportions des débits vers les conteneurs.
//' @param quantité La quantité totale à verser dans les conteneurs.
//' 
//' @return La répartition finale de la quantité dans les conteneurs.
vector<double> distribuer(const vector<double>& conteneurs, 
                const vector<double>& proportion, 
                const double& quantite) {
  
  int N = conteneurs.size();
  vector<double> contenus(N), debit(N);
  
  set<int> disponible;
  
  // Calcul du débit. 
  // 1. On ne tient pas compte du débit des conteneurs nuls.
  // 2. Les débits négatifs sont considérés comme nul.
  for (int i = 0; i < N; ++i) {
    contenus[i] = 0.0;
    if ((conteneurs[i] > 0.0) && (proportion[i] > 0.0)) {
      debit[i] = proportion[i];
      disponible.insert(i);
    } else {
      debit[i] = 0.0;
    }
  }
  
  double debit_tot = 0.0,
         volume_disp = 0.0;
  for (auto& i: disponible) {
    debit_tot += debit[i];
    volume_disp += conteneurs[i];
  }
  // Cas d'absence de débit.
  if (debit_tot == 0) {
    return debit; // les valeurs sont nulles.
  }
  
  // Cas où aucun conteneur n'a de place disponible.
  if (volume_disp == 0) {
    return contenus; // les valeurs sont nulles.
  }
  
  // On vérifie qu'il y a suffisamment de place, ou pas.
  double excedent = quantite - volume_disp;
  
  // Cas simple où ça déborde de partout.
  if (excedent >= 0.0) {
    return conteneurs;
  }
  
  // S'il n'y a pas de quoi tout remplir, on distribue la quantité.
  excedent = quantite;
  set<int>::iterator it, tmp;
  
  while (excedent > 0.0) {
    
    // Calcul du débit restant après coupure des robinets des conteneurs remplis.
    debit_tot = 0.0;
    for (auto& i: disponible) {
      debit_tot += debit[i];
    }
    
    // Remplissage de tout l'excédent, même si ça déborde.
    for (auto& i: disponible) {
      contenus[i] += excedent * debit[i] / debit_tot;
    }
    
    // Récupération de ce qui a débordé.
    excedent = 0.0;
    it = disponible.begin();                   
    
    while (it != disponible.end()) {
      if (contenus[*it] >= conteneurs[*it]) {
        tmp = it;
        ++tmp;
        excedent += contenus[*it] - conteneurs[*it];
        contenus[*it] = conteneurs[*it];
        disponible.erase(it);
        it = tmp;
      } else {
        ++it;
      }
    }
  }
  return contenus;
}

//' La fonction distribuer remplit un ensemble de conteneurs. 
//' Celle-ci est exportée et donne également le reste non distribué.
//' @param conteneurs Un vecteur donnant la capacité de chacun des conteneurs dans l'ordre voulu de remplissage.
//' @param proportion Un vecteur définissant les proportions des débits vers les conteneurs.
//' @param quantité La quantité totale à verser dans les conteneurs.
//' 
//' @return Une liste avec : part = la répartition finale de la quantité dans les conteneurs et reste = la quantité qui n'a pas pu être versée si tout est rempli.
// [[Rcpp::export]]
List distribuer(NumericVector conteneurs, 
                NumericVector proportion, 
                const double& quantite) {
  
  int N = conteneurs.size();
  vector<double> contenus(N), debit(N);
  
  set<int> disponible;
  
  // Calcul du débit. 
  // 1. On ne tient pas compte du débit des conteneurs nuls.
  // 2. Les débits négatifs sont considérés comme nul.
  for (int i = 0; i < N; ++i) {
    contenus[i] = 0.0;
    if ((conteneurs[i] > 0.0) && (proportion[i] > 0.0)) {
      debit[i] = proportion[i];
      disponible.insert(i);
    } else {
      debit[i] = 0.0;
    }
  }
  
  double debit_tot = 0.0,
    volume_disp = 0.0;
  for (auto& i: disponible) {
    debit_tot += debit[i];
    volume_disp += conteneurs[i];
  }
  // Cas d'absence de débit.
  if (debit_tot == 0) {
    return List::create(Named("part", debit), Named("reste", quantite));
  }
  
  // Cas où aucun conteneur n'a de place disponible.
  if (volume_disp == 0) {
    return List::create(Named("part", contenus), Named("reste", quantite));
  }
  
  // On vérifie qu'il y a suffisamment de place, ou pas.
  double excedent = quantite - volume_disp;
  
  // Cas simple où ça déborde de partout.
  if (excedent >= 0.0) {
    return List::create(Named("part", conteneurs), Named("reste", excedent));
  }
  
  // S'il n'y a pas de quoi tout remplir, on distribue la quantité.
  excedent = quantite;
  set<int>::iterator it, tmp;
  
  while (excedent > 0.0) {
    
    // Calcul du débit restant après coupure des robinets des conteneurs remplis.
    debit_tot = 0.0;
    for (auto& i: disponible) {
      debit_tot += debit[i];
    }
    
    // Remplissage de tout l'excédent, même si ça déborde.
    for (auto& i: disponible) {
      contenus[i] += excedent * debit[i] / debit_tot;
    }
    
    // Récupération de ce qui a débordé.
    excedent = 0.0;
    it = disponible.begin();                   
    
    while (it != disponible.end()) {
      if (contenus[*it] >= conteneurs[*it]) {
        tmp = it;
        ++tmp;
        excedent += contenus[*it] - conteneurs[*it];
        contenus[*it] = conteneurs[*it];
        disponible.erase(it);
        it = tmp;
      } else {
        ++it;
      }
    }
  }
  return List::create(Named("part", contenus), Named("reste", 0.0));
}
