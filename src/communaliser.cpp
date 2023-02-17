#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//' La fonction communaliser effectue un regroupement de la matrice des flux (résultant de meaps) selon des groupes origines et destinations.
//' @param flux La matrice des flux de fromid vers toid.
//' @param group_orig Un vecteur d'integers (par ex. code commune) donnant le code du groupe de départ.
//' @param group_dest Un vecteur d'integers (par ex. code commune) donnant le code du groupe d'arrivée.
//' 
//' @return Une matrice contenant les flux agrégés selon les groupes de départ et d'arrivée. 
// [[Rcpp::export]]
NumericMatrix communaliser(NumericMatrix flux,
                           IntegerVector group_orig, 
                           IntegerVector group_dest) {
  
  const int N = flux.nrow(), K = flux.ncol();
  
  if (group_orig.size() != N) { 
    stop("La taille du groupe origine ne correspond pas aux flux");
    }
  if (group_dest.size() != K) { 
    stop("La taille du groupe destination ne correspond pas aux flux");
  }
  
  set<int> communes_orig, communes_dest;
  for (int i = 0; i < N; ++ i) {
    communes_orig.insert(group_orig[i]);
  }
  for (int j = 0; j < K; ++ j) {
    communes_dest.insert(group_dest[j]);
  }
  
  const int Nagg = communes_orig.size(), Kagg = communes_dest.size();
  NumericMatrix agregat(Nagg, Kagg), inter(Nagg, K);
  
  int pos;
  // RowSums.
  for (int i = 0; i < N; ++i) {
    pos = distance(communes_orig.begin(), lower_bound(communes_orig.begin(), communes_orig.end(), group_orig(i)));
    for (int j = 0; j < K; ++j) {
      inter(pos, j) += flux(i,j);
    }
  }
  
  // ColSums.
  for (int j = 0; j < K; ++j) {
    pos = distance(communes_dest.begin(), lower_bound(communes_dest.begin(), communes_dest.end(), group_dest(j)));
    for (int i = 0; i < Nagg; ++i) {
      agregat(i, pos) += inter(i,j);
    }
  }
  agregat.attr("dimnames") = List::create(communes_orig, communes_dest);
  return agregat;
}
