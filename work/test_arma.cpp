#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export]]
inline void convertSparse(S4 mat) {
  
  // obtain dim, i, p. x from S4 object
  IntegerVector dims = mat.slot("Dim");
  arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
  arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
  arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));
  
  int nrow = dims[0], ncol = dims[1];
  
  // use Armadillo sparse matrix constructor
  arma::sp_mat res(i, p, x, nrow, ncol);
}


// meaps_oneboot
void meaps_oneboot(S4 rki, S4 modds) {
arma::sp_mat<int> ranking = convertSparse(S4 rki);
arma::sp_mat<double> odds;
odds = convertSparse(S4 modds);



// Ordre de départs des actifs vers les emplois. Un shuf pour un boot.
std::vector<std::size_t> ishuf(Ns);
for (std::size_t k = 0; k < Ns; ++k) {
  ishuf[k] = shuf(s, k) - 1L;
}

// Le vecteur shuf peut être plus long que le nombre de lignes de rkdist s'il fait repasser plusieurs fois
// la même ligne d'actifs. Dans ce cas, on compte la fréquence de passage de chaque ligne et l'on divise le
// poids de la ligne par cette fréquence.
std::vector<int> freq_actifs(N, 0L);
for (auto i: ishuf) {
  freq_actifs[i]++;
}

// Initialisation des emplois.
std::vector<double> emp(K);
for (std::size_t k = 0; k < K; ++k) {
  emp[k] = emplois[k];
}

for (std::size_t i = ishuf.begin(); i != ishuf.end(); ++i) {
  
  // Récupération des infos depuis la matrice sparse des rangs de passage par chacun des sites.
  // Pour des raisons d'efficience liées au codage en colonne des sparses matrices, une colonne code ce que voit un actif.
  std::vector<int> arrangement(K);
  std::vector<double> od (K);
  std::size_t k_valid = 0;
  for (int k = ranking.begin_col(i); k != ranking.end_col(i); ++k) {
    arrangement[*k - 1L] = k.row();
    od[*k - 1L] = *(odds.begin_col(i) + k);
    k_valid++;
  }
  arrangement.resize(k_valid);
  od.resize(k_valid);
  
  std::vector<double> dispo(k_valid), repartition(k_valid);
  for (std::size_t j = 0; j < k_valid; ++j) {
    dispo[j] = emp[arrangement[j]]; 
  }
  
  Rcout << "arrangement\n" ;
  for (auto k: arrangement) Rcout << *k << " ";
  Rcout << "\nod\n" ;
  for (auto k: od) Rcout << *k << " ";
  // Choix d'une limite basse pour la fuite.
  double fuite = std::max(1e-3, f[i]);
  double actifs_inzone = (1 - fuite) * actifs[i] / freq_actifs[i];
  
  repartition = repartir_actifs(dispo, od, fuite, actifs_inzone);
  }
  