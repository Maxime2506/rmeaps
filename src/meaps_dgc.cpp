#ifdef _OPENMP
# include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <RcppSparse.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "repartir_alt.h"

using namespace Rcpp;
//' La fonction meaps sur plusieurs shufs
//' @param dist La matrice (creuses en colonnes) des distances dans lequel les colonnes j sont passées en revue pour chacune des lignes i.
//' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes).
//' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
//' @param modds La matrice des odds modifiant la chance d'absorption de chacun des sites j pour des résidents en i.
//' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude.
//' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs.
//' @param mode Choix du rôle de l'emploi disponible au cours du processus. Default : continu. Autre choix : discret, subjectif_c ou _d...
//' @param odds_subjectifs Attractivité des sites proches des actifs, pour le mode defini. default : null.
//' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto.
//' @param progress Ajoute une barre de progression. Default : true.
//' @param normalisation Calage des emplois disponibles sur le nombre d'actifs travaillant sur la zone. Default : false.
//' @param fuite_min Seuil minimal pour la fuite d'un actif. Doit être supérieur à 0. Défault = 1e-3.
//' 
//' @return renvoie une matrice avec les estimations du nombre de trajets de i vers j.
// [[Rcpp::export]]
RcppSparse::Matrix meaps_dgc(RcppSparse::Matrix &dist,
                    NumericVector emplois,
                    NumericVector actifs,
                    RcppSparse::Matrix modds,
                    NumericVector f,
                    IntegerMatrix shuf,
                    std::string mode = "continu",
                    Nullable<NumericVector> oddssubjectifs = R_NilValue,
                    int nthreads = 0,
                    bool progress = true,
                    bool normalisation = false,
                    double fuite_min = 1e-3,
                    double seuil_newton = 1e-6) {
  
  const int N = dist.Dim(1L),
            K = dist.Dim(0L),
            Nboot = shuf.nrow(),
            Ns = shuf.ncol();
  
  // Quelques vérifications préalables.
  if (emplois.size() != K) {
    stop("Le vecteur emplois et le nb de colonnes de la matrice rkdist ne correspondent pas.");
  }
  if (f.size() != N) {
    stop("Le vecteur fuite et le nb de lignes de la matrice rkdist ne correspondent pas.");
  }
  if (actifs.size() != N) {
    stop("Le vecteur actifs et le nb de lignes de la matrice rkdist ne correspondent pas.");
  }
  if (modds.ncol() != K) {
    stop("La matrice modds et la matrice rkdist n\'ont pas le même nombre de colonnes.");
  }
  if (modds.nrow() != N) {
    stop("La matrice modds et la matrice rkdist n\'ont pas le même nombre de lignes.");
  }
  
  if ((mode == "subjectif_c" || mode == "subjectif_d")&& oddssubjectifs.isNull()) {
    stop("En mode subjectif, le vecteur oddssubjectifs doit être défini");
  }
  
  std::vector<double> odsub;
  if (mode == "subjectif_c"|| mode == "subjectif_d") {
    odsub = as< std::vector<double> >(oddssubjectifs);
  }
 
#ifdef _OPENMP
  int ntr = nthreads;
  if (ntr == 0) {ntr = omp_get_max_threads();}
  if (ntr > omp_get_max_threads()) { ntr = omp_get_max_threads(); }
  if (progress==TRUE) REprintf("Nombre de threads = %i\n", ntr);
#endif
  
  // Les rangs dans shuf commencent à 1.
  shuf = shuf - 1L;
  
  // Conversion en vecteur de vecteurs C++.
  std::vector<int> shufcpp = as< std::vector<int> >(shuf);
  std::vector< std::vector<int> > ishuf(Nboot, std::vector<int>(Ns));
  for (int i = 0; i < Nboot; ++i) {
    for (int j = 0; j < Ns; ++j) {
      ishuf[i][j] = shufcpp[ i  + j * Nboot];
      if (ishuf[i][j] >= N) { 
        stop("Les rangs de shufs vont au-delà de la taille du vecteur actifs.");
      }
    }
  }
  
  // Choix d'une limite basse pour la fuite.
  f = ifelse(f > fuite_min, f, fuite_min);
  
  // Calage de l'emploi sur les actifs.
  if (normalisation) {
    emplois = emplois * sum(actifs * (1 - f)) / sum(emplois); 
  }
  
  // Conversion en vecteurs c++ 
  std::vector<double> emploisinitial = as< std::vector<double> >(emplois);
  std::vector<double> fcpp = as< std::vector<double> >(f);
  std::vector<double> actifscpp = as< std::vector<double> >(actifs);
  
  // Construction des vecteurs arrangés selon l'ordre subjectifs des actifs.
  std::vecteur<double> les_distances(K);
  std::vecteur<int> les_rangs(K);
  std::vector< std::vector<int> > arrangement(N, std::vector<int>(K));
  std::vector< std::vector<double> > odds(N, std::vector<double>(K));

  int debut, fin, pos, to;

  for (int from = 0; from < N; ++from) {
    debut = dist.p[from];
    fin = dist.p[from + 1L] - 1L;
    if (fin < debut) continue; //Au cas où un from n'accède à aucun to.
    les_distances = dist.x(debut, fin);
    les_rangs = calculer_rang(les_distances);

    pos = modds.p[from];
    std::map<int, double> les_odds(modds.p[from + 1L] - pos);
    for(int k = 0; k < les_odds.size(); ++k) {
        les_odds.insert(make_pair(modds.i[pos + k], modds.x[pos + k]))
    }

    for(int k = 0; k < les_rangs.size(); ++k) {
        to = les_rangs[k] - 1L;
        arrangement[from][to] = dist.i[debut + k];

        if (les_odds.find(dist.i[debut + k]) == les_odds.end()) {
            odds[from][to] = 1; //La valeur par défaut pour modds est un odds ratio fixé à 1.
        } else {
            odds[from][to] = les_odds[dist.i[debut + k]];
        }  
    }
    arrangement[from].resize(fin - debut);
    odds[from].resize(fin - debut);
  }
  // Le vecteur shuf peut être plus long que le nombre de lignes de rkdist s'il fait repasser plusieurs fois
  // la même ligne d'actifs. Dans ce cas, on compte la fréquence de passage de chaque ligne et l'on divise le
  // poids de la ligne par cette fréquence.
  std::vector<int> freq_actifs(N, 0L);
  for (auto i : ishuf[0]) {
    freq_actifs[i]++;
  }
  
  // Initialisation du résultat.
  RcppSparse::Matrix liaisons = dist.clone();
  for (int k =0; k < liaisons.x.size(); ++k) { liaisons.x(k) = 0; }

  // Lancement du bootstrap.
  Progress p(Nboot * Ns, progress);
  #ifdef _OPENMP
  #pragma omp declare reduction(vsum : std::vector<double> :              \
    std::transform(omp_out.begin(), omp_out.end(),                        \
                 omp_in.begin(), omp_out.begin(), std::plus<double>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
  #pragma omp parallel for num_threads(ntr)                                                          \
    shared(Nboot, N, ishuf, emploisinitial, odds, fcpp, actifscpp) \
    reduction (vsum : liaisons.x)
  #endif
    for (int iboot = 0; iboot < Nboot; ++iboot) {
      
      std::vector<int> theshuf = ishuf[iboot]; // deep copy pour un boot.
      
      // Initialisation de l'emploi disponible.
      std::vector<double> emp(emploisinitial); // deep copy.
      
      for (auto from: theshuf) {
        
        // Increment progress_bar.
        p.increment();
        
        std::vector<int> arr = arrangement[from];
        int k_valid = arr.size();
        
        std::vector<double> placeslibres (k_valid), attractivites(k_valid), repartition(k_valid);
        for (int k = 0; k < k_valid; ++k) {
          placeslibres[k] = emp[ arr[k] ]; 
        }
        
        if (mode == "subjectif_d") {
          for (int k = 0; k < k_valid; ++k) {
            attractivites[k] = emploisinitial[ arr[k] ] * odsub[k]; 
          }
        } else if (mode== "subjectif_c") {
          for (int k = 0; k < k_valid; ++k) {
            attractivites[k] = placeslibres[k] * odsub[k]; 
          }
        } else if (mode == "discret") {
          for (int k = 0; k < k_valid; ++k) {
            attractivites[k] = emploisinitial[ arr[k] ]; 
          }
        } else {
          for (int k = 0; k < k_valid; ++k) {
            attractivites[k] = placeslibres[k]; // Cas continu de l'emploi homogène. attractivité liée à la proportion d'emplois disp.
          }
        }
        
        // Nombre d'actifs partant par freq_actif paquets.
        double actifspartant = actifscpp[from] / freq_actifs[from];
        
        repartition = repartir_alt(placeslibres, attractivites, odds[from], fcpp[from], actifspartant, seuil_newton);
       
        // Inscription des résultats locaux dans la perspective globale.
        for(int k = 0; k < k_valid ; ++k) {
          emp[ arr[k] ] -= repartition[k];
          liaisons.x[ liaisons.p[from] + arr[k] ] += repartition[k];
        }
      }
    }
    
  // Division par le nombre de tirage Nboot.
  liaisons.x = liaisons.x / Nboot;
  
  return liaisons;
}



