#include <Rcpp.h>

using namespace Rcpp;
//' @importClassesFrom Matrix dsRMatrix
// [[Rcpp::export]]
 S4 make_dsRMatrix() {
   S4 s("dsRMatrix");
   return s;
 }

// [[Rcpp::depends(Matrix)]]
class RankedRowSM {
public:
    IntegerVector jr, p, dim;
    NumericVector xr;

    // Constructors.
    RankedRowSM(IntegerVector jr1, IntegerVector p1, IntegerVector dim1, NumericVector xr1) : jr(jr1), p(p1), dim(dim1), xr(xr1){};
    RankedRowSM(const S4& mat) {
        if (!mat.hasSlot("xr") || !mat.hasSlot("p") || !mat.hasSlot("jr") || !mat.hasSlot("dim"))
            throw std::invalid_argument("Cannot construct RankedRowSM from this S4 object");
        jr = mat.slot("jr");
        p = mat.slot("p");
        dim = mat.slot("dim");
        xr = mat.slot("xr");
    }
    RankedRowSM(){};

    unsigned int nrow() { return dim[0]; }
    unsigned int ncol() { return dim[1]; }
    unsigned int nvalues() { return xr.size(); }
    unsigned int nzeros() { return dim[0] * dim[1] - xr.size(); };

    // Methods.
    //' Fonction de passage d'une RsparseMatrix (class dsRMatrix du pkg Matrix) à une Ranked Row Sparse Matrix.
    //' Il s'agit d'un reclassement des @x et @j pour réordonner les j des chaques lignes selon le rang des x
    //correspondants.

    // RankedRowSM to_RankedRowSM(IntegerVector j, IntegerVector p, IntegerVector dim, NumericVector x) {
    //     IntegerVector jranked(j.size());
    //     IntegerVector pp(p);
    //     IntegerVector dd(dim);
    //     NumericVector xranked(x.size());
    //     for (auto from = 0; from < dim[0]; ++from) {
    //         std::multimap<double, int> arrangement;
    //         unsigned int debut = p(from);
    //         unsigned int fin = p(from + 1L);
    //         for (auto k = debut; k < fin; ++k) {
    //             arrangement.insert(std::make_pair(x(k), j(k)));
    //         }
    //         int index = 0;
    //         for (auto it = arrangement.begin(); it != arrangement.end(); ++it) {
    //             xranked[debut + index] = (it->first);
    //             jranked[debut + index] = (it->second);
    //             ++index;
    //         }
    //     }
    //     return RankedRowSM(jr = jranked, p = pp, dim = dd, xr = xranked);
    // }
    // 
    S4 to_dsRMatrix(RankedRowSM &mat) {
        std::vector<int> jr = as< std::vector<int> >(mat.jr);
        std::vector<double> xr = as< std::vector<double> >(mat.xr);
        std::vector<int> j(jr);
        std::vector<double> x(xr.size());

        for (auto from = 0; from < mat.nrow(); ++from) {
            int debut = mat.p[from],
                fin = mat.p[from + 1L];
            std::sort(j.begin() + debut, j.begin() + fin);

            for (auto k = debut; k < fin; ++k) {
                auto pos = std::lower_bound(j.begin() + debut, j.begin() + fin, jr[k]);
                auto ind = std::distance(j.begin(), pos);
                x[ind] = xr[k] ;
            }
        }
        
    //     S4 dsR = make_dsRMatrix() {
    //     dsR.slot("p") = mat.p;
    //     dsR.slot("j") = j;
    //     dsR.slot("Dim") = mat.dim;
    //     dsR.slot("x") = x;
    //   
    //     return dsR;
    // };
};


// Expose the classes
RCPP_MODULE(RankedRowsSparseMatrix) {
  using namespace Rcpp;

  class_<RankedRowSM>("RankedRowSM")
      .constructor<IntegerVector, IntegerVector, IntegerVector, NumericVector>("constructor")
      //.method("to_RankedRowSM", &RankedRowSM::to_RankedRowSM)
      //.method("to_dsRMatrix", &RankedRowSM::to_dsRMatrix)
      .field_readonly("jr", &RankedRowSM::jr)
      .field_readonly("p", &RankedRowSM::p)
      .field_readonly("dim", &RankedRowSM::dim)
      .field_readonly("xr", &RankedRowSM::xr);
}