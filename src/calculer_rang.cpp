#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <random>

using namespace std;

std::vector<int> calculer_rang(const std::vector<double> &v)
{
  std::vector<std::size_t> w(v.size());
  std::iota(begin(w), end(w), 0);
  std::sort(begin(w), end(w), 
            [&v](std::size_t i, std::size_t j) { return v[i] < v[j]; });
  
  std::vector<int> r(w.size());

  std::random_device rd;
  std::mt19937 g(rd());

  for (std::size_t n, i = 0; i < w.size(); i += n)
  {
    n = 1;
    while (i + n < w.size() && v[w[i]] == v[w[i+n]]) ++n;
    
    std::vector<int> egaux(n);
    std::iota(begin(egaux), end(egaux), 0L);
    std::shuffle(egaux.begin(), egaux.end(),g);
    for (std::size_t k = 0; k < n; ++k)
    {
      //r[w[i+k]] = i + (n + 1) / 2.0; // average rank of n tied values
      // r[w[i+k]] = i + 1;     // min 
      // r[w[i+k]] = i + n;     // max
      r[ w[i+k] ] = i + egaux[k] + 1; // random order
    }
  }
  return r;
}