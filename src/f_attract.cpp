#include <functional>
#include <cmath>

template <typename T>
using f_attract = T (*)(const std::vector<T>&, const T&);

template <typename T>
T constant(const std::vector<T>& param, const T& dist) { return 1;}

template <typename T>
T marche(const std::vector<T>& param, const T& dist) { 
  return (dist < param[0] ? 1 : param[1]);
  }

template <typename T>
T marche_liss(const std::vector<T>& param, const T& dist) {
  int avant = floor(param[0]);
  if ( dist > avant + 1L ) return 1;
  if ( dist <= avant ) return param[1];
  return param[1] + (1 - param[1]) *  (dist - avant);
}

template <typename T>
T decay(const std::vector<T>& param, const T& dist) {
  return exp(-param[1] * log(dist)) + param[0];
}

template <typename T>
T logistique(const std::vector<T>& param, const T& dist) {
  T ex = exp( (param[0] - dist) / param[1] );
  return param[2] + ex / (ex + 1);
}

