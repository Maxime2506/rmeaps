#ifndef __CLASSES_ATTRACTION__
#define __CLASSES_ATTRACTION__

#include <memory>
#include <vector>

// Définition des fonctions d'attraction.
// Ces fonctions étant renormalisées au sein de meaps, on les définit ici avec les conventions suivantes :
// La limite à l'infini de chaque fonction doit être égale à 1.
// On appelle amplitude la différence entre la valeur en O et celle à l'infini, soit f(0)-1.
class fonction_attraction {
 public:
  virtual ~fonction_attraction(){};
  virtual double operator()(double x) = 0;
};

class constant : public fonction_attraction {
 public:
  virtual double operator()(double x) override { return 1.0; }
};

// Fonction MARCHE.
// marche p1+1 avant rayon (p0), 1 après.
// amplitude = p1.
class marche : public fonction_attraction {
 public:
  std::vector<double> parametres;
  marche(std::vector<double> param) : parametres(param){};
  virtual double operator()(double x) override { return (x < parametres[0] ? 1 + parametres[1] : 1); };
};

// Fonction RAMPE.
// p1 + 1 en 0 puis décroissant linéairement jusqu'à x = p0, vaut 1 à partir de ce point.
// amplitude = p1
class rampe : public fonction_attraction {
 public:
  std::vector<double> parametres;
  rampe(std::vector<double> param) : parametres(param){};
  virtual double operator()(double x) override {
    return (x < parametres[0]
                ? 1 + parametres[1] - parametres[1] / parametres[0] * x
                : 1);
  };
};

// Fonction GRAV_EXP
// f(x) = 1 + p1.exp(-p0.x)
// vitesse de décroissance = p0.
// amplitude = p1.
class grav_exp : public fonction_attraction {
 public:
  std::vector<double> parametres;
  grav_exp(std::vector<double> param) : parametres(param){};
  virtual double operator()(double x) override { return (1 + parametres[1] * exp(-parametres[0] * x)); };
};

// Fonction GRAV_PUISS
// f(x) = 1 + p2.(x + p1)^(-p0)
// amplitude = p2 / p1^p0.
class grav_puiss : public fonction_attraction {
 public:
  std::vector<double> parametres;
  grav_puiss(std::vector<double> param) : parametres(param){};
  virtual double operator()(double x) override {
    return (1 + parametres[2] * (-parametres[0] * log(x + parametres[1])));
  };
};

// Fonction LOGISTIQUE
// logistique qui bascule autour du rayon (p0) avec un smooth (p1) pour tendre vers 1.
// en posant g(x) = exp((p0 - x)/p1), on a f(x) = 1 + p2 * g(x)/(1+g(x))
class logistique : public fonction_attraction {
 public:
  std::vector<double> parametres;
  logistique(std::vector<double> param) : parametres(param){};
  virtual double operator()(double x) override {
    double ex = exp((parametres[0] - x) / parametres[1]);
    return 1 + parametres[2] * ex / (ex + 1);
  };
};

// design pattern i-e factory method
class mode_attraction {
 public:
  virtual ~mode_attraction(){};
  std::shared_ptr<fonction_attraction> create(const std::vector<double>& param, const std::string& type) {
    if (type == "constant") {
      return std::make_shared<constant>();
    } else if (type == "marche") {
      return std::make_shared<marche>(param);
    } else if (type == "rampe") {
      return std::make_shared<rampe>(param);
    } else if (type == "grav_exp") {
      return std::make_shared<grav_exp>(param);
    } else if (type == "grav_puiss") {
      return std::make_shared<grav_puiss>(param);
    } else if (type == "logistique") {
      return std::make_shared<logistique>(param);
    } else {
      return nullptr;
    }
  }
};

#endif  // __CLASSES_ATTRACTION__