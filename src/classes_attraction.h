#ifndef __CLASSES_ATTRACTION__
#define __CLASSES_ATTRACTION__

#include <vector>
#include <memory>

class fonction_attraction {
    public:
    virtual ~fonction_attraction() {};
    virtual double operator()(double x) = 0;
    
};

class constant : public fonction_attraction {
    public:
    virtual double operator()(double x) override { return 1.0; }
};

// 1 avant rayon (p0), plancher (p1) après.
class marche : public fonction_attraction {
    public:
    std::vector<double> parametres;
    marche(std::vector<double> param): parametres(param) {};
    virtual double operator()(double x) override { return (x < parametres[0] ? 1 : parametres[1]);};
};

// 1 avant rayon (p0), plancher (p1) après (pente sur une longueur de 1 entre les deux).
class marche_liss : public fonction_attraction {
public:
    std::vector<double> parametres;
    marche_liss(std::vector<double> param): parametres(param) {};
    virtual double operator()(double x) override { 
        return (x < parametres[0] ? 1 : (x > parametres[0] + 1 ? parametres[1] : 1 + (x - parametres[0]) * (parametres[1] - 1)));
        };
};

// ATTENTION au cas x = 0 qui n'est pas prévu.
// décroissance puissance delta (p0) tendant vers plancher (p1).
class decay : public fonction_attraction {
public:
    std::vector<double> parametres;
    decay(std::vector<double> param): parametres(param) {};
    virtual double operator()(double x) override { return (exp(-parametres[0] * log(x)) + parametres[1]);};
};


// logistique qui bascule autour du rayon (p0) avec une amplitude (p1) pour tendre vers un plancher (p2).
class logistique : public fonction_attraction {
public:
    std::vector<double> parametres;
    logistique(std::vector<double> param): parametres(param) {};
    virtual double operator()(double x) override { 
        double ex = exp( (parametres[0] - x) / parametres[1] );
        return parametres[2] + ex / (ex + 1);
        };
};

// design pattern i-e factory method
class mode_attraction {
public:
    virtual ~mode_attraction(){};
    std::shared_ptr<fonction_attraction> create(const std::vector<double>& param,
                                                const std::string& type) {
        if (type == "constant") {
            return std::make_shared<constant>();
        } else if (type == "marche") {
            return std::make_shared<marche>(param);
        } else if (type == "marche_liss") {
            return std::make_shared<marche_liss>(param);
        } else if (type == "decay") {
            return std::make_shared<decay>(param);
        } else if (type == "logistique") {
            return std::make_shared<logistique>(param);
        } else {
            return nullptr;
        }
    }
};


#endif // __CLASSES_ATTRACTION__