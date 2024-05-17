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
    //constant(){};
    virtual double operator()(double x) override { return x; }
};

class marche : public fonction_attraction {
    public:
    std::vector<double> parametres;
    marche(std::vector<double> param): parametres(param) {};
    virtual double operator()(double x) override { return (x < parametres[0] ? 1 : parametres[1]);};
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
        } else {
            return nullptr;
        }
    }
};


#endif // __CLASSES_ATTRACTION__