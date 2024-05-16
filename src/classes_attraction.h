#ifndef __CLASSES_ATTRACTION__
#define __CLASSES_ATTRACTION__

class fonction_attraction {
    public:
    fonction_attraction() {};
    virtual double operator()(double x) = 0;
    
};

class constant : public fonction_attraction {
    public:
    constant(){};
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
    std::unique_ptr<fonction_attraction> create(const std::vector<double> param,
                                                const std::string& type) {
        if (type == "constant") {
            return std::make_unique<constant>();
        } else if (type == "marche") {
            return std::make_unique<marche>(param);
        } else {
            return nullptr;
        }
    }
};


#endif // __CLASSES_ATTRACTION__