#ifndef __CONSTANTS__
#define __CONSTANTS__
 
constexpr int LIMITE_LOOP = 200; // condition d'arrêt pour les boucles lors de la distribution des résidents vers des emplois.
constexpr double LIMITE_PRECISION_1 = 1e-3; // condition d'arrêt sur le pourcentage de résidents non classés restants.
constexpr double LIMITE_PRECISION_2 = 1e-4; // condition d'arrêt sur la vitesse de reclassement des résidents non occupés.
constexpr double LIMITE_PRECISION_3 = LIMITE_PRECISION_1 / 1e6; // seuil en deçà duquel une ligne d'actifs est supposée vides (tous actifs occupés). 
constexpr double PLANCHER_KL = 1e-6; // gestion de cases nulles dans le calcul de l'entropie relative (KL).

#endif // __CONSTANTS__