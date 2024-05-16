#include "classes.cpp"
using namespace std;

main() {
    vector<int> j = { 1, 2, 4, 3, 2, 0};
    vector<int> p = { 0, 3, 6};
    vector<double> x = { 1, 3, 8, 12, 5, 2};
    vector<double> act = { 4, 5};
    vector<double> job = { 2, 1, 3, 2, 1};
    vector<double> f = { .2, .3};

    Urban urb(j, p , x, f, act, job);
    Residents res();
};