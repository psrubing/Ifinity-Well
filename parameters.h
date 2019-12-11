
#ifndef LAB2_PARAMETERS_H
#define LAB2_PARAMETERS_H

using namespace std;

struct Parameters{

    int N = 100;
    double delta_x = 1.0 /(double)N;
    int S_0 = 400000;
    int S_out = 1000;
    const double Pi = 3.14159265359;
    int n = 4;
    double tau = 0.0;
    double delta_tau = 0.0001;
    int kappa = 10;
    double omega = 3.0 * Pi * Pi / 2.0;


    //tablice
    double psi_R[101];
    double psi_l[101];
    double ham_R[101];
    double ham_l[101];


};

#endif
