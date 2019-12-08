
#ifndef LAB2_PARAMETERS_H
#define LAB2_PARAMETERS_H

using namespace std;

struct Parameters{

    int N = 100;
    double delta_x = 1.0 /(double)N;

    //tablice
    double psi_R[100+1];
    double psi_l[100+1];

    double ham_R[100+1];
    double ham_l[100+1];


};

void count_ham_R(Parameters&, int, int, int);
void count_ham_l(Parameters&, int, int, int);
void count_psi(Parameters&, int, const double);

#endif
