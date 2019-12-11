#include <iostream>
#include <fstream>
#include <cmath>
#include "parameters.h"

using namespace std;

void count_ham_R(Parameters &parameters, int kappa, double omega, double tau){

        for (int k = 1; k < parameters.N; k++){
            double xk = k * parameters.delta_x;
            parameters.ham_R[k] = -0.5*(parameters.psi_R[k+1]+parameters.psi_R[k-1]-2*parameters.psi_R[k])/pow(parameters.delta_x,2)+kappa*(xk-0.5)*parameters.psi_R[k]*sin(omega*tau);
        }
}

void count_ham_l(Parameters &parameters, int kappa, double omega, double tau){

    for (int k = 1; k < parameters.N; k++){
        double xk = k * parameters.delta_x;
        parameters.ham_l[k] = -0.5*(parameters.psi_l[k+1]+parameters.psi_l[k-1]-2*parameters.psi_l[k])/pow(parameters.delta_x,2)+kappa*(xk-0.5)*parameters.psi_l[k]*sin(omega*tau);
    }
}

void count_psi(Parameters &parameters, int n, const double Pi){

    for(int k = 0; k <parameters.N + 1; k++){
        double xk = k * parameters.delta_x;
        parameters.psi_R[k] = sqrt(2) * sin((double)n * Pi * xk);
        parameters.psi_l[k] = 0.0;
    }

}

void savefile(Parameters &parameters, char* oFile,double norm, double x_srd, double epsilon,double t,string wybor,double *gestosc){

    if(wybor == "r"){
        ofstream ofile;
        if (t == 0.0) {
            ofile.open(oFile);
        }
        else {
            ofile.open(oFile,std::ios_base::app);
        }
        if(ofile.is_open()){
            if(t==0){
                ofile << "TIME" << "\t" << "Norma" << '\t'<<"\t" << "x_średnie" << '\t' << '\t' << "epsilon" << endl;
                ofile << t << "\t" << norm << '\t' <<"\t"<< x_srd << '\t' << "\t" << epsilon << endl;
            }
            else{
                ofile << t << "\t" << norm << '\t'<<"\t" << x_srd << '\t' << "\t" << epsilon << endl;
            }
            ofile.close();
        }
        else {
            cout << "Unable to save file!"<<endl;
            exit(1);
        }

    }
    if (wybor == "g"){
        ofstream ofile;
        if (t == 0.0) {
            ofile.open(oFile);
        }
        else {
            ofile.open(oFile,std::ios_base::app);
        }
        if(ofile.is_open()){
                for (int i=0;i<parameters.N + 1;i+=2)ofile<<gestosc[i]<<endl;
                ofile<<"\n"<<endl;
                ofile.close();
        }
        else {
            cout << "Unable to save file!"<<endl;
            exit(1);
        }
    }
}

int main(int argc, char *arg[]){

    Parameters parameters = Parameters();

    //Obliczanie poczatkowego psi_R i psi_l
    count_psi(parameters,parameters.n,parameters.Pi);

    //Warunki brzegowe
    parameters.ham_R[0] = parameters.ham_R[parameters.N] = 0.0;
    parameters.ham_l[0] = parameters.ham_l[parameters.N] = 0.0;

    //Obliczanie poczatkowego ham_R i ham_l
    count_ham_R(parameters,parameters.kappa,parameters.omega,parameters.tau);
    count_ham_l(parameters,parameters.kappa,parameters.omega,parameters.tau);

    parameters.ham_R[0] = parameters.ham_R[parameters.N] = 0.0;
    parameters.ham_l[0] = parameters.ham_l[parameters.N] = 0.0;

    for (int s = 0; s < parameters.S_0; s++){

        double norm = 0.0;
        double x_srd = 0.0;
        double epsilon = 0.0;
        double gestosc[parameters.N+1];

        //Algorytm calkowania
        for (int k = 0 ; k < parameters.N + 1; k++) {
            parameters.psi_R[k] = parameters.psi_R[k] + parameters.ham_l[k] * parameters.delta_tau / 2.0;
        }

        count_ham_R(parameters,parameters.kappa,parameters.omega,(parameters.tau + parameters.delta_tau / 2.0));
        parameters.ham_R[0] = parameters.ham_R[parameters.N] = 0.0;

        for (int k = 0 ; k < parameters.N + 1; k++) {
            parameters.psi_l[k] = parameters.psi_l[k] - parameters.ham_R[k] * parameters.delta_tau;
        }

        count_ham_l(parameters,parameters.kappa,parameters.omega,(parameters.tau + parameters.delta_tau));
        parameters.ham_l[0] = parameters.ham_l[parameters.N] = 0.0;

        for (int k = 0 ; k < parameters.N + 1; k++) {
            parameters.psi_R[k] = parameters.psi_R[k] + parameters.ham_l[k] * parameters.delta_tau / 2.0;
        }

        count_ham_R(parameters,parameters.kappa,parameters.omega,(parameters.tau + parameters.delta_tau / 2.0));
        parameters.ham_R[0] = parameters.ham_R[parameters.N] = 0.0;

        if(s%(int)parameters.S_out == 0){
            for (int k = 0; k < parameters.N + 1; k++){
                double xk = k * parameters.delta_x;
                norm += parameters.delta_x * (pow(parameters.psi_R[k],2) + pow(parameters.psi_l[k],2));
                x_srd += parameters.delta_x * xk * (pow(parameters.psi_R[k],2) + pow(parameters.psi_l[k],2));
                epsilon += parameters.delta_x * (parameters.psi_R[k] * parameters.ham_R[k] + parameters.psi_l[k] * parameters.ham_l[k]);
                gestosc[k] =  pow(parameters.psi_R[k],2) + pow(parameters.psi_l[k],2);
                }

            savefile(parameters,arg[1],norm,x_srd,epsilon,parameters.tau,"r",gestosc);
            savefile(parameters,arg[2],norm,x_srd,epsilon,parameters.tau,"g",gestosc);
        }
        parameters.tau += parameters.delta_tau;
    }

    return 0;
}
