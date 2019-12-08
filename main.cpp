#include <iostream>
#include <fstream>
#include <cmath>
#include "parameters.h"


void count_ham_R(Parameters &parameters, int kappa, int omega, int tau){

        for (int k = 1; k < parameters.N - 1; k++){

            double xk = k * parameters.delta_x;
            parameters.ham_R[k] = ((-0.5) * (parameters.psi_R[k+1] + parameters.psi_R[k-1] - 2.0 * parameters.psi_R[k]) /
                                   (parameters.delta_x * parameters.delta_x)) + kappa * (xk - 0.5) * parameters.psi_R[k]*sin(omega*tau);
        }
}

void count_ham_l(Parameters &parameters, int kappa, int omega, int tau){

    for (int k = 1; k < parameters.N - 1; k++){

        double xk = k * parameters.delta_x;
        parameters.ham_l[k] = ((-0.5) * (parameters.psi_l[k+1] + parameters.psi_l[k-1] - 2.0 * parameters.psi_l[k]) /
                               (parameters.delta_x * parameters.delta_x)) + kappa * (xk - 0.5) * parameters.psi_l[k]*sin(omega*tau);
    }
}

void count_psi(Parameters &parameters, int n, const double Pi){

    for(int k = 0; k <parameters.N; k++){

        double xk = k * parameters.delta_x;
        parameters.psi_R[k] = sqrt(2.0) * sin(n * Pi * xk);
        parameters.psi_l[k] = 0.0;
    }

}

void savefile(Parameters &parameters, char* oFile,double norm, double x_srd, double epsilon,double t,bool b,string wybor,double gestosc){

    if(wybor == "r"){
        ofstream ofile;
        if (b == 1) {
            ofile.open(oFile,std::ios_base::app);
        }
        else {
            ofile.open(oFile);
        }
        if(ofile.is_open()){
            if(b==0){
                ofile << "TIME" << "\t" << "Norma" << '\t'<<"\t" << "x_średnie" << '\t' << "epsilon" << endl;
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
        if (b == 1) {
            ofile.open(oFile,std::ios_base::app);
        }
        else {
            ofile.open(oFile);
        }
        if(ofile.is_open()){
                for (int i=0;i<parameters.N;i++)ofile<<gestosc<<endl;
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

    int n = 1;
    int tau = 1;
    int kappa = 0;
    int omega = 0;
    int n_sym = 18000;
    double delta_tau = 0.0001;
    const double Pi = 3.14;
    double krok_t = 0.0;

    //Warunki brzegowe
    parameters.ham_R[0] = parameters.ham_R[parameters.N] = 0.0;
    parameters.ham_l[0] = parameters.ham_l[parameters.N] = 0.0;

    //Obliczanie poczatkowego psi_R i psi_l
    count_psi(parameters,n,Pi);

    //Obliczanie poczatkowego ham_R i ham_l
    count_ham_R(parameters,kappa,omega,tau);
    count_ham_l(parameters,kappa,omega,tau);

    savefile(parameters,arg[1],0,0,0,krok_t,0,"r",0);
    savefile(parameters,arg[2],0,0,0,krok_t,0,"g",0);

    for (int s = 0; s < n_sym; s++){

        double norm = 0.0;
        double x_srd = 0.0;
        double epsilon = 0.0;
        double gestosc = 0.0;

        //Algorytm calkowania

        for (int k = 0 ; k < parameters.N ; k++) {
            parameters.psi_R[k] = parameters.psi_R[k] + parameters.ham_l[k] * delta_tau / 2.0;
            //cout<<"psi_R z symulacji: "<<parameters.psi_R[k]<<endl;
        }

        count_ham_R(parameters,kappa,omega,tau);

        for (int k = 0 ; k < parameters.N ; k++) {
            parameters.psi_l[k] = parameters.psi_l[k] - parameters.ham_R[k] * delta_tau;
            //cout<<"psi_l z symulacji: "<<parameters.psi_l[k]<<endl;
        }

        count_ham_l(parameters,kappa,omega,tau);

        for (int k = 0 ; k < parameters.N ; k++) {
            parameters.psi_R[k] = parameters.psi_R[k] + parameters.ham_l[k] * delta_tau / 2.0;
            //cout<<"psi_R z symulacji: "<<parameters.psi_R[k]<<endl;
        }

        for (int k = 0; k < parameters.N ; k++){

            double xk = k * parameters.delta_x;
            norm += parameters.delta_x * (parameters.psi_R[k] * parameters.psi_R[k] + parameters.psi_l[k] * parameters.psi_l[k]);
            x_srd += parameters.delta_x * (xk * (parameters.psi_R[k] * parameters.psi_R[k] + parameters.psi_l[k] * parameters.psi_l[k]));
            epsilon += parameters.delta_x * (parameters.psi_R[k] * parameters.ham_R[k] + parameters.psi_l[k] * parameters.ham_l[k]);
            gestosc +=  parameters.psi_R[k] * parameters.psi_R[k] + parameters.psi_l[k] * parameters.psi_l[k];
        }

        if(s%(int)100 == 0){
        savefile(parameters,arg[1],norm,x_srd,epsilon,krok_t,1,"r",gestosc);
        }
        if(s%(int)100 == 0){
            savefile(parameters,arg[2],norm,x_srd,epsilon,krok_t,1,"g",gestosc);
        }

        krok_t += delta_tau;
    }

    return 0;
}