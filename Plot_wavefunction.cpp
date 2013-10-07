#include <iostream>
#include <math.h>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;


int main(int argc, char* argv[]) {



    // Defining variables

    int n_step = atoi(argv[1]);
    double omega = atof(argv[2]);
    double rho_max = atof(argv[3]);
    double h;
    int i, j;
    int n = n_step;
    double rho_i;
    double rho_min;

    mat A = zeros(n,n);
    vec V_i = zeros(n);
    vec d_i = zeros(n);
    vec e_i = zeros(n);

    rho_min = 0.0;
    h = (rho_max - rho_min)/(n_step+1);

    for (i=0; i < n_step; i++) {

        rho_i = rho_min + (i+1)*h;
        V_i(i) = omega*omega*rho_i*rho_i + 1/rho_i;
        d_i(i) = 2.0/(h*h) + V_i(i);
        e_i(i) = -1.0/(h*h);


    }

    for (i=0; i < n; i++) {
        for (j=0; j < n; j++) {

            if (j == i) {
                A(i,j) = d_i(i);
            }

            else if (j == i+1){
                A(i,j) = e_i(i);
            }

            else if (j == i-1){
                A(i,j) = e_i(i);

            }

            else {
                A(i,j) = 0.0;
            }
        }

    }

    // Starting clock
    clock_t start, finish;
    start = clock();

    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, A);

    // Finishing clock
    finish = clock();
    cout << (finish - start)/(double) CLOCKS_PER_SEC << endl;

    //cout << eigval << endl;


    fstream outFile;
        outFile.open("data.dat", ios::out);

        //outFile << n << endl;
        //outFile << rho_max << endl;
        //outFile << omega << endl;
        for (i=0; i < n; i++) {
            outFile << eigvec(i,0) << ' ' << eigvec(i,1) << ' ' << eigvec(i,2) << endl;

        }
        outFile.close();
        return 0;

}
