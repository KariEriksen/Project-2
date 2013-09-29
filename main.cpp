#include <iostream>
#include <math.h>
#include <cstdlib>

using namespace std;


double max_offdiag(double** A, int & k, int & l, int n_step);
void rotation(double **A, int k, int l, int n_step);

int main(int argc, char* argv[]) {



    // Defining variables.
    int n_step;
    //cout <<"Give an integer: ";
    //cin >> n;

    n_step = atoi(argv[1]);
    double h;
    int i, j, k, l;
    int n = n_step;
    double rho_i;
    double rho_min;
    double rho_max;
    double eps, error;

    double** A;
    A = new double* [n];
    double* V_i = new double[n];
    double* d_i = new double[n];
    double* e_i = new double[n];

    for (i=0; i < n; i++) {
        A[i] = new double[n];
    }

    rho_min = 0.0;
    rho_max = n_step*10; // The biggest rho value
    h = (rho_max - rho_min)/(n_step+1);

    for (i=0; i < n_step; i++) {

        rho_i = rho_min + (i+1)*h;
        V_i[i] = rho_i*rho_i;// +1/rho_i[i];
        d_i[i] = 2.0/(h*h) + V_i[i];
        e_i[i] = -1.0/(h*h);


    }

    for (i=0; i < n; i++) {
        for (j=0; j < n; j++) {

            if (j == i) {
                A[i][j] = d_i[i];
            }

            else if (j == i+1){
                A[i][j] = e_i[i];
            }

            else if (j == i-1){
                A[i][j] = e_i[i];

            }

            else {
                A[i][j] = 0.0;
            }
        }

    }

    eps = 10e-10;
    error = max_offdiag(A, k, l, n_step);
    cout << error << endl;
    int counter = 0;

    while ((fabs(error) > eps) && (counter < 10)) {

        rotation(A, k, l, n_step);
        error = max_offdiag(A, k, l, n_step);
        counter ++;
        cout << "hei" << endl;
        cout << error << endl;

    }

    for (i=0; i < n; i++) {
        for (j=0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;

    }
}

// Finding the maximum non-diagonal matrix element in A.
double max_offdiag (double **A, int & k, int & l, int n_step) {

    double max_value;
    max_value = 0.0; // The highest value in A, changeing trough out the loop

    for (int i=0; i < n_step-1; i++) {
        for (int j=i+1; j < n_step-1; j++) {
            if (fabs(A[i][j]) > max_value) {
                max_value = A[i][j];
                l = i;
                k = j;
            }

        }
    }
    return max_value;


}


void rotation (double **A, int k, int l, int n_step) {

    double tau, t_plus, t_minus, s, c;
    double A_ik, A_il, A_kk, A_ll;

    tau = (A[l][l] - A[k][k])/(2*A[k][l]);
    t_plus = -1/(-tau + sqrt(1 + tau*tau));
    t_minus = 1/(tau + sqrt(1 + tau*tau));

    if (t_plus < t_minus) {

        c = 1/(sqrt(1 + t_plus*t_plus));
        s = t_plus*c;
    }

    else {

        c = 1/(sqrt(1 + t_minus*t_minus));
        s = t_minus*c;
    }


    for (int i=0; i < n_step; i++) {

        if ((i != k) && (i != l)) {

            A_ik = A[i][k];
            A_il = A[i][l];
            A[i][k] = A_ik*c - A_il*s;
            A[k][i] = A[i][k];
            A[i][l] = A_il*c + A_ik*s;
            A[l][i] = A[i][l];

        }

    }

    A_kk = A[k][k];
    A_ll = A[l][l];
    A[k][k] = A_kk*c*c - 2*A[k][l]*c*s + A_ll*s;
    A[l][l] = A_ll*c*c + 2*A[k][l]*c*s + A_kk*s*s;
    A[k][l] = 0;
    A[l][k] = 0;

}




