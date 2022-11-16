#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <string>
#include <fstream>
#include <iomanip>
#include <ctime>

#include "steps/smoother.h"
#include "steps/residual.h"
#include "steps/restriction.h"
#include "steps/prolongation.h"
#include "multigrid.h"

#define PI 3.14159265358979323846
using namespace std;

// infinity norm, i.e. largest value in the matrix

double norm_inf(vector<vector<double>> r){
    double r0 = 0;
    for (int i = 0; i < r.size(); i++){
        double temp = *max_element(r[i].begin(), r[i].end());
        if (temp > r0) r0 = temp;
    }
    return r0;
}

vector<vector<double>> mtx_substract(vector<vector<double>> mtx1, vector<vector<double>> mtx2){
    int mtx_size = mtx1.size();
    vector<vector<double>> res(mtx_size, vector<double>(mtx_size, 0.0));
    for (int i = 1; i < mtx_size-1; i++){
        for (int j = 1; j < mtx_size-1; j++) {
            res[i][j] = abs(mtx1[i][j] - mtx2[i][j]);
        }
    }
    return res;
}

int main(){
    time_t start, end;
    start = clock();
    int gamma = 1;                                          // v-cycle if gamma = 1, w-cycle if gamma = 2
    int nu_1 = 1;                                           // # of pre-smoother iteration 
    int nu_2 = 1;                                           // # of post-smoother iteration
    int n = 7;
    int N = pow(2, n) + 1;
    int m = 30;                                             // # of multigrid iteration
    double h = (double)1 / (N-1);
    double x,y;
    vector<vector<double>> f(N, vector<double>(N, 0));
    vector<vector<double>> u(N, vector<double>(N, 0));
    vector<vector<double>> u_sol(N, vector<double>(N, 0));
    
    // Initialize f and u_sol
    for (int i = 1; i < N-1; i++){
        x = i * h;
        for (int j = 1; j < N-1; j++){
            y = j * h;
            f[i][j] = 8 * pow(PI, 2) * sin(2 * PI * x) * sin(2 * PI * y);
            u_sol[i][j] = sin(2 * PI * x) * sin(2 * PI * y);
        }
    }

    // r0, the residual of the first guess
    vector<vector<double>> r0_res = residual(u,f);
    double r0 = norm_inf(r0_res);
    
    vector<double> rel_err;
    vector<double> abs_err;
    rel_err.push_back(r0);
    vector<vector<double>> u_new; 

    for (int i = 0; i < m; i++){
        u_new = multigrid(u, f, gamma, nu_1, nu_2);
        vector<vector<double>> r_res = residual(u_new, f);
        double r = norm_inf(r_res);
        vector<vector<double>> u_minus_u_sol = mtx_substract(u_new, u_sol);                  //exact error
        double err_sol = norm_inf(u_minus_u_sol);
        rel_err.push_back(r/r0);
        abs_err.push_back(err_sol);
        if (r/r0 < 2.5e-13) break;
        cout << "Error of "<< i+1 << "th Iteration: " << rel_err[rel_err.size()-1] << endl;
        u = u_new; 
    }
    end = clock();
    string fileNameBP = "error_";
    string fileExtension = ".dat";
    string fileNameFinal = fileNameBP + to_string(n) + to_string(gamma)
                                + to_string(nu_1) + to_string(nu_2) + to_string(m) + fileExtension;

    cout << "\nWriting data to file... ";
    ofstream fout;
    fout.open(fileNameFinal);
    fout << "Relative Error:" << endl;
    for (int j = 1; j < rel_err.size(); j++) {
        fout << setprecision(3) << scientific << rel_err[j] << endl;
    }
    fout << "Absolute Error" << endl;
    for (int j = 0; j < abs_err.size(); j++){
        fout << setprecision(3) << scientific << abs_err[j] << endl;
    }
    double elapsed_time = (double)(end - start) / CLK_TCK;
    fout << "Elapsed time:" << elapsed_time;
    fout << "\n";
    fout.close();
    cout << "Done!\n";

    return 0;
}
