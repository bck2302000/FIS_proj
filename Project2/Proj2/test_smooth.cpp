#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <string>
#include <fstream>
#include <iomanip>
#include <ctime>

using namespace std;
# define PI 3.14159265358979323846

// infinity norm, i.e. largest value in the matrix
double norm_inf(vector<vector<double>> r){
    double r0 = 0;
    for (int i = 0; i < r.size(); i++){
        double temp = *std::max_element(r[i].begin(), r[i].end());
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


auto smoother_GS(vector<vector<double>>u, vector<vector<double>>f, vector<vector<double>> u_sol){
    struct retVals{
        vector<vector<double>> v1;
        vector<double> v2;
        vector<double> v3;
        int in;
    };
    int u_size = u.size();
    double h = (double)1 / (u_size - 1);
    double hh = pow(h, 2.0);
    vector<vector<double>> u_prev(u_size, vector<double>(u_size, 0.0));
    vector<vector<double>> diff_mtx(u_size, vector<double>(u_size, 0.0));
    vector<vector<double>> diff_con_mtx;
    vector<double> err;
    vector<double> con_max_err;
    u_prev = u;
    int k = 0;
    while (1){
        for (int i = 1; i < u_size - 1; i++){
            for (int j = 1; j < u_size - 1; j++){
                u[i][j] = 0.25 * (hh * f[i][j] + u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
            }
        }
        diff_mtx = mtx_substract(u, u_prev);
        diff_con_mtx = mtx_substract(u, u_sol);
        double diff = norm_inf(diff_mtx);
        double diff_con = norm_inf(diff_con_mtx);
        err.push_back(diff);
        con_max_err.push_back(diff_con);
        k++;
        if (diff < 1e-10){
            break;
        }
        u_prev = u;
    }
    return retVals{u, err, con_max_err, k};    
}

int main() {
    time_t start, end;
    start = clock();
    int N = 10;
    double h = (double)1 / (N-1);
    double x,y;
    //int ite_smoo;
    //vector<double> err_smoo;
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

    vector<double> rel_err;
    vector<double> abs_err;
    //rel_err.push_back(r0);
    //vector<vector<double>> u_new;

    auto [u_new, err_smoo, con_err, ite_smoo] = smoother_GS(u, f, u_sol);

    end = clock();
    std::string fileNameBP = "error_smooth";
    std::string fileExtension = ".dat";
    std::string fileNameFinal = fileNameBP + std::to_string(N) + fileExtension;

    cout << "\nWriting data to file... ";
    ofstream fout;
    fout.open(fileNameFinal);
    fout << "Error of smoother:" << endl;
    for (int j = 0; j < err_smoo.size(); j++) {
        fout << std::setprecision(3) << std::scientific << err_smoo[j] << endl;
    }
    fout << "Use " << ite_smoo << "iteration" << endl;
    
    double elipsed_time = (double)(end - start) / CLK_TCK;
    fout << "Elipsed time:" << elipsed_time;
    fout << "\n";
    fout.close();

    fileNameBP = "converged_max_err";
    fileNameFinal = fileNameBP + std::to_string(N) + fileExtension;
    fout.open(fileNameFinal);
    fout << "Converged maximum error: " << endl;
    for (int i = 0; i < con_err.size(); i++){
        fout << std::setprecision(3) << std::scientific << con_err[i] << endl;
    }
    fout.close();
    cout << "Done!\n";

    return 0;
}