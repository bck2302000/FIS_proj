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

vector<vector<double>> restriction(vector<vector<double>> res){
    int u_size = res.size();
    int N_c = 1 + u_size/2 ;
    vector<vector<double>> new_res(N_c, vector<double>(N_c, 0.0));
    for (int i = 1; i < N_c-1; i++){
        int ii = 2 * i;
        for (int j = 1; j < N_c-1; j++){
            int jj = 2 * j;
            new_res[i][j] = 0.0625 * (res[ii-1][jj-1] + 2*res[ii][jj-1] + res[ii+1][jj-1] + 2*res[ii-1][jj] + 4*res[ii][jj] + 2*res[ii+1][jj]\
                                      + res[ii-1][jj+1] + 2*res[ii][jj+1] + res[ii+1][jj+1]);
        }
    }
    return new_res;
}

vector<vector<double>> prolongation(vector<vector<double>> e){
    int e_size = e.size();
    int ee = 2 * e_size + 1;
    vector<vector<double>> e_p(ee, vector<double>(ee, 0.0));

    for (int i = 1; i < e_size-1; i++){
        int ii = 2*i;
        for (int j = 1; j < e_size-1; j++){
            int jj = 2*j;
            double temp = 0.25 * e[i][j];

            e_p[ii-1][jj-1] += temp; 
            e_p[ii][jj-1] += 2 * temp;
            e_p[ii+1][jj+1] += temp; 
            e_p[ii+1][jj] += 2 * temp;
            e_p[ii+1][jj+1] += temp;
            e_p[ii][jj+1] += 2 * temp;
            e_p[ii-1][jj+1] += temp;
            e_p[ii][jj-1] += 2 * temp;
            e_p[ii][jj] += 4 * temp;
        }
    }
    return e_p;
}

int main() {
    vector<double> err;
    for (int times = 0; times < 2; times++){
        int n = (times == 0) ? 4 : 7;
        int N = pow(2, n) + 1;
        double h = (double)1 / (N-1);
        double x,y;
        vector<vector<double>> u(N, vector<double>(N, 0));
        vector<vector<double>> u_sol(N/2, vector<double>(N/2, 0));
        vector<vector<double>> u_c;
        // Initialize f and u_sol
        for (int i = 1; i < N-1; i++){
            x = i * h;
            for (int j = 1; j < N-1; j++){
                y = j * h;
                u[i][j] = sin(2 * PI * x) * sin(2 * PI * y);
            }
        }
        for (int i = 1; i < N / 2; i++){
            x = i * 2 * h;
            for (int j = 1; j < N / 2; j++){
                y = j * 2 * h;
                u_sol[i][j] = sin(2 * PI * x) * sin(2 * PI * y);
            }
        }
        
        u_c = restriction(u);
        vector<vector<double>> diff_mtx = mtx_substract(u_c, u_sol);
        double diff = norm_inf(diff_mtx);
        err.push_back(diff);
    }

    for (int times = 0; times < 2; times++){
        int n = (times == 0) ? 4 : 7;
        int N = pow(2, n) + 1;
        double h = (double)1 / (N-1);
        vector<vector<double>> u(N/2, vector<double>(N/2, 0));
        vector<vector<double>> u_sol(N, vector<double>(N, 0));
        vector<vector<double>> u_f;
        double x, y;
        // Initialize f and u_sol
        for (int i = 1; i < N-1; i++){
            x = i * h;
            for (int j = 1; j < N-1; j++){
                y = j * h;
                u_sol[i][j] = sin(2 * PI * x) * sin(2 * PI * y);
            }
        }
        for (int i = 1; i < N / 2; i++){
            x = i * 2 * h;
            for (int j = 1; j < N / 2; j++){
                y = j * 2 * h;
                u[i][j] = sin(2 * PI * x) * sin(2 * PI * y);
            }
        }
        
        u_f = prolongation(u);
        vector<vector<double>> diff_mtx = mtx_substract(u_f, u_sol);
        double diff = norm_inf(diff_mtx);
        err.push_back(diff);
    }
    std::string fileNameBP = "error_restriction_and_prolongation";
    std::string fileExtension = ".dat";
    std::string fileNameFinal = fileNameBP + fileExtension;
    ofstream fout;
    fout.open(fileNameFinal);
    fout << "Error of restricion" << endl;
    fout << "n = 4: " << err[0] << endl;
    fout << "n = 7: " << err[1] << endl;
    fout << "Error of prolongation" << endl;
    fout << "n = 4: " << err[2] << endl;
    fout << "n = 7: " << err[3] << endl;
    fout.close();
}