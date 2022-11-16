#include <vector>
#include <math.h>
using namespace std;

vector<vector<double>> multigrid(vector<vector<double>> u, vector<vector<double>> f, int gamma, int nu_1, int nu_2){
    vector<vector<double>> u_bar = smoother_GS(u, f, nu_1); 
    vector<vector<double>> res = residual(u_bar, f);
    vector<vector<double>> res_c = restriction(res);
    vector<vector<double>> e_c(res_c.size(), vector<double>(res_c.size(), 0.0));
    for (int i = 1; i < res_c.size()-1; i++){
            for (int j = 1; j < res_c.size()-1; j++){
                 res_c[i][j] = -res_c[i][j];
            }
        }
    if (res_c.size() == 3){
        double h = 0.5;
        e_c[1][1] = 0.25 * pow(h,2) * res_c[1][1]; 
    }
    else{
        for (int i = 0; i < gamma; i++){
            e_c = multigrid(e_c, res_c, gamma, nu_1, nu_2);
        }
    }
    vector<vector<double>> e = prolongation(e_c);
    vector<vector<double>> u_bar_minus_e(u_bar.size(), vector<double>(u_bar.size(), 0.0));

    for (int i = 1; i < u_bar.size()-1; i++){
        for (int j = 1; j < u_bar.size()-1; j++) {
            u_bar_minus_e[i][j] = u_bar[i][j] - e[i][j];
        }
    }
    vector<vector<double>> u_next = smoother_GS(u_bar_minus_e, f, nu_2);

    return u_next;
}