#include <vector>
#include <math.h>
using namespace std;

vector<vector<double>> smoother_GS(vector<vector<double>>u, vector<vector<double>>f, int nu_1){
    int u_size = u.size();
    double h = (double)1 / (u_size - 1);
    double hh = pow(h, 2.0);
    vector<vector<double>> u_next(u_size, vector<double>(u_size, 0.0));
    for (int k = 0; k < nu_1; k++){
        for (int i = 1; i < u_size - 1; i++){
            for (int j = 1; j < u_size - 1; j++){
                u[i][j] = 0.25 * (hh * f[i][j] + u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
            }
        }
    }
    return u;    
}