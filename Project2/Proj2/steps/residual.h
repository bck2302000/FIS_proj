#include <vector>
#include <math.h>
using namespace std;

vector<vector<double>> residual(vector<vector<double>> u, vector<vector<double>> f){
    int u_size = u.size();
    double h = (double)1 / (u_size - 1);
    double hh_inv = (double)1 / pow(h, 2);
    vector<vector<double>> res(u_size, vector<double> (u_size, 0.0));

    for (int i = 1; i < u_size-1; i++){
        for (int j = 1; j < u_size-1; j++){
            res[i][j] = f[i][j] + hh_inv*(u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] - 4*u[i][j]);
        }
    }
    return res; 
}
