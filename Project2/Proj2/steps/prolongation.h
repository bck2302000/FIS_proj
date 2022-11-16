#include <vector>
#include <math.h>
using namespace std;

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