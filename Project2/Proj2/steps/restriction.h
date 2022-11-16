#include <vector>
#include <math.h>
using namespace std;

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