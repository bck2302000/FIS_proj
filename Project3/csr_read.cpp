#include "csr_read.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
using namespace csr;

CSR_READ::CSR_READ(char* filename){
    COO_read(filename);
    convert_CSR_CSC();
}

CSR_READ::~CSR_READ(){
    //finalise_read();
}

void CSR_READ::COO_read(char* filename){
    ifstream f(filename, ios::in);
    string dummy; 
    

    // Ignore the first line of the mtx file.
    getline(f, dummy);
    f >> m >> n >> nz;

    diag = new double[m];
    int iy, ix;
    double ival;

    for (int i = 0; i < nz; i++){
        f >> iy;
        f >> ix;
        f >> ival;
        if (iy == ix)
            diag[ix - 1] = ival;
        IJV.push_back(iy, ix, ival)
    }
    f.close();
}

void CSR_READ::convert_CSR_CSC() {
    CSR_READ::csc_mtx.col_idx = new int[m];
    CSR_READ::csc_mtx.row = new int[nz];
    CSR_READ::csr_mtx.row_idx = new int[m];
    CSR_READ::csr_mtx.col = new int[nz];

    int last_row = 0;
    int last_col = 0;

    for (int i = 0; i < nz; i++) {
        //CSC
        CSR_READ::csc_mtx.row[i] = coo_mtx.y[i];
        if (coo_mtx.x[i] != last_row){
            CSR_READ::csc_mtx.col_idx[last_row] = coo_mtx.x[i];
        }
        last_row = coo_mtx.x[i];

        //CSR
        CSR_READ::csr_mtx.col[i] = coo_mtx.x[i];
        if (coo_mtx.x[i] != last_col){
            CSR_READ::csr_mtx.row_idx[last_row] = coo_mtx.y[i];
        }
        last_col = coo_mtx.x[i];
    }
}

double* CSR_READ::mtx_vec_product(double *v)
{
    int v_size = sizeof(*v) / sizeof(v[0]);
    double *result = new double[v_size];


    //CSC product
    int cursor = 0;
    for (int i = 0; i < v_size; i++){
        int num_on_row = CSR_READ::csc_mtx.col_idx[i + 1] - CSR_READ::csc_mtx.col_idx[i];
        for (int j = 0; j < num_on_row; j++){
            result[i] += coo_mtx.val[cursor + j] * v[CSR_READ::csc_mtx.row[cursor + j] - 1];
        }
        cursor += num_on_row;
    }
    //CSR product
    cursor = 0;
    for (int i = 0; i < v_size; i++){
        int num_on_row = CSR_READ::csr_mtx.row_idx[i + 1] - CSR_READ::csr_mtx.row_idx[i];
        for (int j = 0; j < num_on_row; j++){
            result[i] += coo_mtx.val[cursor + j] * v[CSR_READ::csr_mtx.col[cursor + j] - 1];
        }
        cursor += num_on_row;
    }
    
    //Subtract the diagnol
    for (int i = 0; i < v_size; i++){
        result[i] -= CSR_READ::diag[i] * v[i];
    }

    return result;
}