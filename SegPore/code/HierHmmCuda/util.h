#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "Parameters.h"
#include "Region.h"
using namespace std;

#define N_B_MAX_MERGE 5
#define N_B_IN_REGION 11
#define MISSING_BLOCK_MU -2000
#define MISSING_BLOCK_SIGMA 2


void readFloatFile(string fileName, vector<float>& floatArray){
    ifstream inFile(fileName);
    string lineStr;
    while (getline(inFile, lineStr)) {
      stringstream ss(lineStr);
      string str;
      while (getline(ss, str, ',')){
        float str_float;
        istringstream istr(str);
        istr >> str_float;
	    floatArray.push_back(str_float);
	 }
  }
}


void read2DFloatFile(string fileName, vector<vector<float>>& _mat){
    ifstream inFile(fileName);
    string line;
    while (getline(inFile, line))
    {
        stringstream ss(line);
        string str;
        vector<float> _row;
        while (getline(ss, str, ',')){
            float str_float;
            istringstream istr(str);
            istr >> str_float;
            _row.push_back(str_float);
        }
        _mat.push_back(_row);
    }
}


void readOnefrom2DFloatFile(string fileName, vector<float>& floatArray, int idx){
    ifstream inFile(fileName);
    string line;
    int i = 0;
    while (getline(inFile, line))
    {
        if(i==idx){
            stringstream ss(line);
            string str;
            while (getline(ss, str, ',')){
               float str_float;
               istringstream istr(str);
               istr >> str_float;
               floatArray.push_back(str_float);
            }
            return;
        }else{
            i++;
        }
    }
}


void readIntFile(string fileName, vector<int>& intArray){
    ifstream inBorderFile(fileName);
    string lineBorderStr;
    while (getline(inBorderFile, lineBorderStr)) {
      stringstream ss(lineBorderStr);
      string strBorder;
      while (getline(ss, strBorder, ',')){
        int str_int;
        istringstream isbordertr(strBorder);
        isbordertr >> str_int;
	    intArray.push_back(str_int);
	 }
  }
}


void read2DIntFile(string fileName, vector<vector<int>>& _mat){
    ifstream inFile(fileName);
    string line;
    while (getline(inFile, line))
    {
        stringstream ss(line);
        string str;
        vector<int> _row;
        while (getline(ss, str, ',')){
            int str_float;
            istringstream istr(str);
            istr >> str_float;
            _row.push_back(str_float);
        }
        _mat.push_back(_row);
    }
}


void readOnefrom2DIntFile(string fileName, vector<int>& intArray, int idx){

    ifstream inFile(fileName);
    string line;
    int i = 0;
    while (getline(inFile, line))
    {
        if(i == idx){
           stringstream ss(line);
           string str;
           while (getline(ss, str, ',')){
               int str_float;
               istringstream istr(str);
               istr >> str_float;
               intArray.push_back(str_float);
           }
           return;
        }else{
            i++;
        }

    }
}


void gene_b_sig_ind_arr(vector<int>& border_arr, vector<int>& b_st_sig_ind_arr, vector<int>& b_en_sig_ind_arr, int N_B){

    assert(N_B == b_st_sig_ind_arr.size());
    assert(N_B == b_en_sig_ind_arr.size());

    for(int i=0;i<N_B;i++){
        b_st_sig_ind_arr[i] = border_arr[i*2];
        b_en_sig_ind_arr[i] = border_arr[i*2+1];
    }
}


void gene_t_sig_ind_arr(vector<int>& border_arr, vector<int>& t_st_sig_ind_arr, vector<int>& t_en_sig_ind_arr, int N_T){

    assert(N_T == t_st_sig_ind_arr.size());
    assert(N_T == t_en_sig_ind_arr.size());

    for(int i=0;i<N_T;i++){
        t_st_sig_ind_arr[i] = border_arr[i*2+1];
        t_en_sig_ind_arr[i] = border_arr[i*2+2];
    }
}


void init_B_mu_sigma(vector<float>& signal_arr, vector<int>& b_st_sig_ind_arr, vector<int>& b_en_sig_ind_arr,
                     vector<vector<float>>& B_mu_mat, vector<vector<float>>& B_sigma_mat){

    int N_B = b_st_sig_ind_arr.size();
    int i, j, k, tmpst, tmpen, tmp_n;
    float tmp_sum, tmp_mu, tmp_sigma;
    float min_sigma = 0.1;

    for(i=0;i<N_B;i++){
        for(j=0;j<N_B_MAX_MERGE;j++){
            if(i + j >= N_B){
                continue;
            }
            tmpst = b_st_sig_ind_arr[i];
            tmpen = b_en_sig_ind_arr[i + j];
            tmp_sum = 0;
            tmp_n = tmpen - tmpst;

            // cal mu
		    for(k=tmpst; k<tmpen; k++)
                tmp_sum += signal_arr[k];
		    tmp_mu = tmp_sum/tmp_n;
            B_mu_mat[i][j] = tmp_mu;

            // cal sigma
            tmp_sum = 0;
            for(k=tmpst; k<tmpen; k++)
                tmp_sum += pow((signal_arr[k] - tmp_mu),2);
            tmp_sigma = sqrt(tmp_sum / tmp_n);
            if(tmp_sigma >= min_sigma){
                B_sigma_mat[i][j] = tmp_sigma;
            }else{
                B_sigma_mat[i][j] = min_sigma;
            }

        }
    }
}


float mean_normal_log_pdf(vector<float>& arr, float m, float s){
    static const float inv_sqrt_2pi = 0.3989422804014327;
    int n = arr.size();
    float res = 0;
    float a;
    for(int i=0;i<n;i++){
        a = (arr[i] - m) / s;
        res += std::log(inv_sqrt_2pi / s * std::exp(-0.5f * a * a));
    }
    return res/n;
}


float calLrLikelihood(vector<float>& y_arr){

    // https://www.amherst.edu/system/files/media/1287/SLR_Leastsquares.pdf
    // sum (yi-k*xi)  / n

    //float k_min = 5.0;
    float sigma_thred = 2.5;


    int n = y_arr.size();
    //if(n<3 || n>=30) return -10;
    if(n<2) return -10;

    vector<float> x_arr;
    for(int i=0;i<n;i++){
        x_arr.push_back(i+1);
    }

    float x_mean=0, y_mean=0;
    for(int i=0; i<n; ++i){
        x_mean += x_arr[i];
        y_mean += y_arr[i];
    }
    x_mean = x_mean / n;
    y_mean = y_mean / n;

    float SS_xy = 0;
    float SS_xx = 0;
    for(int i=0; i<n; ++i){
        SS_xx += (x_arr[i] - x_mean)*(x_arr[i] - x_mean);
        SS_xy += (x_arr[i] - x_mean)*(y_arr[i] - y_mean);
    }

    float k = SS_xy / SS_xx;
    float b = y_mean - k * x_mean;

//    if (abs(k) < abs(k_min)){
//        if(k > 0){
//           k = k_min;
//        }else{
//           k = -k_min;
//        }
//        float tmp = 0;
//        for(int i=0; i<n; ++i) tmp += (y_arr[i] - k *  x_arr[i]);
//        b = tmp / n;
//    }

    vector<float> res_arr;
    for(int i=0;i<n;i++){
        res_arr.push_back(y_arr[i] - k*x_arr[i] - b);
    }

    float tmp_sum = 0;
    int tmp_n = res_arr.size();
    for(int i=0; i<tmp_n; i++){
        tmp_sum += pow((res_arr[i]),2);
    }

	sigma_thred = sqrt(tmp_sum / tmp_n);
    return mean_normal_log_pdf(res_arr, 0, sigma_thred);
}


void cal_T_lr_loglik_cpu(vector<float>& signal_arr, vector<int>& t_st_sig_ind_arr, vector<int>& t_en_sig_ind_arr, vector<float>& T_lr_loglik_arr){

    int i, j;
    vector<float> tmp_arr;
    for(i=0;i<t_st_sig_ind_arr.size();i++){
        for(j=t_st_sig_ind_arr[i]; j<t_en_sig_ind_arr[i]; ++j){
            tmp_arr.push_back(signal_arr[j]);
        }
        T_lr_loglik_arr[i] = calLrLikelihood(tmp_arr);
        tmp_arr.clear();
	    tmp_arr.shrink_to_fit();
    }
}


void cal_prev_next_B_mu_sigma(vector<vector<float>>& B_mu_mat, vector<vector<float>>& B_sigma_mat,
                              vector<vector<float>>& B_prev_mu_mat, vector<vector<float>>& B_prev_sigma_mat,
                              vector<vector<float>>& B_next_mu_mat, vector<vector<float>>& B_next_sigma_mat, int N_B){
    int i, j;
    float prev_mu, prev_sigma, next_mu, next_sigma;
    for(i=0;i<N_B;i++){
        if(i==0){
            prev_mu = MISSING_BLOCK_MU;
            prev_sigma = MISSING_BLOCK_SIGMA;
        }else{
            prev_mu = B_mu_mat[i - 1][0];
            prev_sigma = B_sigma_mat[i - 1][0];
        }
        for(j=0;j<N_B_MAX_MERGE;j++){
            if(i + j >= N_B){
                continue;
            }
            if(i + j == N_B -1){
                next_mu = MISSING_BLOCK_MU;
                next_sigma = MISSING_BLOCK_SIGMA;
            }else{
                next_mu = B_mu_mat[i + j + 1][0];
                next_sigma = B_sigma_mat[i + j + 1][0];
            }

            B_prev_mu_mat[i][j] = prev_mu;
            B_prev_sigma_mat[i][j] = prev_sigma;
            B_next_mu_mat[i][j] = next_mu;
            B_next_sigma_mat[i][j] = next_sigma;
        }
    }
}


void cal_B_hmm_loglik_cpu(Parameters& global_vars,
                          vector<vector<float>>& B_prev_mu_mat, vector<vector<float>>& B_prev_sigma_mat,
                          vector<vector<float>>& B_next_mu_mat, vector<vector<float>>& B_next_sigma_mat, int debug){


    vector<vector<float>> _mu_mat;
    vector<vector<float>> _sigma_mat;
    vector<vector<float>> _loglike_mat;

    read2DFloatFile("debug_tmp_all_reads_data/B_mu_mat.csv", _mu_mat);
    read2DFloatFile("debug_tmp_all_reads_data/B_sigma_mat.csv", _sigma_mat);
    read2DFloatFile("debug_tmp_all_reads_data/B_hmm_loglik_mat.csv", _loglike_mat);

    for(int i=0;i<B_prev_mu_mat.size();i++){
        for(int j=0;j<N_B_MAX_MERGE;j++){
           global_vars.B_mu_mat[i][j] = _mu_mat[i][j];
           global_vars.B_sigma_mat[i][j] = _sigma_mat[i][j];
           global_vars.B_hmm_loglik_mat[i][j] = _loglike_mat[i][j];
        }
    }

    vector<vector<float>>().swap(_mu_mat);
    vector<vector<float>>().swap(_sigma_mat);
    vector<vector<float>>().swap(_loglike_mat);
}


void find_index(vector<int>& arr, vector<int>& key, vector<int>& index){

    for(auto x: key){
        vector<int>::iterator itr = find(arr.begin(), arr.end(), x);
        index.push_back(distance(arr.begin(), itr));
    }
}


void get_final_res(Parameters& global_vars,
                    vector<int>& final_border_arr,
                    vector<float>& final_mu_arr,
                    vector<float>& final_sigma_arr,
                    vector<int>& final_len_arr,
                    vector<int>& final_b_st_sig_ind_arr,
                    vector<int>& final_b_en_sig_ind_arr){

    int M = global_vars.state_arr.size();
    final_border_arr.push_back(0);

    for(int i=0;i<M-1;i++){
        if(global_vars.state_arr[i] - global_vars.state_arr[i+1] == -1 || global_vars.state_arr[i] - global_vars.state_arr[i+1] == 1){
            final_border_arr.push_back(i+1);
        }
    }
    final_border_arr.push_back(M);

    int N_B = final_border_arr.size()/2;

    for(int i=0;i<N_B;i++){
        final_b_st_sig_ind_arr.push_back(final_border_arr[i*2]);
        final_b_en_sig_ind_arr.push_back(final_border_arr[i*2+1]);
    }

    assert(final_b_st_sig_ind_arr.size() == N_B);
    assert(final_b_en_sig_ind_arr.size() == N_B);

    for(int i=0;i<N_B;i++){
        final_len_arr.push_back(final_b_en_sig_ind_arr[i]- final_b_st_sig_ind_arr[i]);
    }

    vector<int> tmp_i_arr;
    vector<int> tmp_j_arr;
    find_index(global_vars.b_st_sig_ind_arr, final_b_st_sig_ind_arr, tmp_i_arr);
    find_index(global_vars.b_en_sig_ind_arr, final_b_en_sig_ind_arr, tmp_j_arr);


    assert(tmp_i_arr.size() == N_B);
    assert(tmp_i_arr.size() == N_B);

    int tmpj, tmpi;
    for(int i=0;i<N_B;i++){
        tmpi = tmp_i_arr[i];
        tmpj = tmp_j_arr[i];
        assert(tmpj - tmpi < 5);
        final_mu_arr.push_back(global_vars.B_mu_mat[tmpi][tmpj - tmpi]);
        final_sigma_arr.push_back(global_vars.B_sigma_mat[tmpi][tmpj - tmpi]);
    }

}


void gen_region_list(int n_total_b_blocks, int offset, Parameters& global_para, vector<Region>& reg_list){

    assert(offset < N_B_IN_REGION);
    if(offset > 0){
        Region reg = Region(0, offset, global_para);
        reg_list.push_back(reg);
    }

    vector<int> st_arr;
    vector<int> en_arr;

    for(int i=offset;i<n_total_b_blocks;i=i+N_B_IN_REGION){
        st_arr.push_back(i);
    }

    for(int i=1;i<st_arr.size();i++){
        en_arr.push_back(st_arr[i]);
    }
    en_arr.push_back(n_total_b_blocks);

    assert(st_arr.size() == en_arr.size());

    for(int i=0;i<st_arr.size();i++){
        Region reg = Region(st_arr[i], en_arr[i], global_para);
        reg_list.push_back(reg);
    }

    vector<int>().swap(st_arr);
    vector<int>().swap(en_arr);
}


void write_final_float_res_to_csv(string filename, vector<float>& final_arr){

    ofstream myFile;
    myFile.open(filename, ios::app);

    int len = final_arr.size();
    for(int i=0;i<len;i++)
    {
        if(i == len -1){
           myFile << final_arr[i] << "\n";
        }else{
           myFile << final_arr[i] << ",";
        }

    }
    // Close the file
    myFile.close();
}


void write_final_int_res_to_csv(string filename, vector<int>& final_arr, int idx){

    ofstream myFile;
    myFile.open(filename, ios::app);

    myFile << idx << ",";

    int len = final_arr.size();
    for(int i=0;i<len;i++)
    {
        if(i == len -1){
           myFile << final_arr[i] << "\n";
        }else{
           myFile << final_arr[i] << ",";
        }

    }
    myFile.close();
}

void write_final_border_to_csv(string filename, vector<int>& final_border_arr){

    ofstream myFile(filename);

    for(auto b: final_border_arr)
    {
        myFile << b << ",";
    }
    myFile << "\n";
    // Close the file
    myFile.close();
}

#endif // UTIL_H_INCLUDED
