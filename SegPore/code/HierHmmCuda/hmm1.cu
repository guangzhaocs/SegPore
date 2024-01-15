#include <iostream>
#include <vector>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "Parameters.h"
#include "util.h"
#include "BTwindow.h"
#include "SortWindow.h"
#include "Region.h"
#include "calEmission.h"
#include "calAlpha.h"
#include "calBeta.h"
#include "calGamma.h"
#include "updateNormalPara.h"
using namespace std;

#define N_B_MAX_MERGE 5
#define N_B_IN_REGION 11
#define MISSING_BLOCK_MU -2000
#define MISSING_BLOCK_SIGMA 2
#define INNER_HMM_MIN_SIGNAL_LEN 4
#define STATE 4
#define NUM_THREADS 32
#define INNER_HMM_N_ITERATION 20
#define OUTER_HMM_N_ITERATION 2

//---------------------------------------------------------------------------//
// Function
//---------------------------------------------------------------------------//

void cal_B_hmm_loglik_gpu(Parameters& global_vars,
                          vector<vector<float>>& B_prev_mu_mat, vector<vector<float>>& B_prev_sigma_mat,
                          vector<vector<float>>& B_next_mu_mat, vector<vector<float>>& B_next_sigma_mat,
                          int debug){

    //---------------------------------------------------------------------------//
    // CPU Parameters
    //---------------------------------------------------------------------------//

    float  *signal_arr_h = nullptr;

    int *b_st_sig_ind_arr_h = nullptr;
    int *b_en_sig_ind_arr_h = nullptr;

    float *B_mu_mat_h = nullptr;
    float *B_sigma_mat_h = nullptr;
    float *B_prev_mu_mat_h = nullptr;
    float *B_prev_sigma_mat_h = nullptr;
    float *B_next_mu_mat_h = nullptr;
    float *B_next_sigma_mat_h = nullptr;

    float *B_hmm_loglik_arr_h = nullptr;

    //---------------------------------------------------------------------------//
    // GPU Parameters
    //---------------------------------------------------------------------------//

    float  *signal_arr_d = nullptr;

    int *b_st_sig_ind_arr_d = nullptr;
    int *b_en_sig_ind_arr_d = nullptr;

    float *B_mu_mat_d = nullptr;
    float *B_sigma_mat_d = nullptr;
    float *B_prev_mu_mat_d = nullptr;
    float *B_prev_sigma_mat_d = nullptr;
    float *B_next_mu_mat_d = nullptr;
    float *B_next_sigma_mat_d = nullptr;

    float *B_hmm_loglik_arr_d = nullptr;


    float *emis_d = nullptr;    //  BATCH, LEN
    float *alpha_d = nullptr;   //  BATCH, STATE*LEN
    float *beta_d = nullptr;    //  BATCH, STATE*LEN
    float *gamma_d = nullptr;   //  BATCH, STATE*LEN
    float *c_d = nullptr;       //  BATCH, LEN
    int *state_d = nullptr;     //  BATCH, LEN

    size_t pitch_emis;
    size_t pitch_alpha;
    size_t pitch_beta;
    size_t pitch_gamma;
    size_t pitch_c;
    size_t pitch_state;

    int M = global_vars.M;

    // generate sorted window regions
    vector<SortWindow> sort_win_arr;
    generateWinRegion(global_vars.N_B, global_vars.b_st_sig_ind_arr, global_vars.b_en_sig_ind_arr,
                      global_vars.B_mu_mat, global_vars.B_sigma_mat,
                      B_prev_mu_mat, B_prev_sigma_mat, B_next_mu_mat, B_next_sigma_mat,sort_win_arr);

//    for(auto win: sort_win_arr){
//        cout << win.b_st_sig_ind << "-" << win.b_en_sig_ind << "    len : " << win.sig_len << endl;
//    }

    int sort_win_arr_start_idx;
    for(sort_win_arr_start_idx=0;sort_win_arr_start_idx<sort_win_arr.size()-1;sort_win_arr_start_idx++){
        if(sort_win_arr[sort_win_arr_start_idx].sig_len < INNER_HMM_MIN_SIGNAL_LEN &&
           sort_win_arr[sort_win_arr_start_idx+1].sig_len >= INNER_HMM_MIN_SIGNAL_LEN){
            break;
        }
    }
    sort_win_arr_start_idx = sort_win_arr_start_idx + 1;

    int win_num = sort_win_arr.size() - sort_win_arr_start_idx;
	int LEN = sort_win_arr[sort_win_arr.size()-1].sig_len;

	if(debug == 1){
        cout << "    - sort_win_arr_start_idx : " << sort_win_arr_start_idx << endl;
        cout << "    - Region arr length : " << win_num << ", the max len : " << LEN << endl;
    }


    // CPU Data
    signal_arr_h = (float*)malloc(sizeof(float)* M);
    b_st_sig_ind_arr_h = (int*)malloc(sizeof(int)* win_num);
    b_en_sig_ind_arr_h = (int*)malloc(sizeof(int)* win_num);
    B_mu_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_sigma_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_prev_mu_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_prev_sigma_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_next_mu_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_next_sigma_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_hmm_loglik_arr_h = (float*)malloc(sizeof(float)* win_num);

    for (int i = 0; i<M; ++i){
        signal_arr_h[i] = global_vars.signal_arr[i];
    }

    for (int i = sort_win_arr_start_idx; i<sort_win_arr.size(); ++i){
        b_st_sig_ind_arr_h[i-sort_win_arr_start_idx] = sort_win_arr[i].b_st_sig_ind;
        b_en_sig_ind_arr_h[i-sort_win_arr_start_idx] = sort_win_arr[i].b_en_sig_ind;
        B_mu_mat_h[i-sort_win_arr_start_idx] = sort_win_arr[i].B_mu;
        B_sigma_mat_h[i-sort_win_arr_start_idx] = sort_win_arr[i].B_sigma;
        B_prev_mu_mat_h[i-sort_win_arr_start_idx] = sort_win_arr[i].B_prev_mu;
        B_prev_sigma_mat_h[i-sort_win_arr_start_idx] = sort_win_arr[i].B_prev_sigma;
        B_next_mu_mat_h[i-sort_win_arr_start_idx] = sort_win_arr[i].B_next_mu;
        B_next_sigma_mat_h[i-sort_win_arr_start_idx] = sort_win_arr[i].B_next_sigma;
        B_hmm_loglik_arr_h[i-sort_win_arr_start_idx] = -10.0;
    }

//    // for debug
//    emis_h = (float*)malloc(sizeof(float)* win_num * STATE * LEN);
//    c_h = (float*)malloc(sizeof(float)* win_num * LEN);
//    state_h = (int*)malloc(sizeof(int)* win_num * LEN);
//    alpha_h = (float*)malloc(sizeof(float)* win_num * STATE * LEN);
//    beta_h  = (float*)malloc(sizeof(float)* win_num * STATE * LEN);
//    gamma_h = (float*)malloc(sizeof(float)* win_num * STATE * LEN);


    // GPU
	// HMM_Param(M, win_num, LEN);
	// 开辟在GPU中的内存
    cudaMalloc((void**)&signal_arr_d, sizeof(float)* M);

    cudaMalloc((void**)&b_st_sig_ind_arr_d, sizeof(int)* win_num);
    cudaMalloc((void**)&b_en_sig_ind_arr_d, sizeof(int)* win_num);

    cudaMalloc((void**)&B_mu_mat_d, sizeof(float)* win_num);
    cudaMalloc((void**)&B_sigma_mat_d, sizeof(float)* win_num);

    cudaMalloc((void**)&B_prev_mu_mat_d, sizeof(float)* win_num);
    cudaMalloc((void**)&B_prev_sigma_mat_d, sizeof(float)* win_num);

    cudaMalloc((void**)&B_next_mu_mat_d, sizeof(float)* win_num);
    cudaMalloc((void**)&B_next_sigma_mat_d, sizeof(float)* win_num);

    cudaMalloc((void**)&B_hmm_loglik_arr_d, sizeof(float)* win_num);

    cudaMallocPitch((void**)&emis_d, &pitch_emis, sizeof(float)*LEN*STATE, win_num);
    cudaMallocPitch((void**)&alpha_d, &pitch_alpha, sizeof(float)*LEN*STATE, win_num);
    cudaMallocPitch((void**)&c_d, &pitch_c, sizeof(float)*LEN, win_num);
    cudaMallocPitch((void**)&state_d, &pitch_state, sizeof(int)*LEN, win_num);
    cudaMallocPitch((void**)&beta_d, &pitch_beta, sizeof(float)*LEN*STATE, win_num);
    cudaMallocPitch((void**)&gamma_d, &pitch_gamma, sizeof(float)*LEN*STATE, win_num);

	// copy data to GPU
    cudaMemcpy(signal_arr_d, signal_arr_h, sizeof(float)* M, cudaMemcpyHostToDevice);
    cudaMemcpy(b_st_sig_ind_arr_d, b_st_sig_ind_arr_h, sizeof(int)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(b_en_sig_ind_arr_d, b_en_sig_ind_arr_h, sizeof(int)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_mu_mat_d, B_mu_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_sigma_mat_d, B_sigma_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_prev_mu_mat_d, B_prev_mu_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_prev_sigma_mat_d, B_prev_sigma_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_next_mu_mat_d, B_next_mu_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_next_sigma_mat_d, B_next_sigma_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);


    // kernel
    dim3 BlockSize(NUM_THREADS);
    dim3 GridSize((win_num + BlockSize.x - 1) / BlockSize.x);

    // start to run hmm
    for(int j=0;j<INNER_HMM_N_ITERATION;j++){

        //cout << "    - Start to iteration hmm " << j << " ... " << endl;

        calEmission << < GridSize, BlockSize >> > (signal_arr_d, b_st_sig_ind_arr_d, b_en_sig_ind_arr_d,
                                                   B_mu_mat_d, B_sigma_mat_d,
                                                   B_prev_mu_mat_d, B_prev_sigma_mat_d,
                                                   B_next_mu_mat_d, B_next_sigma_mat_d,
                                                   emis_d, pitch_emis/sizeof(float));

//        cudaMemcpy2D(emis_h, sizeof(float)*LEN*STATE, emis_d, pitch_emis, sizeof(float)*LEN*STATE, win_num, cudaMemcpyDeviceToHost);

        cudaDeviceSynchronize();



        calAlpha << < GridSize, BlockSize >> > (signal_arr_d, b_st_sig_ind_arr_d, b_en_sig_ind_arr_d,
                                                emis_d, pitch_emis/sizeof(float),
                                                c_d, pitch_c/sizeof(float),
                                                alpha_d, pitch_alpha/sizeof(float),
                                                B_hmm_loglik_arr_d);

//        cudaMemcpy2D(alpha_h, sizeof(float)*LEN*STATE, alpha_d, pitch_alpha, sizeof(float)*LEN*STATE, win_num, cudaMemcpyDeviceToHost);


        cudaDeviceSynchronize();

        calBeta<< < GridSize, BlockSize >> > (b_st_sig_ind_arr_d, b_en_sig_ind_arr_d,
                                                 emis_d, pitch_emis/sizeof(float),
                                                 c_d, pitch_c/sizeof(float),
                                                 beta_d, pitch_beta/sizeof(float));

//        cudaMemcpy2D(beta_h, sizeof(float)*LEN*STATE, beta_d, pitch_beta, sizeof(float)*LEN*STATE, win_num, cudaMemcpyDeviceToHost);


        cudaDeviceSynchronize();


        calGamma<< < GridSize, BlockSize >> > (b_st_sig_ind_arr_d, b_en_sig_ind_arr_d,
                                                  alpha_d, pitch_alpha/sizeof(float),
                                                  beta_d, pitch_beta/sizeof(float),
                                                  gamma_d, pitch_gamma/sizeof(float));

//        cudaMemcpy2D(gamma_h, sizeof(float)*LEN*STATE, gamma_d, pitch_gamma, sizeof(float)*LEN*STATE, win_num, cudaMemcpyDeviceToHost);


        cudaDeviceSynchronize();

        updateNormalPara<< < GridSize, BlockSize >> > (signal_arr_d, b_st_sig_ind_arr_d, b_en_sig_ind_arr_d,
                                                       B_mu_mat_d, B_sigma_mat_d,
                                                       gamma_d, pitch_gamma/sizeof(float),
                                                       state_d, pitch_state/sizeof(int));

        cudaDeviceSynchronize();
    }

    // wirte data to CPU from GPU
    cudaMemcpy((void*)B_mu_mat_h, (void*)B_mu_mat_d, sizeof(float)* win_num, cudaMemcpyDeviceToHost);
    cudaMemcpy((void*)B_sigma_mat_h, (void*)B_sigma_mat_d, sizeof(float)* win_num, cudaMemcpyDeviceToHost);
    cudaMemcpy((void*)B_hmm_loglik_arr_h, (void*)B_hmm_loglik_arr_d, sizeof(float)* win_num, cudaMemcpyDeviceToHost);


	//  Release cuda
    // Release();
    cudaFree(signal_arr_d);
    cudaFree(b_st_sig_ind_arr_d);
    cudaFree(b_en_sig_ind_arr_d);

    cudaFree(B_mu_mat_d);
    cudaFree(B_sigma_mat_d);

    cudaFree(B_prev_mu_mat_d);
    cudaFree(B_prev_sigma_mat_d);

    cudaFree(B_next_mu_mat_d);
    cudaFree(B_next_sigma_mat_d);

    cudaFree(B_hmm_loglik_arr_d);

    cudaFree(emis_d);
    cudaFree(alpha_d);
    cudaFree(beta_d);
    cudaFree(gamma_d);
    cudaFree(c_d);
    cudaFree(state_d);

    // update global_vars
    int row, col;
    for(int i=0;i<sort_win_arr.size();i++){
        row = sort_win_arr[i].old_window_idx / N_B_MAX_MERGE;
        col = sort_win_arr[i].old_window_idx % N_B_MAX_MERGE;
        if(i < sort_win_arr_start_idx){
           global_vars.B_hmm_loglik_mat[row][col] = -10;
        }else{
           global_vars.B_mu_mat[row][col] = B_mu_mat_h[i-sort_win_arr_start_idx];
           global_vars.B_sigma_mat[row][col] = B_sigma_mat_h[i-sort_win_arr_start_idx];
           global_vars.B_hmm_loglik_mat[row][col] = B_hmm_loglik_arr_h[i-sort_win_arr_start_idx];
        }
    }

    vector<SortWindow>().swap(sort_win_arr);
    //FreeCPU();
}


void get_final_state(Parameters& global_vars,
                     vector<float>& final_mu_arr,
                     vector<float>& final_sigma_arr,
                     vector<int>& final_b_st_sig_ind_arr,
                     vector<int>& final_b_en_sig_ind_arr,
                     vector<int>& final_state_arr){

    //---------------------------------------------------------------------------//
    // CPU Parameters
    //---------------------------------------------------------------------------//

    float  *signal_arr_h = nullptr;

    int *b_st_sig_ind_arr_h = nullptr;
    int *b_en_sig_ind_arr_h = nullptr;

    float *B_mu_mat_h = nullptr;
    float *B_sigma_mat_h = nullptr;
    float *B_prev_mu_mat_h = nullptr;
    float *B_prev_sigma_mat_h = nullptr;
    float *B_next_mu_mat_h = nullptr;
    float *B_next_sigma_mat_h = nullptr;

    float *B_hmm_loglik_arr_h = nullptr;
    int *state_h = nullptr;      //  BATCH, LEN

    //---------------------------------------------------------------------------//
    // GPU Parameters
    //---------------------------------------------------------------------------//

    float  *signal_arr_d = nullptr;

    int *b_st_sig_ind_arr_d = nullptr;
    int *b_en_sig_ind_arr_d = nullptr;

    float *B_mu_mat_d = nullptr;
    float *B_sigma_mat_d = nullptr;
    float *B_prev_mu_mat_d = nullptr;
    float *B_prev_sigma_mat_d = nullptr;
    float *B_next_mu_mat_d = nullptr;
    float *B_next_sigma_mat_d = nullptr;

    float *B_hmm_loglik_arr_d = nullptr;


    float *emis_d = nullptr;    //  BATCH, LEN
    float *alpha_d = nullptr;   //  BATCH, STATE*LEN
    float *beta_d = nullptr;    //  BATCH, STATE*LEN
    float *gamma_d = nullptr;   //  BATCH, STATE*LEN
    float *c_d = nullptr;       //  BATCH, LEN
    int *state_d = nullptr;     //  BATCH, LEN

    size_t pitch_emis;
    size_t pitch_alpha;
    size_t pitch_beta;
    size_t pitch_gamma;
    size_t pitch_c;
    size_t pitch_state;


    int M = global_vars.M;
    int win_num = final_mu_arr.size();

    assert(win_num == final_sigma_arr.size());
    assert(win_num == final_b_st_sig_ind_arr.size());
    assert(win_num == final_b_en_sig_ind_arr.size());
    assert(M == final_state_arr.size());

    // process nan value in mu arr
    for(int i=0;i<win_num;i++){
        if(isnan(final_mu_arr[i])){
            float tmp_mu = 0.0;
            for(int j = final_b_st_sig_ind_arr[i];j < final_b_en_sig_ind_arr[i]; j++){
                tmp_mu = tmp_mu + global_vars.signal_arr[j];
            }
            final_mu_arr[i] = tmp_mu/ (final_b_en_sig_ind_arr[i] - final_b_st_sig_ind_arr[i]);
        }
    }

    // generate sorted window regions
    vector<SortWindow> sort_win_arr;
    for(int i=0;i<win_num;i++){
        SortWindow temp;
        temp.b_st_sig_ind = final_b_st_sig_ind_arr[i];
        temp.b_en_sig_ind = final_b_en_sig_ind_arr[i];
        temp.sig_len = final_b_en_sig_ind_arr[i] - final_b_st_sig_ind_arr[i];

        temp.B_mu = final_mu_arr[i];
        temp.B_sigma = final_sigma_arr[i];

        if(i==0 || i==M-1){
            temp.B_prev_mu = final_mu_arr[i];
            temp.B_prev_sigma = final_sigma_arr[i];
            temp.B_next_mu = final_mu_arr[i];
            temp.B_next_sigma = final_sigma_arr[i];
        }else{
            temp.B_prev_mu = final_mu_arr[i-1];
            temp.B_prev_sigma = final_sigma_arr[i-1];
            temp.B_next_mu = final_mu_arr[i+1];
            temp.B_next_sigma = final_sigma_arr[i+1];
        }
        temp.old_window_idx = i;
        sort_win_arr.push_back(temp);
    }
    sort(sort_win_arr.begin(), sort_win_arr.end(), LessSort);
    assert(win_num == sort_win_arr.size());
    int LEN = sort_win_arr[win_num-1].sig_len;

    // CPU Data
    signal_arr_h = (float*)malloc(sizeof(float)* M);
    b_st_sig_ind_arr_h = (int*)malloc(sizeof(int)* win_num);
    b_en_sig_ind_arr_h = (int*)malloc(sizeof(int)* win_num);
    B_mu_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_sigma_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_prev_mu_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_prev_sigma_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_next_mu_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_next_sigma_mat_h = (float*)malloc(sizeof(float)* win_num);
    B_hmm_loglik_arr_h = (float*)malloc(sizeof(float)* win_num);
    state_h = (int*)malloc(sizeof(int)* win_num * LEN);

    for (int i = 0; i<M; ++i){
        signal_arr_h[i] = global_vars.signal_arr[i];
    }

    for (int i=0; i<win_num; ++i){
        b_st_sig_ind_arr_h[i] = sort_win_arr[i].b_st_sig_ind;
        b_en_sig_ind_arr_h[i] = sort_win_arr[i].b_en_sig_ind;
        B_mu_mat_h[i] = sort_win_arr[i].B_mu;
        B_sigma_mat_h[i] = sort_win_arr[i].B_sigma;
        B_prev_mu_mat_h[i] = sort_win_arr[i].B_prev_mu;
        B_prev_sigma_mat_h[i] = sort_win_arr[i].B_prev_sigma;
        B_next_mu_mat_h[i] = sort_win_arr[i].B_next_mu;
        B_next_sigma_mat_h[i] = sort_win_arr[i].B_next_sigma;
        B_hmm_loglik_arr_h[i] = -10.0;
    }

    // GPU
	// HMM_Param(M, win_num, LEN);
	// 开辟在GPU中的内存
    cudaMalloc((void**)&signal_arr_d, sizeof(float)* M);

    cudaMalloc((void**)&b_st_sig_ind_arr_d, sizeof(int)* win_num);
    cudaMalloc((void**)&b_en_sig_ind_arr_d, sizeof(int)* win_num);

    cudaMalloc((void**)&B_mu_mat_d, sizeof(float)* win_num);
    cudaMalloc((void**)&B_sigma_mat_d, sizeof(float)* win_num);

    cudaMalloc((void**)&B_prev_mu_mat_d, sizeof(float)* win_num);
    cudaMalloc((void**)&B_prev_sigma_mat_d, sizeof(float)* win_num);

    cudaMalloc((void**)&B_next_mu_mat_d, sizeof(float)* win_num);
    cudaMalloc((void**)&B_next_sigma_mat_d, sizeof(float)* win_num);

    cudaMalloc((void**)&B_hmm_loglik_arr_d, sizeof(float)* win_num);

    cudaMallocPitch((void**)&emis_d, &pitch_emis, sizeof(float)*LEN*STATE, win_num);
    cudaMallocPitch((void**)&alpha_d, &pitch_alpha, sizeof(float)*LEN*STATE, win_num);
    cudaMallocPitch((void**)&c_d, &pitch_c, sizeof(float)*LEN, win_num);
    cudaMallocPitch((void**)&state_d, &pitch_state, sizeof(int)*LEN, win_num);
    cudaMallocPitch((void**)&beta_d, &pitch_beta, sizeof(float)*LEN*STATE, win_num);
    cudaMallocPitch((void**)&gamma_d, &pitch_gamma, sizeof(float)*LEN*STATE, win_num);

	// Copy data to GPU
    cudaMemcpy(signal_arr_d, signal_arr_h, sizeof(float)* M, cudaMemcpyHostToDevice);
    cudaMemcpy(b_st_sig_ind_arr_d, b_st_sig_ind_arr_h, sizeof(int)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(b_en_sig_ind_arr_d, b_en_sig_ind_arr_h, sizeof(int)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_mu_mat_d, B_mu_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_sigma_mat_d, B_sigma_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_prev_mu_mat_d, B_prev_mu_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_prev_sigma_mat_d, B_prev_sigma_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_next_mu_mat_d, B_next_mu_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);
    cudaMemcpy(B_next_sigma_mat_d, B_next_sigma_mat_h, sizeof(float)* win_num, cudaMemcpyHostToDevice);

    // kernel
    dim3 BlockSize(NUM_THREADS);
    dim3 GridSize((win_num + BlockSize.x - 1) / BlockSize.x);

    calEmission << < GridSize, BlockSize >> > (signal_arr_d, b_st_sig_ind_arr_d, b_en_sig_ind_arr_d,
                                                   B_mu_mat_d, B_sigma_mat_d,
                                                   B_prev_mu_mat_d, B_prev_sigma_mat_d,
                                                   B_next_mu_mat_d, B_next_sigma_mat_d,
                                                   emis_d, pitch_emis/sizeof(float));
    cudaDeviceSynchronize();

    calAlpha << < GridSize, BlockSize >> > (signal_arr_d, b_st_sig_ind_arr_d, b_en_sig_ind_arr_d,
                                                emis_d, pitch_emis/sizeof(float),
                                                c_d, pitch_c/sizeof(float),
                                                alpha_d, pitch_alpha/sizeof(float),
                                                B_hmm_loglik_arr_d);
    cudaDeviceSynchronize();

    calBeta<< < GridSize, BlockSize >> > (b_st_sig_ind_arr_d, b_en_sig_ind_arr_d,
                                                 emis_d, pitch_emis/sizeof(float),
                                                 c_d, pitch_c/sizeof(float),
                                                 beta_d, pitch_beta/sizeof(float));

    cudaDeviceSynchronize();

    calGamma<< < GridSize, BlockSize >> > (b_st_sig_ind_arr_d, b_en_sig_ind_arr_d,
                                                  alpha_d, pitch_alpha/sizeof(float),
                                                  beta_d, pitch_beta/sizeof(float),
                                                  gamma_d, pitch_gamma/sizeof(float));

    cudaDeviceSynchronize();

    updateNormalPara<< < GridSize, BlockSize >> > (signal_arr_d, b_st_sig_ind_arr_d, b_en_sig_ind_arr_d,
                                                       B_mu_mat_d, B_sigma_mat_d,
                                                       gamma_d, pitch_gamma/sizeof(float),
                                                       state_d, pitch_state/sizeof(int));

    cudaMemcpy2D(state_h, sizeof(float)*LEN, state_d, pitch_state, sizeof(float)*LEN, win_num, cudaMemcpyDeviceToHost);

    // Release GPU cuda
    // Release();
    cudaFree(signal_arr_d);
    cudaFree(b_st_sig_ind_arr_d);
    cudaFree(b_en_sig_ind_arr_d);

    cudaFree(B_mu_mat_d);
    cudaFree(B_sigma_mat_d);

    cudaFree(B_prev_mu_mat_d);
    cudaFree(B_prev_sigma_mat_d);

    cudaFree(B_next_mu_mat_d);
    cudaFree(B_next_sigma_mat_d);

    cudaFree(B_hmm_loglik_arr_d);

    cudaFree(emis_d);
    cudaFree(alpha_d);
    cudaFree(beta_d);
    cudaFree(gamma_d);
    cudaFree(c_d);
    cudaFree(state_d);

    // wirte data to final_state_arr from state_h
    int tmp_win_len = 0;
    int tmp_win_start = 0;
    for(int i=0;i<win_num; i++){
        tmp_win_len = sort_win_arr[i].sig_len;
        tmp_win_start = sort_win_arr[i].b_st_sig_ind;
        for(int j=0;j<tmp_win_len;j++){
            final_state_arr[tmp_win_start + j] = state_h[i*LEN + j];
        }
    }

    vector<SortWindow>().swap(sort_win_arr);
    //FreeCPU();
}


void outer_hmm(int debug, vector<float>& signal_arr, vector<int>& border_arr,
              // string save_mu_dir, string save_sigma_dir, string save_len_dir,
               string save_border_dir, string save_state_dir){

    int M = signal_arr.size();
    int N = border_arr.size() - 1;
    int N_B = border_arr.size()/2;
    int N_T = (border_arr.size() - 2)/2;
    if(debug == 1){
        cout << " * Signal size=" << M << ", N=" << N << ", N_B=" << N_B << ", N_T=" << N_T << endl;
    }

    assert(N < M);
    assert(N >= N_B_IN_REGION);
    assert(N%2 != 0 );
    assert(N == N_B + N_T);
    assert(N_B == 1 + N_T);

    vector<int> state_arr(M, 0);

    vector<int> b_st_sig_ind_arr(N_B, 0);
    vector<int> b_en_sig_ind_arr(N_B, 0);

    vector<int> t_st_sig_ind_arr(N_T, 0);
    vector<int> t_en_sig_ind_arr(N_T, 0);

    vector<vector<float>> B_mu_mat(N_B, vector<float> (N_B_MAX_MERGE, 0.0));
    vector<vector<float>> B_sigma_mat(N_B, vector<float> (N_B_MAX_MERGE, 0.0));
    vector<vector<float>> B_hmm_loglik_mat(N_B, vector<float> (N_B_MAX_MERGE, -100));

    vector<float> T_lr_loglik_arr(N_T, -100);

    // 0-B, 1-T
    vector<vector<float>> trans_mat= {{0.99, 0.01}, {0.1, 0.9}};

    // process data
    gene_b_sig_ind_arr(border_arr, b_st_sig_ind_arr, b_en_sig_ind_arr, N_B);
    gene_t_sig_ind_arr(border_arr, t_st_sig_ind_arr, t_en_sig_ind_arr, N_T);

    init_B_mu_sigma(signal_arr, b_st_sig_ind_arr, b_en_sig_ind_arr, B_mu_mat, B_sigma_mat);
    cal_T_lr_loglik_cpu(signal_arr, t_st_sig_ind_arr, t_en_sig_ind_arr, T_lr_loglik_arr);


    Parameters global_vars = Parameters(M, signal_arr, state_arr, trans_mat,
                 N_B, b_st_sig_ind_arr, b_en_sig_ind_arr, B_mu_mat, B_sigma_mat, B_hmm_loglik_mat,
                 N_T, t_st_sig_ind_arr, t_en_sig_ind_arr, T_lr_loglik_arr);


    vector<Region> region_list;
    vector<Region> overlap_reg_list;

    for(int i_iter=0; i_iter<OUTER_HMM_N_ITERATION; i_iter++){

        region_list.clear();
	    region_list.shrink_to_fit();

        overlap_reg_list.clear();
	    overlap_reg_list.shrink_to_fit();

        if(debug == 1){
           cout << " ----------   Outer HMM Iteration : " << i_iter + 1 << "  ----------  " << endl;
        }

        vector<vector<float>> B_prev_mu_mat(N_B, vector<float> (N_B_MAX_MERGE, MISSING_BLOCK_MU));
        vector<vector<float>> B_prev_sigma_mat(N_B, vector<float> (N_B_MAX_MERGE, MISSING_BLOCK_SIGMA));

        vector<vector<float>> B_next_mu_mat(N_B, vector<float> (N_B_MAX_MERGE, MISSING_BLOCK_MU));
        vector<vector<float>> B_next_sigma_mat(N_B, vector<float> (N_B_MAX_MERGE, MISSING_BLOCK_SIGMA));

        cal_prev_next_B_mu_sigma(global_vars.B_mu_mat, global_vars.B_sigma_mat,
                                 B_prev_mu_mat, B_prev_sigma_mat,
                                 B_next_mu_mat, B_next_sigma_mat, global_vars.N_B);
        if(debug == 1){
           cout << "    - Start to running the inner hmm ...  " << endl;
        }



//        cal_B_hmm_loglik_cpu(global_vars,
//                             B_prev_mu_mat, B_prev_sigma_mat,
//                             B_next_mu_mat, B_next_sigma_mat, debug);

        cal_B_hmm_loglik_gpu(global_vars,
                             B_prev_mu_mat, B_prev_sigma_mat,
                             B_next_mu_mat, B_next_sigma_mat,
                             debug);

        gen_region_list(global_vars.N_B, 0, global_vars, region_list);

        if(debug == 1){
           cout << "    - Start to running the first round ...  " << endl;
        }

        for(int reg_i=0; reg_i < region_list.size(); reg_i++){
             region_list[reg_i].find_best_path();
        }

        if(debug == 1){
           cout << "    - Start to running the second round ...  " << endl;
        }

        int tmpi, tmpj;
        for(int ir=0;ir<region_list.size()-1;ir++){
            Region prev_reg = region_list[ir];
            Region next_reg = region_list[ir + 1];

            tmpi = prev_reg.best_path[prev_reg.best_path.size()-1] / N_B_MAX_MERGE;
            tmpj = prev_reg.best_path[prev_reg.best_path.size()-1] % N_B_MAX_MERGE;
            BTwindow prev_win = prev_reg.all_B_win_mat[tmpi][tmpj];

            tmpi = next_reg.best_path[0] / N_B_MAX_MERGE;
            tmpj = next_reg.best_path[0] % N_B_MAX_MERGE;
            BTwindow next_win = next_reg.all_B_win_mat[tmpi][tmpj];

            if(prev_win.get_b_len() > 1 || next_win.get_b_len() > 1){
                Region reg = Region(prev_win.st_b_ind, next_win.en_b_ind, global_vars);
                reg.find_best_path();
                overlap_reg_list.push_back(reg);
            }
        }

        vector<vector<float>>().swap(B_prev_mu_mat);
        vector<vector<float>>().swap(B_prev_sigma_mat);
        vector<vector<float>>().swap(B_next_mu_mat);
        vector<vector<float>>().swap(B_next_sigma_mat);
    }


    if(debug == 1){
           cout << " ---------   Outer HMM Iteration End !  ---------  " << endl;
    }

    for(auto &reg:region_list){
        reg.set_best_path_state(1);
    }
    for(auto &reg:overlap_reg_list){
        reg.set_best_path_state(0);
    }


    vector<int> final_border_arr;
    vector<float> final_mu_arr;
    vector<float> final_sigma_arr;
    vector<int> final_len_arr;
    vector<int> final_b_st_sig_ind_arr;
    vector<int> final_b_en_sig_ind_arr;

    get_final_res(global_vars, final_border_arr, final_mu_arr, final_sigma_arr, final_len_arr, final_b_st_sig_ind_arr, final_b_en_sig_ind_arr);

    if(debug == 1){
        cout << " * Final border arr=" << final_border_arr.size() << ", final N_B=" << final_mu_arr.size() << endl;
    }

    vector<int> final_state_arr(M, 0);
    get_final_state(global_vars, final_mu_arr, final_sigma_arr, final_b_st_sig_ind_arr, final_b_en_sig_ind_arr, final_state_arr);

//    write_final_float_res_to_csv(save_mu_dir, final_mu_arr);
//    write_final_float_res_to_csv(save_sigma_dir, final_sigma_arr);
//    write_final_int_res_to_csv(save_len_dir, final_len_arr);
    write_final_int_res_to_csv(save_border_dir, final_border_arr);
    write_final_int_res_to_csv(save_state_dir, final_state_arr);

    // deallocating the memory
    vector<int>().swap(final_border_arr);
    vector<int>().swap(final_len_arr);
    vector<int>().swap(state_arr);
    vector<int>().swap(b_st_sig_ind_arr);
    vector<int>().swap(b_en_sig_ind_arr);
    vector<int>().swap(t_st_sig_ind_arr);
    vector<int>().swap(t_en_sig_ind_arr);

    vector<int>().swap(final_state_arr);
    vector<int>().swap(final_b_st_sig_ind_arr);
    vector<int>().swap(final_b_en_sig_ind_arr);

    vector<float>().swap(final_mu_arr);
    vector<float>().swap(final_sigma_arr);
    vector<float>().swap(T_lr_loglik_arr);

    vector<vector<float>>().swap(B_mu_mat);
    vector<vector<float>>().swap(B_sigma_mat);
    vector<vector<float>>().swap(B_hmm_loglik_mat);
    vector<vector<float>>().swap(trans_mat);

    vector<Region>().swap(region_list);
    vector<Region>().swap(overlap_reg_list);

}



int main(int argc, char ** argv)
{
    cout << " * Parameters: " << endl;
    for (int i=0;i<argc;i++)
    {
        cout << argv[i] << endl;
    }
    cout << " " << endl;


    string signal_dir = "/signal.csv";
    string border_dir = "/init_border.csv";
    string save_border_dir = "/res_border.csv";
    string save_state_dir = "/res_state.csv";

    signal_dir = argv[1] + signal_dir;
    border_dir = argv[1] + border_dir;
    save_border_dir = argv[2] + save_border_dir;
    save_state_dir = argv[2] + save_state_dir;

    cout << " * Input dir : " << endl;
    cout << signal_dir << endl;
    cout << border_dir << endl;
    cout << " * Onput dir : " << endl;
    cout << save_border_dir << endl;
    cout << save_state_dir << endl;


    // read signal files and border files
    vector<vector<float> > all_signal_arr;
    vector<vector<int> > all_border_arr;

    read2DFloatFile(signal_dir, all_signal_arr);
    read2DIntFile(border_dir, all_border_arr);

    int n_read = all_signal_arr.size();
    assert(n_read == all_border_arr.size());
    cout << " * Read all signal and init border dataset, the size is " << n_read << endl;


    vector<int> error_int_arr{-1};
    vector<float> error_float_arr{-1.0};

    // run hier_hmm
    int debug = 0;
    int N, M;
    for(int i=0;i<n_read;i++){
        N = all_border_arr[i].size();
        M = all_signal_arr[i].size();
        try {
            if(N%2 == 0){
//               outer_hmm(debug, all_signal_arr[i], all_border_arr[i], save_mu_dir, save_sigma_dir, save_len_dir, save_border_dir, save_state_dir);
               outer_hmm(debug, all_signal_arr[i], all_border_arr[i], save_border_dir, save_state_dir);
               cout << " * Read " << i << " finished ! " << endl;
            }else{
              throw(N);
            }
         }catch(int border_len) {
              vector<int> error_state_arr(M, 0);
              cout << " * Read " << i << " failed !  Border length must be even !"<< endl;
//              write_final_float_res_to_csv(save_mu_dir, error_float_arr);
//              write_final_float_res_to_csv(save_sigma_dir, error_float_arr);
//              write_final_int_res_to_csv(save_len_dir, error_int_arr);
              write_final_int_res_to_csv(save_border_dir, error_int_arr);
              write_final_int_res_to_csv(save_state_dir, error_state_arr);
         }

    }

    return 0;
}
