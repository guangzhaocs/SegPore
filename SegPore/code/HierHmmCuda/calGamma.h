#ifndef CALGAMMA_H_INCLUDED
#define CALGAMMA_H_INCLUDED
#define STATE 4
using namespace std;


__global__ void calGamma(int* b_st_sig_ind_arr, int* b_en_sig_ind_arr,
                         float* alpha, size_t alpha_size,
                         float* beta, size_t beta_size,
                         float* gamma, size_t gamma_size){


    // get global index, 取第patch个数据
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int start_idx = b_st_sig_ind_arr[index];
    int end_idx = b_en_sig_ind_arr[index];
    int len = end_idx - start_idx;

    for (int j = 0; j < STATE * len; j++) {
        gamma[index*gamma_size + j] = alpha[index*alpha_size + j] * beta[index*beta_size + j];
    }

    float tmp_wiseSum;

    for (int i = 0; i < len; i++) {

        tmp_wiseSum = 0;
        for (int j = 0; j < STATE; j++) {
            tmp_wiseSum = tmp_wiseSum + gamma[index*gamma_size + i*STATE + j];
        }

        for (int j = 0; j < STATE; j++) {
           gamma[index*gamma_size + i*STATE + j] = gamma[index*gamma_size + i*STATE + j] / tmp_wiseSum;
        }

    }


}


#endif // CALGAMMA_H_INCLUDED
