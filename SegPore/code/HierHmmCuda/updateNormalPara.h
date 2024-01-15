#ifndef UPDATENORMALPARA_H_INCLUDED
#define UPDATENORMALPARA_H_INCLUDED
#define STATE 4


__global__ void updateNormalPara(float* signal, int* start_d, int* end_d,
                                 float* B_mu_mat, float* B_sigma_mat,
                                 float* gamma, size_t gamma_size,
                                 int* state, size_t state_size){

    // get global index, 取第patch个数据
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int start_idx = start_d[index];
    int end_idx = end_d[index];
    int len = end_idx - start_idx;

    // cal state arr
    int tmp_max_index;
    float tmp_max_gamma;
    for (int i = 0; i < len; i++) {
        tmp_max_index = 0;
        tmp_max_gamma = gamma[index*gamma_size + i*STATE + 0];

        for (int j = 1; j < STATE; j++) {
            if (gamma[index*gamma_size + i*STATE + j] > tmp_max_gamma){
                 tmp_max_gamma = gamma[index*gamma_size + i*STATE + j];
                 tmp_max_index = j;
            }
        }
        state[index*state_size + i] = tmp_max_index + 1;
    }


    float sum_signal = 0.0;
    int sub_len = 0;
    float mu = 0.0;
    float sigma = 0.0;
    float sig_min = 1.0;
    float sig_max = 4.0;

    // cal mu
    for (int i = 0; i < len; i++){
        if (state[index*state_size + i] == 1){
           sum_signal += signal[start_idx + i];
           sub_len += 1;
        }
    }

    mu = sum_signal / sub_len;

    // cal sigma
    for (int i = 0; i < len; i++){
        if (state[index*state_size + i] == 1){
            sigma = sigma + (signal[start_idx + i] - mu ) * (signal[start_idx + i] - mu );
        }
    }

    sigma = sigma / sub_len;
    sigma = sqrtf(sigma);


    // update
    B_mu_mat[index] = mu;


    if(sigma > sig_min && sigma < sig_max){
        B_sigma_mat[index] = sigma;
    }else{
        if(sigma <= sig_min){
            B_sigma_mat[index] = sig_min;
        }

        if(sigma >= sig_max){
            B_sigma_mat[index] = sig_max;
        }
    }
}


#endif // UPDATENORMALPARA_H_INCLUDED
