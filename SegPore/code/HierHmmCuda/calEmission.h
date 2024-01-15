#ifndef CALEMISSION_H_INCLUDED
#define CALEMISSION_H_INCLUDED
#define UNI_MAX 130.0
#define UNI_MIN 50.0
#define STATE 4
using namespace std;


__device__ float normal_pdf(float x, float m, float s){
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a); // To be done: attention!!!
}



__global__ void calEmission(float* signal, int* b_st_sig_ind_arr, int* b_en_sig_ind_arr,
                            float* B_mu_mat, float* B_sigma_mat,
                            float* B_prev_mu_mat, float* B_prev_sigma_mat,
                            float* B_next_mu_mat, float* B_next_sigma_mat,
                            float* emis, size_t emis_size){

    // get global index
    int index = threadIdx.x + blockIdx.x * blockDim.x;

    int start_idx = b_st_sig_ind_arr[index];
    int end_idx = b_en_sig_ind_arr[index];
    int len = end_idx - start_idx;

    float mu_arr[3];
    float sigma_arr[3];
    float tmp;
    mu_arr[0] = B_mu_mat[index];
    sigma_arr[0] = B_sigma_mat[index];

    mu_arr[1] = B_prev_mu_mat[index];
    sigma_arr[1] = B_prev_sigma_mat[index];

    mu_arr[2] = B_next_mu_mat[index];
    sigma_arr[2] = B_next_sigma_mat[index];

    for (int i = 0; i < len; i++) {

        for (int j = 0; j < 3; j++) {
            tmp = normal_pdf(signal[start_idx + i], mu_arr[j], sigma_arr[j]);
            if(tmp > 1.0e-38){
                emis[index*emis_size + i*STATE + j] = tmp;
            }else{
                emis[index*emis_size + i*STATE + j] = 1.0e-38;
            }
        }
        emis[index*emis_size + i*STATE + 3] = 1/(UNI_MAX - UNI_MIN);
    }
}


#endif // CALEMISSION_H_INCLUDED
