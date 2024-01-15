#ifndef CALALPHA_H_INCLUDED
#define CALALPHA_H_INCLUDED
#define STATE 4
using namespace std;


__device__ float calArraySum(float Array[], int len){

    float numerator = 0;
    for (int i = 0; i < len; i++) {
         numerator = numerator + Array[i];
    }
    return numerator;
}


__global__ void calAlpha(float* signal, int* b_st_sig_ind_arr, int* b_en_sig_ind_arr,
                         float* emis, size_t emis_size,
                         float* c, size_t c_size,
                         float* alpha, size_t alpha_size,
                         float* loglike_d){


    // get global index
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int start_idx = b_st_sig_ind_arr[index];
    int end_idx = b_en_sig_ind_arr[index];
    int len = end_idx - start_idx;

    static const float piDev[4] = {0.85, 0.05, 0.05, 0.05};

    static float transDev[4][4] = {
               {0.925, 0.025, 0.025, 0.025},
               {0.300, 0.500, 0.100, 0.100},
               {0.300, 0.100, 0.500, 0.100},
               {0.300, 0.100, 0.100, 0.500}};

    // init
    // alpha[:, 0] = self.pi_arr * emis_mat[:, 0]
    float tmp_alpha[STATE];
    for (int i = 0; i < STATE; i++) {
        tmp_alpha[i] = emis[index*emis_size + 0*STATE + i] * piDev[i];
    }

    // c_arr[0] = np.sum(alpha[:, 0])
    float tmp_init_sum = calArraySum(tmp_alpha, STATE);
    c[index*c_size + 0] = tmp_init_sum;

    // alpha[:, 0] = alpha[:, 0] / c_arr[0]
    for (int i = 0; i < STATE; i++) {
        alpha[index*alpha_size + 0*STATE + i] = tmp_alpha[i] / tmp_init_sum;
    }


    float tmp_wiseSum;
    // for i in range(1, n_obs):
    for (int i = 1; i < len; i++) {

        // alpha[:, i] = np.matmul(alpha[:, i - 1], self.trans_mat) * emis_mat[:, i]
        for (int j = 0; j < STATE; j++) {
            tmp_wiseSum = 0;
            for (int z = 0; z < STATE; z++) {
                tmp_wiseSum = tmp_wiseSum + alpha[index*alpha_size + (i-1)*STATE + z] * transDev[z][j];
            }
            tmp_alpha[j] = tmp_wiseSum * emis[index*emis_size + i*STATE + j];
        }

        // c_arr[i] = np.sum(alpha[:, i])
        float tmp_sum = calArraySum(tmp_alpha, STATE);
        c[index*c_size + i] = tmp_sum;

        // alpha[:, i] = alpha[:, i] / c_arr[i]
        for (int j = 0; j < STATE; j++) {
           alpha[index*alpha_size + i*STATE + j] = tmp_alpha[j] / tmp_sum;
        }

    }


    // cal loglike
    float tmp_sum = 0;
    int i;
    for (i = 0; i < len; i++) {
        tmp_sum =  tmp_sum + logf(c[index*c_size + i]);
    }

    //loglike_d[index] = tmp_sum / len;
    loglike_d[index] = tmp_sum / i;
}



#endif // CALALPHA_H_INCLUDED
