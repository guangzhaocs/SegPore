#ifndef CALBETA_H_INCLUDED
#define CALBETA_H_INCLUDED
#define STATE 4
using namespace std;


__global__ void calBeta(int* b_st_sig_ind_arr, int* b_en_sig_ind_arr,
                        float* emis, size_t emis_size,
                        float* c, size_t c_size,
                        float* beta, size_t beta_size){

    static float transDev[4][4] = {
               {0.925, 0.025, 0.025, 0.025},
               {0.300, 0.500, 0.100, 0.100},
               {0.300, 0.100, 0.500, 0.100},
               {0.300, 0.100, 0.100, 0.500}};


    // get global index
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int start_idx = b_st_sig_ind_arr[index];
    int end_idx = b_en_sig_ind_arr[index];
    int len = end_idx - start_idx;

    float tmp_beta[STATE];
    float tmp_wiseSum;

    for (int i = 0; i < STATE; i++) {
        beta[index*beta_size + (len-1)*STATE + i] = 1.0;
    }

    // for i in range(n_obs-2, -1, -1):
    for (int i = len-2; i > -1; i--) {

        // beta[:, i] = np.matmul(self.trans_mat, beta[:, i+1] * emis_mat[:, i+1]) / c_arr[i+1]
        for (int j = 0; j < STATE; j++) {
            tmp_beta[j] = beta[index*beta_size + (i+1)*STATE + j] * emis[index*emis_size + (i+1)*STATE + j];
        }

        for (int j = 0; j < STATE; j++) {
            tmp_wiseSum = 0;
            for (int z = 0; z < STATE; z++) {
                tmp_wiseSum = tmp_wiseSum + transDev[j][z] * tmp_beta[z];
            }
            beta[index*beta_size + i*STATE + j] = tmp_wiseSum / c[index*c_size + (i+1)];
        }

    }
}


#endif // CALBETA_H_INCLUDED
